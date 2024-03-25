r"""The class that bridges the main app with the fit routines.

Basic rundown.
    1. with titration data → analytical concentration
    2. with analytical concentration + betas + other parameters → free concentration
    3. with free concentration + other parameters → calculated magnitudes
    4a. with calculated magnitudes + free concentration → jacobian
    4b. with calculated + measured magnitudes → residual
    5. with jacobian + residual → parameter update
    6. loop back to 1
"""


import collections

import numpy as np

import consts
import libaux
import libeq
import libemf

from calorwidget import CalorWidget
from emfwidget import EmfWidget
from modelwidget import ModelWidget
from nmrwidget import NmrWidget
from otherwidgets import OutputWidget, ExternalDataWidget, IonicWidget, TitrationBaseWidget
from simulationwidgets import SpeciationWidget, TitrationWidget
from specwidget import SpecWidget
from datawidget import DataWidget


class Bridge():
    """Bridge between the GUI and the fitting engine.
    """
    def __init__(self, model, titrationwidgets, datawidgets):
        self.model = model                          # store these values to
        self.titrationwidgets = titrationwidgets    # update their values after
        self.datawidgets = datawidgets              # the fitting

        # related to the variables
        self.parameter = collections.OrderedDict()  # the values needed to calculate everything
        self.parameter_flag = {}                    # the relation parameter/variable
        self.variables = []  # the parameters that are to be refined and passed to the fitting routined
        self.variable_values = []

        # related to the jacobian
        self.jacobian_part = collections.OrderedDict()

        # related to the residual
        self.residual_function = {}
        self.magnitude = collections.OrderedDict()  # the experimental values to fit upon
        self.magnitude_size = {}

        # other variables
        self.constraint = [None, None, None, None, None, None]
        self.stoichiometry = np.array(model.stoich)
        self.stoichiometryx = np.vstack((np.eye(model.number_components, dtype=int),
                                         self.stoichiometry))
        self.titration = {}
        self.free_concentration: dict = {}   # key is the id of the titration associated
                                             # value is the array of free concentration
        self.analyticalc = {}                # key is the id of the titration associated

        # start collecting information

        # betas are always included first
        self.parameter['beta'] = list(model.beta)        # it must be mutable
        self.parameter_flag['beta'] = model.beta_flags   # is should be immutable
        current_jacobian_size: int = 0
        current_residual_size: int = 0

        current_jacobian_size += self._process_flags('beta', model.beta, model.beta_flags)
        self.jacobian_part['beta'] = slice(0, current_jacobian_size)
        self.data_to_titration = {}
        self.data_order = []    # the id of the datawidgets in order of appearance

        for dw in datawidgets:                          # for each datawidget
            self.data_order.append((id(dw), type(dw)))  # store the order of appearance, id and type
            self.titration[id(dw)] = dw._titrationid

            match dw:
                case CalorWidget():
                    raise NotImplementedError
                case EmfWidget():
                    key = (id(dw), 'emf0')
                    self.parameter[key] = list(dw.emf0)      # it must be mutable
                    self.magnitude[id(dw)] = np.array(dw.emf)              # it must be mutable
                    increment_residual = self.magnitude[id(dw)].size
                    self.magnitude_size[id(dw)] = slice(current_residual_size,
                                                        current_residual_size+increment_residual)
                    current_residual_size += increment_residual
                    self.parameter_flag[id(dw)] = dw.emf0_flags
                    increment = self._process_flags(key ,dw.emf0, dw.emf0_flags)
                    if increment > 0:
                        self.jacobian_part[key] = slice(current_jacobian_size,
                                                        current_jacobian_size+increment)
                        current_jacobian_size += increment
                case NmrWidget():
                    raise NotImplementedError
                case SpecWidget():
                    raise NotImplementedError

        self.__refine_titr = {}

        for tw in titrationwidgets:
            rf_init, incr = self._process_titration_aux((id(tw), 'init'), current_jacobian_size,
                                                        tw.initial_amount, tw.init_flags)
            current_jacobian_size += incr
            rf_buret, incr = self._process_titration_aux((id(tw), 'buret'), current_jacobian_size,
                                                         tw.buret, tw.buret_flags)
            current_jacobian_size += incr
            self.__refine_titr[id(tw)] = rf_init or rf_buret
            self.parameter[(id(tw), 'v0')] = tw.starting_volume
            self.parameter[(id(tw), 'titre')] = tuple(tw.titre)

        self.jacobian = np.empty((current_residual_size, current_jacobian_size), dtype=float)
        self.residual = np.empty(current_residual_size, dtype=float)

    def generate_jacobian(self):
        """The jacobian must be an array of dimmensions (number of titration points, number of
        experimental points per titration point, number of parameters to refine).

         ← constants → ← specific parameters  → ← dangerous parameters →
        +----- β -----+--- ε -+----Δ--+-- ΔH --+-- E₀ --+-- t --+-- b --+
        |             |       |       |        |        |       |       |       ↑         ↑
        |    ∂E/∂β    |    0  |    0  |    0   |   1    | ∂E/∂t | ∂E/∂b | potentiometry   |
        |             |       |       |        |        |       |       |       ↓         |
        +-------------+-------+-------+--------+--------+-------+-------+
        |             |       |       |        |        |       |       |       ↑       number
        |    ∂A/∂β    | ∂A/∂ε |    0  |    0   |   0    | ∂A/∂t | ∂A/∂b | spectrometry   of           
        |             |       |       |        |        |       |       |       ↓      magnitudes   
        +-------------+-------+-------+--------+--------+-------+-------+               times        
        |             |       |       |        |        |       |       |       ↑       number       
        |    ∂Q/∂β    |   0   |    0  | ∂Q/∂ΔH |   0    | ∂Q/∂t | ∂Q/∂b |  calorimetry   of           
        |             |       |       |        |        |       |       |       ↓    observations 
        +-------------+-------+-------+--------+--------+-------+-------+
        |             |       |       |        |        |       |       |       ↑         |
        |    ∂δ/∂β    |   0   | ∂δ/∂Δ |    0   |   0    | ∂δ/∂t | ∂δ/∂b |      NMR        |
        |             |       |       |        |        |       |       |       ↓         ↓
        +-------------+-------+-------+--------+--------+-------+-------+
         ← -------------  number of variables to refine -------------- →
        """
        amatrix = {}
        dlc_dlbeta = {}

        def jacobian(values):
            self.update_parameters(values)

            beta = np.array(self.parameter['beta'])
            for titrid, conc in self.free_concentration.items():
                amatrix[titrid] = libeq.jacobian.amatrix(conc, self.stoichiometryx)
                dlc_dlbeta[titrid] = libeq.jacobian.dlogcdlogbeta(amatrix[titrid], conc, self.stoichiometry)

            breakpoint()
            for dataid, datatype in self.data_order:
                titrid = self.titration[dataid]
                row_slice = self.magnitude_size[dataid]

                match datatype:
                    case EmfWidget:
                        for jpart, col_slice in self.jacobian_part.items():
                            if jpart == "beta":
                                insert = libemf.emf_jac_beta(dlc_dlbeta[titrid], beta, slope=1.0)
                                # TODO remove unrefined columns
                            elif data2id == (titrid, 'emf0'):
                                insert = libemf_emf_jac_e0(_size(row_slice, col_slice))
                            elif data2id == (titrid, 'init'):
                                ...
                            elif data2id == (titrid, 'buret'):
                                ...
                            else:       # zeros
                                insert = np.zeros(_size(row_slice, col_slice))
                            self.jacobian[row_slice, col_slice] = insert

        return jacobian

    def generate_fobj(self):
        def fobj(values, free_concentrations):
            for name in self.data_order:
                ...

        return fobj

    def generate_freeconcs(self):
        def fconcs(values):
            self.update_parameters(values)

            beta = np.array(self.parameter['beta'])

            for titrid in set(self.titration.values()):
                if (titrid in self.__refine_titr) or (titrid not in self.analyticalc):
                    init = np.array(self.parameter[(titrid,'init')])
                    buret = np.array(self.parameter[(titrid,'buret')])
                    v0 = self.parameter[(titrid,'v0')]
                    v = np.array(self.parameter[(titrid,'titre')])
                    anc = libaux.build_analyticalc(init, buret, v0, v)
                    self.analyticalc[titrid] = anc
                else:
                    anc = self.analyticalc[titrid]

                if titrid in self.free_concentration:
                    c = libeq.consol.consol(beta, self.stoichiometry, anc, self.free_concentration[titrid])
                else:
                    c = libeq.consol.initial_guess(beta, self.stoichiometry, anc)
                self.free_concentration[titrid] = c

            return self.free_concentration

        return fconcs

    def update_variables(self, values):
        # update the values in the widgets
        #   self.model
        #   self.titrationwidgets
        #   self.datawidgets
        ...

    def update_parameters(self, values):
        for (key, n), value in zip(self.variables, values):
            self.parameter[key][n] = value

        for item in filter(lambda x: x is not None, self.constraint):
            ref = None
            for key, place, value in item:
                if ref is None:
                    ref = value/self.parameter[key][place]
                else:
                    self.parameter[key][place] *= ref

    def weights(self):
        ...

    def _process_flags(self, key, values, flags):
        current_size = 0
        for n, (v, f) in enumerate(zip(values, flags)):
            if f == consts.RF_CONSTANT:
                pass
            elif f == consts.RF_REFINE:
                current_size += 1
                # TODO add value to variables
                self.variables.append((key, n))
            elif consts.RF_CONSTRAINT1 <= f <= consts.RF_CONSTRAINT6:
                nconst = f - consts.RF_CONSTRAINT1
                constraint_signature = (key, n, v)
                if self.constraint[nconst] is None:
                    current_size += 1
                    self.constraint[nconst] = [constraint_signature]
                    # TODO add value to variables
                else:
                    self.constraint[nconst].append(constraint_signature)
        return current_size

    def _process_titration_aux(self, key, init_size, value, flag):
        # key = (_id, what)
        self.parameter[key] = list(value)  # it must be mutable
        self.parameter_flag[key] = flag
        if any(flag):
            increment = self._process_flags(key, value, flag)
            self.jacobian_part[key] = slice(init_size, init_size+increment)
        return any(flag), increment

    def __betas(self, values):
        ...

    def __process_emf(self, widget):
        ...

    def __emf_jacobian(self, free_concentration:np.ndarray, beta:np.ndarray, amatrix: np.ndarray,
                       block: slice, dlc_dt=None, dlc_db=None): 
        """Compose the sub-jacobian related to potentiometry.

        The composition should be the following:

         ← constants → ← specific parameters  → ← dangerous parameters →
        +----- β -----+--- ε -+----Δ--+-- ΔH --+-- E₀ --+-- t --+-- b --+
        |             |       |       |        |        |       |       |
        |    ∂E/∂β    |    0  |    0  |    0   |   1    | ∂E/∂t | ∂E/∂b |
        |             |       |       |        |        |       |       |
        +-------------+-------+-------+--------+--------+-------+-------+

        Parameters:
           free_concentration(:class:`numpy.ndarray`):
           beta(:class:`numpy.ndarray`):
           amatrix(:class:`numpy.ndarray`):
           block(slice):
           dlc_dt(:class:`numpy.ndarray`):
           dlc_db(:class:`numpy.ndarray`):
        """
        dlc_dlbeta = libeq.jacobian.dlogcdlogbeta(amatrix, free_concentration, self.stoichiometry)

        for kw, val in self.jacobian_part:
            view = self.jacobian[block, val]
            hsize = val.end - val.start
            match kw:
                case "beta":
                    view[...] = libemf.emf_jac_beta(dlc_dlbeta, beta)
                case _id_, "emf0":
                    view[...] = libemf.emf_jac_e0(hsize)
                case _id_, "init":
                    view[...] = libemf.emf_jac_init(dlc_dt)
                case _id_, "buret":
                    view[...] = libemf.emf_jac_buret(dlc_db)
                case _id_, other:
                    vsize = block.end - block.start
                    view[...] = np.zeros((vsize, hsize))


def _size(slice1, slice2):
    return (slice1.stop-slice1.start, slice2.stop-slice2.start)
