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
from nmrwidget import NmrWidget
from specwidget import SpecWidget


class Bridge():
    """Bridge between the GUI and the fitting engine.
    """
    def __init__(self, parameters):
        self.parameters = parameters
        self.stoichiometry = parameters.stoichiometry(extended=False)
        self.stoichiometryx = parameters.stoichiometry(extended=True)

        self.titration = {}
        self.free_concentration: dict = {}   # key is the id of the titration associated
                                             # value is the array of free concentration
        self.analyticalc = {}                # key is the id of the titration associated

        self.jacobian = np.empty(parameters.jacobian_shape, dtype=float)
        self.residual = np.empty(parameters.residual_shape, dtype=float)

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
            self.parameters.update_parameters(values)

            beta, beta_refine = self.parameters.beta()
            for titrid, conc in self.free_concentration.items():
                amatrix[titrid] = libeq.jacobian.amatrix(conc, self.stoichiometryx)
                dlc_dlbeta[titrid] = libeq.jacobian.dlogcdlogbeta(amatrix[titrid], conc, self.stoichiometry)

            for dataid, datatype, titrid, row_slice, jpart, col_slice, data in self.parameters.iter_jacobian():
                if datatype is EmfWidget:
                    if jpart == "beta":
                        # TODO replace slope with data['slope']
                        # remove non electroactive elements columns and non-refinable betas
                        _dlcdlbeta = np.squeeze(dlc_dlbeta[titrid][:,data['electroactive'],:][...,beta_refine])
                        insert = libemf.emf_jac_beta(_dlcdlbeta, beta[beta_refine], slope=1.0)
                    elif jpart == (titrid, 'emf0'):
                        insert = libemf.emf_jac_e0(_size(row_slice, col_slice))
                    elif jpart == (titrid, 'init'):
                        ...
                    elif jpart == (titrid, 'buret'):
                        ...
                    else:       # zeros
                        insert = np.zeros(_size(row_slice, col_slice))
                elif datatype is CalorWidget:
                    raise NotImplementedError
                elif datatype is NmrWidget:
                    raise NotImplementedError
                elif datatype is CalorWidget:
                    raise NotImplementedError
                elif datatype is SpecWidget:
                    raise NotImplementedError

                self.jacobian[row_slice, col_slice] = insert
            return self.jacobian

        return jacobian

    def generate_fobj(self):
        def fobj(values):
            for dataid, datatype, titrid, row_slice in self.parameters.iter_ordered_data():
                ...

        return fobj

    def generate_freeconcs(self):
        def fconcs(values):
            self.parameters.update_parameters(values)

            beta, _ = self.parameters.beta()

            for titrid, to_refine in self.parameters.iter_titration():
                if to_refine or (titrid not in self.analyticalc):
                    init = self.parameters.titr_parm(titrid,'init')
                    buret = self.parameters.titr_parm(titrid,'buret')
                    v0 = self.parameters.titr_parm(titrid,'v0')
                    v = self.parameters.titr_parm(titrid,'titre')
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

    def weights(self):
        ...


class Parameters:
    "Class that handles parameters."
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
        self.titration = {}

        # start collecting information

        # betas are always included first
        self.parameter['beta'] = np.array(model.beta)        # it must be mutable
        self.parameter_flag['beta'] = model.beta_flags   # is should be immutable
        jacobian_slice = Slices()
        residual_slice = Slices()

        self.to_refine_beta = self._process_flags('beta', model.beta, model.beta_flags)
        jacobian_slice.step(len(self.to_refine_beta))
        self.jacobian_part['beta'] = jacobian_slice.yield_slice()

        self.data_to_titration = {}
        self.data_order = []    # the id of the datawidgets in order of appearance

        for dw in datawidgets:                          # for each datawidget
            self.data_order.append((id(dw), type(dw)))  # store the order of appearance, id and type
            self.titration[id(dw)] = dw._titrationid

            match dw:
                case CalorWidget():
                    raise NotImplementedError
                case EmfWidget():
                    self._process_emf(dw, self.jacobian_part, jacobian_slice, residual_slice)
                case NmrWidget():
                    raise NotImplementedError
                case SpecWidget():
                    raise NotImplementedError

        self.__refine_titr = {}

        for tw in titrationwidgets:
            rf_init = self._process_titration_aux((id(tw), 'init'), jacobian_slice,
                                                  tw.initial_amount, tw.init_flags)
            rf_buret = self._process_titration_aux((id(tw), 'buret'), jacobian_slice,
                                                   tw.buret, tw.buret_flags)
            self.__refine_titr[id(tw)] = rf_init or rf_buret
            self.parameter[(id(tw), 'v0')] = tw.starting_volume
            self.parameter[(id(tw), 'titre')] = np.array(tw.titre)

        self.jacobian_shape = (residual_slice.stop, jacobian_slice.stop)
        self.residual_shape = residual_slice.stop

    def beta(self) -> np.ndarray:
        "Return an array with betas."
        return self.parameter['beta'], self.to_refine_beta

    def stoichiometry(self, extended=False):
        "Get stoichiometry array."
        if extended:
            return np.vstack((np.eye(self.model.number_components, dtype=int),
                              np.array(self.model.stoich)))
        else:
            return np.array(self.model.stoich)

    def iter_jacobian(self):
        for dataid, datatype, titrid, row_slice in self.iter_ordered_data():
            for jpart, col_slice in self.iter_jacobian_part():
                if datatype is EmfWidget and jpart == "beta":
                    rkeys = ('slope', 'electroactive')

                data = {k: self.parameter[(dataid, k)] for k in rkeys}
                yield dataid, datatype, titrid, row_slice, jpart, col_slice, data

    def iter_jacobian_part(self):
        yield from self.jacobian_part.items()

    def iter_ordered_data(self):
        "For dataset, yield dataid, datatype and titrationid."
        for dataid, datatype in self.data_order:
            yield dataid, datatype, self.titration[dataid], self.magnitude_size[dataid]

    def iter_titration(self):
        yield from ((x,self.__refine_titr[x]) for x in set(self.titration.values()))

    def titr_parm(self, titrid, parm):
        return self.parameter[(titrid, parm)]

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

    def _process_flags(self, key, values, flags):
        to_refine = []
        for n, (v, f) in enumerate(zip(values, flags)):
            if f == consts.RF_CONSTANT:
                pass
            elif f == consts.RF_REFINE:
                to_refine.append(n)
                # TODO add value to variables
                self.variables.append((key, n))
            elif consts.RF_CONSTRAINT1 <= f <= consts.RF_CONSTRAINT6:
                nconst = f - consts.RF_CONSTRAINT1
                constraint_signature = (key, n, v)
                if self.constraint[nconst] is None:
                    to_refine.append(n)
                    self.constraint[nconst] = [constraint_signature]
                    # TODO add value to variables
                else:
                    self.constraint[nconst].append(constraint_signature)
        return to_refine

    def _process_emf(self, widget: EmfWidget, jpart: collections.OrderedDict,
                     jslice, rslice) -> None:
        self.magnitude[id(widget)] = np.array(widget.emf)     # it must be mutable
        rstep = self.magnitude[id(widget)].size
        rslice.step(rstep)
        self.magnitude_size[id(widget)] = rslice.yield_slice()

        key = (id(widget), 'emf0')
        self.parameter[key] = list(widget.emf0)               # it must be mutable
        self.parameter_flag[id(widget)] = widget.emf0_flags
        increment = len(self._process_flags(key, widget.emf0, widget.emf0_flags))
        if increment > 0:
            jslice.step(increment)
            jpart[key] = jslice.yield_slice()
        
        # TODO replace (1,) ... with a proper call to the parameter in EmfWidget
        self.parameter[(id(widget), 'slope')] = (1.0,) *len(widget.active_species)  # it must be immutable
        self.parameter[(id(widget), 'electroactive')] = widget.active_species  # it must be immutable

    def _process_titration_aux(self, key, jslice, value, flag):
        # key = (_id, what)
        self.parameter[key] = np.array(value)  # it must be mutable
        self.parameter_flag[key] = flag
        if any(flag):
            increment = len(self._process_flags(key, value, flag))
            jslice.step(increment)
            if len(jslice):
                self.jacobian_part[key] = jslice.yield_slice()
        return any(flag)


class Slices:
    "Simple counter that returns a slice object."
    def __init__(self, start=0):
        self.start = start
        self.stop = start

    def step(self, step: int) -> None:
        "Increment stop by given amount"
        self.stop += step

    def yield_slice(self) -> slice:
        "Provide slice to be used."
        retval = slice(self.start, self.stop)
        self.start = self.stop
        return retval

    def __len__(self):
        return self.stop - self.start


def _size(slice1, slice2):
    return (slice1.stop-slice1.start, slice2.stop-slice2.start)
