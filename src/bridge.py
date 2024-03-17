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
import enum

import consts
import libaux
import libeq

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

        self.jacobian_function = {}
        self.residual_function = {}
        #self.function_names = []
        self.function_size = []
        self.parameter = collections.OrderedDict()  # the values needed to calculate everything
        self.parameter_flag = {}                    # the relation parameter/variable
        self.variable = collections.OrderedDict()   # the parameters that are to be refined and
                                                    # passed to the fitting routined
        self.magnitude = collections.OrderedDict()  # the experimental values to fit upon
        self.variable_values = []
        self.constraint = [None, None, None, None, None, None]
        self.stoichiometry = np.array(model.stoich)
        self.stoichiometryx = np.vstack((np.eye(model.number_components, dtype=int),
                                         self.stoichiometry))
        self.titrations = []
        self.free_concentration: dict = {}   # key is the id of the titration associated
                                             # value is the array of free concentration
        self.analyticalc = {}                # key is the id of the titration associated
        self.jacobian_part = []
        self.jacobian_size = []

        # start collecting information

        # betas are always included first
        parameter['beta'] = list(model.beta)        # it must be mutable
        parameter_flag['beta'] = model.beta_flags   # is should be immutable
        current_size: int = 0
        for v, f in zip(model.beta, model.beta_flags):
            if f == const.RF_CONSTANT:
                pass
            elif f == const.RF_REFINE:
                current_size += 1
                # TODO add value to variables
            elif const.RF_CONSTRAINT1 <= f <= const.RF_CONSTRAINT6:
                nconst = f - const.RF_CONSTRAINT1
                if self.constraint[nconst] is None:
                    current_size += 1
                    self.constraint[nconst] = [v]
                    # TODO add value to variables
                else:
                    self.constraint[nconst].append(v)
            
        self.data_to_titration = {}
        self.data_order = []    # the id of the datawidgets in order of appearance

        for dw in datawidgets:                          # for each datawidget
            self.data_order.append((id(dw), type(dw)))  # store the order of appearance, id and type

            match dw:
                case CalorWidget():
                    raise NotImplementedError
                case EmfWidget():
                    self.parameter[(id(name), 'emf0')] = list(dw.emf0)      # it must be mutable
                    self.magnitudes[id(dw)] = np.array(dw.emf)              # it must be mutable
                    self.parameter_flag[id(dw)] = dw.emf0_flags
                    if any(dw.emf0_flags):
                        self.jacobian_parts.append((id(dw), 'emf0'))
                    self.data_to_titration[id(dw)] = dw.titration_name
                    self.jacobian_function[id(dw)] = self.__emf_jacobian
                    self.residual_function[id(dw)] = self.__emd_residual
                case NmrWidget():
                    raise NotImplementedError
                case SpecWidget():
                    raise NotImplementedError

        for tw in titrationwidgets:
            self.titrations.append(id(tw))
            self.variables[(id(tw, 'init')] = list(tw.initial_amount)  # it must be mutable
            if any(tw.init_flags):
                self.build_jac_parts.append((id(tw), 'init'))
            self.variables[(id(tw, 'buret')] = list(tw.buret)          # it must be mutable
            if any(tw.buret_flags):
                self.build_jac_parts.append((id(tw), 'buret'))
            self.variables[(id(tw, 'v0')] = tw.starting_volume
            self.variables[(id(tw, 'titre')] = tw.titre
            self.variable_flags[(id(tw), 'init')] = tw.init_flags
            self.variable_flags[(id(tw), 'buret')] = tw.buret_flags


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
        def jacobian(values, free_concentration):
            amatrix = libeq.jacobian.amatrix(free_concentration, self.stoichiometryx)
            dlc_dlbeta = libeq.jacobian.dlogcdlogbeta(amatrix, free_concentration, stoichiometry):
            for dataid, datatype in self.data_order:
                # calculate the beta part
                match datatype:
                    case EmfWidget():           # calculate ∂E/∂β
                        libemf.emf_jac_beta(dlc_dlbeta, beta, slope=1.0):
                for jpart, var in self.jacobian_parts:
                    # calculated the rest as needed
                    ...

        return jacobian

    def generate_fobj(self):
        def fobj(values, free_concentrations):
            for name in self.data_order:
                ...

        return fobj

    def generate_freeconcs(self):
        def fconcs(values, init_concs):
            self.update_variables(values)

            betas = self.variables['beta']

            for titration in self.titrations:
                init = np.array(self.variables[titration+'|init'])
                buret = np.array(self.variables[titration+'|buret'])
                v0 = self.variables[titration+'|v0']
                v = np.array(self.variables[titration+'|titre'])
                anc = libaux.build_analyticalc(init, buret, v0, v0)
                analyticalc[titration] = anc

                if init_concs is None
                    c = libeq.consol.initial_guess(beta, self.stoichiometry, analyticalc)
                else:
                    c = libeq.consol.consol(beta, self.stoichiometry, analyticalc, init_concs)
                self.free_concentrations[titration] = c

            return self.free_concentrations

        return fconcs

    def update_variables(self, values):
        # update the values in the widgets
        #   self.model
        #   self.titrationwidgets
        #   self.datawidgets
        ...

    def weights(self):
        ...

    def __betas(values):
        ...

    def __process_emf(self, widget):
        ...

    def __emf_jacobian(self, name:str, amatrix): 
        """Compose the sub-jacobian related to potentiometry.

        The composition should be the following:

         ← constants → ← specific parameters  → ← dangerous parameters →
        +----- β -----+--- ε -+----Δ--+-- ΔH --+-- E₀ --+-- t --+-- b --+
        |             |       |       |        |        |       |       |
        |    ∂E/∂β    |    0  |    0  |    0   |   1    | ∂E/∂t | ∂E/∂b |
        |             |       |       |        |        |       |       |
        +-------------+-------+-------+--------+--------+-------+-------+
        """
        free_concentration = self.free_concentration[name]
        dlogc_dlogbeta = libeq.jacobian.dlogcdlogbeta(amatrix, free_concentration, self.stoichiometry)
        
