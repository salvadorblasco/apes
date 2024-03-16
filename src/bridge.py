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
    def __init__(self, model, titrationwidgets, datawidgets):
        self.jacobian_function = {}
        self.residual_function = {}
        #self.function_names = []
        self.function_size = []
        self.variables = collections.OrderedDict()
        self.input_variables = []
        self.variable_shift = {}
        self.variable_flags = {}
        self.variable_values = []
        self.constraints = []
        self.stoichiometry = np.array(model.stoich)
        self.stoichiometryx = np.vstack((np.eye(model.number_components, dtype=int),
                                         self.stoichiometry))
        self.model = model
        self.titrations = []
        self.free_concentration = {}
        self.analyticalc = {}
        self.jacobian_parts = []

        # start collecting information

        variables['beta'] = list(model.beta)        # it must be mutable
        variable_flags['beta'] = model.beta_flags   # is should be immutable
        self.data_to_titration = {}
        self.data_order = []
        for dw in datawidgets:
            if dw.name in self.data_order:
                raise ValueError("Two datasets have the same name")
            self.data_order.append(dw.name)

            match dw:
                case CalorWidget():
                    raise NotImplementedError
                case EmfWidget():
                    self.variables[dw.name + '|E0'] = list(dw.emf0)      # it must be mutable
                    self.variable_flags[dw.name] = dw.emf0_flags
                    if any(dw.emf0_flags):
                        self.jacobian_parts.append(dw.name + '|emf0')
                    self.data_to_titration[dw.name] = dw.titration_name
                    self.jacobian_function[dw.name] = self.__emf_jacobian
                    self.residual_function[dw.name] = self.__emd_residual
                case NmrWidget():
                    raise NotImplementedError
                case SpecWidget():
                    raise NotImplementedError

        for tw in titrationwidgets:
            self.titrations.append(tw.name)
            self.variables[tw.name + '|init'] = list(tw.initial_amount)  # it must be mutable
            self.variables[tw.name + '|buret'] = list(tw.buret)          # it must be mutable
            self.variables[tw.name + '|v0'] = tw.starting_volume
            self.variables[tw.name + '|titre'] = tw.titre
            self.variable_flags[tw.name + '|init'] = tw.init_flags
            if any(tw.init_flags):
                self.build_jac_parts
            self.variable_flags[tw.name + '|buret'] = tw.buret_flags


    def generate_jacobian(self):
        """The jacobian must be an array of dimmensions (number of titration points, number of
        experimental points per titration point, number of parameters to refine).

         ← constants → ← specific parameters  → ← dangerous parameters →
        +----- β -----+--- ε -+----Δ--+-- ΔH --+-- E₀ --+-- t --+-- b --+
        |             |       |       |        |        |       |       |       ↑         ↑
        |    ∂E/∂β    |    0  |    0  |    0   |   1    | ∂E/∂t | ∂E/∂b | potentiometry   |
        |             |       |       |        |        |       |       |       ↓         |
        +-------------+-------+-------+--------+--------+-------+-------+       +
        |             |       |       |        |        |       |       |       ↑       number
        |    ∂A/∂β    | ∂A/∂ε |    0  |    0   |   0    | ∂A/∂t | ∂A/∂b | spectrometry   of           
        |             |       |       |        |        |       |       |       ↓      magnitudes   
        +-------------+-------+-------+--------+--------+-------+-------+       +       times        
        |             |       |       |        |        |       |       |       ↑       number       
        |    ∂Q/∂β    |   0   |    0  | ∂Q/∂ΔH |   0    | ∂Q/∂t | ∂Q/∂b |  calorimetry   of           
        |             |       |       |        |        |       |       |       ↓    observations 
        +-------------+-------+-------+--------+--------+-------+-------+       +         
        |             |       |       |        |        |       |       |       ↑         |
        |    ∂δ/∂β    |   0   | ∂δ/∂Δ |    0   |   0    | ∂δ/∂t | ∂δ/∂b |      NMR        |
        |             |       |       |        |        |       |       |       ↓         ↓
        +-------------+-------+-------+--------+--------+-------+-------+
         ← -------------  number of variables to refine ------------- →

        """
        def jacobian(values, free_concentration):
            amatrix = libeq.jacobian.amatrix(free_concentration, self.stoichiometryx)
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
        
