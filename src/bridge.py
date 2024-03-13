import collections

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
        self.function = []
        self.function_names = []
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

        self._fjacobian = {
            EmfWidget: {'beta': libemf.jac_beta, 'init': libemf.jac_init}
        }
        self._fobj = {
            EmfWidget: libemf.fobj
        }

        # start collecting information

        variables['beta'] = list(model.beta)        # it must be mutable
        variable_flags['beta'] = model.beta_flags   # is should be immutable

        for dw in datawidgets:
            match dw:
                case CalorWidget():
                    raise NotImplementedError
                case EmfWidget():
                    self.__process_emf(dw)
                case NmrWidget():
                    raise NotImplementedError
                case SpecWidget():
                    raise NotImplementedError

        for tw in titrationwidgets:
            self.titrations.append(tw.name)
            variables[tw.name + '|init'] = list(tw.initial_amount)  # it must be mutable
            variables[tw.name + '|buret'] = list(tw.buret)          # it must be mutable
            variables[tw.name + '|v0'] = tw.starting_volume
            variables[tw.name + '|titre'] = tw.titre
            variable_flags[tw.name + '|init'] = tw.init_flags
            variable_flags[tw.name + '|buret'] = tw.buret_flags

    def generate_jacobian(self):
        def jacobian(values, free_concentration):
            amatrix = libeq.jacobian.amatrix(free_concentration, self.stoichiometryx)
            ...

        return jacobian

    def generate_fobj(self):
        def fobj(values, free_concentrations):
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
