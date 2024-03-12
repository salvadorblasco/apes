
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
    def __init__(self, model, *datawidgets):
        self.function = []
        self.function_names = []
        self.function_size = []
        self.variables = []
        self.variable_size = []
        self.variable_flags = []
        self.variable_values = []
        self.stoichiometry = np.array(model.stoich)
        self.stoichiometryx = np.vstack((np.eye(model.number_components, dtype=int),
                                         self.stoichiometry))
        self.model = model
        self.free_concentrations = None

        self._fjacobian = {
            EmfWidget: {'beta': libemf.jac_beta, 'init': libemf.jac_init}
        }
        self._fobj = {
            EmfWidget: libemf.fobj
        }

        for dw in datawidgets:
            vnames, vsizes, vflags, values = dw.variables()
            self.variables.extend(vnames)
            self.variable_size.extend(vsizes)
            self.variable_flags.extend(vflags)
            self.variable_values.extend(values)

            fnames, fsize = dw.functions()
            self.functions.extend(
            self.function_names.extend(fnames)
            self.function_size.extend(fsize)

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
            betas = self.__betas(values)
            if init_concs is None
                c = libeq.consol.initial_guess(beta, stoichiometry, analyticalc)
            else:
                c = libeq.consol.consol(beta, stoichiometry, analyticalc, init_concs)
            return c

        return fconcs

    def weights(self):
        ...

    def __betas(values):
        ...
