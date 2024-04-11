r"""The class that bridges the main app with the fit routines.

It basically reads all the data from the widgets that contain the information
to be refined, stores the data in intermediate objects which are used to build
the jacobian and the residual matrices which are required for the refinement
process.

The jacobian must be an array of dimmensions (number of titration points, number of
experimental points per titration point, number of parameters to refine).

                       { JACOBIAN }                                                 { RESIDUALS }
 ← betas → ← specific parameters  → ← dangerous parameters →
+--- β ---+-- ε --+-- Δ --+-- ΔH --+-- E₀ --+-- t --+-- b --+                         +--------+
|         |       |       |        |        |       |       |       ↑         ↑       |        |
|  ∂E/∂β  |    0  |    0  |    0   |   1    | ∂E/∂t | ∂E/∂b | potentiometry   |       | calc E |
|         |       |       |        |        |       |       |       ↓         |       |        |
+---------+-------+-------+--------+--------+-------+-------+                         +--------+
|         |       |       |        |        |       |       |       ↑       number    |        |
|  ∂A/∂β  | ∂A/∂ε |    0  |    0   |   0    | ∂A/∂t | ∂A/∂b | spectrometry   of       | calc A |
|         |       |       |        |        |       |       |       ↓     magnitudes  |        |
+---------+-------+-------+--------+--------+-------+-------+               times     +--------+
|         |       |       |        |        |       |       |       ↑       number    |        |
|  ∂Q/∂β  |   0   |    0  | ∂Q/∂ΔH |   0    | ∂Q/∂t | ∂Q/∂b |  calorimetry   of       | calc Q |
|         |       |       |        |        |       |       |       ↓    observations |        |
+---------+-------+-------+--------+--------+-------+-------+                         +--------+
|         |       |       |        |        |       |       |       ↑         |       |        |
|  ∂δ/∂β  |   0   | ∂δ/∂Δ |    0   |   0    | ∂δ/∂t | ∂δ/∂b |      NMR        |       | calc δ |
|         |       |       |        |        |       |       |       ↓         ↓       |        |
+---------+-------+-------+--------+--------+-------+-------+                         +--------+
 <----------  number of variables to refine --------------->

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
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

import consts
import libaux
import libeq
import libemf

from calorwidget import CalorWidget
from emfwidget import EmfWidget
from nmrwidget import NmrWidget
from specwidget import SpecWidget


class Slices():
    "Simple counter that returns a slice object."
    def __init__(self, start=0):
        self.start = start
        self.stop = start

    def step(self, step: int) -> None:
        "Increment stop by given amount"
        self.stop += step

    def stamp_slice(self, key: str, dict_: dict) -> None:
        if key in dict_:
            raise KeyError(f"key {key} is already in dict")
        if len(self):
            dict_[key] = self.yield_slice()

    def yield_slice(self) -> slice:
        "Provide slice to be used."
        retval = slice(self.start, self.stop)
        self.start = self.stop
        return retval

    def __len__(self):
        return self.stop - self.start


class Bridge():
    """Bridge between the GUI and the fitting engine.
    """
    def __init__(self, parameters):
        self.parameters = parameters
        self.stoichiometry = parameters.stoichiometry(extended=False)
        self.stoichiometryx = parameters.stoichiometry(extended=True)

        self.jacobian = np.empty(parameters.jacobian_shape, dtype=float)
        self.residual = np.empty(parameters.residual_shape, dtype=float)

    def build_matrices(self, values) -> tuple[NDArray[float], NDArray[float]]:
        self.parameters.update_parameters(values)
        betadata = self.parameters.beta
        beta = betadata.beta()
        beta_refine = betadata.to_refine
        betaref = betadata.beta_refine()

        self.update_titrations(beta)

        for dataid, datatype, data in self.parameters.iter_data():
            row_slice = data.vslice
            residual_partial = data.residual()
            self.residual[row_slice] = residual_partial.flat

            for jpart, col_slice in self.parameters.iter_jblock():
                if datatype is EmfWidget:
                    if jpart == "beta":
                        eactiv = data.electroactive
                        full_dl = data.titration.dlcdlbeta
                        # TODO replace slope=1.0 with data['slope']
                        part_dl = np.squeeze(full_dl[:,eactiv,:][...,beta_refine])
                        jac_partial = libemf.emf_jac_beta(part_dl, betaref, slope=1.0)
                    elif jpart == (dataid, 'emf0'):
                        jac_partial = libemf.emf_jac_e0(_size(row_slice, col_slice))
                    elif jpart == (dataid, 'init'):
                        raise NotImplementedError
                    elif jpart == (dataid, 'buret'):
                        raise NotImplementedError
                    else:       # zeros
                        jac_partial = np.zeros(_size(row_slice, col_slice))
                elif datatype is CalorWidget:
                    raise NotImplementedError
                elif datatype is NmrWidget:
                    raise NotImplementedError
                elif datatype is CalorWidget:
                    raise NotImplementedError
                elif datatype is SpecWidget:
                    raise NotImplementedError

                self.jacobian[row_slice, col_slice] = jac_partial
        return self.jacobian, self.residual

    def update_titrations(self, beta) -> None:
        for titration in self.parameters.iter_titrations():
            analc = titration.analc()

            if titration.free_conc is None:
                conc = libeq.consol.initial_guess(beta, self.stoichiometry, analc)
            else:
                conc = libeq.consol.consol(beta, self.stoichiometry, analc, titration.free_conc)
            titration.free_conc = conc
            titration.amatrix = libeq.jacobian.amatrix(conc, self.stoichiometryx)
            titration.dlcdlbeta = libeq.jacobian.dlogcdlogbeta(titration.amatrix, conc,
                                                               self.stoichiometry)

    def weights(self):
        raise NotImplementedError


class Parameters:
    "Class that handles parameters."
    def __init__(self, model, titrationwidgets, datawidgets):
        self.model = model                          # store these values to
        self.titrationwidgets = titrationwidgets    # update their values after
        self.datawidgets = datawidgets              # the fitting

        self.data = {}

        # related to the variables
        # self.parameter = collections.OrderedDict()  # the values needed to calculate everything
        # self.parameter_flag = {}                    # the relation parameter/variable
        self.variables = []  # parameters that are to be refined and passed to the fitting routine

        # related to the jacobian
        self.jacobian_part = collections.OrderedDict()  # this variable is used to reconstruct the
            # jacobian columnwise. The key of this dictionary indicates the identity of the block.
            # It is either a str (like "beta") or a tuple indicatinf the id of the block and the
            # type of data (id, "init"). The value of the dict is the slice of the matrix.

        # related to the residual
        # self.residual_function = {}
        # self.magnitude = collections.OrderedDict()  # the experimental values to fit upon
        # self.magnitude_size = {}

        # other variables
        self.constraint = 6*[None]
        titration_match = {}

        # ~~~~ start collecting information ~~~~
        jacobian_slice = Slices()
        residual_slice = Slices()
        self.data_order = []    # the id of the datawidgets in order of appearance

        # betas are always included first
        self.beta = BetaData(logbeta=np.array(model.beta), beta_flags=model.beta_flags)
        self._process_flags(self.beta, 'logbeta', 'beta_flags', 'to_refine', jacobian_slice)
        jacobian_slice.stamp_slice('beta', self.jacobian_part)

        # datawidgets are processed next
        for dw in datawidgets:                          # for each datawidget
            self.data_order.append((id(dw), type(dw)))  # store the order of appearance, id and type

            match dw:
                case CalorWidget():
                    raise NotImplementedError
                case EmfWidget():
                    data = self._process_emf(dw, jacobian_slice, residual_slice)
                case NmrWidget():
                    raise NotImplementedError
                case SpecWidget():
                    raise NotImplementedError
            self.data[id(dw)] = data
            titration_match[id(dw)] = dw._titrationid

        self.titrations = {id(tw): self._process_titration(tw, jacobian_slice)
                           for tw in titrationwidgets}

        for did, data in self.data.items():
            tid = titration_match[did]
            data.titration = self.titrations[tid]

        self.jacobian_shape = (residual_slice.stop, jacobian_slice.stop)
        self.residual_shape = residual_slice.stop

    def get_temp(self) -> float:
        "Return temperature."
        return self.model.temperature

    def stoichiometry(self, extended=False):
        "Get stoichiometry array."
        if extended:
            return np.vstack((np.eye(self.model.number_components, dtype=int),
                              np.array(self.model.stoich)))
        else:
            return np.array(self.model.stoich)

    def iter_data(self):
        """Iterate over data.

        Yields:
            int: the id of the data begin yielded
            type: the typr of the data being yielded
            data: the data itself
        """
        # return dataid, datatype, data
        for _id, _type in self.data_order:
            yield _id, _type, self.data[_id]

    def iter_jblock(self):
        "Iterate over jacobian column blocks."
        # return jpart, col_slice
        yield from self.jacobian_part.items()

    def update_parameters(self, values):
        "Update values for all variables."
        for variable, value in zip(self.variables, values):
            variable.set_value(value)

    def _process_flags(self, data, attr_values:str, attr_flags:str, attr_toref:str, jslice:Slices):
        to_refine = []
        flags = getattr(data, attr_flags)
        for n, flag in enumerate(flags):
            if flag == consts.RF_CONSTANT:
                pass
            elif flag == consts.RF_REFINE:
                to_refine.append(n)
                variable = Variable(data=data, key=attr_values, position=n)
                self.variables.append(variable)
            elif consts.RF_CONSTRAINT1 <= flag <= consts.RF_CONSTRAINT6:
                nconst = flag - consts.RF_CONSTRAINT1
                if self.constraint[nconst] is None:
                    constraint = Constraint()
                    to_refine.append(n)
                    self.variables.append(constraint)
                    self.constraint[nconst] = constraint
                else:
                    constraint = self.constraint[nconst]
                constraint.append(data, attr_values, n)
        if to_refine:
            jslice.step(len(to_refine))
        setattr(data, attr_toref, to_refine)

    def _process_emf(self, widget, jslice, rslice):
        data = EmfData(
            emf0=np.array(widget.emf0),
            emf0_flags = widget.emf0_flags,
            emf = np.array(widget.emf),
            slope = np.array(widget.slope),
            electroactive = widget.active_species,  # it must be immutable
            temperature = self.get_temp()
        )
        rstep = data.emf.size
        rslice.step(rstep)
        data.vslice = rslice.yield_slice()
        self._process_flags(data, 'emf0', "emf0_flags", "emf0_torefine", jslice)
        key = (id(widget), type(widget))
        jslice.stamp_slice(key, self.jacobian_part)
        return data

    def _process_titration(self, twidget, jslice):
        data = TitrationData(init=np.array(twidget.initial_amount),
                             init_flags=twidget.init_flags,
                             buret=np.array(twidget.buret),
                             buret_flags=twidget.buret_flags,
                             starting_volume=twidget.starting_volume,
                             titre=np.array(twidget.titre))

        self._process_flags(data, 'init', 'init_flags', 'rf_init', jslice)
        key = (id(twidget), 'init')
        jslice.stamp_slice(key, self.jacobian_part)

        self._process_flags(data, 'buret', 'buret_flags', 'rf_buret', jslice)
        key = (id(twidget), 'buret')
        jslice.stamp_slice(key, self.jacobian_part)

        return data

@dataclass
class BetaData:
    "Placeholder for equilibrium constants."
    logbeta: np.ndarray     # the value of log10(β)
    beta_flags: tuple[int]
    to_refine: tuple[int] = field(init=False)

    def beta(self) -> np.ndarray:
        "Return beta values."
        return 10**self.logbeta

    def beta_refine(self):
        return 10**self.logbeta[self.to_refine]


@dataclass
class TitrationData():
    "Placeholder for titration data."
    init: np.ndarray           # values of the initial amount
    buret: np.ndarray          # values of the buret
    titre: np.ndarray
    starting_volume: float
    init_flags: tuple[int]
    buret_flags: tuple[int]
    free_conc: NDArray[float] | None = field(init=False, default=None)
    anal_conc: NDArray[float] = field(init=False)
    amatrix: NDArray[float] = field(init=False)
    dlcdlbeta: NDArray[float] = field(init=False)
    rf_init: tuple[int] = field(init=False)  # indices of elements to be refined
    rf_buret: tuple[int] = field(init=False )# indices of elements to be refined
    refine: bool = field(init=False, default=False)
    init_slice: Slices = field(init=False)
    buret_slice: Slices = field(init=False)

    def __post_init__(self):
        self.refine = any(self.init_flags) or any(self.buret_flags)
        self.__calc_analc()

    def analc(self) -> np.ndarray:
        "Return the analytical concentration."
        if self.refine:
            self.__calc_analc()
        return self.anal_conc

    def __calc_analc(self):
        self.anal_conc = libaux.build_analyticalc(self.init, self.buret,
                                                  self.starting_volume, self.titre)


@dataclass
class EmfData():
    "Placeholder for emf data."
    emf0: NDArray[float]
    emf0_flags: tuple[int]
    emf: NDArray[float]
    slope: NDArray[float]
    titration: TitrationData = field(init=False)
    electroactive: tuple[int]
    temperature: float
    vslice: slice = field(init=False)
    emf0_torefine: tuple[int] = field(init=False)

    def emf_calc(self) -> np.ndarray:
        "Return calculated emf values."
        hconc = libemf.hselect(self.titration.free_conc, self.electroactive)
        return libemf.nernst(hconc, self.emf0, self.slope, 0.0, self.temperature)

    def residual(self) -> np.ndarray:
        "Return residual."
        return self.emf - self.emf_calc()


class Variable:
    "Handling of variables."
    def __init__(self, data, key, position):
        if not hasattr(data, key):
            raise ValueError(f"data {type(data)} does not contain property {key}")
        self.data = data
        self.key = key
        self.position = position

    def get_value(self) -> float:
        "Get the value of the variable."
        dataholder = getattr(self.data, self.key)
        return dataholder[self.position]

    def set_value(self, value: float) -> None:
        """Set value for all parameters in this Constraint.

        Args:
            new_value (float): the new value to include.
        """
        dataholder = getattr(self.data, self.key)
        dataholder[self.position] = value


class Constraint:
    "Transparent handling of constraints."
    def __init__(self):
        self.__values = []
        self.__data = []

    def append(self, data, parname: str, position: int) -> None:
        """Append variable to constraint.

        Args:
            data (EmfData): the data object reference that contains the constraint.
            parname (str): the name of the parameter
            position (int): the index of the parameter
        """
        value = getattr(data, parname)
        self.__data.append((data, parname, position))
        self.__values.append(value)

    def get_value(self) -> float:
        "Get the value of the variable."
        return self.__values[0]

    def set_value(self, new_value: float) -> None:
        """Set value for all parameters in this Constraint.

        Args:
            new_value (float): the new value to include.
        """
        new_values = (new_value*self.__values[0]/val for val in self.__values)
        for (data, parname, position), value in zip(self.__data, new_values):
            dataholder = getattr(data, parname)
            dataholder[position] = value

    def __len__(self):
        return len(self.__values)


def _size(*slices):
    return tuple(len(s) for s in slices)
