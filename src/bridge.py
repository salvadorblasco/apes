r"""The class that bridges the main app with the fit routines.

It basically reads all the data from the widgets that contain the information
to be refined, stores the data in intermediate objects which are used to build
the jacobian and the residual matrices which are required for the refinement
process.

The jacobian must be an array of dimmensions (number of titration points, number of
experimental points per titration point, number of parameters to refine).

                       { JACOBIAN }                                                 { RESIDUALS }
 ← betas → ← specific parameters  → ← dangerous parameters →
+--- β ---+-- Δ --+-- ΔH --+-- ε --+-- E₀ --+-- t --+-- b --+                         +--------+
|         |       |        |       |        |       |       |       ↑         ↑       |        |
|  ∂E/∂β  |    0  |    0   |    0  |   1    | ∂E/∂t | ∂E/∂b | potentiometry   |       | calc E |
|         |       |        |       |        |       |       |       ↓         |       |        |
+---------+-------+--------+-------+--------+-------+-------+                         +--------+
|         |       |        |       |        |       |       |       ↑       number    |        |
|  ∂A/∂β  |    0  |    0   | ∂A/∂ε |   0    | ∂A/∂t | ∂A/∂b | spectrometry   of       | calc A |
|         |       |        |       |        |       |       |       ↓     magnitudes  |        |
+---------+-------+--------+-------+--------+-------+-------+               times     +--------+
|         |       |        |       |        |       |       |       ↑       number    |        |
|  ∂Q/∂β  |    0  | ∂Q/∂ΔH |   0   |   0    | ∂Q/∂t | ∂Q/∂b |  calorimetry   of       | calc Q |
|         |       |        |       |        |       |       |       ↓    observations |        |
+---------+-------+--------+-------+--------+-------+-------+                         +--------+
|         |       |        |       |        |       |       |       ↑         |       |        |
|  ∂δ/∂β  | ∂δ/∂Δ |    0   |   0   |   0    | ∂δ/∂t | ∂δ/∂b |      NMR        |       | calc δ |
|         |       |        |       |        |       |       |       ↓         ↓       |        |
+---------+-------+--------+-------+--------+-------+-------+                         +--------+
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
import math
import typing
from functools import reduce, partial

import numpy as np
from numpy.typing import NDArray
import scipy

import consts
import libaux
import libeq
import libemf
import libspec

from modelwidget import ModelWidget
from calorwidget import CalorWidget
from emfwidget import EmfWidget
from nmrwidget import NmrWidget
from specwidget import SpecWidget
from otherwidgets import TitrationBaseWidget


prod_tuple = partial(reduce, lambda x, y: x*y)


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
    def __init__(self, parameters, report_buffer=None):
        self.parameters = parameters
        self.stoichiometry = parameters.stoichiometry(extended=False)
        self.stoichiometryx = parameters.stoichiometry(extended=True)

        self.report_buffer = report_buffer

        self.jacobian = np.empty(parameters.jacobian_shape, dtype=float)
        self.residual = np.empty(parameters.residual_shape, dtype=float)

    def accept_values(self):
        self.parameters.accept_values()

    def size(self) -> tuple[int]:
        return self.parameters.jacobian_shape

    def step_values(self, increments):
        self.parameters.increment_variables(increments)

    def build_matrices(self) -> tuple[NDArray[float], NDArray[float]]:
        temp = self.parameters.get_temp()
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
                        part_dl = full_dl[:,eactiv,:][...,beta_refine]
                        flatten_shape = (prod_tuple(part_dl.shape[:-1]), part_dl.shape[-1])
                        flat_dl = np.reshape(part_dl, flatten_shape)
                        # TODO replace slope=1.0 with data['slope']
                        jac_partial = libemf.emf_jac_beta(flat_dl, slope=1.0, temperature=temp)
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

    def report_step(self, **kwargs):
        if self.report_buffer is None:
            return
        write = self.report_buffer.write

        if 'iteration' in kwargs:
            write(f"  \niteration {kwargs['iteration']}  \n")
            write(f" chi-squared={kwargs['chisq']:.3e}; sigma={kwargs['sigma']:.3e}  \n")

        write(f" NUM  VALUE    STEP  \n")
        for n, v in enumerate(self.parameters.beta.variables):
            vvalue = v.get_value() / consts.LOGK
            vstep = v.last_increment / consts.LOGK
            write(f" {n:4}  {vvalue:7.4f}  {vstep:7.4f}  \n")

        # print refined betas


    def update_titrations(self, beta) -> None:
        for titration in self.parameters.titrations.values():
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
        # TODO change this to real values
        return np.ones_like(self.residual)


class Parameters:
    "Class that handles parameters."
    def __init__(self, model, titrationwidgets, datawidgets):
        self.model = model                          # store these values to
        self.titrationwidgets = titrationwidgets    # update their values after
        self.datawidgets = datawidgets              # the fitting

        self.data = {}

        # related to the variables
        self.variables: list[FreeVariable] = []       # parameters that are to be refined and
                                                      # passed to the fitting routine

        # related to the jacobian
        self.jacobian_part = collections.OrderedDict()  # this variable is used to reconstruct the
            # jacobian columnwise. The key of this dictionary indicates the identity of the block.
            # It is either a str (like "beta") or a tuple indicatinf the id of the block and the
            # type of data (id, "init"). The value of the dict is the slice of the matrix.

        # other variables
        self.constraint = 6*[None]
        titration_match = {}
        self.spectraldata = None

        # ~~~~ start collecting information ~~~~
        jacobian_slice = Slices()
        residual_slice = Slices()
        self.data_order = []    # the id of the datawidgets in order of appearance

        # betas are always included first
        # BEWARE! the widget provides the values of the constants as LOG10
        # from this point on all calculations are done in the LN of the constants
        self.beta = BetaData(logbeta=consts.LOGK*np.array(model.beta_raw), beta_flags=model.beta_flags)
        self._process_flags(self.beta, 'logbeta', 'beta_flags', 'to_refine', jacobian_slice)
        jacobian_slice.stamp_slice('beta', self.jacobian_part)
        self.beta.variables = self.variables.copy()

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
                    if self.spectraldata is None:
                        self.spectraldata = SpectralData()
                    data = self._process_spectrum(dw, jacobian_slice, residual_slice)
                case _:
                    raise ValueError(f"Widget type {type(dw)} not recognised")

            self.data[id(dw)] = data
            titration_match[id(dw)] = id(dw.titration)

        self.titrations = {id(tw): self._process_titration(tw, jacobian_slice)
                           for tw in titrationwidgets}

        for did, data in self.data.items():
            tid = titration_match[did]
            data.titration = self.titrations[tid]

        self.jacobian_shape = (residual_slice.stop, jacobian_slice.stop)
        self.residual_shape = residual_slice.stop

    def dump_to_widgets(self):
        "Dump data used in refinement to the corresponding widgets."
        self.beta.dump(self.model)

        for titwidget in self.titrationwidgets:
            titdata = self.titrations[id(titwidget)]
            titdata.dump(titwidget)

        for datwidget in self.datawidgets:
            data = self.data[id(datwidget)]
            data.dump(datwidget)

    def get_temp(self) -> float:
        "Return temperature."
        return self.model.temperature

    def initial_values(self):
        for var in self.variables:
            yield var.get_value()

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

    def accept_values(self):
        for variable in self.variables:
            variable.accept_value()

    def increment_variables(self, increments):
        "Increment variable values."
        for variable, value in zip(self.variables, increments):
            variable.increment_value(value)

    def set_errors(self, values):
        "Set the variable errors"
        for variable, value in zip(self.variables, values):
            variable.set_error(value)

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

    def _process_spectrum(self, dw, jacobian_slice, residual_slice):
        data = SpectrumData(absorbance=np.array(dw.spectrum),
                            wavelength=np.fromiter(dw.wavelengths, dtype=float),
                            optical_path=dw.optical_path)
        self.spectraldata.add_spectrum(dw.optically_active, data)
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
    amatrix: NDArray[float] = field(init=False, repr=False)
    dlcdlbeta: NDArray[float] = field(init=False, repr=False)
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

    def dump(self, widget: TitrationBaseWidget) -> None:
        "Dump data into the widget to update the GUI."
        widget.free_conc = self.free_conc
        if self.refine:
            return
        if any(self.init_flags):
            widget.initial_amount = self.init
        if any(self.buret_flags):
            widget.buret = self.buret

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
    titration: TitrationData = field(init=False, repr=False)
    electroactive: tuple[int]
    temperature: float
    vslice: slice = field(init=False)
    emf0_torefine: tuple[int] = field(init=False)

    def dump(self, widget: EmfWidget) -> None:
        "Dump data into the widget to update the GUI."
        widget.emffitted = self.emf_calc()
        if any(self.emf0_flags):
            widget.emf0 = self.emf0

    def emf_calc(self) -> np.ndarray:
        "Return calculated emf values."
        hconc = libemf.hselect(self.titration.free_conc, self.electroactive)
        return libemf.nernst(hconc, self.emf0, self.slope, 0.0, self.temperature)

    def residual(self) -> np.ndarray:
        "Return residual."
        return self.emf - self.emf_calc()


@dataclass
class BetaData:
    "Placeholder for equilibrium constants."
    logbeta: np.ndarray     # the value of log(β)
    beta_flags: tuple[int]
    to_refine: tuple[int] = field(init=False)
    errors: NDArray[float] = field(init=False)
    variables: list["FreeVariable"] = field(init=False)

    def __post_init__(self):
        self.errors = np.zeros_like(self.logbeta)

    def beta(self) -> np.ndarray:
        "Return beta values."
        return np.exp(self.logbeta)

    def beta_refine(self):
        return np.exp(self.logbeta[self.to_refine])

    def dump(self, widget: ModelWidget) -> None:
        "Dump data into the widget to update the GUI."
        if any(self.beta_flags):
            widget.beta_raw = self.logbeta/consts.LOGK      # dump value as LOG10
            widget.beta_error = self.errors/consts.LOGK


class SpectralData:
    def __init__(self):
        self.optically_active = None
        self.xdata = np.empty(0)
        self.ydata = np.empty(0)
        self.spanning = 5.0  # nanometres
        self.cubicspline = None

    def add_spectrum(self, optically_active: tuple[bool], spectrum: "SpectrumData") -> None:
        self._assert_active_consistency(optically_active)
        new_xdata = spectrum.wavelength
        new_ydata = self._getestimation(spectrum)
        self.xdata = np.append(self.xdata, new_xdata)
        self.ydata = np.append(self.ydata, new_ydata)

    def setup(self):
        self.cubicspline = scipy.interpolate.CubicSpline(self.xdata, self.ydata, axis=1)
        working_wavelength = np.arange(np.min(self.xdata), np.max(self.xdata), self.spanning)
        self.data = self.interpolate(working_wavelength)

    def interpolate(self, wavelength):
        return self.cubicspline(wavelength)

    def _assert_active_consistency(self, optically_active: tuple[bool]) -> None:
        if self.optically_active is None:
            self.optically_active = optically_active
        else:
            if not all(i==j for i, j in zip(self.optically_active, optically_active)):
                raise ValueError("optically active values are not consistent")

    def _getestimation(self, spectrum):
        activ_free_conc = spectrum.titration.free_conc[:,self.optically_active]
        return spectrum.optical_path * spectrum.absorbance @ np.pinv(activ_free_conc.T)


@dataclass
class SpectrumData:
    absorbance: NDArray[float]  # axis0->wavelength, axis1->titration point
    wavelength: NDArray[float]  # 1D array
    optical_path: float
    titration: TitrationData = field(init=False)
    epsilon: SpectralData = field(init=False)

    def dump(self, widget: SpecWidget) -> None:
        "Dump data into the widget to update the GUI."
        raise NotImplementedError

    def residual(self) -> np.ndarray:
        "Return residual."
        free_conc = self.titration.free_conc
        optical_activity = self.epsilon.interpolate(self.wavelength)
        return self.absorbance - libspec.spec_function(free_conc, optical_activity, self.optical_path)


class FreeVariable(typing.Protocol):
    max_increment: float 
    stored_value: float | None
    error: float

    def accept_value(self):
        ...

    def increment_value(self, increment):
        ...

    def get_value(self) -> float:
        ...

    def set_error(self, value: float) -> None:
        ...

    def set_value(self, value: float) -> None:
        ...


class Variable:
    "Handling of variables."
    def __init__(self, data, key, position):
        if not hasattr(data, key):
            raise ValueError(f"data {type(data)} does not contain property {key}")
        self.data: float = data
        self.key = key
        self.keyerror = "errors"
        self.position: int = position
        self.max_increment: float = math.inf
        self.stored_value: float | None = None
        self.last_increment: float

    def accept_value(self):
        self.stored_value = None

    def increment_value(self, increment):
        if self.stored_value is None:
            self.stored_value = self.get_value() 
        if abs(increment) > self.max_increment:
            increment = math.copysign(self.max_increment, increment)
        self.last_increment = increment
        new_value = self.stored_value + increment
        self.set_value(new_value)

    def get_value(self) -> float:
        "Get the value of the variable."
        dataholder = getattr(self.data, self.key)
        return dataholder[self.position]

    def get_error(self) -> float:
        dataholder = getattr(self.data, self.keyerror)
        return dataholder[self.position]

    def set_error(self, value: float) -> None:
        dataholder = getattr(self.data, self.keyerror)
        dataholder[self.position] = value

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
        self.max_increment = math.inf
        self.stored_value: float | None = None
        self.last_increment: float

    def accept_value(self):
        "accept and store current value."
        self.stored_value = None

    def increment_value(self, increment):
        "increment variable value."
        if self.stored_value is None:
            self.stored_value = self.get_value()
        if abs(increment) > self.max_increment:
            increment = math.copysign(self.max_increment, increment)
        self.last_increment = increment
        new_value = self.stored_value + increment
        self.set_value(new_value)

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

    def set_error(self, value: float) -> None:
        new_values = (value*self.__values[0]/val for val in self.__values)
        for (data, parname, position), value in zip(self.__data, new_values):
            dataholder = getattr(data, "errors")
            dataholder[position] = value

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


def trivial_capping(x, dx):
    "Capping function where there is no capping"
    return x + dx


def max_ratio_capping(x, dx, ratio):
    "Capping to a fraction of change"
    if (aux := np.abs(dx)/x) < ratio:
        return x+dx
    else:
        return x*(1+aux)


def abs_capping(x, dx, maximum):
    if abs(dx) < maximum:
        return x + dx
    else:
        return x + maximum
