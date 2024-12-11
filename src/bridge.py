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
from typing import Sequence, Dict
from functools import reduce, partial, cache
import itertools

import numpy as np
from numpy.typing import NDArray
import scipy

import consts
import libaux
import libeq
import libemf
import libspec

from datawidget import DataWidget
from modelwidget import ModelWidget
from calorwidget import CalorWidget
from emfwidget import EmfWidget
from nmrwidget import NmrWidget
from specwidget import SpecWidget
from otherwidgets import TitrationBaseWidget


prod_tuple = partial(reduce, lambda x, y: x*y)


class Slices():
    """Utility class to generate and manage slice objects.

    This class provides an incrementing slice, allowing for controlled and organized
    slicing of arrays, especially useful in matrix construction (e.g., Jacobian assembly).

    Attributes:
        start (int): Starting index for the slice.
        stop (int): Stopping index, dynamically updated.

    Methods:
        step(step: int): Increments the stop index by a given amount.
        stamp_slice(key: str, dict_: dict): Stores the current slice in the provided dictionary.
        yield_slice() -> slice: Returns a slice object from start to stop.
        __len__() -> int: Returns the length of the current slice.

    Example:
        >>> slices = Slices(0)
        >>> slices.step(5)
        >>> print(slices.yield_slice())  # slice(0, 5)
        >>> slices.step(3)
        >>> print(slices.yield_slice())  # slice(5, 8)
    """

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

    This class serves as a central point to connect GUI data inputs (e.g., titration,
    widget data) with the underlying fitting engine. It manages the intermediate
    calculations necessary for fitting, including building Jacobian and residual matrices,
    updating parameters, and refining estimates.
    """

    def __init__(self, parameters: "Parameters", report_buffer: "StringIO|None"=None):
        self.parameters = parameters
        self.stoichiometry = parameters.stoichiometry(extended=False)
        self.stoichiometryx = parameters.stoichiometry(extended=True)

        self.report_buffer = report_buffer

        self.jacobian = np.empty(parameters.jacobian_shape, dtype=float)
        self.residual = np.empty(parameters.residual_shape, dtype=float)

        self.chisq_hist: list[float] = []   # record fitting parameter history
        self.sigma_hist: list[float] = []

    def accept_values(self):
        "Change definetly the values of the variables."
        self.parameters.accept_values()

    def size(self) -> tuple[int]:
        return self.parameters.jacobian_shape

    def iteration_history(self, chisq:float, sigma:float) -> None:
        self.chisq_hist.append(chisq)
        self.sigma_hist.append(sigma)

    def step_values(self, increments):
        "Change provisionally the values of the variables."
        self.parameters.increment_variables(increments)

    def build_matrices(self) -> tuple[NDArray[float], NDArray[float]]:
        """Builds the Jacobian and residual matrices required for refinement.

        This method iterates over available data to calculate partial derivatives and residuals
        for each data type (e.g., EmfWidget, CalorWidget), storing results in pre-allocated
        matrices. Each data type's specific processing function is called to compute the 
        respective Jacobian and residual components.

        Returns:
            tuple[NDArray[float], NDArray[float]]: The computed Jacobian and residual matrices.
        """
        temp = self.parameters.get_temp()
        betadata = self.parameters.beta
        beta = betadata.beta()
        beta_refine = betadata.to_refine
        betaref = betadata.beta_refine()

        self.update_titrations(beta)

        data_methods = {
            EmfWidget: self._process_emf_data,
            CalorWidget: self.__not_implemented,
            NmrWidget: self.__not_implemented,
            SpecWidget: self.__not_implemented
        }

        for dataid, datatype, data in self.parameters.iter_data():
            row_slice = data.vslice
            residual_partial = data.residual()
            self.residual[row_slice] = residual_partial.flat

            for jpart, col_slice in self.parameters.iter_jblock():
                if datatype in data_methods:
                    jac_partial = data_methods[datatype](dataid, data, jpart, beta_refine, row_slice, col_slice, temp)
                else:
                    jac_partial = np.zeros(_size(row_slice, col_slice))
                
                self.jacobian[row_slice, col_slice] = jac_partial
        return self.jacobian, self.residual

    def _process_emf_data(self, dataid, data, jpart, beta_refine, row_slice, col_slice, temp):
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
        return jac_partial

    def __not_implemented(self):
        raise NotImplementedError

    def report_step(self, **kwargs):
        self.iteration_history(chisq=kwargs['chisq'], sigma=kwargs['sigma'])
        if self.report_buffer is None:
            return
        write = self.report_buffer.write

        if 'iteration' in kwargs:
            write(f"  \niteration {kwargs['iteration']}  \n")
            write(f" chi-squared={kwargs['chisq']:.3e}; sigma={kwargs['sigma']:.3e}  \n")

        write(" NUM  OLD VALUE    STEP   NEW VALUE  \n")
        write(" ---  ---------  -------  ---------  \n")
        for n, v in enumerate(self.parameters.beta.variables):
            ovalue = v.previous_value / consts.LOGK
            nvalue = v.get_value() / consts.LOGK
            vstep = v.last_increment / consts.LOGK
            write(f"{n+1: 4d}  {ovalue:8.4f}  {vstep:8.4f}  {nvalue:8.4f}  \n")

    def report_raw(self, txt):
        if self.report_buffer is None:
            return
        self.report_buffer.write(txt)


    def update_titrations(self, beta: np.ndarray) -> None:
        """Update titration data by recalculating free concentrations and analytical matrices.

        For each titration instance, calculates new free concentrations and updates the Jacobian
        with respect to the provided `beta` values. If `free_conc` is uninitialized, an initial
        guess is computed, otherwise, it refines the existing concentration values.

        Args:
            beta (np.ndarray): Array of current beta values used to update titration calculations.

        Notes:
            This function assumes that the stoichiometry has been initialized and that `free_conc`
            has either initial guesses or prior calculated values.
        """
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

    @cache
    def weights(self):
        """Calculate the weight matrix.

        .. note :: It is not necessary to call :func:`update_titrations` first because no free
        concentrations are used in this function.
        """
        # TODO change this to real values
        #return np.ones_like(self.residual)
        w = np.zeros_like(self.residual)
        wdiag = np.zeros_like(self.residual)
        for dataid, datatype, data in self.parameters.iter_data():
            row_slice = data.vslice
            weight_partial = data.weight()
            w[row_slice] = weight_partial.flat
            wdiag[row_slice] = data.variance().flat

        covar = np.diag(wdiag) + w[:,None] @ w[None, :]
        return 1/covar


class Parameters:
    """Class that handles parameters and creates a bond between the dataclasses and regular Python
    data structures.
    """
    def __init__(self,
                 model: ModelWidget,
                 titrationwidgets: Sequence[TitrationBaseWidget],
                 datawidgets: Sequence[DataWidget]):
        self.model = model                          # store these values to
        self.titrationwidgets = titrationwidgets    # update their values after
        self.datawidgets = datawidgets              # the fitting

        self.data: Dict[int, DataWidget] = {}

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
        titration_match = {}    # keys are id of datawidgets and values are id of the matching TitrationBaseWidget
        self.spectraldata = None

        # ~~~~ start collecting information ~~~~
        jacobian_slice = Slices()
        residual_slice = Slices()
        self.data_order: list = []    # the id of the datawidgets in order of appearance

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

    def iter_unique_pairs_data(self):
        """Iterate over unique pair of data with repetition.

        Yields:
            int: the id of the first data begin yielded
            type: the typr of the first data being yielded
            data: the first data itself
            int: the id of the second data begin yielded
            type: the typr of the second data being yielded
            data: the second data itself
        """
        for data1, data2 in itertools.combinations_with_replacement(self.data_order, 2):
            id1, type1 = data1
            id2, type2 = data2
            yield ((id1, type1, self.data[id1]),(id2, type2, self.data[id2]))

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
        """Processes flags for refinement, setting up constraints and variables.

        This method inspects the flags within a specified attribute of `data` to determine
        which parameters are refinable, which are constants, and which are constrained. It then
        updates the provided slice to reflect any variables identified for refinement.

        Args:
            data (Any): Object containing the data and associated flags.
            attr_values (str): Name of the attribute containing the variable values.
            attr_flags (str): Name of the attribute containing the flags.
            attr_toref (str): Name of the attribute to store indices of refinable variables.
            jslice (Slices): Slice object used to incrementally build Jacobian slices.

        Notes:
            Constraints are added if the flag falls within a predefined range, and new `Variable`
            objects are created for each refinable parameter.
        """
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
            temperature = self.get_temp(),
            error_emf = np.array(widget.emf0_error)
        )
        rstep = data.emf.size
        rslice.step(rstep)
        data.vslice = rslice.yield_slice()
        self._process_flags(data, 'emf0', "emf0_flags", "emf0_torefine", jslice)
        key = (id(widget), type(widget))
        jslice.stamp_slice(key, self.jacobian_part)
        return data

    def _process_spectrum(self, datawidget, jacobian_slice, residual_slice):
        data = SpectrumData(absorbance=np.array(datawidget.spectrum),
                            wavelength=np.fromiter(datawidget.wavelengths, dtype=float),
                            optical_path=datawidget.optical_path)
        self.spectraldata.add_spectrum(datawidget.optically_active, data)
        return data

    def _process_titration(self, twidget, jslice):
        data = TitrationData(init=np.array(twidget.initial_amount),
                             init_flags=twidget.init_flags,
                             buret=np.array(twidget.buret),
                             buret_flags=twidget.buret_flags,
                             starting_volume=twidget.starting_volume,
                             titre=np.array(twidget.titre),
                             error_volume=twidget.volume_error)

        self._process_flags(data, 'init', 'init_flags', 'rf_init', jslice)
        key = (id(twidget), 'init')
        jslice.stamp_slice(key, self.jacobian_part)

        self._process_flags(data, 'buret', 'buret_flags', 'rf_buret', jslice)
        key = (id(twidget), 'buret')
        jslice.stamp_slice(key, self.jacobian_part)

        return data


@dataclass
class TitrationData():
    """Container for titration data required for analysis and refinement.

    Attributes:
        init (np.ndarray): Initial amounts.
        buret (np.ndarray): Buret measurements.
        titre (np.ndarray): Titre values.
        starting_volume (float): Starting volume for titration.
        error_volume (float): Error associated with volume measurements.
        init_flags (tuple[int]): Flags for initial amounts.
        buret_flags (tuple[int]): Flags for buret measurements.
        free_conc (NDArray[float] | None): Free concentration, calculated post-initialization.
        anal_conc (NDArray[float]): Analytical concentration.
        amatrix (NDArray[float]): Analytical matrix for Jacobian computation.
        dlcdlbeta (NDArray[float]): Derivative matrix with respect to beta values.
        rf_init (tuple[int]): Indices of refinable initial amounts.
        rf_buret (tuple[int]): Indices of refinable buret measurements.
        refine (bool): Indicates if this titration is being refined.

    Example:
        >>> init_values = np.array([0.1, 0.2])
        >>> buret_values = np.array([0.05, 0.1])
        >>> data = TitrationData(init=init_values, buret=buret_values, ...)
        >>> data.analc()
    """
    init: np.ndarray           # values of the initial amount
    buret: np.ndarray          # values of the buret
    titre: np.ndarray
    starting_volume: float
    error_volume: float
    init_flags: tuple[int]
    buret_flags: tuple[int]
    free_conc: NDArray[float] | None = field(init=False, default=None)
    anal_conc: NDArray[float] = field(init=False)
    amatrix: NDArray[float] = field(init=False, repr=False)
    dlcdlbeta: NDArray[float] = field(init=False, repr=False)
    rf_init: tuple[int] = field(init=False)  # indices of elements to be refined
    rf_buret: tuple[int] = field(init=False )# indices of elements to be refined
    refine: bool = field(init=False, default=False)
    init_slice: Slices = field(init=False, repr=False)
    buret_slice: Slices = field(init=False, repr=False)

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
    error_emf: float

    def dump(self, widget: EmfWidget) -> None:
        "Dump data into the widget to update the GUI."
        widget.emffitted = self.emf_calc()
        if any(self.emf0_flags):
            widget.emf0 = self.emf0

    def emf_calc(self) -> np.ndarray:
        "Return calculated emf values."
        hconc = libemf.hselect(self.titration.free_conc, self.electroactive)
        return libemf.nernst(hconc, self.emf0, self.slope, 0.0, self.temperature)

    def variance(self) -> NDArray[float]:
        return np.full_like(self.emf, self.error_emf**2)

    def weight(self) -> NDArray[float]:
        return np.gradient(self.emf.T, self.titration.titre, axis=1) * self.titration.error_volume

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
    """Handles a single variable's data and operations for fitting routines.

    This class provides methods to retrieve, update, and increment a variable's value,
    as well as to manage associated error data. Each `Variable` instance represents a
    single parameter that can be refined in a fitting process, with constraints on the
    maximum allowable change per iteration.

    Attributes:
        data (Any): The object that holds the actual value of the variable.
        key (str): The attribute name in `data` that contains the variable value.
        keyerror (str): The attribute name in `data` that holds the variable error.
        position (int): The index of this variable within the data array.
        max_increment (float): The maximum allowed change in value per increment step.
        stored_value (float | None): Stores the current value when incrementing;
                                     reset after accepting the new value.
        last_increment (float): The last increment applied to the variable's value.

    Methods:
        accept_value(): Resets `stored_value`, indicating acceptance of the current value.
        increment_value(increment: float): Adjusts the variable's value by a specified increment,
                                           respecting `max_increment`.
        get_value() -> float: Retrieves the current value of the variable.
        get_error() -> float: Retrieves the current error associated with the variable.
        set_error(value: float): Sets the error for the variable.
        set_value(value: float): Updates the variable's value directly.

    Example:
        >>> variable = Variable(data=some_data_obj, key="attribute_name", position=0)
        >>> variable.get_value()   # Retrieve the current value
        >>> variable.increment_value(0.05)  # Increment by 0.05, respecting max_increment
        >>> variable.accept_value()  # Accept the current increment, resetting stored_value
    """
    def __init__(self, data, key, position):
        if not hasattr(data, key):
            raise ValueError(f"data {type(data)} does not contain property {key}")
        self.data: float = data
        self.key = key
        self.keyerror = "errors"
        self.position: int = position
        self.max_increment: float = math.inf
        self.stored_value: float | None = None
        self.increment: float = 0.0
        self.previous_value: float = 0.0
        self.last_increment: float

    def accept_value(self) -> None:
        self.previous_value = self.stored_value
        self.last_increment = self.increment
        self.set_value(self.stored_value + self.increment)
        self.stored_value = None
        self.increment = 0.0

    def increment_value(self, increment: float) -> None:
        if abs(increment) > self.max_increment:
            increment = math.copysign(self.max_increment, increment)

        if self.stored_value is None:
            self.stored_value = self.get_value()

        self.increment = increment

    def get_value(self) -> float:
        "Get the value of the variable."
        if self.stored_value is None:
            dataholder = getattr(self.data, self.key)
            return dataholder[self.position]
        else:
            return self.stored_value + self.increment

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
        # setattr(self.data, self.key, dataholder)


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
        for (data, _, position), val in zip(self.__data, new_values):
            dataholder = getattr(data, "errors")
            dataholder[position] = val

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


def trivial_capping(x: NDArray[float], dx: NDArray[float]) -> NDArray[float]:
    "Capping function where there is no capping"
    return x + dx


def max_ratio_capping(x: NDArray[float], dx: NDArray[float], ratio: float) -> NDArray[float]:
    "Capping to a fraction of change"
    if (aux := np.abs(dx)/x) < ratio:
        return x+dx
    return x*(1+aux)


def abs_capping(x, dx, maximum):
    if abs(dx) < maximum:
        return x + dx
    return x + maximum
