"""Data encapsulation classes.

.. module:: data.py
.. moduleauthor:: Salvador Blasco <salvador.blasco@gmail.com>

Data encapsulation that serves as middleclass between the GUI and the backend
that does the actual calculations.

- :class:`Data`
- :class:`ModelData`
- :class:`Dataset`
- :class:`CaloriData`
- :class:`EMFData`
- :class:`SpecData`
- :class:`NMRData`
- :class:`SimulationData`
- :class:`SpeciationData`
- :class:`TitrationData`
"""

# pylint: disable=invalid-name
# pylint: disable=no-member
# pylint: disable=attribute-defined-outside-init

import itertools
import xml.etree.ElementTree as ET

import numpy as np

import consts
import libaux
import libeq
import libplot


class Data:
    """This class handles all potentiometric data together in order to better
    manage data transfer between different routines and help to keep
    consistency with the data.
    """
    def __init__(self):
        self.clear()

    # -------------
    # ↓↓ methods ↓↓
    # -------------

    def addEquilibrium(self, b=0.0, pos=-1, p=None):
        """Add a new equilibrium to the current model.

        Parameters:
            b (float): The value of the constant in decimal logarithmic units
            pos (int): The index of the position in which the new equilibrium
                will be inserted. Defaults to -1 (the last position).
            p (:class:`numpy.ndarray`): A 1D array of *S* ints representing
                the stoichiometry associated to the new equilibrium and will
                be inserted as a row in the :term:`stoichiometry array`.  If
                not provided, zeros will be used.
        """
        self.__currentmodel.addEquilibrium(b, pos, p)
        self.ui.model_widget.refreshWidgets()

    def addReagent(self, label='', pos=-1):
        """Add new reagent.

        Adds one additional reagent to the :term:`independent components`
        to all models and datasets.

        Parameters:
            label (str): A label for this reagent
            pos (int): The position in which will be inserted
        """
        if pos == -1:
            pos = len(self.labels)

        self.__labels.insert(pos, label)
        for m in itertools.chain(self.__models,
                                 *[iter(i) for i in self.__alldata]):
            m.addReagent(pos)

    def clear(self):
        """Clear the content of this object and restores the initial values."""

        self.__currentmodel = ModelData(parent=self, name="model #1")
        self.__models = [self.__currentmodel]
        self.__labels = ['L', 'H']
        self.__datasets = []
        self.__speciations = []
        self.__titrations = []
        self.__calordatasets = []
        self.__emfdatasets = []
        self.__specdatasets = []
        self.__alldata = (self.__speciations, self.__titrations,
                          self.__calordatasets, self.__emfdatasets,
                          self.__specdatasets)
        self.__title = "no title"
        self.__temperature = 25.00
        self.__chisq = []
        self.__htmlout = None

    def delete(self, what):
        """Delete the specified object.

        Parameters:
            what (:class:`Dataset` or :class:`SimulationData`): The object to
                be deleted
        """
        for s in self.__alldata:
            if what in s:
                n = s.index(what)    # may rise ValueError
                del s[n]
                break

    def fit(self, method=consts.METHOD_LM, method_pars=None):
        """Fit data.

        Method for fitting the free variables defined in the flags to the
        experimental data. It takes no arguments and returns no value. The data
        is taken from this object and updates the values of the refined
        variables.

        Parameters:
            method (int): Indicates which algorithm will be used for fitting.
                Accepted values are 0 for Levenverg-Marquardt algorithm
                (non-linear least squares, default) and 1 for Nelder-Mead
                algorithm (simplex). Any other value will raise
                :class:`ValueError`.
            method_pars (dict): A dict with the parameters needed for either
                method. They will be unpacked and passed upstream. See
                documentation for them.

        """
        self._callfit(method, method_pars)

    def indexEMF(self, d):
        "Return the index of a given EMFData."
        return self.__emfdatasets.index(d)

    def indexSpeciation(self, d):
        "Return the index of a given Speciation."
        return self.__speciations.index(d)

    def indexTitration(self, d):
        "Return the index of a given Titration."
        return self.__titrations.index(d)

    def iterate(self):
        """Perform one iteration.

        This routine does the same as :func:`fit()` but it does not perform
        the full fitting until convergence but only one iteration is done
        and the the results are updated"""

        self._callfit({'one_iter': True})

    def newEmfDs(self):
        """Create an empty :class:`EMFData` and returns a reference to it.

        Returns:
            :class:`EMFData`: A reference to the newly created dataset
        """

        ds = EMFData(parent=self)
        self.__emfdatasets.append(ds)
        return ds

    def newCalorDs(self):
        """Creates an empty :class:`CaloriData` and returns a reference to it.

        Returns:
            :class:`CaloriData`: A reference to the newly created dataset
        """

        ds = CaloriData(parent=self)
        self.__calordatasets.append(ds)
        self.__currentmodel.addCalorimetry()
        return ds

    def newSpecDs(self):
        """Create an empty :class:`SpecData` and returns a reference to it.

        Returns:
            :class:`SpecData`: A reference to the newly created dataset
        """
        ds = SpecData(parent=self)
        self.__specdatasets.append(ds)
        return ds

    def newModel(self):
        "Create an empty model and sets as current."
        m = ModelData(parent=self, name="new model")
        self.__models.append(m)
        self.__currentmodel = m

    def newSpeciation(self):
        """Create an empty :class:`SpeciationData`.

        Returns:
            :class:`SpeciationData`: A reference to the newly created
                speciation.
        """
        ds = SpeciationData(parent=self)
        self.__speciations.append(ds)
        return ds

    def newTitration(self):
        """Add new titration.

        Creates an empty :class:`TitrationData` and returns a reference to it.

        Returns:
            :class:`TitrationData`: A reference to the newly created
                speciation.
        """
        ds = TitrationData(parent=self)
        self.__titrations.append(ds)
        return ds

    def removeEquilibrium(self, pos=-1):
        # TODO if there is calorimetry, remove also the corresponding enthalpy
        self.model.removeEquilibrium(pos)

    def removeReagent(self, n):
        "Remove reagent from the list of reagents."
        libaux.assert_type(int, n)
        del self.labels[n]
        for m in itertools.chain(self.__models,
                                 *[iter(i) for i in self.__alldata]):
            m.removeReagent(n)

    def resetChisq(self):
        "Reset fitting history."
        self.__chisq = []

    def _callfit(self, method=consts.METHOD_LM, method_pars=None, xkwarg=None):
        """Method for fitting the free variables defined in the flags to the
        experimental data. It takes no arguments and returns no value. The data
        is taken from this object and updates the values of the refined
        variables.

        Parameters:
            method (int): Indicates which algorithm will be used for fitting.
                Accepted values are 0 for Levenverg-Marquardt algorithm
                (non-linear least squares, default) and 1 for Nelder-Mead
                algorithm (simplex). Any other value will raise
                :class:`ValueError`.
            method_pars (dict): A dict with the parameters needed for either
                method. They will be unpacked and passed upstream. See
                documentation for them.
            xkwarg (dict): Extra parameters to be unpacked
        """

        # TODO split cases into functions
        if xkwarg is None:
            xkwarg = {}
        else:
            libaux.assert_type(dict, xkwarg)

        if self.__herald is not None:
            self.__herald('I am fitting')

        # potentiometry
        if self.NEmf > 0 and self.NCal == 0:
            import libemf
            use = True
            kwa = {}
            for d in self.iterEmfDs(use):
                d.setWeightingScheme('auto')
            kwa['weights'] = tuple(d.weights for d in self.iterEmfDs(use))

            for d in self.iterEmfDs(use):
                if d.C is None:
                    d.calcC()

            c0 = tuple(d.C for d in self.iterEmfDs(use))
            if all(c is not None for c in c0):
                kwa['c0'] = c0
            kwa['out_chisq'] = self.__chisq
            if self.__herald is not None:
                kwa['qt_handle'] = self.__herald
            if self.__htmlout is not None:
                kwa['htmlout'] = self.__htmlout

            for k in xkwarg:
                kwa[k] = xkwarg[k]

            args = (self.model.reflogB,
                    self.model.refBflags,
                    self.model.refP,
                    tuple(d.emf_read for d in self.iterEmfDs(use)),
                    tuple(d.getTitration() for d in self.iterEmfDs(use)),
                    tuple(d.getElectrode() for d in self.iterEmfDs(use)))

            nlnB, C, alt = libemf.emffit(*args, method=method, **method_pars,
                                         **kwa)
            self.model.reflogB = nlnB
            if 'error' in alt:
                self.model.referrlogB = alt['error']

            assert isinstance(C, list)
            for c, d in zip(C, self.iterEmfDs()):
                d.C = c

        # Calorimetry
        if self.NEmf == 0 and self.NCal > 0:
            import libcal
            use = True
            kwa = {}
            kwa['weights'] = tuple(d.getW(1) for d in self.iterEmfDs(use))
            kwa['T0_flags'] = tuple(d.Tflags for d in self.iterEmfDs(use))
            args = (self.model.reflogB,
                    self.model.refBflags,
                    self.model.enthalpy,
                    self.model.enthalpy_flags,
                    self.model.refP,
                    tuple(d.Q for d in self.iterEmfDs(use)),
                    tuple(d.getTitration() for d in self.iterEmfDs(use)))
            nlnB, nH, C, alt = libcal.calfit(*args, method=method,
                                             **method_pars, **kwa)
            self.model.reflogB = nlnB
            self.model.enthalpy = nH
            if 'error' in alt:
                self.model.referrlogB = alt['error']

            assert isinstance(C, list)
            for c, d in zip(C, self.iterEmfDs()):
                d.C = c

    # ----------------------------------
    # ↓↓ properties, getters, setters ↓↓
    # ----------------------------------

    @property
    def B(self):
        "The values of the :term:`equilibrium constants array`."
        return 10**self.__currentmodel.reflogB

    @property
    def Bflags(self):
        "The values of the flags for constants refinement"
        return self.__currentmodel.refBflags

    # DELETE
    @Bflags.setter
    def Bflags(self, value):
        self.__currentmodel.Bflags = value

    @property
    def chiSq(self):
        "a list containing the χ₂ values of each iteration"
        return self.__chisq

    @property
    def errlogB(self):
        """The value of the error associated to the log:sub:`10` values of the
        constants of the current model"""
        return self.__currentmodel.referrlogB

    @errlogB.setter
    def errlogB(self, value):
        self.__currentmodel.errlogB = value

    @property
    def fullLabels(self):
        "full list of labels"
        return libplot.make_labels(self.labels, self.P)

    @property
    def herald(self):
        "This is a handlet for sending messages into the GUI."
        return self.__herald

    @herald.setter
    def herald(self, h):
        self.__herald = h

    @property
    def htmlout(self):
        "This is a handlet for sending messages into the GUI."
        return self.__htmlout

    @htmlout.setter
    def htmlout(self, h):
        self.__htmlout = h

    @property
    def labels(self):
        "Returns the labels for the independent components."
        return self.__labels

    @labels.setter
    def labels(self, value):
        self.__labels = value

    @property
    def logB(self):
        "The value of the log:sub:`10` of the constants of the current model"
        return self.model.reflogB

    @logB.setter
    def logB(self, value):
        self.model.reflogB = value

    @property
    def model(self):
        "The current model"
        return self.__currentmodel

    @property
    def models(self):
        "A list of :class:`ModelData` being used"
        return self.__models

    @property
    def NCal(self):
        "The number of calorimetry datasets in use"
        return sum(1 if d.use else 0 for d in self.__calordatasets)

    @property
    def NDatasets(self):
        "The number of datasets currently being used."
        raise NotImplementedError

    @property
    def NEmf(self):
        "The number of potentiometry datasets in use"
        return sum(1 if d.use else 0 for d in self.__emfdatasets)

    @property
    def NSpec(self):
        "The number of spectrometry datasets in use"
        return sum(1 if d.use else 0 for d in self.__specdatasets)

    @property
    def Nt(self):
        "return total number of experimental points"
        raise NotImplementedError

    @property
    def Nvt(self):
        "return total number of not excluded experimental points"
        raise NotImplementedError

    @property
    def P(self):
        "The :term:`stoichiometry array` of the current model."
        return self.__currentmodel.refP

    @P.setter
    def P(self, value):
        self.__currentmodel.refP = value

    @property
    def temperature(self):
        "The temperature of the system."
        return self.__temperature

    @temperature.setter
    def temperature(self, T):
        if not isinstance(float, T):
            raise TypeError("T must be float")
        self.__temperature = T

    @property
    def title(self):
        "The title of the project"
        return self.__title

    @title.setter
    def title(self, t):
        self.__title = t

    @property
    def E(self):
        "(int) the number of equilibria"
        return self.model.refP.shape[0]

    @property
    def S(self):
        "The number of independent components in the current model"
        return self.model.P.shape[1]

    def setCurrentModel(self, n):
        if not isinstance(n, int):
            raise TypeError
        self.__currentmodel = self.__models[n]

    def getUse(self, ds):
        return self.__mask[ds]

    def setWeightingScheme(self, scheme='auto'):
        if scheme == 'auto' or scheme == 'unit':
            for d in self.__emfdatasets:
                d.setWeightingScheme(scheme)
        else:
            raise ValueError("Invalid weighting scheme")

    # -----------------
    # ↓↓  iterators  ↓↓
    # -----------------

    def iterTitr(self):
        "Iterate over titrations"
        return iter(self.__titrations)

    def iterSimu(self):
        "Iterate over simulations"
        return iter(self.__speciations)

    def iterCalorDs(self, use=True):
        """Iterate over calorimetry data.

        Parameters:
            use (bool): If True, iterate over the datasets that have the flag
                'use' set. If False, iterate over all datasets.
        """
        return (d for d in self.__calordatasets if (not use) or d.use)

    def iterEmfDs(self, use=False):
        """Iterate over potentiometric data.

        Parameters:
            use (bool): If True, iterate over the datasets that have the flag
                'use' set. If False, iterate over all datasets.
        """
        return (d for d in self.__emfdatasets if (not use) or d.use)

    def iterLabels(self):
        "Iterate over labels."
        return iter(self.__labels)

    def iterSpecDs(self, use=True):
        """Iterate over spectrometric data.

        Parameters:
            use (bool): If True, iterate over the datasets that have the flag
                'use' set. If False, iterate over all datasets.
        """
        return (d for d in self.__specdatasets if (not use) or d.use)

    def plotCalor(self):
        """Provide data for plotting calorimetry data.

        Yields:
            #. number of datasets to plot
            #. titre used
            #. titre ignored or None
            #. heat data used
            #. heat data ignored or None
            #. volume for fit (all points)
            #. free concentration array
            #. calculated heat
            #. residuals for used points
            #. residuals for ignored points
        """

        use = True
        ct = True
        yield len(self.__calordatasets)
        yield (d.getV(1, caltweak=ct) for d in self.iterCalorDs(use))
        yield (d.getV(-1, caltweak=ct) for d in self.iterCalorDs(use))
        yield (d.getQ(1) for d in self.iterCalorDs(use))
        yield (d.getQ(-1) for d in self.iterCalorDs(use))
        yield (d.V for d in self.iterCalorDs(use))
        yield (d.C for d in self.iterCalorDs(use))
        yield (d.getCalcQ() for d in self.iterCalorDs(use))
        yield (d.getResiduals(1) for d in self.iterCalorDs(use))
        yield (d.getResiduals(-1) for d in self.iterCalorDs(use))

    def plotEMF(self, plot_ignored=True):
        """Returns a generator with contains all the data needed to plot the
        potentiometry data.

        Yields:
            #. The number of datasets to plot
            #. The use flag
            #. Volume of titre
            #. Experimental data (masked)
            #. Fitted emf (masked)
            #. Free concentrations array
            #. Residuals (masked)
        """

        use = True
        if plot_ignored:
            yield len(self.__emfdatasets)
        else:
            yield sum(1 for i in self.iterEmfDs(use))

        yield (d.use for d in self.iterEmfDs(use))
        yield (d.V for d in self.iterEmfDs(use))
        yield (d.emf_read for d in self.iterEmfDs(use))
        yield (d.emf_fit for d in self.iterEmfDs(use))
        yield (d.C for d in self.iterEmfDs(use))
        yield (d.residuals for d in self.iterEmfDs(use))

    def plotSpeciation(self, d, ref=None, pX=False):
        """Yield data needed to draw a spectiation plot.

        Returns a generator with all the data needed to plot the species
        distribution.

        Parameters:
            d (:class:`Dataset`): The index for the species simulation to plot
            ref (int or None): Optional. The index of the species for the
                relative concentration to be referred to. If it is none
                the absolute concentrations are plotted.
            pX (bool): Optional. Default is False. If True, the x-axis is
                in logarithmic units.

        Yields:
            #. the x-axis array
            #. the y-values
            #. the errors (or None)
            #. the labels
        """
        full_labels = self.fullLabels

        if pX:
            new_x = -np.log10(d.X)
        else:
            new_x = d.X
        yield new_x

        if ref is None:
            new_y = d.C
            err_y = d.eC
            labels = full_labels
        else:
            if not isinstance(ref, int):
                raise TypeError

            # find reference species
            # !ref_label = ptype[1:]
            # !ref = self.data.labels.index(ref_label)
            rP = self.P[:, ref].tolist()

            # find species to plot (all those that contain reference
            # species)
            to_plot = [n+self.S for n, i in enumerate(rP) if i != 0]
            to_plot.insert(0, ref)
            to_plot.sort()
            q_plot = [i for i in rP if i != 0]
            q_plot.insert(0, 1)

            # Tp = [[a, b] for a, b in zip(d.T0, d.T1)]
            # T = libaux.build_T(Tp, d.pX, d.N)
            T = d.T

            # remove species not to plot
            pcC = libeq.percent_distribution(d.C, self.P, T, ref)

            # arrange labels
            # !full_labels = libplot.make_labels(labels, P)
            labels = [full_labels[i] for i in to_plot]
            err_y = libeq.percent_distribution(d.eC, self.P, T, ref)
            new_y = pcC
            # titl = '% Formation Relative to ' + ref_label

        yield new_y
        yield err_y
        yield labels

    def plotSpectr(self):
        """Yield data needed to plot a spectrophotometric titration.

        Returns a generator with contains all the data needed to plot the
        spectrophotometric data.

        Yields:
            #. The number of datasets to plot
            #. The total number of spectra sets to plot
            #. A generator that return all spectra in turn
            #. A generator with all residuals in turn
            #. fit spectra
        """
        # yield number of datasets
        NDs = len(self.__specdatasets)
        yield NDs

        # yield total number of spectra
        yield sum([d.Nds for d in self.__specdatasets])

        # yield a generator with all spectra in turn
        yield itertools.chain(*(i.spectra for i in self.__specdatasets))

        # yield a generator with all residuals in turn
        yield itertools.chain(*(i.residuals for i in self.__specdatasets))

        # yield model
        yield NDs * [SpecData.MODEL_LINEAR]

        # yield fit spectra
        yield (i.epsilon for i in self.__specdatasets)

    # --------------------------------
    # ↓↓ importers, savers, loaders ↓↓
    # --------------------------------

    def loadXML(self, f):
        """Load an XML file into self.

        Parameters:
            f (file or str): The file to read data from
        """
        tree = ET.parse(f)
        root = tree.getroot()

        if root.tag != 'apes':
            raise ValueError("XML file does not contain APES data")

        self.clear()

        if root.attrib['version'] == '1.0':
            # process 1.0 version of the xml file

            self.__labels = root.find('labels').text.split()
            models = root.find('models')
            self.__models = []
            for model in models.iter('model'):
                m = ModelData(self)
                m.loadXML(model)
                self.__models.append(m)
            self.__currentmodel = self.__models[0]

            # self.__speciations = []
            for esim in root.iter('distri'):
                spe = self.newSpeciation()
                spe.loadXML(esim)

            for etit in root.iter('simu'):
                tit = self.newTitration()
                tit.loadXML(etit)

            tags = ('potentiometricdata', 'calordata')
            calls = (self.newEmfDs, self.newCalorDs)
            for t, c in zip(tags, calls):
                dats = root.find(t)
                if dats:
                    for curve in dats.iter('curve'):
                        d = c()
                        d.loadXML(curve)

        else:
            raise ValueError("error reading file")

    def saveXML(self, f):
        """Save the data into an XML file.

        Parameters:
            f (str): A file to save the data into. If file exists it will be
                overwritten.
        """

        root = ET.Element('apes', attrib={'version': '1.0'})
        tree = ET.ElementTree(root)

        title = ET.SubElement(root, 'title')
        title.text = self.title

        # TODO add metadata support
        # metadata = ET.SubElement(root, 'metadata')

        ET.SubElement(root, 'labels').text = " ".join(self.labels)

        models = ET.SubElement(root, 'models')
        # iterate over models
        for model in iter(self.__models):
            models.append(model.saveXML())

        for d in itertools.chain(self.iterTitr(), self.iterSimu()):
            root.append(d.saveXML())

        A = ('potentiometricdata', 'specdata', 'calordata')
        B = (len(i) for i in (self.__emfdatasets, self.__specdatasets,
                              self.__calordatasets))
        C = (self.iterEmfDs(use=False), self.iterSpecDs(use=False),
             self.iterCalorDs(use=False))

        for a, b, c in zip(A, B, C):
            if b:
                dtag = ET.SubElement(root, a)
                for t in c:
                    dtag.append(t.saveXML())

        tree.write(f)

    def importFromHyperquad(self, f):
        """Import data from a hyperquad file and loads it into this object.

        Parameters:
            fn (str): A valid file name
        """
        import libio
        d = libio.importHyperquad(f)
        self.clear()

        self.__labels = d['labels']
        self.logB = np.array(d['logB'])
        self.Bflags = d['Bkeys']
        self.P = np.array(d['P'])
        self.title = d['title']
        E = self.P.shape[0]
        self.errlogB = E*[0.0]

        if 'titration' in d:
            t = self.newTitration()
            t.V0 = d['titration']['V0']
            t.V1 = d['titration']['V1']
            t.T0 = d['titration']['T0']
            t.buret = np.array(d['titration']['buret'])

        if 'simulation' in d:
            # V0, T0, T1
            t = self.newSpeciation()
            V = d['simulation']['V0']
            t.T0 = np.array(d['simulation']['T0'])/V
            t.T1 = np.array(d['simulation']['T1'])/V

        # TODO complete with datasets

    def importFromSuperquad(self, fn):
        """Import data from a superquad file and loads it into this object.

        Parameters:
            fn (str): A valid file name
        """
        import libio
        data = libio.importSuperquad(fn)
        self.clear()
        self.__temperature = data['temperature']
        self.__labels = data['labels']
        self.logB = data['logB']
        self.errlogB = np.zeros_like(self.logB)
        self.Bflags = data['flags']
        self.P = data['P']

        self.__emfdatasets = []
        z = zip(data['V'], data['emf'], data['E0'], data['n'], data['E0flags'],
                data['error_E0'], data['V0'], data['error_V0'], data['T0'],
                data['T0flags'], data['buret'], data['hindex'], data['fRTnF'])
        for v, e, e0, n, ef, ee0, v0, ev0, t0, tf, b, h, f in z:
            d = EMFData(parent=self)
            d.V = v
            d.emf_read = e
            d.E0 = [e0]
            d.n = n
            d.E0flags = [ef]
            d.errE0 = [ee0]
            d.V0 = v0
            d.errV0 = ev0
            d.T0 = t0
            d.Tflags = tf
            d.buret = b
            d.hindex = h
            d.fRTnF = f
            self.__emfdatasets.append(d)


class ModelData:
    """This class contains and manages the data corresponding to a model,
    namely the equilibrium constants, the stoichiometric coefficients and
    other.

    Here the following data are encapsulated:

    - the equilibrium constants which is an :class:`numpy.ndarray` that
        contains the values of the equilibrium constants in decimal
        logarithmic units.
    - the error of the equilibrium constants is an :class:`numpy.ndarray`
        that represents the error of the equilibrium constants in decimal
        logarithmic units.
    - the refinement flags for the constants which is a list of ints
        representing what to do with the constants. Accepted values are
        -1 for ignoring this constant, 0 for not refine this constant, 1 for
        refining this constant, 2 and more for constraining the constant.
    - the stoichiometry array.
    """

    def __init__(self, parent, name="noname"):
        self._parent = parent

        # 1D array with the equilibrium constants in log10 units
        self.__logB = np.array([8.65, -13.73])
        # 1D array with the errors of equilibrium constants
        self.__errlogB = np.array([0.0, 0.0], dtype=np.float)
        # list of refimenet flags for the equilibrium constants
        self.__Bflags = 2*[consts.RF_REFINE]
        self.__P = np.array([[1, 1], [0, -1]])
        self.__name = name
        self.__iscalori = False  # whether model contains calorimetry data
        self.__isspec = False    # whether model contains spectroscopy data
        self.__specdat = None
        self.__isnmr = False     # whether model contains nmr data
        self.__nmdat = None

    def isCalori(self):
        return self._parent.NCal > 0

    def addCalorimetry(self):
        "Add calorimetry data to the current model."
        if self.__iscalori:
            return
        self.__iscalori = True
        self.__enthalpy = np.array(self.E * [0.0])
        self.__enthalpyflags = self.E * [0]
        self.__errorenthalpy = np.array(self.E * [0.0])

    def addReagent(self, pos=-1, p=None):
        """Adds a new reagent (:term:`independent component`) to the model

        Parameters:
            pos (int): The index of the position in which the new reagent will
                be inserted. Defaults to -1 (the last position).
            p (:class:`numpy.ndarray`): An int or 1D array of *E* ints
                representing the stoichiometry associated to the new reagent
                with respect to all the equilibria and will be inserted as a
                column in the :term:`stoichiometry array`. If not provided,
                zeros will be used. If it is an int, that number will be used
                for all.
        """

        if p is None:
            p = 0
        else:
            libaux.assert_array_dim(1, p)
            if len(p) != self.E:
                raise ValueError("p must have %d elements" % self.E)

        self.__P = np.insert(self.__P, p, pos, axis=1)

    def removeReagent(self, pos=-1):
        "Removes one reagent"
        self.__P = np.delete(self.__P, pos, axis=1)

    def addEquilibrium(self, b=0.0, pos=-1, p=None):
        """Adds a new equilibrium to the model

        Parameters:
            b (float): The value of the constant in decimal logarithmic units
            pos (int): The index of the position in which the new equilibrium
                will be inserted. Defaults to -1 (the last position).
            p (:class:`numpy.ndarray`): A 1D array of *S* ints representing
                the stoichiometry associated to the new equilibrium and will
                be inserted as a row in the :term:`stoichiometry array`.  If
                not provided, zeros will be used.
        """

        if p is None:
            p = np.zeros(self.S, np.int)
        else:
            libaux.assert_array_dim(1, p)
            if len(p) != self.E:
                raise ValueError("p must have %d elements" % self.E)
        if pos == -1:
            pos = self.E

        self.__P = np.insert(self.__P, pos, p, axis=0)
        self.__logB = np.insert(self.__logB, pos, b)
        self.__errlogB = np.insert(self.__errlogB, pos, 0.0)
        self.__Bflags.insert(pos, 0)
        if self.__iscalori:
            self.__enthalpy = np.insert(self.__enthalpy, pos, 0.0)
            self.__enthalpyflags.insert(pos, 0)
            self.__errorenthalpy = np.insert(self.__errorenthalpy, pos, 0.0)

    def removeEquilibrium(self, pos=-1):
        self.__P = np.delete(self.__P, pos, axis=0)

    def loadXML(self, xmle):
        """Loads the contents of an :class:`xml.etree.ElementTree.Element`
        into this class.

        Parameters:
            :class:`xml.etree.ElementTree.Element`: An Element object
                containing the information about the model to be loaded.

        Raises:
            ValueError: if **xmle** has the wrong type or is not the right
                xml substructure.
        """

        if not isinstance(xmle, ET.Element):
            raise ValueError("xmle must be an ElementTree.Element instance")

        if xmle.tag != 'model':
            raise ValueError("Wrong XML piece was loaded")

        # t = " ".join(map(str, range(1000)))
        # In [17]: %timeit a = np.fromiter(iter(t.split()), np.int)
        # 10000 loops, best of 3: 146 µs per loop
        # In [18]: %timeit a = np.array([int(i) for i in t.split()])
        # 1000 loops, best of 3: 297 µs per loop

        self.logB = np.fromiter(iter(xmle.find('beta').text.split()), np.float)
        E = len(self.logB)
        eB = xmle.find('error')
        if eB is not None:
            self.errlogB = np.fromiter(iter(eB.text.split()), np.float)
        else:
            self.errlogB = np.array(E*[0.0])
        lP = np.fromiter(iter(xmle.find('p').text.split()), np.int)
        self.P = lP.reshape((E, len(lP)//E))
        self.Bflags = [int(i) for i in xmle.find('key').text.split()]
        c = xmle.find('enthalpy')
        if c is not None:
            self.addCalorimetry()
            self.enthalpy = np.fromiter(iter(c.text.split()), np.float)

        c = xmle.find('enthalpy_error')
        if c is not None:
            self.__errorenthalpy = np.fromiter(iter(c.text.split()), np.float)

    def saveXML(self):
        """Returns an :class:`xml.etree.ElementTree.Element` containing the
        information from this class in order to write in the xml file.

        Returns:
            :class:`xml.etree.ElementTree`: An ElementTree containing the
                information about the model.
        """

        def j(a):
            return " ".join(map(str, a))

        tm = ET.Element('model')
        ET.SubElement(tm, 'beta').text = j(self.logB.tolist())
        ET.SubElement(tm, 'error').text = j(self.errlogB.tolist())
        ET.SubElement(tm, 'p').text = j(self.P.flat)
        ET.SubElement(tm, 'key').text = j(self.Bflags)

        if self.__iscalori:
            ET.SubElement(tm, 'enthalpy').text = j(self.enthalpy.tolist())
            ET.SubElement(tm, 'enthalpy_error').text = \
                j(self.error_enthalpy.tolist())
            ET.SubElement(tm, 'enthalpy_key').text = j(self.enthalpy_flags)
            ET.SubElement(tm, 'enthropy').text = j(self.enthropy.tolist())
            ET.SubElement(tm, 'enthropy_error').text = \
                j(self.error_enthropy.tolist())

        return tm

    # ------------------
    # ↓↓  properties  ↓↓
    # ------------------

    @property
    def B(self):
        """The values of the :term:`equilibrium constants array`.

        .. note:: This is the **full** array, including ignored values.
        """
        return np.power(10.0, self.__logB)

    @property
    def Bflags(self):
        """The flags for the :term:`equilibrium constants array`.

        .. note:: This is the **full** array, including ignored values.
        """
        return self.__Bflags

    @Bflags.setter
    def Bflags(self, value):
        libaux.assert_type(list, value)
        self.__Bflags = value

    @property
    def E(self):
        return self.__P.shape[0]

    @property
    def enthalpy(self):
        return self._chcalo(self.__enthalpy)

    @enthalpy.setter
    def enthalpy(self, H):
        self.__enthalpy = H

    @property
    def error_enthalpy(self):
        return self._chcalo(self.__enthalpy)

    @error_enthalpy.setter
    def error_enthalpy(self, eH):
        self.__errorenthalpy = eH

    @property
    def enthalpy_flags(self):
        return self._chcalo(self.__enthalpyflags)

    @property
    def enthropy(self):
        "Returns the enthropy of the system in ΔS units"
        H = self.enthalpy
        T = self._parent.temperature
        # R = 8.314472  # gas constant(J ∕ (K⋅mol))
        # TΔS = ΔH + RTlogβ
        assert isinstance(H, np.ndarray)
        assert H.ndim == 1
        assert len(H) == self.E
        return H/T + consts.R*consts.LOGK*self.logB

    @property
    def error_enthropy(self):
        "Returns the error of enthropy of the system in ΔS units"
        eH = self.error_enthalpy
        T = self._parent.temperature
        # R = 8.314472  # gas constant(J ∕ (K⋅mol))
        elnB = consts.LOGK*self.errlogB
        # e²(ΔS) = e²(ΔH)/T + R e²(logβ)
        return np.sqrt(eH**2/T + consts.R*elnB**2)

    def _chcalo(self, prop):
        if self.__iscalori:
            return prop
        else:
            raise RuntimeError('No calorimetry data available')

    @property
    def labels(self):
        return self._parent.labels

    @labels.setter
    def labels(self, value):
        # TODO check input
        self._parent.labels = value

    @property
    def logB(self):
        """The values of the :term:`equilibrium constants array` in
        logarithmic units.

        .. note:: This is the **full** array, including ignored values.
        """
        return self.__logB

    @logB.setter
    def logB(self, value):
        self.__logB = value

    @property
    def name(self):
        "The name for this object"
        return self.__name

    @name.setter
    def name(self, value):
        self.__name = value

    @property
    def errlogB(self):
        return self.__errlogB

    @errlogB.setter
    def errlogB(self, value):
        self.__errlogB = value

    @property
    def P(self):
        "The :term:`stoichiometry array` of the current model."
        return self.__P

    @P.setter
    def P(self, value):
        self.__P = value

    @property
    def parent(self):
        return self._parent

    @property
    def refBflags(self):
        i = self._nridx()
        assert isinstance(i, list)
        assert isinstance(self.__logB, np.ndarray)
        return np.array([self.__Bflags[j] for j in i])

    @property
    def reflogB(self):
        """A view of the :class:`numpy.ndarray` of the non-ignored constants
        only"""

        return self._refget(self.__logB)

    @reflogB.setter
    def reflogB(self, value):
        self._refset(self.__logB, value)

    @property
    def referrlogB(self):
        return self._refget(self.__errlogB)

    @referrlogB.setter
    def referrlogB(self, value):
        self._refset(self.__errlogB, value)

    @property
    def refP(self):
        """The stoichiometry array including only the entries that are
        not ignored"""
        # TODO
        i = self._nridx()
        assert isinstance(i, list)
        assert isinstance(self.__logB, np.ndarray)
        return self.__P[i, :]

    @property
    def refenthalpy(self):
        return self._refget(self.__enthalpy)

    @refenthalpy.setter
    def refenthalpy(self, H):
        self._refset(self.__enthalpy, H)

    @property
    def refenthalpyflags(self):
        return self._refget(self.__enthalpyflags)

    @refenthalpyflags.setter
    def refenthalpyflags(self, hf):
        self._refset(self.__enthalpyflags, hf)

    @property
    def residuals(self):
        return self.__residuals

    @property
    def S(self):
        "The number of independent components"
        return self.__P.shape[1]

    # ---------------------------
    # ↓↓  auxiliary functions  ↓↓
    # ---------------------------

    def _nridx(self):
        """Returns a list with the indices of the equilibria not to be
        ignored.

        >>> self.Bflags = [0, -1, 0, 1, 2, 2, 1, 0, -1, 0, -1]
        >>> self._nridx()
        [0, 2, 3, 4, 5, 6, 7, 9]
        """

        return [i for i, j in enumerate(self.__Bflags) if j != -1]

    def _refget(self, prop):
        i = self._nridx()
        assert isinstance(i, list)
        assert isinstance(prop, np.ndarray)
        return prop[i]

    def _refset(self, prop, value):
        i = self._nridx()
        assert isinstance(i, list)
        libaux.assert_same_len(i, value)
        prop[i] = value


class Dataset:
    """This is the base class for all types of datasets"""

    def __init__(self, parent):
        if not isinstance(parent, Data):
            raise TypeError
        self._parent = parent
        self._name = "no name"
        self._use = True           # bool or list of bools
        self._C = None              # array or list of arrays
        self._V = np.array([0.0])  # volume of titre in mL
        self._V0 = 0.0
        self._error_V0 = 0.003
        self._T0 = self.S*[0.01]
        self._Tflags = self.S*[0]
        self._buret = self.S*[0.0]
        # self._vmask = None
        self._istitr = True

        # a masked array or None. This represents the experimental data
        # in which titration points go along axis=0 and may have any number
        # of dimmensions >= 2
        self._expdata = None

        # An array with the data fitted. Must be of the same size than
        # _expdata. Shares mask with it.
        self._fitdata = None

        # An array that works as auxiliary and contains the residuals.
        # Shares mask with _expdata.
        self._residuals = None

        # Weights array for the fitting.
        self._weights = None

    def addReagent(self, n):
        "Remove reagent from the list of reagents."
        libaux.assert_type(int, n)
        if n == -1:
            n = len(self._T0)
        if self.__istitr:
            self._T0.insert(n, 0.0)
            self._Tflags.insert(n, 0)
            self._buret.insert(n, 0.0)
        else:
            self._T = np.insert(self._T, n, 0.0, axis=1)
        self._C = None

    def calcC(self):
        "Updates the free concentrations array"
        if self._C is not None:
            kwargs = {'x0': self.C[:, :self.S]}
        else:
            kwargs = {}
        self._C = libeq.consol(self.B, self.P, self.T, **kwargs)
        libaux.assert_same_shapeN(0, self.C, self.T)

    def loadCurveXML(self, xmle):
        "Loads the titration data from an XML file into self."

        if 'title' in xmle.attrib:
            self.__name = xmle.attrib['title']
        if 'use' in xmle.attrib:
            self.use = bool(xmle.attrib['title'])

        def fif(x):
            "auxiliary function"
            return np.fromiter((float(i) for i in xmle.find(x).text.split()),
                               np.float)
        self.T0 = fif('T0')
        self.Tflags = list(map(int, xmle.find('T0key').text.split()))
        self.buret = fif('buret')
        self.V0, self.__error_V0 = map(float, xmle.find('V0').text.split())
        self.V = fif('titre')

        c = xmle.find('c')
        if c is not None:
            self._C = np.fromiter(map(float, c.text.split()),
                                  np.float).reshape(self.N, self.E+self.S)

    def saveCurveXML(self, tcurve):
        "Puts self data into XML object"
        def j(x):
            return " ".join(map(str, x))

        if self._istitr:
            ET.SubElement(tcurve, 'T0', {'unit': 'mmol'}).text = j(self.T0)
            ET.SubElement(tcurve, 'buret', {'unit': 'mmol/mL'}).text = \
                j(self.buret)
            ET.SubElement(tcurve, 'T0key').text = j(self.Tflags)
        else:
            ET.SubElement(tcurve, 'T', {'unit': 'mmol/mL'}).text = j(self.T0)

        ET.SubElement(tcurve, 'V0', {'unit': 'mL'}).text = \
            "%f %f" % (self.V0, self.errV0)

        ET.SubElement(tcurve, 'titre').text = " ".join(map(str, self.V.flat))

        if self.C is not None:
            ET.SubElement(tcurve, 'c').text = " ".join(map(str, self.C.flat))

    def getTitration(self):
        "Returns 'titration' dict which can be fed to low level routines"
        if self._istitr:
            return {'T': self.T, 'V0': self.V0, 'V': self.V}
        else:
            return {'T0': self.T0, 'Tflags': self.Tflags,
                    'buret': self.buret, 'V0': self.V0, 'V': self.V}

    def maskCol(self, index):
        self._expdata[:, index, ...] = np.ma.masked

    def unmaskCol(self, index):
        self._expdata[:, index, ...] = np.ma.nomask

    def maskPoint(self, index):
        self._expdata[index, ...] = np.ma.masked

    def unmaskPoint(self, index):
        self._expdata[index, ...] = np.ma.nomask

    def maskValues(self, indices):
        self._expdata[indices] = np.ma.masked

    def unmaskValues(self, indices):
        self._expdata[indices] = np.ma.nomask

    def removeReagent(self, n):
        "Remove reagent from the list of reagents."
        libaux.assert_type(int, n)
        if self.__istitr:
            self._T0 = np.delete(self._T0, n)
            self._Tflags.pop(n)
            self._buret = np.delete(self._buret, n)
        else:
            self.__T = np.delete(self._T, n, 1)
        self._C = None

    # DELETE. Use mask methods instead
    def usePoint(self, luse):
        """Flags one or some titration point as 'not ignored' to be taken
        into account in fitting.

        Parameters:
            ign (int or iterable): The indices of the points that will be
                flagged.
        """
        self._ignuse(luse, True)

    # ------------------
    # ↓↓  properties  ↓↓
    # ------------------

    @property
    def B(self):
        "The equilibrium constants array"
        return 10**self._parent.logB

    @property
    def buret(self):
        "an array containing the titre in mL for each point"
        return self._buret

    @buret.setter
    def buret(self, buret):
        self._buret = buret

    @property
    def C(self):
        "The :term:`free concentrations array`"
        return self._C

    @C.setter
    def C(self, C):
        a, b = C.shape
        if a != self.N or b != (self.E + self.S):
            raise ValueError
        self._C = C

    @property
    def expdata(self):
        return self._expdata

    @property
    def E(self):
        "the number of equilibrium constants"
        return self._parent.E

    @property
    def labels(self):
        return self._parent.labels

    @property
    def mask(self):
        return self._expdata.mask

    @property
    def N(self):
        "The number of experimental datapoints"
        return len(self._V)

    @property
    def name(self):
        "The name of this dataset"
        return self._name

    @name.setter
    def name(self, n):
        self._name = n

    @property
    def P(self):
        "The :term:`stoichiometry array` of the current model."
        return self._parent.P

    @property
    def S(self):
        "the number of principal components"
        return self._parent.S

    @property
    def use(self):
        "A flag indicating whether this dataset will be used or not"
        return self._use

    @use.setter
    def use(self, u):
        self._use = bool(u)

    @property
    def errV0(self):
        return self._error_V0

    @errV0.setter
    def errV0(self, errV0):
        self._error_V0 = errV0

    @property
    def T(self):
        """The :term:`total concentrations array`. Notice that this important
        parameter can be set directly by assigning a valis array of values
        to this variable, or indirectly by setting the parameters *T0*, *V*,
        *V0* and *buret*."""
        if self._istitr:
            return libaux.build_T_titr2(np.array(self.T0),
                                        np.array(self.buret),
                                        self.V0, self.V)
        else:
            return self._T

    @T.setter
    def T(self, T):
        self._istitr = False
        self._T = T

    @property
    def T0(self):
        "The initial total amount for each component, in mmol"
        return self._T0

    @T0.setter
    def T0(self, T0):
        self._T0 = T0

    @property
    def Tflags(self):
        "The refinement flags for T0"
        return self._Tflags

    @Tflags.setter
    def Tflags(self, Tflags):
        self._Tflags = Tflags

    @property
    def V(self):
        "the volume of titre in mL"
        return self._V

    @V.setter
    def V(self, V):
        libaux.assert_array_dim(1, V)
        self._V = V

    @property
    def V0(self):
        "The initial volume in mL"
        return self._V0

    @V0.setter
    def V0(self, V0):
        self._V0 = V0

    @property
    def weights(self):
        return self._weights

    @weights.setter
    def weights(self, w):
        self._weights = np.ma.array(w, mask=self.mask)

    def getV(self, mask=0, caltweak=False):
        # ?if caltweak:
        # ?    return getPrMask(self._V[1:], self._vmask[1:], mask)
        # ?else:
        # ?    return getPrMask(self._V, self._vmask, mask)
        m = self._expdata.mask
        assert m.shape[0] + caltweak == len(self.V)

        if caltweak:
            n = 1
        else:
            n = 0
        v = np.moveaxis(np.tile(self.V[n:], (*m.shape[1:], 1)), -1, 0)
        if mask == 0:
            ret = v
        elif mask == -1:
            ret = np.ma.MaskedArray(v, ~m)
        elif mask == 1:
            ret = np.ma.MaskedArray(v, m)
        else:
            raise ValueError

        return ret

    def _setExpData(self, dat, mask=False):
        self._expdata = np.ma.masked_array(dat, mask)
        dmask = self._expdata.mask
        self._weights = np.ma.masked_array(np.ones_like(dat), dmask)
        self._fitdata = np.ma.masked_array(np.zeros_like(dat), dmask)
        self._residuals = self._expdata - self._fitdata


class CaloriData(Dataset):
    """Abstraction for calorimetric data"""

    def __init__(self, parent):
        super(CaloriData, self).__init__(parent)
        self.__Q = None
        self._vmask = None

    def getCalcQ(self):
        # TODO duplicated code in libcal.py
        H = self._parent.model.enthalpy
        if self._C is None:
            self.calcC()
        # V = self.getV(mask, caltweak=True)
        # C = self.getC(mask)
        # if C is None:
        #     C = getPrMask(np.zeros((self.N, self.E+self.S)),
        #                   self._vmask,
        #                   mask)
        aux1 = self.V[:, np.newaxis] * self.C[:, self.S:]
        aux2 = aux1[:-1, :] - aux1[1:, :]
        return np.sum(aux2*H, axis=1)

    def getQ(self, mask=0):
        m = self._expdata.mask
        p = self._expdata.data
        if mask == 0:
            ret = p
        elif mask == -1:
            ret = np.ma.MaskedArray(p, ~m)
        elif mask == 1:
            ret = np.ma.MaskedArray(p, m)
        else:
            raise ValueError

        return ret

    def getResiduals(self, mask=0):
        return self.Q - self.getCalcQ()

    @property
    def Q(self):
        return self.__Q

    @Q.setter
    def Q(self, q):
        libaux.assert_array_dim(1, q)
        # self.__Q = q
        # self._vmask = np.array(len(q)*[True])
        self._setExpData(q)
        self.__Q = self._expdata

    def loadXML(self, xmle):
        """Loads the contents of an :class:`xml.etree.ElementTree.Element`
        into this class.

        Parameters:
            xmle (:class:`xml.etree.ElementTree`): An ElementTree containing
                the information to be loaded.

        Raises:
            ValueError: if **xmle** has the wrong type or is not the right
                xml substructure.
        """

        if xmle.tag != 'curve':
            raise ValueError("Wrong XML piece was loaded")

        # TODO check data
        self.loadCurveXML(xmle)
        self.Q = np.fromiter(iter(xmle.find('heat').text.split()), np.float)

    def saveXML(self):
        """Returns an :class:`xml.etree.ElementTree.Element` containing the
        information from this class in order to write in the xml file.

        Returns:
            :class:`xml.etree.ElementTree`: An ElementTree containing the
                information about the calorimetry data.
        """

        def j(x):
            return " ".join(map(str, x))

        tcurve = ET.Element('curve', {'title': self.name,
                                      'use': str(self.use)})
        self.saveCurveXML(tcurve)
        ET.SubElement(tcurve, 'heat').text = j(self.Q)
        if not np.all(self.vmask):
            ET.SubElement(tcurve, 'key').text = j(self.vmask)
        return tcurve


class EMFData(Dataset):
    """This class encapsulates the data for a potentiometric titration."""

    def __init__(self, parent):
        super(EMFData, self).__init__(parent)

        # array with the emf data read. Titration points must be on axis=1 and
        # each electrode read (if more than one) on axis=0.
        self._setExpData(np.array([0.0]))
        assert isinstance(self.emf_read, np.ma.core.MaskedArray)

        # The standard potential per each electrode, in mV
        self.__E0 = [0.0]

        # The number of electrons for each electrode
        self.__n = [1]

        # Nernst's tweak parameter
        self.__f = None

        # The refinement flag for each electrode.
        self.__E0flags = [0]

        # The error associated with E0, for each electrode.
        self.__error_E0 = [0.03]

        # The index for the species to which each electrode acts upon.
        self.__hindex = [-1]

        # The sum of the residuals squared. Each member of the list represents
        # one iteration.
        self.__chisq = []

        # Proportionality constant for Nernst's equation. One per electrode.
        self.__fRTnF = [consts.fRTnF]  # mV

    def calcC(self):
        super().calcC()
        # H = self.C[:, self.__hindex, ...]
        # fitd = np.squeeze(self.__E0 + self.__fRTnF*np.log(H))
        # assert fitd.shape == self._expdata.shape
        # self._fitdata = np.ma.MaskedArray(fitd, self._expdata.mask)
        # self._residuals = self._expdata - self._fitdata
        self._recC()

    def getElectrode(self):
        """Returns a dict with the data of the electrodes which will be used
        to feed the fitting routines."""
        return {'E0': self.E0, 'E0flags': self.E0flags, 'hindex': self.hindex,
                'fRTnF': self.fRTnF}

    def loadXML(self, xmle, temperature=(25+273.15)):
        """Loads the contents of an :class:`xml.etree.ElementTree.Element`
        into this class.

        Parameters:
            xmle (:class:`xml.etree.ElementTree`): An ElementTree containing
                the information about the model to be loaded.
            temperature (float): The temperature of the system, in kelvin

        Raises:
            ValueError: if **xmle** has the wrong type or is not the right
                xml substructure.
        """

        if xmle.tag != 'curve':
            raise ValueError("Wrong XML piece was loaded")

        # TODO check data
        self.__E0 = [float(i) for i in xmle.find('emf0').text.split()]
        nelec = len(self.__E0)
        self.__error_E0 = [float(i)
                           for i in xmle.find('erroremf0').text.split()]
        self.__hindex = [int(i) for i in xmle.find('active').text.split()]
        self.__E0flags = [int(i) for i in xmle.find('emf0flags').text.split()]
        self.__n = [int(i) for i in xmle.find('n').text.split()]

        f = xmle.find('f')
        if not f:
            self.__f = nelec*[1.0]
        else:
            self.__f = [str(i) for i in f.split()]

        fRTnF = xmle.find('fRTnF')
        if not fRTnF:
            # 0.086173424 mV/K
            self.__fRTnF = [f*temperature*consts.RoverF/n
                            for f, n in zip(self.__f, self.__n)]
        else:
            self.__fRTnF = [float(i) for i in fRTnF.text.split()]

        def fif(x):
            return np.fromiter((float(i) for i in xmle.find(x).text.split()),
                               np.float)

        self.loadCurveXML(xmle)

        # RESHAPE
        emf = np.squeeze(fif('emfread').reshape((self.N, self.NEl)))
        t = xmle.find('key')
        if t is not None:
            mask = np.fromiter((not bool(int(i)) for i in t.text.split()),
                               np.bool)
        else:
            mask = False  # np.array(self.N*[True])
        # self.vmask = _mask
        self._setExpData(emf, mask)
        # self.__emf_read = self._expdata

    def saveXML(self):
        """Returns an :class:`xml.etree.ElementTree.Element` containing the
        information from this class in order to write in the xml file.

        Returns:
            :class:`xml.etree.ElementTree`: An ElementTree containing the
                information about the model.
        """

        def j(x):
            if type(x) is list or type(x) is tuple:
                return " ".join(map(str, x))
            else:
                return str(x)

        tcurve = ET.Element('curve', {'title': self.name,
                                      'use': str(self.use)})
        ET.SubElement(tcurve, 'emf0', {'unit': 'mV'}).text = j(self.E0)
        ET.SubElement(tcurve, 'erroremf0', {'unit': 'mV'}).text = j(self.errE0)
        ET.SubElement(tcurve, 'active').text = j(self.hindex)
        ET.SubElement(tcurve, 'n').text = j(self.n)
        ET.SubElement(tcurve, 'emf0flags').text = j(self.E0flags)
        if self.f is not None:
            ET.SubElement(tcurve, 'f').text = j(self.f)
        ET.SubElement(tcurve, 'fRTnF').text = j(self.fRTnF)

        self.saveCurveXML(tcurve)

        ET.SubElement(tcurve, 'emfread').text = \
            " ".join(str(i) for i in self.emf_read.data.flat)
        ET.SubElement(tcurve, 'key').text = \
            " ".join("0" if i else "1" for i in self.mask.flat)
        return tcurve

    def setWeightingScheme(self, scheme='auto'):
        if scheme == 'unit':
            self.weights = np.ones(self.N)
        elif scheme == 'auto':
            import libemf
            self.weights = libemf.weighting(self.V, self.emf_read, self.errV0,
                                            self.errE0)
        else:
            raise ValueError("Invalid weighting scheme")

    # ------------------
    # ↓↓  properties  ↓↓
    # ------------------

    @property
    def C(self):
        return super().C

    def _recC(self):
        H = self.C[:, self.__hindex, ...]
        fitd = np.squeeze(self.__E0 + self.__fRTnF*np.log(H))
        assert fitd.shape == self._expdata.shape
        self._fitdata = np.ma.array(fitd, mask=self._expdata.mask)
        self._residuals = self._expdata - self._fitdata

    @C.setter
    def C(self, C):
        super(EMFData, self.__class__).C.__set__(self, C)
        self._recC()

    @property
    def E0(self):
        "The standard potential in mV"
        return self.__E0

    @E0.setter
    def E0(self, E0):
        self.__E0 = E0

    @property
    def E0flags(self):
        "Refinement flags for standard potentials"
        return self.__E0flags

    @E0flags.setter
    def E0flags(self, E0f):
        self.__E0flags = E0f

    @property
    def emf_fit(self):
        return self._fitdata

    @property
    def emf_read(self):
        return self._expdata

    @emf_read.setter
    def emf_read(self, emf):
        libaux.assert_array(emf)
        self._setExpData(emf)

    @property
    def residuals(self):
        return self._residuals

    @property
    def errE0(self):
        "Errors in standard potential"
        return self.__error_E0

    @errE0.setter
    def errE0(self, errE0):
        self.__error_E0 = errE0

    @property
    def f(self):
        return self.__f

    @f.setter
    def f(self, f):
        self.__f = f

    @property
    def fRTnF(self):
        return self.__fRTnF

    @fRTnF.setter
    def fRTnF(self, fRTnF):
        self.__fRTnF = fRTnF

    @property
    def hindex(self):
        return self.__hindex

    @hindex.setter
    def hindex(self, hindex):
        self.__hindex = hindex

    @property
    def n(self):
        "the number of electrons reacted with the electrode"
        return self.__n

    @n.setter
    def n(self, nn):
        self.__n = nn

    @property
    def NEl(self):
        return len(self.__E0)


class SpecData(Dataset):
    """This class handles all spectrometric data that belongs to the same
    titration.  """

    MODEL_LINEAR = 0

    def __init__(self, parent):
        super(SpecData, self).__init__(parent)

        self.__spectra = []      # spectra
        self.__calcspectra = []  # calculated spectra
        self.__spectratype = []  # list of type of spectra
        self.__residuals = []
        self.__model = 0         # model for data fitting
        self.__coloured = self.S*[True]
        self.__epsilon = []
        self.__wavelengths = None
        self.__ignoredwl = []    # list of indices of ignored wavelengths
        self.__ignoredsp = []    # list of indices of ignored spectra
        self.__ignoredpt = []    # list of indices of ignored individual points

    # ----------------------------
    # ↓↓        methods         ↓↓
    # ----------------------------

    def addSpectra(self, spectra, sptype=0):
        """Adds a spectra set to the collection. A spectra set is usually a
        titration with N experimental points each with a spectrum of Nl
        wavelengths + 1 column for wavelength values.

        Parameters:
            spectra (:class:`numpy.ndarray`): A 2D array containing the data.
            stype (int): Optional. The type of spectrum.
        """

        libaux.assert_array_dim(2, spectra)
        if len(self.__spectra) == 0:
            self.__spectra.append(spectra)
        else:
            if self.__spectra.shape[0] == spectra.shape[0]:
                self.__spectra.append(spectra)
            else:
                txt = "wrong dimmensions for spectra: " + \
                      "{} provided, ".format(self.__spectra.shape[0]) + \
                      "{} expected".format(spectra.shape[0])
                raise ValueError(txt)

        self.__spectratype.append(sptype)
        # TODO ↓ uncomment when resolved
        # !self.updateParms()

        self.__calcspectra = []

    def addSpectraFromFile(self, f, sptype=0):
        "Load spectra from a data file using numpy routines"
        s = np.loadtxt(f)
        # TODO check s
        self.addSpectra(s, sptype)

    def calcSpectra(self):
        """Returns the calculated spectra."""

        C = self._C

        if self.__model == self.MODEL_LINEAR:
            self.__calcspectra = [C*self.__epsilon[t]
                                  for t in self.__spectratype]
        else:
            raise NotImplementedError

    def getSpectrum(self, n):
        self._ignorecmm(n)
        return self.__spectra[n]

    def _ignorecmm(self, nset):
        if not (0 < nset < len(self.__spectra)):
            raise IndexError

    def ignoreWl(self, nset, wl):
        """Ignores a specific wavelength from a given spectrum set.

        Parameters:
            nset (int): Index of the spectra set to which the wavelength is
                referred.
            wl (int): Index of the wavelength that will be ignored.

        Raises:
            IndexError: if any of the indices are incorrect
        """

        assert len(self.__ignoredwl) == len(self.__spectra)
        self._ignorecmm(nset)
        self.__ignoredwl[nset].append(wl)

    def ignoreSp(self, nset, sp):
        """Ignores a specific spectrum from a given spectrum set.

        Parameters:
            nset (int): Index of the spectra set to which the wavelength is
                referred.
            sp (int): Index of the spectrum that will be ignored.

        Raises:
            IndexError: if any of the indices are incorrect
        """
        self._ignorecmm(nset)
        self.__ignoredsp[nset].extend(sp)

    def ignorePt(self, nset, wl, sp):
        """Ignores a specific point from a given spectrum set.

        Parameters:
            nset (int): Index of the spectra set to which the wavelength is
                referred.
            sp (int): Index of the spectrum that will be ignored.
            wl (int): Index of the wavelength that will be ignored.

        Raises:
            IndexError: if any of the indices are incorrect
        """
        self._ignorecmm(nset)
        self.__ignoredpt[nset].extend(zip(sp, wl))

    def removeSpectra(self, n):
        """Removes one spectra set from this dataset.

        Parameters:
            n (int): The index of the spectra array that will be deleted.

        Raises:
            IndexError: if the index provided does not exist"""

        del self.__spectra[n]
        del self.__spectratype[n]
        if len(self.__calcspectra) != 0:
            del self.__calcspectra[n]

    def useWl(self, nset, wl):
        """Unignores a specific wavelength from a given spectrum set.

        Parameters:
            nset (int): Index of the spectra set to which the wavelength is
                referred.
            wl (int): Index of the wavelength that will be un-ignored.

        Raises:
            IndexError: if any of the indices are incorrect
        """

        assert len(self.__ignoredwl) == len(self.__spectra)
        self._ignorecmm(nset)
        del self.__ignoredwl[nset][wl]

    def useSp(self, nset, sp):
        """Ignores a specific spectrum from a given spectrum set.

        Parameters:
            nset (int): Index of the spectra set to which the wavelength is
                referred.
            sp (int): Index of the spectrum that will be ignored.

        Raises:
            IndexError: if any of the indices are incorrect
        """
        self._ignorecmm(nset)
        self.__ignoredsp[nset].append(sp)

    def usePt(self, nset, wl, sp):
        """Ignores a specific point from a given spectrum set.

        Parameters:
            nset (int): Index of the spectra set to which the wavelength is
                referred.
            sp (int): Index of the spectrum that will be ignored.
            wl (int): Index of the wavelength that will be ignored.

        Raises:
            IndexError: if any of the indices are incorrect
        """
        self._ignorecmm(nset)
        self.__ignoredpt[nset].append((sp, wl))

    # ---------------
    # ↓↓ iterators ↓↓
    # ---------------

    def iterResidual(self):
        return iter(self.__residuals)

    def iterSpectra(self):
        return iter(self.__spectra)

    # ----------------------------------
    # ↓↓ properties, getters, setters ↓↓
    # ----------------------------------

    @property
    def model(self):
        return self.__model

    # !@property
    # !def N(self):
    # !    "The number of experimental points available"
    # !    return self.__spectra.shape[0]

    @property
    def Nds(self):
        "The number of spectra sets"
        return len(self.__spectra)

    @property
    def Nl(self):
        "The number of wavelength available"
        return self.__spectra.shape[1]

    @property
    def w(self):
        return self.__weights

    @property
    def coloured(self):
        return self.__coloured

    @property
    def epsilon(self):
        return self.__epsilon

    @property
    def fitSpectra(self):
        if self.epsilon is None:
            return None
        else:
            return np.dot(self.epsilon, self._C)

    @property
    def residuals(self):
        if self.epsilon is None:
            return self.Nds * [None]
        for a, b in zip(self.__spectra, self.fitSpectra):
            assert all(a.shape == b.shape)
            yield a-b

    @property
    def spectra(self):
        return self.__spectra

    def getUse(self, ds):
        return self.__mask[ds]

    def setUse(self, ds, value):
        self.__mask[ds] = value
        return


class NMRData(Dataset):
    def __init__(self, parent):
        super(NMRData, self).__init__(parent)
        self.__nucleilabels = None
        self.__shifts = None


class SimulationData:
    """This class provides encapsulation for the data for a titration
    simulation."""

    def __init__(self, parent):
        self._parent = parent
        self.__N = 50
        self.__erlvl = 0.95
        # ↓ an array with x values column 0 and y vals in the rest
        self.__externaldata = None
        # ↓ a string with the title to be written in the twin y axis
        self.__externaldatatitle = None
        self.__externaldatalabels = None    # a list of strings
        self.__C = None
        self.__X = None
        self.__eC = None

    # ------------------
    # ↓↓   methods    ↓↓
    # ------------------

    def clearExternalData(self):
        "Clears external data."
        self.__externaldata = None
        self.__externaldatatitle = None
        self.__externaldatalabels = None

    def getExternalData(self):
        return self.__externaldatatitle, \
               self.__externaldatalabels, self.__externaldata

    def setExternalData(self, labels, data):
        libaux.assert_array(data)
        if not isinstance(labels, (list, tuple)):
            raise ValueError("'labels' must be a list or a tuple")
        if data.shape[1] != 2*len(labels):
            raise ValueError("'labels' and 'data' are of different size")
        self.__externaldata = data
        self.__externaldatalabels = labels

    # ------------------
    # ↓↓  properties  ↓↓
    # ------------------

    @property
    def B(self):
        "The values of the :term:`equilibrium constants array`."
        return self._parent.B

    @property
    def C(self):
        "the free concentrations array."
        return self.__C

    @C.setter
    def C(self, c):
        # TODO input check
        self.__C = c

    @property
    def eC(self):
        "the errors associated with the free equilibrium concentrations array"
        return self.__eC

    @eC.setter
    def eC(self, ec):
        # TODO input check
        self.__eC = ec

    @property
    def erlvl(self):
        """The error level at which the uncertanty of the concentrations is
        calculated"""
        return self.__erlvl

    @erlvl.setter
    def erlvl(self, n):
        self.__erlvl = n

    @property
    def externalDataTitle(self):
        return self.__externaldatatitle

    @externalDataTitle.setter
    def externalDataTitle(self, t):
        self.__externaldatatitle = t

    @property
    def logB(self):
        """The values of the :term:`equilibrium constants array` in
        logarithmic units"""
        return self._parent.logB

    @property
    def errlogB(self):
        return self._parent.errlogB

    @property
    def N(self):
        "The number of points to be simulated"
        return self.__N

    @N.setter
    def N(self, N):
        self.__N = N

    @property
    def P(self):
        "The :term:`stoichiometry array` of the current model."
        return self._parent.P

    @property
    def S(self):
        "The number of :term:`independent components`."
        return self._parent.S

    @property
    def labels(self):
        return self._parent.labels

    @property
    def X(self):
        return self.__X

    @X.setter
    def X(self, x):
        # TODO input check
        self.__X = x

    # ------------------
    # ↓↓  iterators   ↓↓
    # ------------------

    def iterExternalData(self):
        """Returns an iterator that contains the external data to be plotted
        or None if there is none. The iterator returns in turn a triplet
        containing the label, x and y for each set of data."""
        d = self.__externaldata
        if d is None:
            return None
        for n, t in enumerate(self.__externaldatalabels):
            yield t, d[:, 2*n], d[:, 2*n+1]


class SpeciationData(SimulationData):
    "This class provides encapsulation for the data for a speciation."

    def __init__(self, parent):
        super(SpeciationData, self).__init__(parent)
        self.__T0 = self.S*[0.10]
        self.__T1 = self.S*[0.10]
        self.__pX = self.S*[0]
        self.__T0[-1] = 2
        self.__T1[-1] = 11
        self.__pX[-1] = 1

    def addReagent(self, n):
        "Add new reagent to the list of reagents."
        libaux.assert_type(int, n)
        if n == -1:
            n = len(self.__T0)
        self.__T0.insert(n, 0.0)
        self.__T1.insert(n, 0)
        self.__pX.insert(n, False)

    def calcC(self, x, erlvl=None):
        """This function calculates the free concentrations and the errors.

        Parameters:
            x (int): a valid index indicating which of the independent
                components will be used as independent variable.
            erlvl (float): The error level at which the error in the free
                concentrations will be calculated. If None, no error will
                be calculated."""

        # import pdb
        # pdb.set_trace()
        if erlvl is not None:
            self.erlvl = erlvl

        Tp = [[a, b] for a, b in zip(self.T0, self.T1)]
        T = libaux.build_T(Tp, self.pX, self.N)
        X, C = libeq.simulation(self.B, self.P, T, x=x)
        fullC = np.insert(C, x, X, axis=1)

        if self.erlvl is not None:
            err_C = libeq.beta_uncertainty(fullC, np.array(self.errlogB),
                                           self.P, self.erlvl)
        else:
            err_C = None

        self.C = fullC
        self.X = X
        self.eC = err_C

    def loadXML(self, xmle):
        """Loads XML information into self.

        Parameters:
            xmle (:class:`xml.etree.ElementTree`): object with the XML info.
        """
        # speciation
        def f(x, t):
            return list(map(t, xmle.find(x).text.split()))
        self.__T0 = f('initial', float)
        self.__T1 = f('final', float)
        self.__pX = f('pX', bool)
        # TODO add support for external data

    def removeReagent(self, n):
        "Remove reagent from the list of reagents."
        libaux.assert_type(int, n)
        # TODO
        del self.__T0[n]
        del self.__T1[n]
        del self.__pX[n]

    def saveXML(self):
        tdistri = ET.Element('distri', {'active': 'no'})
        ET.SubElement(tdistri, 'initial', {'unit': 'mol/L'}).text = \
            " ".join(map(str, self.T0))
        ET.SubElement(tdistri, 'final', {'unit': 'mol/L'}).text = \
            " ".join(map(str, self.T1))
        ET.SubElement(tdistri, 'pX').text = " ".join(map(str, self.pX))
        # TODO add external data
        return tdistri

    @property
    def T(self):
        Tp = [[a, b] for a, b in zip(self.T0, self.T1)]
        return libaux.build_T(Tp, self.pX, self.N)

    @property
    def T0(self):
        return self.__T0

    @T0.setter
    def T0(self, T0):
        self.__T0 = T0

    @property
    def T1(self):
        return self.__T1

    @T1.setter
    def T1(self, T1):
        self.__T1 = T1

    @property
    def pX(self):
        return self.__pX

    @pX.setter
    def pX(self, pX):
        self.__pX = pX


class TitrationData(SimulationData):
    def __init__(self, parent):
        super(TitrationData, self).__init__(parent)
        self.__T0 = self.S*[0.01]
        self.__buret = self.S*[0.0]
        self.__V0 = 20.0
        self.__V1 = 21.0

    def addReagent(self, n):
        "Add new reagent to the list of reagents."
        libaux.assert_type(int, n)
        if n == -1:
            n = len(self.__T0)
        self.__T0.insert(n, 0.0)
        self.__buret.insert(n, 0)

    def calcC(self, erlvl=None):
        """This function calculates the free concentrations and the errors"""
        if erlvl is not None:
            self.erlvl = erlvl

        fullC = libeq.consol(self.B, self.P, self.T)
        X = np.linspace(self.V0, self.V0+self.incrV, self.N)

        elB = np.array(self.errlogB)
        if self.erlvl is not None:
            err_C = libeq.beta_uncertainty(fullC, elB, self.P, self.erlvl)
        else:
            err_C = None

        self.C = fullC
        self.X = X
        self.eC = err_C

    def loadXML(self, xmle):
        """Loads XML information into self.

        Parameters:
            xmle (:class:`xml.etree.ElementTree`): object with the XML info.
        """
        # TODO check input
        T0 = list(map(float, xmle.find('initial').text.split()))
        self.__T0 = T0
        b = list(map(float, xmle.find('buret').text.split()))
        self.__buret = b
        V0, V1 = map(float, xmle.find('V').text.split())
        self.__V0 = V0
        self.__V1 = V1
        # TODO add support for external data

    def removeReagent(self, n):
        "Remove reagent from the list of reagents."
        libaux.assert_type(int, n)
        del self.__T0[n]
        del self.__buret[n]

    def saveXML(self):
        """returns an :class:`xml.etree.ElementTree` with the information of
        this :class:`data.SpeciationData` to be saved.
        """
        ttitr = ET.Element('simu', {'active': 'no'})
        ET.SubElement(ttitr, 'initial', {'unit': 'mmol'}).text = \
            " ".join(map(str, self.T0))
        ET.SubElement(ttitr, 'buret', {'unit': 'mol/L'}).text = \
            " ".join(map(str, self.buret))
        ET.SubElement(ttitr, 'V').text = "%f %f" % (self.V0, self.V1)
        # TODO add external data
        return ttitr

    @property
    def buret(self):
        "The concentration in buret in mmol/mL"
        return self.__buret

    @buret.setter
    def buret(self, buret):
        self.__buret = buret

    @property
    def incrV(self):
        "Volume increment for each titration point in mL"
        return (self.V1-self.V0)/self.N

    @property
    def T(self):
        "The :term:`total concentrations array`"
        return libaux.build_T_titr(np.array(self.T0), np.array(self.buret),
                                   self.V0, self.incrV, self.N)

    @property
    def T0(self):
        "Starting amounts of principal components in mmol"
        return self.__T0

    @T0.setter
    def T0(self, T0):
        self.__T0 = T0

    @property
    def V0(self):
        "Initial volume in mL"
        return self.__V0

    @V0.setter
    def V0(self, V0):
        self.__V0 = V0

    @property
    def V1(self):
        "Final volume in mL"
        return self.__V1

    @V1.setter
    def V1(self, V1):
        self.__V1 = V1
