import datetime
import itertools
import xml.etree.ElementTree as ET

import numpy as np

import consts
import libqt
import modelwidget
import calorwidget
import emfwidget
import mainwindow
import nmrwidget
import otherwidgets
import specwidget
import simulationwidgets
import tabwidget
from excepts import DataFileException


def loadXML(app: "mainwindow.MainWindow", f: str) -> None:
    """Load an XML file into application.

    Parameters:
        app (:class:`MainWindow`): a main window instance
        f (file or str): The file to read data from.

    .. warning:: This routine assumes that MainWindow in empty.
    """
    tree = ET.parse(f)
    root = tree.getroot()

    if root.tag != 'apes':
        raise ValueError("XML file does not contain APES data")

    assert root.attrib['version'] == '1.0'  # process 1.0 version of the xml file

    if (xdata := root.find('metadata')) is not None:
        if (_title := xdata.find('title')) is not None:
            app.project.title = _title.text
        if (_author := xdata.find('author')) is not None:
            app.project.author = _author.text
        if (_comments := xdata.find('comments')) is not None:
            app.project.comments = _comments.text
        if (_lastm := xdata.find('lastmodified')) is not None:
            app.project.last_modified = _lastm.text
        if (_created := xdata.find('created')) is not None:
            app.project.created = _created.text

    labels = root.find('labels').text.split()

    _temp = root.find('temperature')
    if _temp is None:
        temperature = 298.15
    else:
        temperature = float(_temp.text)
        if _temp.attrib['units'] == 'celsius':
            temperature += 273.15
        elif _temp.attrib['units'] == 'kelvin':
            pass
        else:
            raise KeyError('unit not recognized')

    app.project.temperature = temperature

    kwargs = {'models': {'labels': labels}}

    tags = ('models', 'distri', 'simu', 'fittinggroup', 'externaldata')
    loaders = (loadModelsXML, loadSpeciationXML, loadTitrationXML, loadFittingXML, loadExternalXML)
    callables = (app.new_model, app.new_speciation, app.new_titration, app.new_fitting_group,
                 app.new_external_data)

    for tag, call, loader in zip(tags, callables, loaders):
        for section in root.iter(tag):
            widget = call()
            if tag in kwargs:
                loader(widget, section, **kwargs[tag])
            else:
                loader(widget, section)


def loadFittingXML(widget, xmle):
    assert isinstance(widget, tabwidget.TabWidget)

    for titration in xmle.findall('titration'):
        twidget = widget.add_titration()
        index = widget.indexOf(twidget)
        widget.setTabText(index, titration.attrib['name'])
        loadTitrationBaseXML(twidget, titration)

    titration_names = tuple(widget.get_titration_names())

    for spectrum in xmle.findall('specdata'):
        twidget = widget.add_spectrumuv()
        index = widget.indexOf(twidget)
        widget.setTabText(index, spectrum.attrib['name'])
        twidget.populate_cb_titration(titration_names)
        loadSpectrXML(twidget, spectrum)

    for potentio in xmle.findall('potentiometricdata'):
        twidget = widget.add_emf()
        index = widget.indexOf(twidget)
        widget.setTabText(index, potentio.attrib['name'])
        twidget.populate_cb_titration(titration_names)
        loadEmfXML(twidget, potentio)

    for nmrdata in xmle.findall('nmrdata'):
        twidget = widget.add_nmr()
        index = widget.indexOf(twidget)
        widget.setTabText(index, nmrdata.attrib['name'])
        twidget.populate_cb_titration(titration_names)
        loadNmrXML(twidget, nmrdat)


def loadExternalXML(widget, xmle):
    """
        <externaldata>
         <title>y2-axis title</title>
         <!-- order="0" applies to all -->
         <x unit="mL">0.1 0.2 0.3 (...) 1.4</x>
         <data  order="1" unit="chemical shift" label="δ / ppm">
          <y>123.3 (...)</y>
          <ey>0.1 (...)</ey>
         </data>
        </externaldata>
    """
    assert isinstance(widget, otherwidgets.ExternalDataWidget)
    title = xmle.find('title').text
    widget.set_title(title)
    # breakpoint()
    x = tuple(map(float, xmle.find('x').text.split()))
    # data = [x]
    max_depth = len(x)

    table = widget.ui.table
    table.clear()
    table.setRowCount(len(x)+1)

    # breakpoint()
    alldata = xmle.findall('data')
    labels = [_.attrib.get('label', '#') for _ in alldata]
    y = [map(float, dt.find('y').text.split()) for dt in alldata] 
    _ey = [dt.find('ey') for dt in alldata] 
    if with_errors := any(_ is not None for _ in _ey):
        ey = [map(float, _.text.split()) if _ is not None  else max_depth*(0.0) for _ in _ey]
    #     table.setColumnCount(1+len(y)+len(ey))
    #     for n, (a, b) in enumerate(zip(y, ey)):
    #         libqt.fill_column(table, n+1, a, row0=1)
    #         libqt.fill_column(table, n+2, b, row0=1)
    # else:
    #     table.setColumnCount(1+len(y))
    #     for n, a in enumerate(y):
    #         libqt.fill_column(table, n+1, a, row0=1)
    # libqt.replace_nones(table)

    #     # ordr = dt.attrib.get('order', '0')
    #     # type_ = dt.attrib['type']
    #     # sordr = dt.attrib.get('suborder', None)
    #     # TODO 'unit' is ignored for the time being
    #     labels.append(dt.attrib.get('label', ''))
    #     _dat = tuple(map(float, dt.text.split()))
    #     if len(_dat) > max_depth:
    #         max_depth = len(_dat)
    #     data.append((ordr, type_, sordr, labels, _dat))

    widget.feed_data(x, y, (ey if with_errors else None))
    # widget.feed_data(data, max_depth)
    widget.set_labels(labels)
    return widget


def loadModelsOnlyXML(modelwidget, filename):
    tree = ET.parse(filename)
    root = tree.getroot()

    if root.tag != 'apes':
        raise ValueError("XML file does not contain APES data")

    assert root.attrib['version'] == '1.0'
    # process 1.0 version of the xml file
    labels = root.find('labels').text.split()
    loadModelsXML(modelwidget, root.find('models'), labels)


def loadModelsXML(widget, xmle, labels):
    """Load the contents of an XML into :class:`apes.ModelWidget`.

    Parameters:
        :class:`xml.etree.ElementTree.Element`: An Element object
            containing the information about the model to be loaded.
    """
    _checkXML(xmle, 'models')
    widget.clear()
    widget.labels = labels
    for model in xmle.iter('model'):
        widget.append(loadModelXML(model))

    # if 'active' in xmle.attrib:
    #     widget.setCurrentModel(int(xmle.attrib['active']))
    # else:
    #     widget.setCurrentModel(0)
    current_model = int(xmle.attrib.get('active', '0'))
    widget.setCurrentModel(current_model)


def loadModelXML(xmle):
    """Load the contents of an XML and return a :class:`consts.Model`.

    Parameters:
        :class:`xml.etree.ElementTree.Element`: An Element object
            containing the information about the model to be loaded.
    Returns:
        :class:`consts.Model`: The data read
    Raises:
        ValueError: if **xmle** has the wrong type or is not the right
            xml substructure.
    """
    # from modelwidget import ModelData
    _checkXML(xmle, 'model')
    model = modelwidget.ModelData()
    # TODO include 'title'

    model.name = xmle.attrib.get('title', 'no title')

    model.const = [float(number) for number in xmle.find('beta').text.split()]
    E = len(model.const)
    eB = xmle.find('error')
    if eB is not None:
        errlogB = [float(number) for number in eB.text.split()]
    else:
        errlogB = E*[0.0]
    model.const_error = errlogB

    p = list(map(int, xmle.find('p').text.split()))
    S = len(p) // E
    model.stoich = [p[S*n:S*(n+1)] for n in range(E)]
    # P = lP.reshape((E, len(lP)//E))
    model.const_flags = [int(i) for i in xmle.find('key').text.split()]

    c = xmle.find('enthalpy')
    if c is not None:
        model.enthalpy = [float(i) for i in c.text.split()]

    c = xmle.find('enthalpy_key')
    if c is not None:
        model.enthalpy_flags = [int(i) for i in c.text.split()]

    c = xmle.find('enthalpy_error')
    if c is not None:
        model.enthalpy_error = [float(i) for i in c.text.split()]

    return model


def loadEmfXML(widget, xmle):
    """Load the contents of XML into this class.

    Parameters:
        widget (:class:`apes.EMDDSWidget`): widget to load data into.
        xmle (:class:`xml.etree.ElementTree`): An ElementTree containing
            the information about the model to be loaded.

    Raises:
        ValueError: if **xmle** has the wrong type or is not the right
            xml substructure.
    """
    _checkXML(xmle, 'potentiometricdata')

    with libqt.signals_blocked(widget):
        widget.ui.table_params.clearContents()
        # IMPORTANT: set always first emf0 in order fot th widget to reshape the table
        widget.emf0 = tuple(_read_seq(xmle, 'emf0'))
        widget.slope = tuple(_read_seq(xmle, 'slope'))
        widget.slope_flags = tuple(_read_seq(xmle, 'slopeflags', dtype=int))
        widget.nelectrons = tuple(_read_seq(xmle, 'n', dtype=int))
        widget.emf0_error = _read_seq(xmle, 'erroremf0')
        widget.active_species = _read_seq(xmle, 'active', dtype=int)
        widget.emf0_flags = tuple(_read_seq(xmle, 'emf0flags', dtype=int))

        nelectrod = widget.nelectrodes

        widget.titration = xmle.attrib['titration']

        exppoints = widget.titration.n_points

        # flat_emf = tuple(_read_seq(xmle, 'emfread'))
        # emf = np.array(flat_emf).reshape(exppoints, nelectrod)
        emf = np.fromiter(_read_seq(xmle, 'emfread'), dtype=float).reshape(exppoints, nelectrod)
        widget.nelectrodes = nelectrod
        widget.npoints = exppoints

        libqt.fill_column(widget.ui.table_data, 0, widget.titration.titre)
        libqt.array2tab(widget.ui.table_data, emf, col0=1)


def loadNmrXML(widget, xmle):
    raise NotImplementedError


# deprecate
def loadCurveXML(widget, xmle):
    "Loads the titration data from an XML file into self."
    raise DeprecationWarning

    if 'title' in xmle.attrib:
        widget.name = xmle.attrib['title']
    if 'use' in xmle.attrib:
        widget.use = bool(xmle.attrib['title'])

    def fif(x):
        "auxiliary function"
        return np.fromiter((float(i) for i in xmle.find(x).text.split()),
                           np.float)
    widget.initial_amount = fif('T0')
    widget.amount_flags = list(map(int, xmle.find('T0key').text.split()))
    widget.buret = fif('buret')
    widget.starting_volume, widget.volume_error = \
        map(float, xmle.find('V0').text.split())
    widget.titre = fif('titre')

    c = xmle.find('c')
    if c is not None:
        conc = np.fromiter(map(float, c.text.split()),
                           np.float).reshape(widget.N, widget.E + widget.S)
        xmle.set_free_concentration(conc)


def loadCalorXML(self, xmle):
    """Loads the contents of an :class:`xml.etree.ElementTree.Element`
    into this class.

    Parameters:
        xmle (:class:`xml.etree.ElementTree`): An ElementTree containing
            the information to be loaded.

    Raises:
        ValueError: if **xmle** has the wrong type or is not the right
            xml substructure.
    """
    # TODO check data
    self.loadCurveXML(xmle)
    self.Q = np.fromiter(iter(xmle.find('heat').text.split()), np.float)


def loadSpeciationXML(widget, xmle):
    """Load speciation XML information.

    Parameters:
        xmle (:class:`xml.etree.ElementTree`): object with the XML info.
    """
    # speciation
    def f(x, t):
        return list(map(t, xmle.find(x).text.split()))

    def g(x):
        return [t == 'True' for t in xmle.find(x).text.split()]

    widget.set_n_points(int(xmle.attrib.get('points', '50')))
    widget.set_initial_concentration(f('initial', float))
    widget.set_final_concentration(f('final', float))
    widget.set_pX(g('pX'))

    xreftag = xmle.find('xref')
    if xreftag is not None:
        isp = _bool(xreftag.get('p', 'no'))
        widget.set_referencex(int(xreftag.text), isp)
    yreftag = xmle.find('xref')
    if yreftag is not None:
        widget.set_referencey(int(yreftag.text)-1)

    # TODO add support for external data


def loadSpectrXML(widget, xmle):
    raise NotImplementedError


def loadTitrationBaseXML(widget, xmle) -> None:
    """Load XML information into self.

    Parameters:
        xmle (:class:`xml.etree.ElementTree`): object with the XML info.
    """
    name = xmle.attrib['name']
    widget.name = name

    verr = float(xmle.find('volumeerror').text)
    widget.volume_error = verr

    init = __read_floats(xmle, 'init')
    widget.initial_amount = init

    initk = __read_ints(xmle, 'initkey')
    widget.init_flags = initk

    buret = __read_floats(xmle, 'buret')
    widget.buret = buret

    buretk = __read_ints(xmle,  'buretkey')
    widget.buret_flags = buretk

    vinit = float(xmle.find('startingvolume').text)
    widget.starting_volume = vinit

    if (titre := xmle.find('titre')) is not None:
        widget.titre = tuple(map(float, titre.text.split()))
        widget.set_volume_implicit(False)
    else:
        vfinal = xmle.find('finalvolume')
        npoints = xmle.find('totalpoints')
        widget.final_volume = float(vfinal.text)
        widget.n_points = int(npoints.text)
        widget.set_volume_implicit(True)


def loadTitrationXML(widget, xmle):
    """Load XML information into self.

    Parameters:
        xmle (:class:`xml.etree.ElementTree`): object with the XML info.
    """
    # TODO check input
    T0 = list(map(float, xmle.find('initial').text.split()))
    widget.set_initial_amount(T0)
    b = list(map(float, xmle.find('buret').text.split()))
    widget.set_buret(b)
    V0, V1 = map(float, xmle.find('V').text.split())
    widget.set_starting_volume(V0)
    widget.set_final_volume(V1)
    # TODO add support for external data


def _bool(x):
    return x.strip().lower() in ('true', 'yes', '1')


def _checkXML(xmle, tag):
    if not isinstance(xmle, ET.Element):
        raise ValueError("xmle must be an ElementTree.Element instance")
    if xmle.tag != tag:
        raise ValueError("Wrong XML piece was loaded")
    return
