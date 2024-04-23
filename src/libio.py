"""Routines for input/output and import/export data."""

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


__version__ = '0.2'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  SAVE ROUTINES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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
        loadTitrationBaseXML(twidget, titration)

    for spectrum in xmle.findall('specdata'):
        twidget = widget.add_spectrumuv()
        loadSpectrXML(twidget, spectrum)

    for potentio in xmle.findall('potentiometricdata'):
        twidget = widget.add_emf()
        loadEmfXML(twidget, potentio)

    for nmrdat in xmle.findall('nmrdata'):
        twidget = widget.add_nmr()
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
    widget.emf0 = _read_seq(xmle, 'emf0')
    widget.nelectrons = _read_seq(xmle, 'n', dtype=int)
    widget.emf0_error = _read_seq(xmle, 'erroremf0')
    widget.active_species = _read_seq(xmle, 'active', dtype=int)
    widget.flags_emf0 = _read_seq(xmle, 'emf0flags', dtype=int)
    # widget.nelectrons = _read_seq(xmle, 'n', dtype=int)
    # loadCurveXML(widget, xmle)

    nelectrod = widget.nelectrodes
    flat_emf = tuple(_read_seq(xmle, 'emfread'))

    try:
        titre = tuple(_read_seq(xmle, 'titre'))
        exppoints = len(titre)
    except AttributeError:
        exppoints = len(flat_emf) // nelectrod

    emf = np.array(flat_emf).reshape((exppoints, nelectrod))
    widget.nelectrodes = nelectrod
    widget.npoints = exppoints
    # table.setColumnCount(nelectrod+1)
    # table.setRowCount(exppoints)
    with libqt.table_locked(widget.ui.table_data) as table:
        libqt.array2tab(table, emf, col0=1)
    # flat_emf = tuple(_read_seq(xmle, 'emfread'))
    # data_feed = itertools.chain((titre,
    # *[flat_emf[(exppoints*i):(exppoints*(i+1))]for i in range(nelectrod)]))
    # data_feed = itertools.chain((titre, *zip(*[iter(flat_emf)]*exppoints)))
    # widget.reshape_data_table(exppoints, 1+nelectrod)
    # widget.feed_data_table(data_feed)

    # t = xmle.find('key')
    # if t is not None:
    #     mask = tuple(not bool(int(i)) for i in t.text.split())
    # else:
    #     mask = False
    # widget.mask = mask


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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  SAVE ROUTINES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def saveXML(app, f):
    """Save the data into an XML file.

    Parameters:
        main (:class:`apes.MainWindow`): The main application widget.
        f (file or str): A file to save the data into. If file exists it will
            be overwritten.
    """
    root = ET.Element('apes', attrib={'version': '1.0'})
    tree = ET.ElementTree(root)

    metadata = ET.SubElement(root, 'metadata')
    ET.SubElement(metadata, 'title').text = app.project.title
    ET.SubElement(metadata, 'author').text = app.project.author
    ET.SubElement(metadata, 'comments').text = app.project.comments
    ET.SubElement(metadata, 'lastmodified').text = datetime.datetime.now().ctime()
    ET.SubElement(metadata, 'created').text = app.project.created

    ET.SubElement(root, 'labels').text = " ".join(app.modelwidget.labels)

    types = (modelwidget.ModelWidget, simulationwidgets.SpeciationWidget,
             simulationwidgets.TitrationWidget,
             emfwidget.EmfWidget, nmrwidget.NmrWidget,
             specwidget.SpecWidget, calorwidget.CalorWidget,
             otherwidgets.ExternalDataWidget,
             otherwidgets.TitrationBaseWidget)
    calls = (saveModelWidgetXML, saveSpeciationXML, saveTitrationXML,
             saveEmfXML, saveNmrXML, saveSpectrXML, saveCalorXML,
             saveExternalXML, saveTitrationBaseXML)

    for window in app.ui.mdiArea.subWindowList():
        widget = window.widget()
        if type(widget) not in types:
            continue
        which = types.index(type(widget))
        call = calls[which]
        dtag = call(widget)
        root.append(dtag)

    tree.write(f)


def saveFittingXML(twidget):
    if not isinstance(twidget, tabwidget.TabWidget):
        raise TypeError
    xmle = ET.Element('fittinggroup', attrib={'name': twidget.name})

    types = (emfwidget.EmfWidget, nmrwidget.NmrWidget, specwidget.SpecWidget,
             calorwidget.CalorWidget, otherwidgets.TitrationBaseWidget)
    calls = (saveEmfXML, saveNmrXML, saveSpectrXML, saveCalorXML, saveTitrationBaseXML)

    for widget in twidget.widgets_to_save():
        if type(twidget) not in types:
            continue

        which = types.index(type(widget))
        call = calls[which]
        dtag = call(widget)
        xmle.append(dtag)

    return xmle


def saveModelWidgetXML(widget):
    if not isinstance(widget, modelwidget.ModelWidget):
        raise TypeError
    xmle = ET.Element('models')
    xmle.extend([saveModelXML(model) for model in widget])
    return xmle


def saveModelXML(model):
    """Returns an :class:`xml.etree.ElementTree.Element` containing the
    information from this class in order to write in the xml file.

    Parameters:
        model (:class:`modelwidget.ModelData`):

    Returns:
        :class:`xml.etree.ElementTree`: An ElementTree containing the
            information about the model.
    """
    assert isinstance(model, modelwidget.ModelData)
    tm = ET.Element('model')
    tm.attrib['title'] = model.name

    def faux(tag, dat):
        ET.SubElement(tm, tag).text = j(dat)

    faux('beta', model.const)
    faux('error', model.const_error)
    faux('p', itertools.chain(*model.stoich))
    faux('key', model.const_flags)

    # if model.enthalpy is not None:
    #     faux('enthalpy', model.enthalpy)
    #     faux('enthalpy_error', model.enthalpy_error)
    #     faux('enthalpy_key', model.enthalpy_flags)

    return tm


def saveNmrXML(data):
    raise NotImplementedError


def saveEmfXML(widget):
    """Returns an :class:`xml.etree.ElementTree.Element` containing the
    information from this class in order to write in the xml file.

    Returns:
        :class:`xml.etree.ElementTree`: An ElementTree containing the
            information about the model.
    """
    prop = {'unit': 'mV'}

    def j2(x):
        if isinstance(x, (list, tuple)):
            return " ".join(map(str, x))
        else:
            return str(x)

    tcurve = ET.Element('potentiometricdata', {'title': widget.name,
                                               'use': str(widget.use),
                                               'titration': widget.ui.cb_titration.currentText()})
    ET.SubElement(tcurve, 'emf0', prop).text = j2(widget.emf0)
    ET.SubElement(tcurve, 'erroremf0', prop).text = j2(widget.emf0_error)
    ET.SubElement(tcurve, 'active').text = j2(widget.active_species)
    ET.SubElement(tcurve, 'n').text = j2(widget.nelectrons)
    ET.SubElement(tcurve, 'emf0flags').text = j2(widget.emf0_flags)
    # if widget.f is not None:
    #     ET.SubElement(tcurve, 'f').text = j2(widget.f)
    ET.SubElement(tcurve, 'fRTnF').text = j2(widget.fRTnF)

    # saveCurveXML(widget, tcurve)

    ET.SubElement(tcurve, 'emfread').text = j(itertools.chain(*widget.emf))
    ET.SubElement(tcurve, 'key').text = " ".join("0" if i else "1" for i in widget.mask)
    return tcurve


def saveExternalXML(widget):
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
    exdat = ET.Element('externaldata')
    ET.SubElement(exdat, 'title').text = widget.title()
    # TODO add unit if needed
    label = widget.labels()
    ET.SubElement(exdat, 'x').text = j(widget.xdata())
    for n, (data, edata) in enumerate(widget.iter_data()):
        print('>>>', n)
        emldata = ET.SubElement(exdat, 'data', attrib={'order': str(n+1), 'label': next(label)})
        ET.SubElement(emldata, 'y').text = j(data)
        if edata is not None:
            ET.SubElement(emldata, 'ey').text = j(edata)
    return exdat


def saveCurveXML(widget, xmle):
    "Puts self data into XML object"
    def _sub(tag, kw, what):
        ET.SubElement(xmle, tag, kw).text = j(what)

    conc_kw = {'unit': 'mmol'}
    prop2 = {'unit': 'mmol/mL'}

    if hasattr(widget, 'initial_amount'):
        _sub('T0', conc_kw, widget.initial_amount)
        ET.SubElement(xmle, 'buret', prop2).text = \
            j(widget.buret)
        ET.SubElement(xmle, 'T0key').text = j(widget.amount_flags)
    else:
        ET.SubElement(xmle, 'T', prop2).text = j(widget.analyticalc)

    ET.SubElement(xmle, 'V0', {'unit': 'mL'}).text = \
        "%f %f" % (widget.starting_volume, widget.volume_error)

    ET.SubElement(xmle, 'titre').text = j(widget.titre)


def saveCalorXML(calorwidget):
    """Returns an :class:`xml.etree.ElementTree.Element` containing the
    information from this class in order to write in the xml file.

    Returns:
        :class:`xml.etree.ElementTree`: An ElementTree containing the
            information about the calorimetry data.
    """
    tcurve = ET.Element('curve', {'title': calorwidget.name,
                                  'use': str(calorwidget.use)})
    calorwidget.saveCurveXML(tcurve)
    ET.SubElement(tcurve, 'heat').text = j(calorwidget.Q)
    if not np.all(calorwidget.vmask):
        ET.SubElement(tcurve, 'key').text = j(calorwidget.vmask)
    return tcurve


def saveSpeciationXML(widget):
    tdistri = ET.Element('distri', {'active': 'no'})
    props = {'unit': 'mol/L'}
    ET.SubElement(tdistri, 'initial', props).text = \
        j(widget.initial_concentration())
    ET.SubElement(tdistri, 'final', props).text = \
        j(widget.final_concentration())
    ET.SubElement(tdistri, 'pX').text = j(widget.pX())

    # isp = _abool(widget.xscalelog())
    # xref = str(widget.referencex())
    xref, isp = widget.referencex()
    ET.SubElement(tdistri, 'xref', p=_abool(isp)).text = str(xref)

    yref = widget.referencey()
    if yref is not None:
        _yref = str(yref)
    else:
        _yref = '-1'
    ET.SubElement(tdistri, 'yref').text = _yref

    # TODO add external data
    return tdistri


def saveSpectrXML(widget, xmle):
    raise NotImplementedError


def saveTitrationBaseXML(widget):
    """Save titration data into an XML object

    returns an :class:`xml.etree.ElementTree` with the information of
    this :class:`data.SpeciationData` to be saved.

    Parameters:
        widget (:class:`otherwidgets.TitrationBaseWidget`): The widget to
            get the data from.
    Returns:
        :class:`xml.etree.ElementTree.Element`: the XML object.
    """
    assert isinstance(widget, otherwidgets.TitrationBaseWidget)
    xmle = ET.Element('titration', {'name': widget.name})
    ET.SubElement(xmle, 'init', {'unit': 'mmol'}).text = j(widget.initial_amount)
    ET.SubElement(xmle, 'initkey').text = j(widget.init_flags)
    ET.SubElement(xmle, 'buret', {'unit': 'mmol/mL'}).text = j(widget.buret)
    ET.SubElement(xmle, 'buretkey').text = j(widget.buret_flags)
    ET.SubElement(xmle, 'volumeerror').text = str(widget.volume_error)
    ET.SubElement(xmle, 'startingvolume', {'unit':'mL'}).text = str(widget.starting_volume)
    if widget.is_titre_implicit():
        ET.SubElement(xmle, 'finalvolume', {'unit':'mL'}).text = str(widget.final_volume)
        ET.SubElement(xmle, 'totalpoints').text = str(widget.n_points())
    return xmle


def saveTitrationXML(widget):
    """Save titration data into an XML object

    returns an :class:`xml.etree.ElementTree` with the information of
    this :class:`data.SpeciationData` to be saved.

    Parameters:
        widget (:class:`simulationwidgets.TitrationWidget`): The widget to
            get the data from.
    Returns:
        :class:`xml.etree.ElementTree.Element`: the XML object.
    """
    ttitr = ET.Element('simu', {'active': 'no'})
    ET.SubElement(ttitr, 'initial',
                  {'unit': 'mmol'}).text = j(widget.initial_amount())
    ET.SubElement(ttitr, 'buret', {'unit': 'mol/L'}).text = j(widget.buret())
    ET.SubElement(ttitr, 'V').text = "%f %f" % (widget.starting_volume(),
                                                widget.final_volume())
    # TODO add external data
    return ttitr


def _checkXML(xmle, tag):
    if not isinstance(xmle, ET.Element):
        raise ValueError("xmle must be an ElementTree.Element instance")
    if xmle.tag != tag:
        raise ValueError("Wrong XML piece was loaded")
    return


def _read_seq(xmle, x, dtype=float):
    return (dtype(i) for i in xmle.find(x).text.split())


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  IMPORT OTHER FORMATS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def import_pasat(filename):
    """Import data from a file which complies with PASAT file format.

    Parameters:
        filename (str): A readable source from where the data is read.
    Returns:
        dict: Containing all the information read from the file
    """
    titre = []
    emf = []
    with open(filename, 'r', encoding='iso-8859-1') as fhandler:
        for i in range(15):
            fhandler.readline()
        for l in fhandler:
            if l.strip() == '':
                break

            a = l.split()
            titre.append(float(a[0]))
            emf.append(float(a[1]))

    return titre, emf


def import_spectrum_text(filename):
    data = np.loadtxt(filename)
    # first column iw wavelength, next block is data
    return data[:, 0], data[:, 1:]


def import_tiamo(filename):
    """Import data from a file which complies with TIAMO file format.

    Parameters:
        filename (str): A readable source from where the data is read.
    Returns:
        tuple: Containing two lists, being the titre and the potential
    """
    with open(filename, 'r', encoding='iso-8859-1') as fhandler:
        data = np.loadtxt(fhandler, skiprows=5, usecols=(0, 1))
    return data[:, 0].tolist(), data[:, 1].tolist()


def importSuperquadApp(app, filename):
    """Import data from a superquad file.

    Parameters:
        app (mainwindow): instance to the main window.
        filename (str): The file to read data from.
    """
    # from modelwidget import ModelData
    data = importSuperquad(filename)
    app.title = next(data)
    _ = next(data)    # control numbers (unused)
    model = app.ui.tab_main.add_model()
    model.labels = list(next(data))
    app.temperature = next(data)
    logB = next(data)

    model.clear()
    modeldata = model.newModel()
    modeldata.name = 'model #0'
    modeldata.const = logB
    modeldata.stoich = next(data)
    modeldata.const_flags = next(data)
    modeldata.const_error = len(logB)*[0.0]

    # model.append(model)
    model.setCurrentModel(0)

    for emfd in data:
        titr_widget = app.ui.tab_main.add_titrationbase()
        titr_widget.set_volume_implicit(False)
        titr_widget.set_labels(model.labels)
        data_widget = app.ui.tab_main.add_emf()

        # cascade unpacking
        amounts, electr, dataset = emfd
        plot_keys, order, t0, buret, tflag = amounts
        V0, errV, n, hindex, emf0, erremf0 = electr
        V, emf = dataset

        titr_widget.starting_volume = V0
        titr_widget.volume_error = errV
        titr_widget.initial_amount = t0
        titr_widget.buret = buret
        data_widget.emf0 = emf0
        data_widget.emf0_error = erremf0
        data_widget.nelectrons = (n,)
        data_widget.active_species = order.index(hindex)
        data_widget.emf = emf
        data_widget.titre = V


def importSuperquad(filename):
    """Import data from Superquad file.

    This function imports data from a file which complies with SUPERQUAD
    file format.

    Parameters:
        filename (string or file): A readable source from where the data is read.

    Yields:
        * title (str): title of the project
        * control numbers (sequence)
        * labels (list of str): the labels of the principal components
        * temperature (float): the temperature
        * logB (:class:`numpy.ndarray`): the constants in log10 units
        * P (:class:`numpy.ndarray`): the stoichiometric coefficients
        * flags (:class:`numpy.ndarray`): the refinement flags
        * emf (generator of :class:`numpy.ndarray`): the potential read
        * V (generator of :class:`numpy.ndarray`): the volume of titre
        * E0 (generator of floats): the standard potential
        * n (generator of int): the number of electrons involved
        * E0flags (generator of int): the refinement flags for E0
        * error_E0 (generator of float): the error associated to E0
        * V0 (generator of float): the initial volume
        * error_V0 (generator of float): the error associated to volume
        * T0 (generator of :class:`numpy.ndarray`): the initial amounts for the
          principal components
        * T0flags (generator of :class:`numpy.ndarray`): the refinement flags
          for T0.
        * buret (generator of :class:`numpy.ndarray`): the concentration in the
          buret.
        * hindex (generator of int): The index of the electroactive component.
        * fRTnF (generator of float): Nernst's propocionality number.
    """
    def read_amounts(handler):
        plot_keys = []
        order = []
        t0 = []
        buret = []
        tflag = []
        for line in handler:
            if line.strip() == '':
                if len(t0) == 0:
                      return None
                else:
                      return plot_keys, order, t0, buret, tflag
                # break
                # return plot_keys, order, t0, buret, tflag

            aux = line.split()
            plot_keys.append(int(aux[0]))
            order.append(int(aux[1]))
            t0.append(float(aux[2]))
            buret.append(float(aux[3]))
            tflag.append(int(aux[4]))

    def read_electrodes(handler):
        volume, err_volume = map(float, handler.readline().split())
        aux = handler.readline().split()
        n, hindex = map(int, aux[0:2])
        emf0, erremf0 = map(float, aux[2:4])
        assert handler.readline().strip() == ''
        return volume, err_volume, n, hindex, emf0, erremf0

    def read_data(handler):
        aux = []
        for line in handler:
            if line.strip() == '':
                break
            aux.append(tuple(map(float, line.split())))
        return tuple(zip(*aux))

    def read_titration(handler):
        while True:
            amm = read_amounts(handler)
            if amm is None:
                return
            elc = read_electrodes(handler)
            dat = read_data(handler)

            yield (amm, elc, dat)

    with open(filename, "r") as f:
        yield f.readline()      # title
        numbers = tuple(int(i) for i in f.readline().split())
        yield numbers   # control numbers
        num_species = numbers[2]
        yield tuple(f.readline().strip() for i in range(num_species))  # labels
        yield float(f.readline())  # temperature

        B = []
        P = []
        keys = []
        for line in f:
            if line.strip() == '':
                break
            b_, *p_, k_ = line.split()
            B.append(float(b_))
            P.append([int(_) for _ in p_])
            keys.append(int(k_))

        yield B
        yield P
        yield keys
        yield from read_titration(f)


def importHyperquad(filename):
    """Import data from Hyperquad file.

    This function imports data from a Hyperquad file. Those are XML files.

    Parameters:
        filename (str): The file to import data from
    """
    tree = ET.parse(filename)

    ret = {}

    ret['title'] = tree.findtext('model/title')
    S = int(tree.findtext('model/reagents/number'))
    ret['labels'] = [tree.findtext('model/reagents/name%d' % (r+1))
              for r in range(S)]

    E = int(tree.findtext('model/betas/number'))
    ret['betas'] = [float(tree.findtext('model/betas/b%d/value' % (r+1)))
         for r in range(E)]
    ret['stoich'] = [[int(tree.findtext('model/betas/b%d/a%d' % (r+1, c+1)))
         for c in range(S)] for r in range(E)]
    ret['betas_key'] = [int(tree.findtext('model/betas/b%d/refine' % (r+1)))
             for r in range(E)]

    if tree.find('HySS').text: 
        t_V0 = float(tree.findtext('HySS/titration/volume'))
        t_V1 = float(tree.findtext('HySS/titration/final'))
        t_T0 = [float(tree.findtext('HySS/titration/r%d/total' % (r+1)))
                for r in range(S)]
        t_buret = [float(tree.findtext('HySS/titration/r%d/burette' % (r+1)))
                   for r in range(S)]
        titration = {'V0': t_V0, 'V1': t_V1, 'T0': t_T0, 'buret': t_buret}

        s_V0 = float(tree.findtext('HySS/simulation/volume'))
        s_T0 = [float(tree.findtext('HySS/simulation/r%d/initial' % (r+1)))
                for r in range(S)]
        s_T1 = [float(tree.findtext('HySS/simulation/r%d/final' % (r+1)))
                for r in range(S)]
        simulation = {'V0': s_V0, 'T0': s_T0, 'T1': s_T1}
    else:
        titration = None
        simulation = None

    emfdata = tree.find('potentiometric_data')
    if emfdata is not None:
        retemf = []
        prev_point = 0
        numcurves = int(emfdata.findtext('curves'))
        for et in (emfdata.find('curve%d' % (_+1)) for _ in range(numcurves)):
            thiscurve = {}
            num_points = int(et.findtext('points'))
            thiscurve['startingvolume'] = float(et.findtext('vinit'))
            thiscurve['sigma_v'] = float(et.findtext('sigmav'))
            thiscurve['initial amount'] = [float(et.findtext('r%d/totmm' % (r+1))) for r in range(S)]
            thiscurve['buret'] = [float(et.findtext('r%d/addc' % (r+1))) for r in range(S)]
            thiscurve['emf0'] = float(et.findtext('e1/ezero'))
            thiscurve['err_emf0'] = float(et.findtext('e1/sigmae'))
            thiscurve['n_electrode'] = int(et.findtext('e1/ifc'))

            thiscurve['titre'] = [float(et.findtext('p%d/titre' % (r+1))) for r in range(prev_point, num_points)]
            thiscurve['emf'] = [float(et.findtext('p%d/emf1' % (r+1))) for r in range(prev_point, num_points)]
            retemf.append(thiscurve)
            prev_point = num_points

    retval = {'title': ret['title'], 'labels': ret['labels'], 
            'logB': ret['betas'], 'P': ret['stoich'],
            'Bkeys': ret['betas_key'], 'emfdata': retemf}
    if simulation is not None:
        retval['simulation'] = simulation
    if titration is not None:
        retval['titration'] = titration

    return retval


def importK88(filename):
    """Import data from a file which complies with K88 file format.

    K88 files are text files with the format
    10.096     0.000     0.000    40.344  1211.900A1
    10.096     0.000     0.000    37.899  1226.900      0.330

    Parameters:
        filename (str): The file to import data from
    Returns:
        tuple: of lists containing the following information information:

        - list: the amount of titre at the beginning of the point
        - list: the amount of titre at the end of the point
        - list: the id of the point
        - list: the heat measured
    """
    A, B, idp, V1, V2, Q = [], [], [], [], [], []
    flag = True
    with open(filename, 'r') as fhandler:
        for l in fhandler.readlines():
            if flag:
                A.append(list(map(float, l[0:40].split())))
                V1.append(float(l[41:50]))
                idp.append(l[50:].strip())
            else:
                B.append(list(map(float, l[0:40].split())))
                V2.append(float(l[41:50]))
                Q.append(float(l[50:]))
            flag = not flag

    return A, B, V1, V2, idp, Q


def importHyperquadApp(app, f):
    '''Import data from Hyperquad data into application.

    .. note:: This function assumes that the main application is empty.
    '''
    # from othermodules import SpeciationWidget
    data = importHyperquad(f)
    n_equil = len(data['logB'])
    n_species = len(data['labels'])
    model = app.modelwidget
    model.clear(n_equil, n_species)
    model.beta_raw=data['logB']
    model.beta_error=n_equil*[0.0]
    model.stoich=data['P']
    model.beta_flags=data['Bkeys']
    # app.modelwidget[0] = model
    # E, S = data['logB'].shape # model.stoich.shape

    app.title = data['title']
    app.modelwidget.labels = data['labels']

    if 'titration' in data:
        widget = app.newTitration()
        titration = data['titration']
        widget.set_starting_volume(titration['V0'])
        widget.set_final_volume(titration['V1'])
        widget.set_initial_amount(titration['T0'])
        widget.set_buret(titration['buret'])

    if 'simulation' in data:
        widget = app.newSpeciation()
        simulation = data['simulation']
        widget.set_initial_concentration(simulation['T0'])
        widget.set_final_concentration(simulation['T1'])

    # for i in app.data.iterSimu():
    #     sdw = simulationwidgets.SpeciationWidget(i)
    #     app.ui.tab_main.addTab(sdw, "Speciation")

    # for i in app.iterEmfDs():
    #     app.newTabEmfDs()

    # # TODO include other dataset types here.

    # app.set_mode(consts.FM_EMF)
    # # app._fmode = consts.FM_EMF
    # # app.__notifymode()


def importK88App(app):
    """Import calorimetry data from a K88 file."""
    filters = "k88 Files (*.k88);;All Files (*.*)"
    fhandler, ok = app.open_dialog(filters)
    if not ok:
        return

    with open(fhandler, 'r') as fi:
        C1, __, V1, V2, __, Q = importK88(fi)

    n = len(C1[0])
    if app.S < n:
        msg = ("The data file read contains {} components, while the "
               "current model contains only {}. Do you want to insert "
               "new {} component(s)?").format(n, app.S, n - app.S)
        if libqt.confirm_delete(app, msg):
            t = "ABCDEFGHIJKLMNOPQRSTUVWZYZ"
            q = 0
            for i in range(n-app.data.S):
                while t[q] in app.data.labels:
                    q += 1
                    if q > len(t):
                        raise RuntimeError
                app.data.addReagent(t[q])

    T0 = C1[0]
    aux1 = np.gradient(np.array(C1), axis=0)
    DV = np.array(V2) - np.array(V1)
    buret = (aux1/np.tile(DV/1000.0, (4, 1)).T)[0] / 1000.0
    V0 = V1[0]/1000.0
    V = np.array([0.0] + [v/1000.0 - V0 for v in V2])

    w = app.ui.tab_main.currentWidget()
    d = w.data
    d.Q = np.array(Q)
    d.T0 = np.array(T0)
    d.V0 = V0
    d.V = V
    d.buret = buret

    app.model_widget.refreshWidgets()
    w.refreshWidgets()
    app.canvas.plotCalorFit()
    app.set_mode(consts.FM_CALOR)
    # app._fmode = consts.FM_CALOR
    # app.__notifymode()


def importPasatApp(app):
    "import data from a PASAT file into current dataset and update."

    widget = app.ui.tab_main.currentWidget()
    assert isinstance(widget, emfwidget.EmfWidget)

    filters = "PASAT Files (*.ptr);;All Files (*.*)"
    fhandler, ok = app.open_dialog(filters)
    if not ok:
        return

    titre, emf = import_pasat(fhandler)
    widget.data.V = titre
    widget.data.emf_read = emf
    # widget.refreshWidgets()
    app.set_mode(consts.FM_EMF)


def importTiamoApp(app):
    "import data from a TIAMO file into current dataset and update."

    widget = app.ui.tab_main.currentWidget()
    assert isinstance(widget, emfwidget.EmfWidget)
    filters = "TIAMO Files (*.csv);;All Files (*.*)"
    fhandler, ok = app.open_dialog(filters)
    if not ok:
        return

    titre, emf = import_tiamo(fhandler)
    widget.data.V = titre
    widget.data.emf_read = emf
    # widget.refreshWidgets()
    app.set_mode(consts.FM_EMF)


def importTxtSpApp(app):
    """Import spectra from data file and load into current widget."""

    assert app.currentTabType() == 'spectrometry'
    widget = app.ui.tab_main.currentWidget()
    filters = "Text Files (*.txt);;Data files (*.dat);;All Files (*.*)"
    fhandler, ok = app.open_dialog(filters)
    if not ok:
        return
    data = widget.data
    data.addSpectraFromFile(fhandler)
    # widget.refreshWidgets()
    app.set_mode(consts.FM_SPEC)


def j(a):
    'A shortcut.'
    return " ".join(map(str, a))


def _bool(x):
    return x.strip().lower() in ('true', 'yes', '1')


def _abool(x):
    return 'yes' if x else 'no'


def __read_floats(xmle, tag):
    return map(float, xmle.find(tag).text.split())


def __read_ints(xmle, tag):
    return map(int, xmle.find(tag).text.split())
