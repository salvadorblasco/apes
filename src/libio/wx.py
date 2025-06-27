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


def _abool(x):
    return 'yes' if x else 'no'


def _checkXML(xmle, tag):
    if not isinstance(xmle, ET.Element):
        raise ValueError("xmle must be an ElementTree.Element instance")
    if xmle.tag != tag:
        raise ValueError("Wrong XML piece was loaded")
    return


def _read_seq(xmle, x, dtype=float):
    return (dtype(i) for i in xmle.find(x).text.split())


def j(a):
    'A shortcut.'
    return " ".join(map(str, a))
