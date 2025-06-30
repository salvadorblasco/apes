#!/usr/bin/python3

import random
import unittest

import sys
from PyQt5 import QtWidgets

sys.path.append('../src')
import libio
import libio.wx
import libio.rx


def _xml(xml_text):
    from xml.etree import ElementTree
    xmle = ElementTree.fromstring(xml_text)
    return xmle


def _gen_test_model(n_species=3, n_equils=4):
    import modelwidget

    def _gen_stoich():
        return [[random.randint(0, 3) for _ in range(n_species)]
                for _ in range(n_equils)]

    model1 = modelwidget.ModelData(n_equils, n_species)
    model1.name = 'model 1'
    model1.stoich = _gen_stoich()
    model1.const = [1.0, 2.0, 3.0, 4.0]
    model1.const_error = [0.1, 0.2, 0.3, 0.4]

    model2 = modelwidget.ModelData(n_equils=4, n_species=3)
    model2.name = 'model 2'
    model2.stoich = _gen_stoich()
    model2.const = [1.0, 2.0, 3.0, 4.0]
    model2.const_error = [0.1, 0.2, 0.3, 0.4]

    widget = modelwidget.ModelWidget()
    widget.clear()
    widget.append(model1)
    widget.append(model2)
    return widget


def _gen_rand_float(size=10):
    yield from (random.random() for i in range(size))


def _validate_units(*args):
    valid_units = ('mV', 'chemical shift', 'cal')
    for arg in args:
        if arg not in valid_units:
            raise ValueError


class TestLoad(unittest.TestCase):
    def test_load(self):
        import mainwindow
        fname = '_load_all1.xml'
        # fname = '../data/pylh.xml'
        app = QtWidgets.QApplication(sys.argv)
        mainw = mainwindow.MainWindow()
        libio.loadXML(mainw, fname)

        self.assertListEqual(mainw.modelwidget.labels, ['L', 'H'])
        self.assertEqual(mainw.modelwidget.temperature, 298.15)
        metadata = {
            'author': 'author 1',
            'title': 'my title',
            'comments': 'comment 1',
            'created': '07/02/2016 15:25',
            'last_modified': '07/02/2016 15:25'
        }
        for key, val in metadata.items():
            self.assertEqual(getattr(mainw.project, key), val)

    def test_load_model(self):
        text = """<model title="title example" energy_units='kcal/mol'>
           <beta>9.25 10.19</beta>
           <error>0.06 0.01</error>
           <p>1 1 0 -1</p>
           <key>1 0</key>
           <enthalpy>-11.20 0.28</enthalpy>
           <enthalpy_error>0.06 0.01</enthalpy_error>
           <enthalpy_key>1 0</enthalpy_key>
           <enthropy>-0.22 0.19</enthropy>
           <enthropy_error>0.06 0.09</enthropy_error>
          </model>"""
        xmle = _xml(text)
        model = libio.rx.loadModelXML(xmle)

        self.assertEqual(model.name, "title example")
        self.assertListEqual(model.const, [9.25, 10.19])
        self.assertListEqual(model.stoich, [[1, 1],[0, -1]])
        self.assertListEqual(model.const_flags, [1, 0])
        self.assertListEqual(model.enthalpy, [-11.20, 0.28])
        self.assertListEqual(model.enthalpy_error, [0.06, 0.01])
        self.assertListEqual(model.enthalpy_flags, [1, 0])

    @unittest.skip('DEBUG!!!')
    def test_load_models(self):
        app = QtWidgets.QApplication(sys.argv)
        import modelwidget
        text = """<models active="1">
           <model title="title example" energy_units='kcal/mol'>
            <beta>9.25 10.19</beta>
            <error>0.06 0.01</error>
            <p>1 1 0 -1</p>
            <key>1 0</key>
          </model>
           <model title="title example" energy_units='kcal/mol'>
            <beta>8.25 11.19</beta>
            <error>0.02 0.03</error>
            <p>1 1 2 2</p>
            <key>0 1</key>
           </model>
          </models>"""
        xmle = _xml(text)
        widget = modelwidget.ModelWidget()
        labels = ['A', 'B']
        libio.loadModelsXML(widget, xmle, labels)

        self.assertTrue(len(widget) == 2)
        self.assertTupleEqual(tuple(widget.beta), (8.25, 11.19))

    @unittest.skip('DEBUG!!!')
    def test_load_emf(self):
        import emfwidget
        with open('_load_emf.xml', 'r') as f:
            text = f.read()
        xmle = _xml(text)
        app = QtWidgets.QApplication(sys.argv)
        modelw = _gen_test_model()
        modelw.setCurrentModel(0)
        emfw = emfwidget.EmfWidget(modelw)

        libio.loadEmfXML(emfw, xmle)

        self.assertTupleEqual(emfw.emf0, (359.32,))
        self.assertTupleEqual(emfw.nelectrons, (1,))
        self.assertTupleEqual(emfw.emf0_error, (0.03,))
        self.assertTupleEqual(emfw.emf, tuple((10.0,) for _ in range(11)))
        self.assertTupleEqual(tuple(emfw.mask), 11*(False,))

        self.assertEqual(emfw.npoints, 11)
        self.assertEqual(emfw.starting_volume, 30.0)
        self.assertEqual(emfw.volume_error, 0.01)
        self.assertTupleEqual(tuple(emfw.titre), tuple(_/10 for _ in range(11)))
        # TODO Fails
        self.assertTupleEqual(emfw.active_species, (1,))
        # TODO f fRTnF T0 T0key buret

    @unittest.skip('DEBUG!!')
    def test_load_external(self):
        app = QtWidgets.QApplication(sys.argv)
        import otherwidgets

        DATA_MAX = 20
        DATA_MIN = 1

        for i in range(3):
            with self.subTest(i=i):
                filename = '_external%d.xml' % i
                edwidget = otherwidgets.ExternalDataWidget()
                with open(filename, 'r') as f:
                    n_orders = int(f.readline())
                    keys = [f.readline().split() for line in range(n_orders)]
                    raw_text = f.read()

                lst = []
                for lvl1 in keys:
                    lng = random.randint(DATA_MIN, DATA_MAX)
                    for lvl2 in lvl1:
                        lst.append(tuple(_gen_rand_float(lng)))

                feed = (" ".join(l) for l in lst)
                text = raw_text.format(*feed)
                xmle = _xml(text)
                libio.loadExternalXML(edwidget, xmle)

                data = edwidget.get_data()
                data_test1 = {'0': [('x', x), ('y', y), ('ey', ey)]}
                self.assertDictEqual(data, data_test1)

    @unittest.skip('DEBUG!!')
    def test_load_external2(self):
        # TODO delete when the other is complete
        # import pudb
        # pudb.set_trace()
        app = QtWidgets.QApplication(sys.argv)
        import otherwidgets
        edwidget = otherwidgets.ExternalDataWidget()
        x = tuple(_gen_rand_float())
        y = tuple(_gen_rand_float())
        ey = tuple(_gen_rand_float())
        _x = " ".join(map(str, x))
        _y = " ".join(map(str, y))
        _ey = " ".join(map(str, ey))
        text = """
         <externaldata>
          <title>y2-axis title</title>
          <!-- order="0" applies to all -->
          <data type="x"  unit="mL">{}</data>
          <data type="y"  suborder="1" unit="mV" label="δ / ppm">{}</data>
          <data type="ey" suborder="1" unit="mV" label="δ / ppm">{}</data>
         </externaldata>""".format(_x, _y, _ey)
        xmle = _xml(text)
        libio.loadExternalXML(edwidget, xmle)

        data = edwidget.get_data()
        data_test1 = {'0': [('x', x), ('y', y), ('ey', ey)]}
        self.assertDictEqual(data, data_test1)

    @unittest.skip('not implemented')
    def test_load_nmr(self):
        pass

    @unittest.skip('not implemented')
    def test_load_calor(self):
        pass

    @unittest.skip('DEBUG')
    def test_load_speciation(self):
        import otherwidgets
        app = QtWidgets.QApplication(sys.argv)
        text = """<distri name="noname" points="100">
          <initial unit="mol/L">0.001 0.001 2</initial>
          <final unit="mol/L">0.001 0.001 11</final>
          <pX>False False True</pX>
          <xref p="yes">2</xref>
          <yref>1</yref>
         </distri>"""
        xmle = _xml(text)
        modelw = _gen_test_model()
        modelw.setCurrentModel(0)
        widget = otherwidgets.SpeciationWidget(modelw)

        libio.loadSpeciationXML(widget, xmle)

        self.assertEqual(widget.n_points(), 100)
        analc = widget.analyticalc()
        self.assertTupleEqual(analc[0], (0.001, 0.001, 1e-2))
        self.assertTupleEqual(analc[-1], (0.001, 0.001, 1e-11))
        self.assertTupleEqual(widget.initial_concentration(),
                              (0.001, 0.001, 2))
        self.assertTupleEqual(widget.final_concentration(),
                              (0.001, 0.001, 11))
        self.assertTupleEqual(widget.pX(), (False, False, True))
        self.assertEqual(widget.referencex(), 2)
        self.assertEqual(widget.referencey(), 1)

    @unittest.skip('not implemented')
    def test_load_spectr(self):
        pass

    @unittest.skip('not implemented')
    def test_load_titration(self):
        pass


class TestSave(unittest.TestCase):
    @unittest.skip('not implemented')
    def test_save(self):
        pass

    @unittest.skip('not implemented')
    def test_save_model_widget(self):
        pass

    def test_save_model(self):
        from modelwidget import ModelData
        n_equils = 4
        n_species = 3
        model = ModelData(n_equils, n_species)
        model.name = 'test name'
        model.const = [1.0, 2.0, 3.0, 4.0]
        model.stoich = [[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4]]
        model.const_error = [0.1, 0.2, 0.3, 0.4]

        xmle = libio.wx.saveModelXML(model)
        self.assertEqual('test name', xmle.attrib['title'])
        _read_beta = [float(number) for number in xmle.find('beta').text.split()]
        self.assertListEqual(model.const, _read_beta)
        _read_berror = [float(number) for number in xmle.find('error').text.split()]
        self.assertListEqual(model.const_error, _read_berror)

        p = list(map(int, xmle.find('p').text.split()))
        self.assertListEqual(p, [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4])

        flags = [int(i) for i in xmle.find('key').text.split()]
        self.assertListEqual(flags, [0, 0, 0, 0])

    @unittest.skip('not implemented')
    def test_save_nmr(self):
        pass

    @unittest.skip('not implemented')
    def test_save_emf(self):
        pass

    @unittest.skip('not implemented')
    def test_save_external(self):
        pass

    @unittest.skip('not implemented')
    def test_save_curve(self):
        pass

    @unittest.skip('not implemented')
    def test_save_calor(self):
        pass

    @unittest.skip('not implemented')
    def test_save_speciation(self):
        pass

    @unittest.skip('not implemented')
    def test_save_spectr(self):
        pass

    @unittest.skip('not implemented')
    def test_save_titration(self):
        pass


if __name__ == '__main__':
    unittest.main()
