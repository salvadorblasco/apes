#!/usr/bin/python3

import random
import sys
import unittest

from PyQt5.QtWidgets import QApplication

sys.path.append('../src/')
import apes

app = QApplication(sys.argv)


class CanvasTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        import canvas
        self.main = apes.MainWindow()
        self.canvas = canvas.MyCanvas()

    def tearDown(self):
        self.main.close()

    def test_set_style(self):
        for i in range(-1, 6):
            with self.subTest(i=i):
                if i < 0 or i > 4:
                    with self.assertRaises(ValueError):
                        self.canvas.setStyle(i)
                else:
                    self.canvas.setStyle(i)
                    self.assertEqual(self.canvas.current_style, i)

    def test_set_color(self):
        for i in range(-1, 6):
            with self.subTest(i=i):
                if i < 0 or i > 4:
                    with self.assertRaises(ValueError):
                        self.canvas.setColor(i)
                else:
                    self.canvas.setColor(i)
                    self.assertEqual(self.canvas.current_color, i)


class MainWindowTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.main = apes.MainWindow()

    def tearDown(self):
        self.main.close()

    def test_properties(self):
        properties = ('modelwidget', 'canvas', 'project', 'ui')
        msg = "property '{}' not found"
        for prop in properties:
            self.assertTrue(hasattr(self.main, prop), msg.format(prop))

        # with self.assertRaises(AttributeError):
        #     self.main.modelwidget = None

        self.main.modelwidget.temperature = 25
        self.assertTrue(self.main.modelwidget.temperature == 25)

        with self.assertRaises(ValueError):
            self.main.modelwidget.temperature = -20

        title_test = "title test"
        self.main.title = title_test
        self.assertEqual(self.main.title, title_test)


class ModelTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.mainw = apes.MainWindow()
        self.model = self.mainw.modelwidget

    def tearDown(self):
        self.mainw.close()

    def test_add_component(self):
        nc1 = self.model.number_components
        i = random.randint(0, nc1)
        self.model.addComponent('X', position=i)
        self.assertEqual(self.model.labels[i], 'X', 
                         msg=f"i={i}, labels={self.model.labels}")
        self.assertEqual(self.model.number_components, nc1+1)

    def test_add_equilibrium(self):
        neqs = self.model.number_equilibria
        i = random.randint(0, neqs)
        self.model.addEquilibrium(position=i)
        self.assertEqual(self.model.number_equilibria, neqs+1)

    def test_clear(self):
        # TODO self.model.clear()
        pass

    def test_equilibrium(self):
        "test addEquilibrium and removeEquilibrium."
        num_equil = self.model.number_equilibria
        self.model.addEquilibrium()
        self.assertEqual(self.model.number_equilibria, num_equil + 1)
        self.model.addEquilibrium()
        self.assertEqual(self.model.number_equilibria, num_equil + 2)
        self.model.removeEquilibrium(confirm=False)
        self.assertEqual(self.model.number_equilibria, num_equil + 1)

    def test_new_model(self):
        # TODO self.model.newModel()
        pass

    def test_properties(self):
        # TODO check property returns
        properties = ('beta', 'beta_error', 'beta_flags', 'number_equilibria',
                      'number_components',  'labels', 'modelname',
                      'stoich', 'temperature')
        msg = "property '{}' not found"
        for prop in properties:
            self.assertTrue(hasattr(self.model, prop), msg.format(prop))

    @unittest.skip('broken')
    def test_tablerw(self):
        import numpy as np
        import consts
        compare_arrays = np.testing.assert_array_almost_equal

        def generate_model(iscal, num_comps, num_samples):
            for E in np.random.randint(2, high=15, size=(num_resamples,)):
                # self.assertEqual(n, len(self.model))
                if iscal:
                    _h = np.random.rand(E)
                    _eh = np.random.rand(E)
                    _fh = np.random.randint(low=0, high=7+1, size=(E,))
                    assert(np.all(np.logical_and(_fh <= 7, _fh >= 0)))
                    self.model.addCalorimetry()
                else:
                    _h = None
                    _eh = None
                    _fh = None
                    self.model.removeCalorimetry()
                _p = np.random.randint(low=-5, high=11, size=(E, S))
                _fb = np.random.randint(low=-1, high=7+1, size=(E,))
                assert(np.all(np.logical_and(_fb <= 7, _fb >= -1)))
                newmodel = consts.Model(name='No name',
                                        const=np.random.rand(E),
                                        stoich=_p,
                                        const_flags=_fb,
                                        const_error=np.random.rand(E),
                                        enthalpy=_h,
                                        enthalpy_error=_eh,
                                        enthalpy_flags=_fh)
                yield newmodel

        num_samples = 20
        num_resamples = 10
        labels = [letter for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ']
        for iscal in (False, True):
            for S in np.random.randint(low=2, high=15, size=(num_samples,)):
                self.model.clear()
                self.model.labels = labels[:S]
                modelgen = generate_model(iscal, S, num_resamples)
                for n, model in enumerate(modelgen):
                    self.model.mainwend(model)
                    self.assertEqual(n+1, len(self.model))
                    self.model.setCurrentModel(-1)
                    self.assertEqual(iscal, self.model.isCalori())

                    flagged = model.const_flags != -1
                    compare_arrays(self.model.beta, model.const[flagged])
                    compare_arrays(self.model.beta_error,
                                   model.const_error[flagged])
                    compare_arrays(self.model.stoich,
                                   model.stoich[np.nonzero(flagged)[0], :])
                    # if np.any(self.model.beta_flags != model.const_flags[flagged]):
                    #     import pdb
                    #     pdb.set_trace()
                    np.testing.assert_equal(self.model.beta_flags,
                                            model.const_flags[flagged])
                    if iscal:
                        np.testing.assert_equal(self.model.enthalpy_flags,
                                                model.enthalpy_flags[flagged])
                        compare_arrays(self.model.enthalpy, model.enthalpy[flagged])
                        compare_arrays(self.model.enthalpy_error,
                                       model.enthalpy_error[flagged])
                    else:
                        self.assertIsNone(self.model.enthalpy)
                        self.assertIsNone(self.model.enthalpy_error)
                        self.assertIsNone(self.model.enthalpy_flags)
                        self.assertIsNone(self.model.enthropy)
                        self.assertIsNone(self.model.enthropy_error)


class EmfTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        self.main = apes.MainWindow()
        # self.main.show()
        # QTest.qWaitForWindowShown(self.main)
        self.emf = self.main.newEmf()
        # n_data = 20
        # import random
        # self.emf.emf = [random.random() for r in range(n_data)]
        # self.emf.emf0 = random.random()
        # self.emf.emf0_error = random.random()

    def test_dummy(self):
        return True

    def test_emf(self):
        n_data = 20
        import random
        made_emf = [random.random() for r in range(n_data)]
        self.emf.emf = made_emf
        read_emf = self.emf.emf
        # import pudb
        # pudb.set_trace()
        # np.testing.assert_array_almost_equal(made_emf, read_emf)
        for a, b in zip(made_emf, read_emf):
            self.assertAlmostEqual(a, b[0], places=3)

    def test_emf0(self):
        import random
        made_emf0 = [100*random.random()]
        self.emf.emf0 = made_emf0
        read_emf0 = tuple(self.emf.emf0)
        self.assertAlmostEqual(made_emf0[0], read_emf0[0])

        self.emf._nelectrs_changed(2)
        made_emf0 = 100*random.random(), 100*random.random()
        self.emf.emf0 = made_emf0
        read_emf0 = self.emf.emf0
        self.assertTupleEqual(made_emf0, read_emf0)

    @unittest.skip('Fails')
    def test_properties(self):
        # TODO check property returns
        properties = ('active_species', 'emf', 'emf0_error', 'emf0',
                      'initial_amount', 'emf0_flags', 'nelectrodes',
                      'nelectrons', 'starting_volume', 'volume_error',
                      'titre')
        msg = "property '{}' not found"
        for prop in properties:
            self.assertTrue(hasattr(self.emf, prop), msg.format(prop))

    def tearDown(self):
        self.main.close()


class ExternalDataTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self):
        from otherwidgets import ExternalDataWidget
        self.mainw = apes.MainWindow()
        self.widget = ExternalDataWidget

    def tearDown(self):
        self.mainw.close()

    def test_general(self):
        pass


if __name__ == '__main__':
    unittest.main()
