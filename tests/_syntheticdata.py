import sys
sys.path.append('../src')

import consts


def load_lmh():
    import data_lmh

    from modelwidget import ModelWidget, ModelData
    model = ModelWidget()
    model.clear(n_equils=data_lmh.N_EQUILS, n_species=data_lmh.N_COMPON)
    mdata = ModelData(n_equils=data_lmh.N_EQUILS, n_species=data_lmh.N_COMPON)
    mdata.stoich = data_lmh.stoich
    mdata.const = data_lmh.logbeta
    mdata.const_flags = 6*[consts.RF_CONSTANT] + 3*[consts.RF_REFINE] + [consts.RF_CONSTANT]
    model.append(mdata)
    model.setCurrentModel(-1)

    from otherwidgets import TitrationBaseWidget
    titr1 = TitrationBaseWidget(model)
    titr1.name = 'Titration 1'
    titr1.initial_amount = data_lmh.t1_init
    titr1.buret = data_lmh.t1_buret
    titr1.starting_volume = data_lmh.t1_startingvol
    titr1.final_volume = data_lmh.t1_endvol
    titr1.n_points = 101

    titr2 = TitrationBaseWidget(model)
    titr2.name = 'Titration 2'
    titr2.initial_amount = data_lmh.t2_init
    titr2.buret = data_lmh.t2_buret
    titr2.starting_volume = data_lmh.t2_startingvol
    titr2.final_volume = data_lmh.t2_endvol
    titr2.n_points = 101

    from emfwidget import EmfWidget
    emfw1 = EmfWidget(model)
    emfw1.emf0 = data_lmh.t1_emf0
    emfw1.slope = (1.0, 1.0)
    emfw1.slope_flags = (0, 0)
    emfw1.nelectrons = (1, 1)
    emfw1.emf0_error = (0.01,0.01)
    emfw1.active_species = (2, 1)
    emfw1.emf0_flags = (0,0)
    emfw1.titration = titr1
    emfw1.titre = titr1.titre
    emfw1.emf = data_lmh.t1_emf

    emfw2 = EmfWidget(model)
    emfw2.emf0 = data_lmh.t2_emf0
    emfw2.slope = (1.0, 1.0)
    emfw2.slope_flags = (0, 0)
    emfw2.nelectrons = (1, 1)
    emfw2.emf0_error = (0.01,0.01)
    emfw2.active_species = (2, 1)
    emfw2.emf0_flags = (0,0)
    emfw2.titration = titr2
    emfw2.emf = data_lmh.t2_emf

    from bridge import Parameters
    b = Parameters(model, [titr1, titr2], [emfw1, emfw2])
    return data_lmh, b


def load_hexaprotic():
    import hexaprotic

    from modelwidget import ModelWidget, ModelData
    model = ModelWidget()
    mdata = ModelData(n_equils=7, n_species=2)
    mdata.stoich = hexaprotic.stoich
    mdata.const = hexaprotic.logbeta
    mdata.const_flags = 6*[consts.RF_REFINE] + [consts.RF_CONSTANT]
    model.append(mdata)
    model.setCurrentModel(-1)

    from otherwidgets import TitrationBaseWidget
    titr = TitrationBaseWidget(model)
    titr.name = 'Titration'
    titr.initial_amount = hexaprotic.init
    titr.buret = hexaprotic.buret
    titr.starting_volume = hexaprotic.v0
    titr.titre = hexaprotic.titre.tolist()

    from emfwidget import EmfWidget
    emfw = EmfWidget(model)
    emfw.emf0 = hexaprotic.emf0
    emfw.emf = hexaprotic.emf
    emfw.titration = titr

    from bridge import Parameters
    b = Parameters(model, [titr], [emfw])
    return hexaprotic, b

