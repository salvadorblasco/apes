import xml.etree.ElementTree as ET


def import_hyperquad_data(filename):
    """Import data from Hyperquad file.

    This function imports data from a Hyperquad file. Those are XML files.

    HQD files are usually malformed because they contain two root elements. An exception
    in triggered. This implementation removes the <HySS> tag to avoid this error.

    Parameters:
        filename (str): The file to import data from
    """
    try:
        tree = ET.parse(filename)
    except ET.ParseError:
        with open(filename, 'r') as f:
            text = f.read().replace('<HySS></HySS>', '')
        tree = ET.fromstring(text)

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

    if (t_hyss := tree.find('HySS') is not None) and t_hyss.text:
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
            thiscurve['starting volume'] = float(et.findtext('vinit'))
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


def import_hyperquad_app(app, f):
    '''Import data from Hyperquad data into application.

    .. note:: This function assumes that the main application is empty.
    '''
    # from othermodules import SpeciationWidget
    data = import_hyperquad_data(f)
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

