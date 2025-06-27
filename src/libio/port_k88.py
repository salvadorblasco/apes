import numpy as np

import libqt


def import_K88_data(filename):
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


def import_K88_app(app):
    """Import calorimetry data from a K88 file."""
    filters = "k88 Files (*.k88);;All Files (*.*)"
    fhandler, ok = app.open_dialog(filters)
    if not ok:
        return

    with open(fhandler, 'r') as fi:
        C1, __, V1, V2, __, Q = import_K88_data(fi)

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

