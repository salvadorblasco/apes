def import_tiamo_data(filename):
    """Import data from a file which complies with TIAMO file format.

    Parameters:
        filename (str): A readable source from where the data is read.
    Returns:
        tuple: Containing two lists, being the titre and the potential
    """
    with open(filename, 'r', encoding='iso-8859-1') as fhandler:
        data = np.loadtxt(fhandler, skiprows=5, usecols=(0, 1))
    return data[:, 0].tolist(), data[:, 1].tolist()


def import_tiamo_app(app):
    "import data from a TIAMO file into current dataset and update."

    widget = app.ui.tab_main.currentWidget()
    assert isinstance(widget, emfwidget.EmfWidget)
    filters = "TIAMO Files (*.csv);;All Files (*.*)"
    fhandler, ok = app.open_dialog(filters)
    if not ok:
        return

    titre, emf = import_tiamo_data(fhandler)
    widget.data.V = titre
    widget.data.emf_read = emf
    # widget.refreshWidgets()
    #app.set_mode(consts.FM_EMF)

