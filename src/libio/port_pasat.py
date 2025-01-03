def import_pasat_data(filename):
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


def import_pasat_app(app):
    "import data from a PASAT file into current dataset and update."

    widget = app.ui.tab_main.currentWidget()
    assert isinstance(widget, emfwidget.EmfWidget)

    filters = "PASAT Files (*.ptr);;All Files (*.*)"
    fhandler, ok = app.open_dialog(filters)
    if not ok:
        return

    titre, emf = import_pasat_data(fhandler)
    widget.data.V = titre
    widget.data.emf_read = emf
    # widget.refreshWidgets()
    # app.set_mode(consts.FM_EMF)
