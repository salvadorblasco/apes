import libqt


def import_superquad_app(app, filename: str):
    """Import data from a superquad file.

    Parameters:
        app (mainwindow): instance to the main window.
        filename (str): The file to read data from.
    """
    app.blockSignals(True)
    # from modelwidget import ModelData
    # breakpoint()
    data = import_superquad_data(filename)
    app.title = next(data)
    _ = next(data)    # control numbers (unused)
    model = app.new_model()
    labels = list(next(data))
    app.temperature = next(data)
    logB = next(data)
    model.clear(n_equils = len(logB), n_species=len(labels))
    model.labels = labels

    model.ui.table_model.setRowCount(len(logB))
    model.beta_raw = logB
    model.stoich = next(data)
    model.beta_flags = next(data)
    model.beta_error = len(logB) * [0.0]
    model.modelname = 'model #0'

    fgroup = app.new_fitting_group()
    data_titr_pairs = []

    for emfd in data:
        # cascade unpacking
        amounts, electr, dataset = emfd
        plot_keys, order, t0, buret, tflag = amounts
        V0, errV, n, hindex, emf0, erremf0 = electr
        V, emf = dataset

        titr_widget = fgroup.add_titration()
        with libqt.signals_blocked(titr_widget):
            titr_widget.set_volume_implicit(False)
            titr_widget.titre = V

        titr_widget.set_labels(model.labels)
        data_widget = fgroup.add_emf()
        titr_widget.final_volume = V0 + V[-1]
        titr_widget.starting_volume = V0
        titr_widget.volume_error = errV
        titr_widget.initial_amount = t0
        titr_widget.buret = buret
        data_widget.titre = titr_widget.titre
        data_widget.emf0 = (emf0,)
        data_widget.emf0_error = (erremf0,)
        data_widget.nelectrons = (n,)
        data_widget.active_species = (order.index(hindex),)
        data_widget.emf = emf
        data_titr_pairs.append((data_widget, titr_widget))

    for data_widget, titr_widget in data_titr_pairs:
        data_widget.titration = titr_widget



    app.blockSignals(False)


def import_superquad_data(filename):
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
