#!/usr/bin/python3
# -*- encoding: utf-8 -*-

# FILE: manual-tk.py

import sys
import tkinter as tk

import numpy

import matplotlib


import libeq
import libplot
import libaux
import ui_manualfit as mainui

class MainWindow(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

        # default values. Modify this in the future. TODO
        self.nreagents = 2      
        self.nbetas = 1

        #self.__simu_checkboxes.setExclusive(False)
        #self.__fill_checkboxes()

        #self.make_connections()

        #l = QtWidgets.QVBoxLayout(self.ui.widget_plot)
        #self.canvas = MyCanvas(self.ui.widget_plot, width=5, height=4, dpi=100)
        #l.addStretch(0)
        #l.addWidget(self.canvas)
        #l.addStretch(0)

        self.data = None
        self.__import_Superquad()

    def createWidgets(self):
        pass

    def make_connections(self):
        #connect = QtCore.QObject.connect
        #SIGNAL = QtCore.SIGNAL
        ui = self.ui

        ui.actionNew.triggered.connect(self.__menu_new)
        ui.actionSave.triggered.connect(self.__menu_save_file)

        ui.pb_addbeta.clicked.connect(self.add_beta)
        ui.pb_removebeta.clicked.connect(self.remove_beta)
        ui.pb_addreag.clicked.connect(self.add_reagent)
        ui.pb_removereag.clicked.connect(self.remove_reagent)
        ui.pb_go.clicked.connect(self.go)

        #ui.table_species.cellChanged.connect(self.__table_species_cell_changed)
        ui.actionOpen.triggered.connect(self.__openfile)


    def add_beta(self):
        """This function handles the fact that a beta is added. This implies
        Adding one row to table_model."""

        tab = self.ui.table_model 
        cr = tab.currentRow()
        if cr == -1:
            cr = tab.rowCount() - 1
        tab.insertRow(cr+1)
        tab.setVerticalHeaderLabels(
            QtCore.QStringList(
                map(str, range(1, tab.rowCount()+1))
            )
        )

    def collect_data(self):
        '''This function checks data inserted in tables and makes sure
        they are correct before submitting and the return True. If some 
        data are not correct, this function returns False before 
        highlighting what was wrong.'''
    
        def check(tab, row, col, cell_type):
            n = tab.item(row, col)
            _error = "not valid"
            _go = True
            if n is None:
                tab.setItem(row, col, QtWidgets.QTableWidgetItem(_error))
                _go = False
                newc = None
            else:
                try:
                    newc = cell_type(n.text())
                except ValueError:
                    tab.setItem(row, col, QtWidgets.QTableWidgetItem(_error))
                    newc = None
            return _go, newc

        def check_rows(tab, col, cell_type):
            _go = True
            _rlist = list()
            for r in range(tab.rowCount()):
                _go2, it = check(tab, r, col, cell_type)
                _go = _go and _go2
                _rlist.append(it)
    
            return _go, _rlist

        def check_cols(tab, row, cell_type):
            _go = True
            _rlist = list()

            for c in range(tab.columnCount()):
                _go2, it = check(tab, row, c, cell_type)
                _go = _go and _go2
                _rlist.append(it)
    
            return _go, _rlist

        tab1 = self.ui.table_model 
        c1 = tab1.currentColumn()

        if self.current_mode():
            tab2 = self.ui.table_species
        else:
            tab2 = self.ui.table_titration

        labels = [ str(tab2.item(r, 0).text()) for r in range(tab2.rowCount()) ]

        go = True

        # Check betas are float numbers
        beta_column = 0
        #go2, lB = check_rows(tab1, beta_column, float)
        #go = go and go2

        _re_beta = re.compile(r'^([+-]?\d+\.\d+)(\((\d)\))?$')
        lB = list()
        elB = list()
        for r in range(tab1.rowCount()):
            n = str(tab1.item(r, beta_column).text())
            m = _re_beta.match(n)
            if not m:
                self.message("Error in beta #%d" % r)
                #print n
                go = False
                break
            else:
                b = m.group(1)
                lB.append(float(b))
                eb = m.group(3)
                if eb is None:
                    elB.append(0.0)
                else:
                    nd = len(b)-b.find('.')-1   # number of decimal positions
                    elB.append(float(eb)*10**(-nd))

        # Check stoich. are integers
        lP = list()
        for c in range(1,tab1.columnCount()):
            go2, lP2 = check_rows(tab1, c, int)
            go = go and go2
            lP.append(lP2)
            
        # Check bounduary conditions are float
        if self.current_mode():
            tab2 = self.ui.table_species
            c_initial = 1
            c_final = 2
            c_pX = 3
            pX = list()

            go2, Tinit = check_rows(tab2, c_initial, float)
            go = go and go2

            go2, Tfinal = check_rows(tab2, c_final, float)
            go = go and go2
    
            #TODO checkboxes
            pX = [ False ] * tab2.rowCount()
            for b in self.__simu_checkboxes.buttons():
                n = self.__simu_checkboxes.id(b)
                #print "n = ", n, " of ", tab1.rowCount()
                pX[n] = ( b.checkState() == QtCore.Qt.Checked )
            
            ret3 = [Tinit, Tfinal, pX ]
        else:
            tab2 = self.ui.table_titration
            c_initial = 1
            c_buret = 2

            go2, Tinit = check_rows(tab2, c_initial, float)
            go = go and go2

            go2, Tburet = check_rows(tab2, c_buret, float)
            go = go and go2

            V = self.ui.dsb_V0.value(), self.ui.dsb_Vf.value()

            ret3 = [ Tinit, Tburet, V ]

        #TODO
        # Check that in cb_Yaxis and cb_Xaxis are not the same variables
        x = self.ui.cb_Xaxis.currentIndex() / 2
        y = self.ui.cb_Yaxis.currentIndex() - 1
        #print "x = ", x, \
        #    "(",self.ui.cb_Xaxis.currentIndex() ,"), y = ", \
        #    y, "(", self.ui.cb_Yaxis.currentIndex(),")"
        if x == y:
            self.message('X and Y for plot are the same variable')
            go = False

        if go:
            self.message('Data are consistent')

        return go, lB, elB, lP, labels, ret3

    def get_B(self):
        tab1 = self.ui.table_model 
        B = list()
        for r in range(tab1.rowCount()):
            it = tab1.item(r,0)
            assert it is not None
            B.append(float(it.text()))

        return numpy.array(B)

    def get_P(self):
        tab1 = self.ui.table_model 
        P = list()
        for r in range(tab1.rowCount()):
            rt = list()
            for c in range(1, tab1.columnCount()):
                it = tab1.item(r,c)
                assert it is not None
                rt.append(int(it.text()))
            P.append(rt)

        return numpy.array(P)

    def get_labels(self):
        tab2 = self.ui.table_species
        label_col = 0
        labels = []
        for r in range(tab2.rowCount()):
            labels.append(str(tab2.item(r, label_col).text()))
    
        return labels

    def remove_beta(self):
        tab = self.ui.table_model 
        if tab.rowCount() < 2:
            return
        cr = tab.currentRow()
        if cr == -1:
            self.statusBar().showMessage("Select a row")
        else:
            tab.removeRow(cr)

        tab.setVerticalHeaderLabels(
            QtCore.QStringList(
                map(str, range(1, tab.rowCount()+1))
            )
        )

    def add_reagent(self):
        tab1 = self.ui.table_model 
        c1 = tab1.currentColumn()

        if self.current_mode():
            tab2 = self.ui.table_species
            tab3 = self.ui.table_titration
        else:
            tab2 = self.ui.table_titration
            tab3 = self.ui.table_species

        c2 = tab2.currentRow()

        if c2 == -1:
            if c1 != -1:
                new_row = c1+1
                #tab1.insertColumn(c1+2)
                #tab2.insertRow(c1+1)
                #tab3.insertRow(c1+1)
            else:
                rc = tab2.rowCount()
                #print rc, tab1.columnCount() 
                assert rc == tab3.rowCount() == tab1.columnCount() - 1
                new_row = rc + 1
                #tab1.insertColumn(rc+1)
                #tab2.insertRow(rc)
                #tab3.insertRow(rc)
        else:
            new_row= c2 + 1 
            #tab1.insertColumn(c2+2)
            #tab2.insertRow(c2+1)
            #tab3.insertRow(c2+1)

        tab1.insertColumn(new_row + 1)
        tab2.insertRow(new_row)
        tab3.insertRow(new_row)
        #self.__pX.insert(new_row, False)
        self.__fill_checkboxes()

    def remove_reagent(self):
        tab1 = self.ui.table_model 
        c1 = tab1.currentColumn()

        if self.current_mode():
            tab2 = self.ui.table_species
            tab3 = self.ui.table_titration
        else:
            tab2 = self.ui.table_titration
            tab3 = self.ui.table_species
            c2 = tab2.currentRow()

        if tab1.rowCount() <= 3:
            return
        cc = tab1.currentColumn()
        if cc == -1:
            self.message("Select a column")
        elif cc == 0:
            self.message("Can't delete this column")
        else:
            tab1.removeColumn(cc)
    
    def current_mode(self):
        '''Return 1 if mode is species distribution and 0 for titration
        simulation'''

        if self.ui.rb_distri.isChecked():
            return 1
        elif self.ui.rb_simu.isChecked():
            return 0
        else:
            #This should never happen. But just in case check rb_distri and come back
            self.ui.rb_distri.click()
            return 0

    def list_of_reagents(self):
        'returns a list with the labels of the reagents'
        # find out type of graphic 
        if self.current_mode():
            table = self.ui.table_species
        else:
            table = self.ui.table_titration
            
        # get and return column values
        column_labels = 0
        return [ table.item(row, column_labels).text() \
            for row in range(table.rowCount()) ]

    def message(self, text):
        self.statusBar().showMessage(text)

    def go(self):
        'Start calculations'
        
        # Default number of points. Modify this in the future
        N = 50
        go, lB, elB, lP, labels, ret3 = self.collect_data()
        if not go: 
            print("bad data")
            return

        self.canvas.gca().clear()

        assert isinstance(labels, list)

        logB = numpy.array(lB)
        B = 10**logB
        P = numpy.array(lP).T

        E, S = P.shape

        x = floor(self.ui.cb_Xaxis.currentIndex() / 2)

        #print logB, P
        if self.current_mode():     # species distribution
            # TODO
            Tinit, Tfinal, pX = ret3
            E = len(B)
            T0 = [[a, b] for a, b in zip(Tinit, Tfinal)]
            X, C = libeq.simulation(B.reshape((E,1)), P, T0, pX, N=N, x=x)
            if pX[x]:
                new_x = -numpy.log10(X)
            else:
                new_x = X

            T = libaux.build_T(T0, pX, N)
            if self.ui.cb_Yaxis.currentText() == 'conc.':
                self.canvas.plot_conc(new_x, C, libeq.make_labels(labels, P), \
                    self.ui.cb_Xaxis.currentText())
            if self.ui.cb_Yaxis.currentText()[0] == '%':
                # find reference species
                ref_label = self.ui.cb_Yaxis.currentText()[1:]
                ref = labels.index(ref_label)
                rP = P[:,ref].tolist()
                
                # find species to plot (all those that contain reference species)
                to_plot = [ n+S for n, i in enumerate(rP) if i != 0 ]
                to_plot.insert(0, ref)
                to_plot.sort()
                q_plot = [ i for i in rP if i != 0 ]
                q_plot.insert(0, 1)
                
                # remove species not to plot
                new_C = C[:, to_plot]
                
                # divide by T
                T_ref = T[:, ref]
                T_ref_ext = T_ref.reshape((N,1)) * numpy.ones((1, len(to_plot)))
                new_C = new_C / T_ref_ext
                
                # divide by stoichiometric coef.
                new_C = new_C / (numpy.ones((N,1)) * numpy.array(q_plot))
                
                # multiply by 100
                new_C = 100* new_C
                
                # arrange labels
                full_labels = libplot.make_labels(labels, P) 
                plot_labels = [ full_labels[i] for i in to_plot ]
                
                self.canvas.plot_relative(new_x, new_C, 
                    plot_labels, 
                    self.ui.cb_Xaxis.currentText(), 
                    '% Formation Relative to ' + ref_label)

        else:                       # titration simulation
            # TODO
            Tinit, Tburet, V = ret3

    def __assertrb(self):
        pass

    def __fill_checkboxes(self):
        column_checkboxes = 3
        t = self.ui.table_species
        #assert len(self.__pX) == t.rowCount()
        for r in range(t.rowCount()):
            if t.item(r, column_checkboxes) is None:
                c = QtWidgets.QCheckBox("no", t)
                self.__simu_checkboxes.addButton(c, r)
                t.setCellWidget(r, column_checkboxes, c)

    def __checkboxchanged(self, n):
        cb = self.__simu_checkboxes.button(n)
        if cb.checkState() == QtCore.Qt.Checked:
            cb.setText("yes")
        if cb.checkState() == QtCore.Qt.Unchecked:
            cb.setText("no")

    def __menu_new(self):
        '''Clear all the forms.'''

        # TODO if modified, ask for saving file

        ui = self.ui

        QtCore.QObject.disconnect(ui.table_species, QtCore.SIGNAL("cellChanged(int, int)"), 
            self.__table_species_cell_changed)

        # clear tables
        ui.table_species.clearContents()
        ui.table_model.clearContents()
        ui.table_titration.clearContents()

        ui.table_titration.setRowCount(self.nreagents)
        ui.table_species.setRowCount(self.nreagents)
        ui.table_model.setRowCount(self.nbetas)

        ui.dsb_V0.setValue(0.0)
        ui.dsb_Vf.setValue(0.0)

        ui.cb_Xaxis.clear()
        ui.cb_Yaxis.clear()
        self.__fill_checkboxes()

        QtCore.QObject.connect(ui.table_species, QtCore.SIGNAL("cellChanged(int, int)"), 
            self.__table_species_cell_changed)

    def __menu_save_file(self):
        self.__savefile()

    def __openfile(self):
        import xml.etree.ElementTree as ET
        
        default_dir = '.'
        filters = "XML Files (*.xml);;All Files (*.*)"
        #filename = QtWidgets.QFileDialog.getOpenFileName(
        #    parent = self, 
        #    caption = 'Choose file to open', 
        #    directory = default_dir, 
        #    filter = filters) 
        filename = 'data.xml'
        
        tree = ET.parse(filename)
        root = tree.getroot()
        version = root.attrib['version']
        version_major, version_minor = map(int, version.split('.'))
        model = root.find('model')
        title = root.find('title')
        distri = root.find('distri')
        simu = root.find('simu')

        self.canvas.set_title(title.text)

        labels = model.find('labels').text.split()
        Blist = list()
        eBlist = list()
        Plist = list()
        for line in model.findall('line'):
            b = line.find('beta')
            eb = line.find('error')
            Blist.append(b.text)
            eBlist.append(b.text)       # might return None. This is normal but must be handled.

            p = line.find('p').text.split()
            Plist.append([ i for i in p ])

        distri_Tinitial = [ i for i in distri.find('initial').text.split() ]
        distri_Tfinal   = [ i for i in distri.find('final').text.split() ]
        distri_pX = [ i == 'True' for i in distri.find('pX').text.split() ]

        simu_Tinitial = simu.find('initial').text.split() 
        simu_Tburet   = simu.find('buret').text.split() 
        simu_V = [ float(i) for i in simu.find('V').text.split() ]

        #TODO: Verify data

        # Put data on the table
        S = len(labels)
        E = len(Blist)

        tab1 = self.ui.table_model 
        tab2 = self.ui.table_species
        tab3 = self.ui.table_titration

        QtCore.QObject.disconnect(self.ui.table_species, QtCore.SIGNAL("cellChanged(int, int)"), self.__table_species_cell_changed)

        #def add_rows(t, n):
        #    for i in range(t.rowCount(), n):
        #        t.insertRow(i)

        #def remove_rows(t, n):
        #    for i in range(t.rowCount(), n, -1):
        #        t.removeRow(i)

        tab1.clearContents()            #tab1 is E x (S+1)
        tab1.setRowCount(E)
        tab1.setColumnCount(S+1)
        tab1.setHorizontalHeaderLabels( [ "log B" ] + labels )

        tab2.clearContents()            #tab2 is S x 4
        tab2.setRowCount(S)

        tab3.clearContents()            #tab3 is S x 3
        tab3.setRowCount(S)

        #TODO remove this 'if' in final version
        if version_major >= 0 and  version_minor >= 1:
            #TODO read <styles>
            #TODO read <plots>
            pass

        def str2twi(txt):
            assert isinstance(txt, str)
            return QtWidgets.QTableWidgetItem(txt)

        for s in range(S):
            tab2.setItem(s, 0, str2twi(labels[s]))
            tab2.setItem(s, 1, str2twi(distri_Tinitial[s]))
            tab2.setItem(s, 2, str2twi(distri_Tfinal[s]))
            tab3.setItem(s, 0, str2twi(labels[s]))
            tab3.setItem(s, 1, str2twi(simu_Tinitial[s]))
            tab3.setItem(s, 2, str2twi(simu_Tburet[s]))
            for e in range(E):
                tab1.setItem(e,s+1, str2twi(Plist[e][s]))

        for e in range(E):
            #TODO rebuild new label from Blist and eBlist
            currB = Blist[e]
            curreB = eBlist[e]
            #TODO 1. decompose curreB in number of decimal places and significant digit
            #TODO 2. reformat currB 


            tab1.setItem(e,0, str2twi(Blist[e]))
            
        self.ui.dsb_V0.setValue(simu_V[0])
        self.ui.dsb_Vf.setValue(simu_V[1])

        QtCore.QObject.connect(self.ui.table_species, QtCore.SIGNAL("cellChanged(int, int)"), self.__table_species_cell_changed)
    
        self.__fill_checkboxes()
        column_checkboxes = 3
        for r in range(tab2.rowCount()):
            it = tab2.cellWidget(r, column_checkboxes)
            assert isinstance(it, QtWidgets.QCheckBox)
            it.setChecked(distri_pX[r])
            self.__checkboxchanged(r)
    
        self.__populate_cb_Xaxis()
        self.__populate_cb_Yaxis()

    def __import_Superquad(self):
        from importsuperquad import importsuperquad
        default_dir = '.'
        filters = "Superquad Files (*.sup);;All Files (*.*)"
        #filename = QtWidgets.QFileDialog.getOpenFileName(
        #    parent = self, 
        #    caption = 'Choose file to open', 
        #    directory = default_dir, 
        #    filter = filters) 
        filename = ("/home/salvador/proyectos/eslib/hpal.sup", 0)

        #TODO add try block here
        self.data = importsuperquad(filename[0])
        print(self.data)
        self.__fill_all_data() 

        data = { 'V': self.data['V'], 'emf_read': self.data['emf']}
        self.canvas.setDataGraphs(data)
        self.canvas.draw()

    def __fill_all_data(self):
        #TODO fill table_model        
        # B P labels flags
        E, S = self.data["P"].shape
        columns = S + 2
        rows = E
        t = self.ui.table_model
        t.setRowCount(rows)
        t.setColumnCount(columns)

        # column 0 -> QDoubleSpinBox AND column E+1 -> qcombobox
        self.__tabmodelwidgets_beta = []
        self.__tabmodelwidgets_flags = []
        for r in range(rows):
            self.__tabmodelwidgets_beta.append(QtWidgets.QDoubleSpinBox())
            self.__tabmodelwidgets_beta[-1].setMinimum(-1000.0)
            self.__tabmodelwidgets_beta[-1].setMaximum(1000.0)
            self.__tabmodelwidgets_beta[-1].setValue(self.data["logB"][r])
            t.setCellWidget(r, 0, self.__tabmodelwidgets_beta[-1])

            self.__tabmodelwidgets_flags.append(QtWidgets.QComboBox())
            self.__tabmodelwidgets_flags[-1].addItems(["refine","constant","ignore"])
            t.setCellWidget(r, S+1, self.__tabmodelwidgets_flags[-1])

        # column 1 to E+1 -> int
        t.setHorizontalHeaderLabels([ 'Value'] + self.data['labels'] + ['Flag'])
        for c in range(1,columns-1):
            # TODO add labels to column headers
            for r in range(rows):
                t.setItem(r, c, QtWidgets.QTableWidgetItem(str(self.data["P"][r,c-1])))

        # populate cbdataset
        NDS = len(self.data["emf"])
        self.__visible_dataset = 0
        self.ui.cbdataset.addItems([ "dataset #%s" % (i+1) for i in range(NDS)])

        # fill table_titration
        # labels, T, buret 
        self.ui.table_titration.setRowCount(S)
        for r in range(S):
            self.ui.table_titration.setItem(r, 0, QtWidgets.QTableWidgetItem(self.data["labels"][r]))
            self.ui.table_titration.setItem(r, 1, QtWidgets.QTableWidgetItem(str(self.data["T"][self.__visible_dataset][r])))
            self.ui.table_titration.setItem(r, 2, QtWidgets.QTableWidgetItem(str(self.data["buret"][self.__visible_dataset][r])))
        
        #TODO fill table_data
        # V, emf
        N = len(self.data['emf'][self.__visible_dataset])
        self.ui.table_data.setRowCount(N)
        for r in range(N):
            self.ui.table_data.setItem(r, 0, QtWidgets.QTableWidgetItem(str(self.data["V"][self.__visible_dataset][r])))
            self.ui.table_data.setItem(r, 1, QtWidgets.QTableWidgetItem(str(self.data["emf"][self.__visible_dataset][r])))

        #TODO update plot

    def __savefile(self):
        import xml.etree.ElementTree as ET

        default_dir = '.'
        filters = "XML Files (*.xml);;All Files (*.*)"
        #filename = QtWidgets.QFileDialog.getOpenFileName(
        #    parent = self, 
        #    caption = 'Choose file to open', 
        #    directory = default_dir, 
        #    filter = filters) 
        #filename = 'data.xml'

        data_version = '0.0'
        xml_version = '1.0'
        encoding = 'UTF-8'

        #TODO modify title
        title = 'untitled data'

        text = '''<?xml version="%s" encoding="%s" ?>
            <simdata version="%s">\n''' % (xml_version, encoding, data_version)
    
        text += ' <title>' + title + '</title>\n <model>\n'
        text += '<labels>' + " ".join(self.get_labels()) + '</labels>\n'
    
        B = self.get_B()
        P = self.get_P()
        assert isinstance(P, numpy.ndarray)
        E, S = P.shape
    
        for n in range(E):
            text += '<line> <beta>'
            text += str(B[n])
            text += '</beta> <p>' 
            text += " ".join(str(P[n]))
            text += '</p> <key>0</key> </line>'
    
        text += ' </model>\n <distri active="no">\n'
    
        def get_distri_pars():
            tab = self.ui.table_species
            dp = []
            dp.append([ str(tab.item(r, 1).text()) for r in range(tab.rowCount()) ])
            dp.append([ str(tab.item(r, 2).text()) for r in range(tab.rowCount()) ])
            dp.append([ str(tab.cellWidget(r, 3).isChecked()) for r in range(tab.rowCount()) ])
            return dp    

        dp = get_distri_pars()
        text += '  <initial units="mol/L">' + " ".join(dp[0]) + '</initial>\n'
        text += '  <final units="mol/L">' + " ".join(dp[1]) + '</final>\n'
        text += '<pX>' + " ".join(dp[2]) + '</pX>\n</distri>'

        def get_simu_pars():
            tab = self.ui.table_titration
            sp = []
            sp.append([ str(tab.item(r, 1).text()) for r in range(tab.rowCount()) ])
            sp.append([ str(tab.item(r, 2).text()) for r in range(tab.rowCount()) ])
            sp.append([ str(self.ui.dsb_V0.value()), str(self.ui.dsb_Vf.value()) ])
            return sp

        sp = get_simu_pars()
        text += ' <simu active="yes">\n'
        text += '  <initial units="mmol">' + " ".join(sp[0]) + '</initial>\n'
        text += '  <buret units="mol/L">' + " ".join(sp[1]) + '</buret>\n'
        text += '  <V units="mL">' + " ".join(sp[2]) + '</V>\n </simu>\n</simdata>'

        #print text
        #TODO save in file


    def __table_species_cell_changed(self, row, column):
        #print "changed %d, %d" % (row, column)
        if column == 0:
            new_label = self.ui.table_species.item(row, column)
            self.ui.table_model.setHorizontalHeaderItem(row+1, new_label)

        #print "new content is ", self.ui.table_species.item(row, column).text()

        self.__populate_cb_Xaxis()
        self.__populate_cb_Yaxis()

    def __ignorelowerthan(self, v):
        self.canvas.ignorelowerthan = v

    def __unlabellowerthan(self, v):
        self.canvas.unlabellowerthan = v 

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

class MyCanvas(FigureCanvas):
    '''This class contains the canvas where all the graphs are displayed. Also
    handles events and such.'''

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        # Default parameters
        self.__pickradius = 2
        self.ignorelowerthan = 4
        self.unlabellowerthan = 9
        self.title = None

        fig = Figure(figsize=(width, height), dpi=dpi)
        #self.axes = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        #self.axes.hold(False)

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

        self.data = None
        self.fit_data = None

        #self.plot_objects = list()
        #self.title_object = None

        #self.__setup_popupmenu()
        #cid = self.mpl_connect('button_press_event', self.on_click)


    def setDataGraphs(self, data):
        self.data = data

    def setFitGraphs(self, data):
        self.fit_data = data

    def draw(self):
        ax = self.figure.add_subplot(111)
        ax.plot(numpy.arange(0, 3.2, 0.1), numpy.sin(numpy.arange(0, 3.2, 0.1)),'-')

        #assert self.data is not None
        #if self.fit_data is None:
        #    spln1 = 2
        #else:
        #    spln1 = 1

        #spln2 = len(self.data) 

        #self.fig.clear()
        #assert 'V' in self.data.keys()
        #assert 'emf_read' in self.data.keys()
        #axes = []
        #for n, (v, e) in enumerate(zip(self.data['V'],self.data['emf_read'])):
        #    axes.append(self.fig.add_subplot(spln1, spln2, n))
        #    axes[-1].scatter(v, e) 

    def __setup_popupmenu(self):
        self.menu = QtWidgets.QMenu(self)

        self.__action_set_title = QtWidgets.QAction("set title", self)
        self.menu.addAction(self.__action_set_title)

    def gca(self):
        return self.axes

    def on_click(self, event):
        #print "Mouse clicked on (%d, %d)" % (event.x, event.y)

        # Check click on curves
        for obj in self.plot_objects:
            if isinstance(obj, Line2D):
                c, a = obj.contains(event)
                if c: 
                    #print "clicked on object ", obj
                    cl = obj.get_color()
                    ls = obj.get_linestyle()
                    lw = obj.get_linewidth()

                    d = LinePropertiesDialog(self)
                    d.set_initial(cl, ls, lw)
                    d.exec_()
                    if d.result() == QtWidgets.QDialog.Accepted:
                        ans = d.get_result()
                        print(ans)
                        obj.set_linestyle(ans['linestyle'])
                        obj.set_color(ans['color'])
                        obj.set_linewidth(ans['linewidth'])
                        self.draw()
                    return

        # TODO check click on titles

        if self.title_object is not None:
            c, d =  self.title_object.contains(event)
            if c:
                print("clicked on title: ", d)
                #TODO open dialog

        xaxis_title = self.axes.get_xaxis().get_label()
        c, d = xaxis_title.contains(event)
        if c:
            print("clicked on xaxes title")
            #TODO open dialog

        yaxis_title = self.axes.get_yaxis().get_label()
        c, d = yaxis_title.contains(event)
        if c:
            print("clicked on yaxes title")
            #TODO open dialog

        # TODO check click on axis

        # TODO Correct coordinates
        point = QtCore.QPoint(event.x, event.y)

        # TODO If right click, open popup menu
        if event.button == 2:
            self.menu.popup( self.mapToGlobal(point) )

        # TODO If hold click, drag element if dragable

        # TODO If double click, open element properties if aplicable
        if event.button == 0 and event.dblclick:
            pass

    def place_labels(self, label_text, label_xy):
        assert len(label_text) == len(label_xy)
        label_objs = list()
        self.drlbls = list()
        for xy, l in zip(label_xy, label_text):
            assert isinstance(xy, tuple)
            assert len(xy) == 2
            assert isinstance(xy[0], float)
            assert isinstance(xy[1], float)
            label_objs.append(
                self.axes.text(
                    xy[0], xy[1], l,
                    verticalalignment='bottom',
                    horizontalalignment='center',
                    color='k',
                    family="serif"
                )
            )
            drlbl = DragableLabel(label_objs[-1])
            drlbl.connect()
            self.drlbls.append(drlbl)

    def plot(self, X, Y, labels, xlabel=None, ylabel=None):
        self.axes.clear()

        if xlabel is not None:
            self.axes.set_xlabel(xlabel)
        if ylabel is not None:
            self.axes.set_ylabel(ylabel)

        new_plot_objects = self.axes.plot(X, Y, pickradius=self.__pickradius)
        self.plot_objects.extend(new_plot_objects)

        window = (self.axes.get_xlim(), self.axes.get_ylim())
        label_xy = libplot.position_labels(X, Y, window)
        self.place_labels(labels, label_xy)

        if self.title is not None:
            self.title_object = self.axes.set_title(self.title)

        self.draw()

    def plot_conc(self, X, C, labels, xlabel):
        self.plot(X, C, labels, xlabel, 'Concentration / M')

    def plot_relative(self, X, rC, labels, xlabel, ylabel):
        to_plot = []
        for i, s in enumerate(rC.T):
            assert s.shape == X.shape
            if max(s) >= self.ignorelowerthan:
                to_plot.append(i)

        Y = rC.T[to_plot].T
        labels_plot = [ l for i, l in enumerate(labels) if i in to_plot ]
        self.plot(X, Y, labels_plot, xlabel, ylabel)

    def set_title(self, title):
        assert isinstance(title, str)
        self.title = title

    def __line_clicked(self, event):
        print("line clicked!!")

    def __prepare_menus(self):
        m = QtWidgets.QMenu(parent = self)
        action_aed = QtWidgets.QAction("add external data")
        QtCore.QObject.connect(action_aed, QtCore.SIGNAL("triggered()"), self.add_external_data)
        m.addAction(action_aed)

from ui_linepropertiesdialog import Ui_LinePropertiesDialog
class LinePropertiesDialog(QtWidgets.QDialog):
    def __init__(self, parent):
        QtWidgets.QDialog.__init__(self, parent=parent)
        self.setModal(True)
        self.ui = Ui_LinePropertiesDialog()
        self.ui.setupUi(self)
        self.show()
        #print "Hello, dialog!"

        self.__linetypes = {
            'solid line style'      : '-',
            'dashed line style'     : '--',
            'dash-dot line style'   : '-.',
            'dotted line style'     : ':',
            'point marker'          : '.',
            'pixel marker'          : ',',
            'circle marker'         : 'o',
            'triangle_down marker'  : 'v',
            'triangle_up marker'    : '^',
            'triangle_left marker'  : '<',
            'triangle_right marker' : '>',
            'tri_down marker'       : '1',
            'tri_up marker'         : '2',
            'tri_left marker'       : '3',
            'tri_right marker'      : '4',
            'square marker'         : 's',
            'pentagon marker'       : 'p',
            'star marker'           : '*',
            'hexagon1 marker'       : 'h',
            'hexagon2 marker'       : 'H',
            'plus marker'           : '+',
            'x marker'              : 'x',
            'diamond marker'        : 'D',
            'thin_diamond marker'   : 'd',
            'vline marker'          : '|',
            'hline marker'          : '_' }
        self.__reverselinetypes = dict(zip(self.__linetypes.values(), self.__linetypes.keys()))

        self.ui.cb_linetype.addItems([str(x) for x in self.__linetypes.keys()])

        self.__colors = {
            'blue'   : 'b',	
            'green'  : 'g',	
            'red'    : 'r',	
            'cyan'   : 'c',	
            'magenta': 'm',	
            'yellow' : 'y',	
            'black'  : 'k',	
            'white'  : 'w'	}
        self.__reversecolors = dict(zip(self.__colors.values(), self.__colors.keys()))

        self.ui.cb_color.addItems([ str(x) for x in self.__colors.keys()])

    def set_initial(self, color, linestyle, linewidth):
        assert color in self.__colors.values()
        assert linestyle in self.__linetypes.values()
        assert isinstance(linewidth, float)

        self.ui.sb_thickness.setValue(linewidth)

        self.ui.cb_color.setCurrentIndex(
            self.ui.cb_color.findText(
                self.__reversecolors[color]))

        self.ui.cb_linetype.setCurrentIndex(
            self.ui.cb_linetype.findText(
                self.__reverselinetypes[linestyle]))

    def get_result(self):
        color = self.__colors[str(self.ui.cb_color.currentText())]
        ls = self.__linetypes[str(self.ui.cb_linetype.currentText())]
        lw = self.ui.sb_thickness.value()

        return { 'color': color, 'linestyle': ls, 'linewidth': lw }


# TODO Duplicated code in distri_drag.py
class DragableLabel:
    def __init__(self, label):
        self.label = label
        self.press = None

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.label.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.label.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.label.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.label.axes: return

        contains, attrd = self.label.contains(event)
        if not contains: return

        x0, y0 = self.label.get_position()
        self.press = x0, y0, event.xdata, event.ydata

    def on_motion(self, event):
        'on motion we will move the rect if the mouse is over us'
        if self.press is None: return
        if event.inaxes != self.label.axes: return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        #print 'x0=%f, xpress=%f, event.xdata=%f, dx=%f, x0+dx=%f'%(x0, xpress, event.xdata, dx, x0+dx)
        self.label.set_x(x0+dx)
        self.label.set_y(y0+dy)

        self.label.figure.canvas.draw()

    def on_release(self, event):
        'on release we reset the press data'
        self.press = None
        self.label.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.label.figure.canvas.mpl_disconnect(self.cidpress)
        self.label.figure.canvas.mpl_disconnect(self.cidrelease)
        self.label.figure.canvas.mpl_disconnect(self.cidmotion)


if __name__ == '__main__':
    root = tk.Tk()
    app = MainWindow(master=root)
    app.mainloop()
    sys.exit(0)


