PYUIC=pyuic5
SRCDIR=./src
UIDIR=./forms
TESTDIR=tests
PY=python3

all: deps gui

deps:
	pip3 install -r requirements.txt

doc:
	cd doc
	make html

# UIFILES=$(wildcard $(UIDIR)/*.ui)
UIFILES=$(notdir $(wildcard $(UIDIR)/*.ui))
COMPILEDUI = $(UIFILES:%.ui=$(SRCDIR)/ui_%.py)
#COMPILEDUI = $(notdir $(UIFILES:%.ui=%.py))

gui : $(COMPILEDUI)

#gui : $(SRCDIR)/ui_*.py
 
$(SRCDIR)/ui_%.py : $(UIDIR)/%.ui
	$(PYUIC) $< -o $@

debug:
	@echo $(UIFILES)
	@echo $(COMPILEDUI)

# gui: $(SRCDIR)/ui_%.py
# 
# $(SRCDIR)/ui_%.py: $(UIDIR)/%.ui
# 	$(PYUIC) $< -o $@

#	$(PYUIC) $(UIDIR)/about.ui                > $(SRCDIR)/ui_about.py
#	$(PYUIC) $(UIDIR)/docpros.ui              > $(SRCDIR)/ui_docpros.py
#	$(PYUIC) $(UIDIR)/calords.ui              > $(SRCDIR)/ui_calords.py
#	$(PYUIC) $(UIDIR)/emfds.ui                > $(SRCDIR)/ui_emfds.py
#	$(PYUIC) $(UIDIR)/form.ui                 > $(SRCDIR)/ui_form.py
#	$(PYUIC) $(UIDIR)/linepropertiesdialog.ui > $(SRCDIR)/ui_linepropertiesdialog.py
#	$(PYUIC) $(UIDIR)/model.ui                > $(SRCDIR)/ui_model.py
#	$(PYUIC) $(UIDIR)/modelextra.ui           > $(SRCDIR)/ui_modelextra.py
#	$(PYUIC) $(UIDIR)/nmrds.ui                > $(SRCDIR)/ui_nmrds.py
#	$(PYUIC) $(UIDIR)/optionsdialog.ui        > $(SRCDIR)/ui_optionsdialog.py
#	$(PYUIC) $(UIDIR)/output.ui               > $(SRCDIR)/ui_output.py
#	$(PYUIC) $(UIDIR)/qsi.ui                  > $(SRCDIR)/ui_qsi.py
#	$(PYUIC) $(UIDIR)/sdistri.ui              > $(SRCDIR)/ui_sdistri.py
#	$(PYUIC) $(UIDIR)/specds.ui               > $(SRCDIR)/ui_specds.py
#	$(PYUIC) $(UIDIR)/titrsimu.ui             > $(SRCDIR)/ui_titrsimu.py
#	$(PYUIC) $(UIDIR)/ui_mainwindow.ui        > $(SRCDIR)/ui_mainwindow.py
#	$(PYUIC) $(UIDIR)/csolver.ui              > $(SRCDIR)/ui_csolver.py

tests:
	env --chdir="$(TESTDIR)" $(PY) ./alltests.py

.PHONY: deps tests
