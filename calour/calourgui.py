#!/usr/bin/env python

# Calour GUI - a full GUI wrapping for calour functions

# ----------------------------------------------------------------------------
# Copyright (c) 2016--,  Calour development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import os
from logging import getLogger
from pkg_resources import resource_filename

from PyQt5 import QtWidgets, QtCore, uic
from PyQt5.QtWidgets import (QHBoxLayout, QVBoxLayout,
                             QWidget, QPushButton, QLabel,
                             QComboBox, QLineEdit, QCheckBox, QSpinBox,
                             QDialog, QDialogButtonBox)

import calour as ca
import calour.cahelper as cah
import calour.analysis

logger = getLogger(__name__)


class AppWindow(QtWidgets.QMainWindow):
    # the experiments loaded for analysis
    _explist = {}

    def __init__(self):
        super().__init__()
        # load the gui
        uic.loadUi(get_ui_file_name('CalourGUI.ui'), self)

        # handle button clicks
        self.wLoad.clicked.connect(self.load)
        self.wPlot.clicked.connect(self.plot)

        # the experiment list right mouse menu
        self.wExperiments.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.wExperiments.customContextMenuRequested.connect(self.listItemRightClicked)

        # add functions
        # init the action group list
        action_groups = ['sample', 'feature', 'analysis']
        self.actions = {}
        for caction in action_groups:
            self.actions[caction] = {}

        self.add_action_button('sample', 'barvaz', self.test)
        self.add_action_button('sample', 'Sort', self.sample_sort_samples)
        self.add_action_button('sample', 'Filter', self.sample_filter_samples)
        self.add_action_button('sample', 'Cluster', self.sample_cluster_samples)
        self.add_action_button('sample', 'Join fields', self.sample_join_fields)
        self.add_action_button('sample', 'Filter orig. reads', self.sample_filter_orig_reads)

        self.add_action_button('feature', 'Cluster', self.feature_cluster)
        self.add_action_button('feature', 'Filter min reads', self.feature_filter_min_reads)
        self.add_action_button('feature', 'Filter taxonomy', self.feature_filter_taxonomy)
        self.add_action_button('feature', 'Filter fasta', self.feature_filter_fasta)
        self.add_action_button('feature', 'Sort Abundance', self.feature_sort_abundance)

        self.add_action_button('analysis', 'Diff. abundance', self.analysis_diff_abundance)

        # load sample dataset for debugging
        exp = ca.read_taxa('/Users/amnon/Projects/centenarians/final.withtax.biom', '/Users/amnon/Projects/centenarians/map.txt')
        exp._studyname = 'centenarians'
        self.addexp(exp)
        self.show()

    def get_exp_from_selection(self):
        '''Get the experiment from the selection in wExperiments

        Parameters
        ----------

        Returns
        -------
        expdat : Experiment
            the first selected experiment in the wExperiments list
        '''
        item = self.wExperiments.selectedItems()[0]
        cname = str(item.text())
        if cname not in self._explist:
            logger.warn('experiment not found. name=%s' % cname)
            return None
        expdat = self._explist[cname]
        return expdat

    def plot(self):
        '''
        Plot the experiment
        '''
        expdat = self.get_exp_from_selection()
        sort_field_vals = ['<none>']+list(expdat.sample_metadata.columns)
        res = dialog([{'type': 'label', 'label': 'Plot experiment %s' % expdat._studyname},
                      {'type': 'combo', 'label': 'Field', 'items': sort_field_vals},
                      {'type': 'bool', 'label': 'sort'}], expdat=expdat)
        if res is None:
            return
        if res['Field'] == '<none>':
            field = None
        else:
            field = res['Field']
        if res['sort'] and field is not None:
            logger.debug('sort')
            newexp = expdat.sort_by_metadata(field, axis=0)
        else:
            newexp = expdat
        newexp.plot(gui='qt5', sample_field=field)

    def sample_sort_samples(self):
        expdat = self.get_exp_from_selection()
        res = dialog([{'type': 'label', 'label': 'Sort Samples'},
                      {'type': 'field', 'label': 'Field'},
                      {'type': 'string', 'label': 'new name'}], expdat=expdat)
        if res is None:
            return
        if res['new name'] == '':
            res['new name'] = '%s-sort-%s' % (expdat._studyname, res['field'])
        newexp = expdat.sort_by_metadata(res['field'], axis=0)
        newexp._studyname = res['new name']
        self.addexp(newexp)

    def sample_cluster_samples(self):
        expdat = self.get_exp_from_selection()
        newexp = expdat.cluster_data(axis=1)
        newexp._studyname = newexp._studyname + '-cluster-samples'
        self.addexp(newexp)

    def sample_filter_samples(self):
        expdat = self.get_exp_from_selection()
        logger.debug('filter samples for study: %s' % expdat._studyname)
        res = dialog([{'type': 'label', 'label': 'Filter Samples'},
                      {'type': 'field', 'label': 'Field'},
                      {'type': 'value', 'label': 'value'},
                      {'type': 'bool', 'label': 'negate'},
                      {'type': 'string', 'label': 'new name'}], expdat=expdat)
        if res is None:
            return
        if res['new name'] == '':
            if res['negate']:
                res['new name'] = '%s-%s-not-%s' % (expdat._studyname, res['field'], res['value'])
            else:
                res['new name'] = '%s-%s-%s' % (expdat._studyname, res['field'], res['value'])
        newexp = expdat.filter_by_metadata(res['field'], res['value'], negate=res['negate'])
        newexp._studyname = res['new name']
        self.addexp(newexp)

    def sample_join_fields(self):
        expdat = self.get_exp_from_selection()
        res = dialog([{'type': 'label', 'label': 'Join Fields'},
                      {'type': 'combo', 'label': 'Field1', 'items': expdat.sample_metadata.columns},
                      {'type': 'combo', 'label': 'Field2', 'items': expdat.sample_metadata.columns},
                      {'type': 'string', 'label': 'new name'}], expdat=expdat)
        if res is None:
            return
        if res['new name'] == '':
            res['new name'] = '%s-join-%s-%s' % (expdat._studyname, res['Field1'], res['Field2'])
        newexp = cah.join_fields(expdat, field1=res['Field1'], field2=res['Field2'])
        newexp._studyname = res['new name']
        self.addexp(newexp)

    def sample_filter_orig_reads(self):
        expdat = self.get_exp_from_selection()
        res = dialog([{'type': 'label', 'label': 'Filter Original Reads'},
                      {'type': 'int', 'label': 'Orig Reads', 'max': 100000, 'default': 10000},
                      {'type': 'string', 'label': 'new name'}], expdat=expdat)
        if res is None:
            return
        if res['new name'] == '':
            res['new name'] = '%s-min-%d' % (expdat._studyname, res['Orig Reads'])
        newexp = cah.filter_orig_reads(expdat, minreads=res['Orig Reads'])
        newexp._studyname = res['new name']
        self.addexp(newexp)

    def feature_filter_min_reads(self):
        expdat = self.get_exp_from_selection()
        res = dialog([{'type': 'label', 'label': 'Filter minimal reads per feature'},
                      {'type': 'int', 'label': 'min reads', 'max': 50000, 'default': 10},
                      {'type': 'string', 'label': 'new name'}], expdat=expdat)
        if res is None:
            return
        if res['new name'] == '':
            res['new name'] = '%s-minreads-%d' % (expdat._studyname, res['min reads'])
        newexp = cah.filter_min_reads(expdat, minreads=res['min reads'])
        newexp._studyname = res['new name']
        self.addexp(newexp)

    def feature_filter_taxonomy(self):
        expdat = self.get_exp_from_selection()
        res = dialog([{'type': 'label', 'label': 'Filter Taxonomy'},
                      {'type': 'string', 'label': 'Taxonomy'},
                      {'type': 'bool', 'label': 'Exact'},
                      {'type': 'bool', 'label': 'Negate'},
                      {'type': 'string', 'label': 'new name'}], expdat=expdat)
        if res is None:
            return
        if res['new name'] == '':
            res['new name'] = '%s-tax-%s' % (expdat._studyname, res['Taxonomy'])
        newexp = cah.filter_taxonomy(expdat, res['Taxonomy'], negate=res['Negate'], exact=res['Exact'])
        newexp._studyname = res['new name']
        self.addexp(newexp)

    def feature_cluster(self):
        expdat = self.get_exp_from_selection()
        res = dialog([{'type': 'label', 'label': 'Cluster Features'},
                      {'type': 'int', 'label': 'min reads', 'max': 50000, 'default': 10},
                      {'type': 'string', 'label': 'new name'}], expdat=expdat)
        if res is None:
            return
        if res['new name'] == '':
            res['new name'] = '%s-cluster-features-min-%d' % (expdat._studyname, res['min reads'])
        newexp = cah.cluster_features(expdat, minreads=res['min reads'])
        newexp._studyname = res['new name']
        self.addexp(newexp)

    def feature_filter_fasta(self):
        expdat = self.get_exp_from_selection()
        res = dialog([{'type': 'label', 'label': 'Filter Fasta'},
                      {'type': 'filename', 'label': 'Fasta File'},
                      {'type': 'bool', 'label': 'Negate'},
                      {'type': 'string', 'label': 'new name'}], expdat=expdat)
        if res is None:
            return
        if res['new name'] == '':
            res['new name'] = '%s-cluster-fasta-%s' % (expdat._studyname, res['Fasta File'])
        newexp = cah.filter_fasta(expdat, filename=res['Fasta File'], negate=res['Negate'])
        newexp._studyname = res['new name']
        self.addexp(newexp)

    def feature_sort_abundance(self):
        expdat = self.get_exp_from_selection()
        res = dialog([{'type': 'label', 'label': 'Sort features by abundance'},
                      {'type': 'field', 'label': 'Field', 'withnone': True},
                      {'type': 'value', 'label': 'value'},
                      {'type': 'string', 'label': 'new name'}], expdat=expdat)
        if res is None:
            return
        if res['new name'] == '':
            res['new name'] = '%s-sort-abundance' % expdat._studyname
        newexp = cah.sort_freq(expdat, field=res['field'], value=res['value'])
        newexp._studyname = res['new name']
        self.addexp(newexp)

    def analysis_diff_abundance(self):
        expdat = self.get_exp_from_selection()
        res = dialog([{'type': 'label', 'label': 'Differential abundance'},
                      {'type': 'field', 'label': 'Field', 'withnone': True},
                      {'type': 'value', 'label': 'Value group 1'},
                      {'type': 'value', 'label': 'Value group 2'},
                      {'type': 'string', 'label': 'new name'}], expdat=expdat)
        if res is None:
            return
        if res['new name'] == '':
            res['new name'] = '%s-diff-%s' % (expdat._studyname, res['field'])
        newexp = calour.analysis.getpbfdr(expdat, field=res['field'], val1=res['Value group 1'], val2=res['Value group 2'])
        if newexp is None:
                QtWidgets.QMessageBox.information(self, "No enriched terms found", "No enriched annotations found")
                return
        newexp._studyname = res['new name']
        self.addexp(newexp)

    def test(self):
        # fname,_ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open table file')
        expdat = self.get_exp_from_selection()
        logger.debug(expdat._studyname)

    def add_action_button(self, group, name, function):
        self.actions[group][name] = QPushButton(text=name)
        if group == 'sample':
            self.wSample.addWidget(self.actions[group][name])
        elif group == 'feature':
            self.wFeature.addWidget(self.actions[group][name])
        elif group == 'analysis':
            self.wAnalysis.addWidget(self.actions[group][name])

        self.actions[group][name].clicked.connect(function)

    def listItemRightClicked(self, QPos):
        self.listMenu = QtWidgets.QMenu()
        menurename = self.listMenu.addAction("Rename")
        menurename.triggered.connect(self.menuRename)
        menuremove = self.listMenu.addAction("Delete")
        menuremove.triggered.connect(self.menuRemove)
        menusave = self.listMenu.addAction("Save biom")
        menusave.triggered.connect(self.menuSave)
        menuinfo = self.listMenu.addAction("Info")
        menuinfo.triggered.connect(self.expinfo)
        menusavecommands = self.listMenu.addAction("Save commands")
        menusavecommands.triggered.connect(self.menuSaveCommands)
        parentPosition = self.wExperiments.mapToGlobal(QtCore.QPoint(0, 0))
        self.listMenu.move(parentPosition + QPos)
        self.listMenu.show()

    def expinfo(self):
        pass
        # items = self.bMainList.selectedItems()
        # if len(items) != 1:
        #     print("Need 1 item")
        #     return
        # for citem in items:
        #     cname = str(citem.text())
        #     cexp = self.explist[cname]
        #     # listwin = ListWindow(cexp.filters, cexp.studyname)
        #     res = listwin.exec_()

    def menuRename(self):
        expdat = self.get_exp_from_selection()
        val, ok = QtWidgets.QInputDialog.getText(self, 'Rename experiment', 'old name=%s' % expdat._studyname)
        if ok:
            self.removeexp(expdat)
            expdat._studyname = val
            self.addexp(expdat)

    def menuRemove(self):
        if len(self.wExperiments.selectedItems()) > 1:
            if QtWidgets.QMessageBox.warning(self, "Remove experiments?", "Remove %d experiments?" % len(self.wExperiments.selectedItems()),
                                             QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No) == QtWidgets.QMessageBox.No:
                return
        for currentItemName in self.wExperiments.selectedItems():
            currentItemName = str(currentItemName.text())
            cexp = self._explist[currentItemName]
            self.removeexp(cexp)

    def menuSave(self):
        expdat = self.get_exp_from_selection()
        fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save experiment')
        fname = str(fname)
        if fname == '':
            return
        print('saving')
        # TODO: change to hdf5 once biom bug is solved
        expdat.save_biom(fname, fmt='json')
        print('done')

    def menuSaveCommands(self):
        pass
        # cname = str(self.bMainList.currentItem().text())
        # fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save experiment python commands', '.py')
        # fname = str(fname)
        # hs.savecommands(self.explist[cname], fname)
        # QtWidgets.QMessageBox.information(self, 'Analysis', 'experiment %s commands saved to file:\n%s' % (cname, fname))

    def addexp(self, expdat):
        '''Add a new experiment to the list of experiments

        Parameters
        ----------
        expdat : Experiment
            the experiment to add (note it needs also the _studyname field)
        '''

        # make sure the experiment is not already in the list
        # if so, give a new unique name
        expname = expdat._studyname
        cnum = 2
        expnames = [cexp._studyname for cexp in self._explist.values()]
        while expname in expnames:
            expname = expdat._studyname + '(' + str(cnum) + ')'
            cnum += 1
        expdat._studyname = expname
        expdname = '%s (%s-S, %s-F)' % (expname, expdat.get_num_samples(), expdat.get_num_features())
        expdat._displayname = expdname
        self._explist[expdname] = expdat
        self.wExperiments.addItem(expdname)
        self.wExperiments.clearSelection()
        self.wExperiments.setCurrentRow(self.wExperiments.count()-1)
        logger.debug('experiment %s added' % expname)

    def removeexp(self, exp):
        """
        remove an experiment from the list (and clear)
        """
        expdname = exp._displayname
        del self._explist[expdname]
        items = self.wExperiments.findItems(expdname, QtCore.Qt.MatchExactly)
        for item in items:
            self.wExperiments.takeItem(self.wExperiments.row(item))

    def load(self):
        win = LoadWindow()
        res = win.exec_()
        if res == QtWidgets.QDialog.Accepted:
            tablefname = str(win.wTableFile.text())
            mapfname = str(win.wMapFile.text())
            if mapfname == '':
                mapfname = None
            expname = str(win.wNewName.text())
            exptype = str(win.wType.currentText())
            if exptype == 'Amplicon':
                try:
                    expdat = ca.read_taxa(tablefname, mapfname)
                except:
                    logger.warn('Load for table %s map %s failed' % (tablefname, mapfname))
                    return
                expdat._studyname = expname
            self.addexp(expdat)
            # for biom table show the number of reads`


class LoadWindow(QtWidgets.QDialog):
    def __init__(self):
        super(LoadWindow, self).__init__()
        uic.loadUi(get_ui_file_name('CalourGUILoad.ui'), self)
        self.wTableFileList.clicked.connect(self.browsetable)
        self.wMapFileList.clicked.connect(self.browsemap)

    def browsemap(self):
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open map file')
        fname = str(fname)
        self.wMapFile.setText(fname)

    def browsetable(self):
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open table file')
        fname = str(fname)
        self.wTableFile.setText(fname)
        mname = os.path.dirname(fname) + '/map.txt'
        self.wMapFile.setText(mname)
        pname = os.path.basename(fname)
        self.wNewName.setText(pname)


def dialog(items, expdat=None,  title=None):
    '''Create a dialog with the given items for then experiment

    Parameters
    ----------
    items : list of dict
        Entry for each item to display. fields in the dict are:
            'type' :
                'string' : a string (single line)
                'field' : a field from the sample_metadata
                'value' : a value input for the field in the field item
                'bool' : a boolean
                'label' : a label to display (text in 'label' field)
            'default' : the value to initialize the item to
            'label' : str
                label of the item (also the name in the output dict)
    expdat : Experiment (optional)
        the experiment to use to get the field/values items (needed if item is 'field'/'value')
    title : str (optional)
        title of the dialog

    Returns
    -------
    output : dict or None
        if cancel was selected, return None
        otherwise, a dict with label as key, value as val
    '''
    class DialogWindow(QDialog):
        def __init__(self, items, title=None, expdat=None):
            super().__init__()
            self._expdat = expdat
            if title:
                self.setWindowTitle(title)

            self.main_widget = QWidget(self)
            self.layout = QVBoxLayout(self)

            self.widgets = {}
            for citem in items:
                if citem['type'] == 'label':
                    widget = QLabel(text=citem.get('label'))
                    self.add(widget)
                elif citem['type'] == 'string':
                    widget = QLineEdit(citem.get('default'))
                    self.add(widget, label=citem.get('label'), name=citem.get('label'))
                elif citem['type'] == 'int':
                    widget = QSpinBox()
                    if 'max' in citem:
                        widget.setMaximum(citem['max'])
                    if 'default' in citem:
                        widget.setValue(citem.get('default'))
                    self.add(widget, label=citem.get('label'), name=citem.get('label'))
                elif citem['type'] == 'combo':
                    widget = QComboBox()
                    widget.addItems(citem.get('items'))
                    self.add(widget, label=citem.get('label'), name=citem.get('label'))
                elif citem['type'] == 'field':
                    if expdat is None:
                        logger.warn('Experiment is empty for dialog %s' % title)
                        return None
                    widget = QComboBox()
                    if citem.get('withnone', False):
                        items = ['<none>'] + list(expdat.sample_metadata.columns.values)
                    else:
                        items = expdat.sample_metadata.columns.values
                    widget.addItems(items)
                    self.add(widget, label=citem.get('label'), name='field')
                elif citem['type'] == 'value':
                    if expdat is None:
                        logger.warn('Experiment is empty for dialog %s' % title)
                        return None
                    widget = QLineEdit()
                    self.add(widget, label=citem.get('label'), name=citem.get('label'), addbutton=True)
                elif citem['type'] == 'filename':
                    widget = QLineEdit()
                    self.add(widget, label=citem.get('label'), name='value', addfilebutton=True)
                elif citem['type'] == 'bool':
                    widget = QCheckBox()
                    self.add(widget, label=citem.get('label'), name=citem.get('label'))

            buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)

            buttonBox.accepted.connect(self.accept)
            buttonBox.rejected.connect(self.reject)

            self.layout.addWidget(buttonBox)

        def add(self, widget, name=None, label=None, addbutton=False, addfilebutton=False):
            hlayout = QHBoxLayout()
            if label is not None:
                label_widget = QLabel(label)
                hlayout.addWidget(label_widget)
            hlayout.addWidget(widget)
            if addbutton:
                bwidget = QPushButton(text='...')
                bwidget.clicked.connect(lambda: self.field_vals_click(widget))
                hlayout.addWidget(bwidget)
            if addfilebutton:
                bwidget = QPushButton(text='...')
                bwidget.clicked.connect(lambda: self.file_button_click(widget))
                hlayout.addWidget(bwidget)
            self.layout.addLayout(hlayout)
            self.widgets[name] = widget

        def field_vals_click(self, widget):
            cfield = str(self.widgets['field'].currentText())
            if cfield not in self._expdat.sample_metadata.columns:
                return
            val, ok = QtWidgets.QInputDialog.getItem(self, 'Select value', 'Field=%s' % cfield, list(set(self._expdat.sample_metadata[cfield].astype(str))))
            if ok:
                widget.setText(val)

        def file_button_click(self, widget):
            fname, _x = QtWidgets.QFileDialog.getOpenFileName(self, 'Open fasta file')
            fname = str(fname)
            if fname != '':
                widget.setText(fname)

        def get_output(self, items):
            output = {}
            for citem in items:
                cname = citem.get('label')
                if citem['type'] == 'string':
                    output[cname] = str(self.widgets[cname].text())
                if citem['type'] == 'int':
                    output[cname] = self.widgets[cname].value()
                elif citem['type'] == 'combo':
                    output[cname] = str(self.widgets[cname].currentText())
                elif citem['type'] == 'field':
                    output['field'] = str(self.widgets['field'].currentText())
                    if output['field'] == '<none>':
                        output['field'] = None
                elif citem['type'] == 'value':
                    output[cname] = str(self.widgets[cname].text())
                elif citem['type'] == 'file':
                    output[cname] = str(self.widgets[cname].text())
                elif citem['type'] == 'bool':
                    output[cname] = self.widgets[cname].checkState() > 0
            return output

    aw = DialogWindow(items, expdat=expdat)
    aw.show()
    # if app_created:
    #     app.references.add(self.aw)
    aw.adjustSize()
    res = aw.exec_()
    # if cancel pressed - return None
    if not res:
        return None
    output = aw.get_output(items)
    return output


def get_ui_file_name(filename):
    '''Get the full path to a ui file name filename

    Parameters
    ----------
    filename : str
        file name of the ui file

    Returns
    -------
    uifile : str
        full path to the ui file filename
    '''
    uifile = resource_filename(__name__, 'ui/%s' % filename)
    logger.debug('full path for ui file %s is %s' % (filename, uifile))
    return uifile


def main():
    logger.info('starting Calour GUI')
    app = QtWidgets.QApplication(sys.argv)
    window = AppWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
