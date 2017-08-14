from PyQt4 import QtGui, QtCore, QtNetwork
import FileDialog

import sys
import webbrowser
import time
import os

import general_helpers
import database_helpers
import db_wizard
import sifi_pipeline
import imageviewer

from Resources.ui_sifi2015 import Ui_MainWindow

# pyrcc4 sifi_2015.qrc > sifi_2015_rc.py
# pyuic4 sifi2015.ui > ui_sifi2015.py
# pyuic4 wizard.ui > wizard_ui.py


class MyPopup(QtGui.QWidget):
    def __init__(self, images_location):
        QtGui.QWidget.__init__(self)
        self.setWindowTitle('sifi21_1.2.3-0008')
        self.images_location = images_location
        label = QtGui.QLabel(self)
        pixmap = QtGui.QPixmap(self.images_location + "about.tif")
        label.setPixmap(pixmap)

        label.mousePressEvent = self.open_url

    def open_url(self, event):
        webbrowser.open("http://www.snowformatics.com/")

    def center(self):
        screen = QtGui.QDesktopWidget().screenGeometry()
        size = self.frameSize()
        self.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)


class MyMainWindow(QtGui.QMainWindow):
    def __init__(self, *args):
        """Main window of si-Fi."""
        QtGui.QMainWindow.__init__(self, *args)
        self.ui = Ui_MainWindow()

        # resolution
        #res = QtGui.QDesktopWidget().screenGeometry()
        #self.resize(100, 200)

        # Some important locations
        # Does not work reliable on all OS, that's why we split of "local" and add application folder name manually
        self.data_location = str(QtGui.QDesktopServices.storageLocation(QtGui.QDesktopServices.DataLocation))
        self.data_location = self.data_location.split('Local')[0] + '/Local/'
        self.app_location = self.data_location + '/siFi2015/'
        self.home_location = str(QtGui.QDesktopServices.storageLocation(QtGui.QDesktopServices.HomeLocation))
        self.temp_location = str(QtGui.QDesktopServices.storageLocation(QtGui.QDesktopServices.TempLocation))
        self.images_location = self.app_location + '/Images/'
        self.bowtie_location = self.app_location + '/Bowtie/'
        self.rnaplfold_location = self.app_location + '/RNAplfold/'
        self.db_location = self.app_location + '/Databases/'


        # Create required folders and copy required files
        try:
            general_helpers.create_folders([self.app_location, self.images_location, self.bowtie_location,
                                            self.rnaplfold_location, self.db_location],
                                           [self.app_location + 'logging_sifi.txt', self.app_location + 'last_db.txt'])
            general_helpers.copying_files(self.app_location)
            self.log_file = self.app_location + '/logging_sifi.txt'
            self.lastdb_file = self.app_location + '/last_db.txt'

        except (IOError, OSError):
            self.show_info_message('Could not create folder or file.')

        self.setWindowIcon(QtGui.QIcon(self.images_location + "siFi21_icon64.tif"))

        # Setup UI

        self.ui.setupUi(self)
        self.ui.plainTextEdit_seq.setFocus()
        self.ui.tabWidget.tabBar().setTabButton(0, QtGui.QTabBar.RightSide, None)
        self.ui.tabWidget.tabBar().setTabButton(1, QtGui.QTabBar.RightSide, None)
        #self.ui.statusbar.showMessage("System Status | Normal | " + self.app_location)

        # Aks user which mode to use
        self.show_mode_question()

        # Connect widgets
        self.create_connects()

        # Inspect all bowtie databases created so far by the user and update combo
        self.available_databases()


    def create_connects(self):
        """Connect widgets to functions."""
        self.ui.commandLinkButton_remove.clicked.connect(self.delete_dbs)
        self.ui.commandLinkButton_add.clicked.connect(self.start_db_wizard)
        self.ui.actionRNAi_design.triggered.connect(self.set_design_settings)
        self.ui.actionOff_target_prediction_2.triggered.connect(self.set_offtarget_settings)
        self.ui.pushButton_run.clicked.connect(self.start_sifi_pipeline)
        self.ui.pushButton_default.clicked.connect(self.default_settings)
        self.ui.actionDocumentation.triggered.connect(self.show_help)
        self.ui.actionAbout.triggered.connect(self.show_about_message)
        self.ui.pushButton_openDB.clicked.connect(self.open_file_and_insert_seq)
        self.connect(self.ui.comboBox_miss, QtCore.SIGNAL("currentIndexChanged(const QString&)"), self.test)

    def test(self):
        widget_lst = [self.ui.checkBox_end,  self.ui.doubleSpinBox_end, self.ui.checkBox_access,
                      self.ui.doubleSpinBox_access,  self.ui.spinBox_xmer, self.ui.label_5,
                      self.ui.checkBox_strand, self.ui.checkBox_terminal]
        if self.ui.comboBox_miss.currentIndex() > 0:
            for widget in widget_lst:
                widget.setEnabled(False)
        else:
            self.default_settings()


    def start_sifi_pipeline(self):
        """
        Starts the siFi pipeline.
        Two modes available, off-target prediction and RNAi design.
        Parameters:
        bowtie_db, db_location, query_sequences,
        sirna_size, mismatches,
        accessibility_check, accessibility_window,
        rnaplfold_location, bowtie_location, mode,
        strand_check, end_check,
        end_stability_treshold,
        target_site_accessibility_treshold, temp_location
        """

        no_efficience = True
        if self.ui.checkBox_access.isEnabled() and self.ui.checkBox_strand.isEnabled() and self.ui.checkBox_end.isEnabled() and self.ui.checkBox_terminal.isEnabled():
            no_efficience = False

        # Get the sequence and validate it.
        seq_file = self.validate_seq()
        # If sequence is IUPAC, start si-Fi
        if seq_file:
            # Catch no main target selection
            sifi = sifi_pipeline.SifiPipeline(str(self.ui.comboBox_db.currentText()), self.db_location, seq_file,
                                              int(self.ui.spinBox_size.value()), int(self.ui.comboBox_miss.currentText()),
                                              int(self.ui.checkBox_access.checkState()), int(self.ui.spinBox_xmer.value()),
                                              self.rnaplfold_location, self.bowtie_location, self.mode,
                                              int(self.ui.checkBox_strand.checkState()), int(self.ui.checkBox_end.checkState()),
                                              float(self.ui.doubleSpinBox_end.value()),
                                              float(self.ui.doubleSpinBox_access.value()), self.temp_location,
                                              int(self.ui.checkBox_terminal.checkState()), no_efficience)
            temp_img_file, json_data, table_data, main_target_range, off_target_range, msg = sifi.run_pipeline
            self.save_last_db(str(self.ui.comboBox_db.currentText()))

            if temp_img_file == None and json_data == None:
                self.show_info_message(msg)

            # # If run was successful, show the results plot in a new window.
            # if temp_img_file != None and json_data != None:
            #     if os.path.exists(temp_img_file):
            #         self.save_last_db(str(self.ui.comboBox_db.currentText()))
            #         image_window = imageviewer.ImageViewer(self, temp_img_file, self.mode, table_data,
            #                                                main_target_range, off_target_range, seq_file)
            #         image_window.closed.connect(self.show)
            #         image_window.show()
            #         image_window.open(temp_img_file)
            # else:
            #     self.show_info_message(msg)

    def validate_seq(self):
        """
        Validate the input sequence.
        Nucleotides should contain only IUPAC characters.
        Single fasta file or plain text allowed.
        """

        sequences = str(self.ui.plainTextEdit_seq.toPlainText())

        plain_text = general_helpers.validate_seq(sequences)
        fasta = general_helpers.validate_fasta_seq(sequences)

        if plain_text:
            sequences = '>' + 'my_sequence' + '\n' + sequences
            sequence_temp_file = general_helpers.save_seq_file(sequences)
        elif fasta == 1:
            sequence_temp_file = general_helpers.save_seq_file(sequences)
        elif fasta > 1:
            sequence_temp_file = False
            self.show_info_message("Please enter only one sequence or use the batch mode.")
        else:
            sequence_temp_file = False
            self.show_info_message('Please enter a valid nucleic acid sequence!')
        return sequence_temp_file

    def start_db_wizard(self):
        """Starts DB wizard."""
        self.wizard = db_wizard.DBWizard(self, self.data_location, self.app_location,
                                         self.home_location, self.db_location, self.bowtie_location)
        self.wizard.show()

    def available_databases(self):
        """Get all created databases and update widgtes"""
        database_dict = database_helpers.all_dbs(str(self.db_location))
        # Update Combobox with database names
        self.ui.comboBox_db.clear()
        # Get last DB and show it on first position
        last_db = self.get_last_db()

        all_dbs = database_dict['Database name']
        try:
            index = all_dbs.index(last_db)
            all_dbs.insert(0, all_dbs.pop(index))
        except:
            pass
        self.ui.comboBox_db.addItems(all_dbs)
        database_dict = database_helpers.all_dbs(str(self.db_location))
        # Update Table with database information's
        self.ui.tableWidget.clear()
        self.ui.tableWidget.setRowCount(len(database_dict['Database name']))
        self.ui.tableWidget.setColumnCount(len(database_dict.keys()))
        table_headers = database_dict.keys()
        for row, key in enumerate(database_dict.keys()):
            for column, item in enumerate(database_dict[key]):
                new_item = QtGui.QTableWidgetItem(str(item))
                if row == 0:
                    new_item.setCheckState(QtCore.Qt.Unchecked)
                self.ui.tableWidget.setItem(column, row, new_item)
        self.ui.tableWidget.setHorizontalHeaderLabels(table_headers)

    def delete_dbs(self):
        """Delete all databases the user has selected."""
        db_to_delete = []
        for row in range(self.ui.tableWidget.rowCount()):
            for column in range(self.ui.tableWidget.rowCount()):
                item = self.ui.tableWidget.item(row, column)
                if item is not None:
                    if item.checkState() == QtCore.Qt.Checked:
                        checked_name = self.ui.tableWidget.item(row, column).text()
                        if not checked_name in db_to_delete:
                            db_to_delete.append(str(checked_name))

        if db_to_delete:
            delete = QtGui.QMessageBox.question(self, u"Delete database?",
                                                u"Are you sure that you want to delete the database(s)?",
                                                QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if delete == QtGui.QMessageBox.Yes:
                info_message = database_helpers.delete_databases(db_to_delete, self.db_location)
                time.sleep(1)
                self.show_info_message(info_message[0])
        else:
            self.show_info_message("Please select a database!")
        # Update widgets
        self.available_databases()

    def save_last_db(self, db):
        """Saves the last database choice."""
        f_in = open(self.lastdb_file, 'w')
        f_in.write(str(self.ui.comboBox_db.currentText()))
        f_in.close()

    def get_last_db(self):
        """Gets the last database choice."""
        last_db = open(self.lastdb_file, 'r').read()
        return last_db

    def show_mode_question(self):
        """Pop up a mode message questions.
           Adjust the settings corresponding to choice."""
        msg_box = QtGui.QMessageBox()
        font = QtGui.QFont()
        font.setPointSize(12)
        msg_box.setFont(font)
        msg_box.setText('Which mode would you like to use?')
        msg_box.setWindowTitle('Choose mode')
        btn_design = QtGui.QPushButton(' RNAi design ')
        msg_box.addButton(btn_design, QtGui.QMessageBox.YesRole)
        btn_check = QtGui.QPushButton(' Off-target prediction ')
        msg_box.addButton(btn_check, QtGui.QMessageBox.NoRole)
        self.mode = msg_box.exec_()

        # For off-target prediction, disable efficiency settings
        if self.mode == 1:
            self.set_offtarget_settings()

    def open_file_and_insert_seq(self):
        file_format = u"Sequence (*.fasta *.fas *.txt)"
        # Clear plainTextEdit
        self.ui.plainTextEdit_seq.clear()
        # Get path for file path

        path = str(self.open_sequence_file(file_format))
        if path != "None":
            self.ui.label_1.setText(path)
            self.insert_seq(path, self.ui.plainTextEdit_seq)

    def insert_seq(self, file_location, plain_textedit):
        """Insert a sequence file content into a text edit widget."""
        f_seq = open(file_location, 'r')
        sequence = f_seq.read()
        plain_textedit.clear()
        plain_textedit.insertPlainText(sequence)
        f_seq.close()

    def open_sequence_file(self, file_format):
        """Open a sequence file and return the path."""
        sequence_file_location = QtGui.QFileDialog.getOpenFileName(
                                    self,
                                    u"Open a sequence file",
                                    self.home_location,
                                    file_format)
        if not sequence_file_location.isNull():
            if os.path.exists(sequence_file_location):
                return str(sequence_file_location)

    def default_settings(self):
        """Sets default settings."""
        if self.mode == 1:
            self.set_offtarget_settings()
        else:
            self.set_design_settings()

    def set_offtarget_settings(self):
        """Sets the off-target default settings."""
        widget_lst = [self.ui.checkBox_end,  self.ui.doubleSpinBox_end, self.ui.checkBox_access,
                      self.ui.doubleSpinBox_access,  self.ui.spinBox_xmer, self.ui.label_5,
                      self.ui.checkBox_strand, self.ui.checkBox_terminal]
        for widget in widget_lst:
            widget.setEnabled(True)
        self.ui.doubleSpinBox_end.setValue(1)
        self.ui.doubleSpinBox_access.setValue(0.1)
        self.ui.spinBox_xmer.setValue(8)
        self.ui.tabWidget.setTabText(0, 'Off-target prediction')
        self.ui.groupBox_1.setTitle("Paste RNAi trigger sequence or import a file")
        self.ui.spinBox_size.setValue(21)
        self.ui.comboBox_miss.setCurrentIndex(0)
        self.ui.checkBox_end.setChecked(False)
        self.ui.checkBox_access.setChecked(False)
        self.ui.checkBox_strand.setChecked(True)
        self.ui.checkBox_terminal.setChecked(True)
        self.ui.plainTextEdit_seq.setFocus()
        self.mode = 1


    def set_design_settings(self):
        """Sets the RNAi design default settings."""
        widget_lst = [self.ui.checkBox_end,  self.ui.doubleSpinBox_end, self.ui.checkBox_access,
                      self.ui.doubleSpinBox_access,  self.ui.spinBox_xmer, self.ui.label_5,
                      self.ui.checkBox_strand, self.ui.checkBox_terminal]
        for widget in widget_lst:
            widget.setEnabled(True)
        self.ui.checkBox_end.setChecked(True)
        self.ui.checkBox_access.setChecked(True)
        self.ui.checkBox_strand.setChecked(True)
        self.ui.checkBox_terminal.setChecked(True)
        self.ui.tabWidget.setTabText(0, 'RNAi design')
        self.ui.groupBox_1.setTitle("")
        self.ui.spinBox_size.setValue(21)
        self.ui.comboBox_miss.setCurrentIndex(0)
        self.ui.doubleSpinBox_end.setValue(1)
        self.ui.doubleSpinBox_access.setValue(0.1)
        self.ui.spinBox_xmer.setValue(8)
        self.ui.plainTextEdit_seq.setFocus()
        self.mode = 0

    def show_info_message(self, message):
        """Pop up an info message."""
        QtGui.QMessageBox.information(self,
                                      u"Information",
                                      message)

    def show_about_message(self):
        """Shows an about box."""
        self.w = MyPopup(self.images_location)
        self.w.setGeometry(QtCore.QRect(100, 100, 481, 380))
        self.w.center()
        self.w.show()


    def show_help(self):
        """Show online documentation in a browser."""
        webbrowser.open("http://labtools.ipk-gatersleben.de/si-Fi/Quick_help.pdf")


def main():
    app = QtGui.QApplication(sys.argv)
    screen_size = app.desktop().screenGeometry()
    width, height = screen_size.width(), screen_size.height()

    QtGui.QApplication.setStyle(QtGui.QStyleFactory.create("gtk"))
    QtGui.QApplication.setPalette(QtGui.QApplication.style().standardPalette())
    myapp = MyMainWindow()

    if height < 830:
        #Resize width and height
        myapp.resize(640, height-50)

    myapp.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
