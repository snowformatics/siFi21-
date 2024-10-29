from PyQt5 import QtCore
import database_helpers


class CreateDBThread(QtCore.QThread):
    def __init__(self, database_file_location, db_location, db_name, bowtie_location):
        QtCore.QThread.__init__(self)

        self.database_file_location = database_file_location
        self.db_location = db_location
        self.db_name = db_name
        self.bowtie_location = bowtie_location

    def run(self):
        """Create Bowtie DB."""
        info_message, bowtie_path = database_helpers.create_bowtie_database(self.db_name,
                                                                            self.database_file_location,
                                                                            self.db_location,
                                                                            self.bowtie_location)
        self.finished.connect(print(f"{info_message[0]}, Bowtie ok"))
        # self.emit(QtCore.SIGNAL("threadDone(QString, QString)"), info_message[0], 'Bowtie ok')
