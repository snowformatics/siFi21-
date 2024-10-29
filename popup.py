from PyQt5 import QtCore, QtWidgets


class ListSelection(QtWidgets.QDialog):
    def __init__(self, item_ls, parent=None):
        """Pop up dialog to choose the main targets for design plot."""
        super(ListSelection, self).__init__(parent)
        self.setWindowTitle('Select true target')
        layout = QtWidgets.QGridLayout()
        row=0
        self.listWidget = QtWidgets.QListWidget()
        layout.addWidget(self.listWidget, row, 0, 1, 3)
        if item_ls:
            w_item = QtWidgets.QListWidgetItem('siRNA Hits' + '\t' + 'Hit name')
            self.listWidget.addItem(w_item)

            for item in item_ls:
                w_item = QtWidgets.QListWidgetItem(str(item[1]) + '\t' + str(item[0]))
                w_item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
                w_item.setCheckState(QtCore.Qt.Unchecked)
                self.listWidget.addItem(w_item)

        else:
            label = QtWidgets.QLabel("Sorry, no targets found!\nPress Ok to plot the efficiency of the query sequence.")
            label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
            layout.addWidget(label, row, 0, 1, 3)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel,
            QtCore.Qt.Horizontal, self)


        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.setLayout(layout)
        self.setGeometry(300, 200, 460, 350)
        self.centerOnScreen()

    def get_seletced(self, item_ls):
        """Get the selected hits."""

        selected = []
        for index in xrange(self.listWidget.count()):
            item = self.listWidget.item(index)
            if item.checkState() == 2:
                selected.append(item_ls[index-1][0])
        return selected


    def centerOnScreen (self):
        '''Centers the window on the screen.'''
        resolution = QtWidgets.QDesktopWidget().screenGeometry()
        self.move((resolution.width() / 2) - (self.frameSize().width() / 2),
                  (resolution.height() / 2) - (self.frameSize().height() / 2))

    def show_info_message(self, message):
        """Pop up an info message."""
        QtWidgets.QMessageBox.information(self,
                                      u"Information",
                                      message)





