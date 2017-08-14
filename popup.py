from PyQt4 import QtCore, QtGui


class ListSelection(QtGui.QDialog):
    def __init__(self, item_ls, parent=None):
        """Pop up dialog to choose the main targets for design plot."""
        super(ListSelection, self).__init__(parent)
        self.setWindowTitle('Select true target')
        layout = QtGui.QGridLayout()
        row=0
        self.listWidget = QtGui.QListWidget()
        layout.addWidget(self.listWidget, row, 0, 1, 3)
        if item_ls:
            w_item = QtGui.QListWidgetItem('siRNA Hits' + '\t' + 'Hit name')
            self.listWidget.addItem(w_item)

            for item in item_ls:
                w_item = QtGui.QListWidgetItem(str(item[1]) + '\t' + str(item[0]))
                w_item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
                w_item.setCheckState(QtCore.Qt.Unchecked)
                self.listWidget.addItem(w_item)

        else:
            label = QtGui.QLabel("Sorry, no targets found!\nPress Ok to plot the efficiency of the query sequence.")
            label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
            layout.addWidget(label, row, 0, 1, 3)

        buttons = QtGui.QDialogButtonBox(
            QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
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
        resolution = QtGui.QDesktopWidget().screenGeometry()
        self.move((resolution.width() / 2) - (self.frameSize().width() / 2),
                  (resolution.height() / 2) - (self.frameSize().height() / 2))

    def show_info_message(self, message):
        """Pop up an info message."""
        QtGui.QMessageBox.information(self,
                                      u"Information",
                                      message)





