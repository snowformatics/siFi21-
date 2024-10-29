# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'db_wizard.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtWidgets, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtWidgets.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtWidgets.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtWidgets.QApplication.translate(context, text, disambig)

class Ui_wizard(object):
    def setupUi(self, wizard):
        wizard.setObjectName(_fromUtf8("wizard"))
        wizard.resize(676, 387)
        self.wizardPage_2 = QtWidgets.QWizardPage()
        self.wizardPage_2.setObjectName(_fromUtf8("wizardPage_2"))
        self.layoutWidget = QtWidgets.QWidget(self.wizardPage_2)
        self.layoutWidget.setGeometry(QtCore.QRect(20, 20, 601, 141))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.layoutWidget)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.commandLinkButton_seq = QtWidgets.QCommandLinkButton(self.layoutWidget)
        self.commandLinkButton_seq.setMaximumSize(QtCore.QSize(200, 16777215))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/Images/Images/open.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.commandLinkButton_seq.setIcon(icon)
        self.commandLinkButton_seq.setIconSize(QtCore.QSize(48, 48))
        self.commandLinkButton_seq.setObjectName(_fromUtf8("commandLinkButton_seq"))
        self.verticalLayout_3.addWidget(self.commandLinkButton_seq)
        self.label_5 = QtWidgets.QLabel(self.layoutWidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_5.setFont(font)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.verticalLayout_3.addWidget(self.label_5)
        wizard.addPage(self.wizardPage_2)
        self.wizardPage_3 = QtWidgets.QWizardPage()
        self.wizardPage_3.setObjectName(_fromUtf8("wizardPage_3"))
        self.layoutWidget1 = QtWidgets.QWidget(self.wizardPage_3)
        self.layoutWidget1.setGeometry(QtCore.QRect(20, 20, 401, 91))
        self.layoutWidget1.setObjectName(_fromUtf8("layoutWidget1"))
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.layoutWidget1)
        self.verticalLayout_4.setObjectName(_fromUtf8("verticalLayout_4"))
        self.label_6 = QtWidgets.QLabel(self.layoutWidget1)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_6.setFont(font)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.verticalLayout_4.addWidget(self.label_6)
        self.lineEdit = QtWidgets.QLineEdit(self.layoutWidget1)
        self.lineEdit.setMaximumSize(QtCore.QSize(200, 16777215))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.lineEdit.setFont(font)
        self.lineEdit.setObjectName(_fromUtf8("lineEdit"))
        self.verticalLayout_4.addWidget(self.lineEdit)
        wizard.addPage(self.wizardPage_3)
        self.wizardPage2 = QtWidgets.QWizardPage()
        self.wizardPage2.setObjectName(_fromUtf8("wizardPage2"))
        self.layoutWidget2 = QtWidgets.QWidget(self.wizardPage2)
        self.layoutWidget2.setGeometry(QtCore.QRect(20, 10, 351, 131))
        self.layoutWidget2.setObjectName(_fromUtf8("layoutWidget2"))
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.layoutWidget2)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.label_2 = QtWidgets.QLabel(self.layoutWidget2)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_2.setFont(font)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout_2.addWidget(self.label_2)
        self.progressBar = QtWidgets.QProgressBar(self.layoutWidget2)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.progressBar.setFont(font)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName(_fromUtf8("progressBar"))
        self.verticalLayout_2.addWidget(self.progressBar)
        self.label_3 = QtWidgets.QLabel(self.layoutWidget2)
        self.label_3.setText(_fromUtf8(""))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout_2.addWidget(self.label_3)
        wizard.addPage(self.wizardPage2)

        self.retranslateUi(wizard)
        QtCore.QMetaObject.connectSlotsByName(wizard)

    def retranslateUi(self, wizard):
        wizard.setWindowTitle(_translate("wizard", "Create new database", None))
        self.commandLinkButton_seq.setText(_translate("wizard", "Choose source file", None))
        self.label_5.setText(_translate("wizard", "No source file selected", None))
        self.label_6.setText(_translate("wizard", "Enter a database name:", None))
        self.label_2.setText(_translate("wizard", "Progress:", None))

import sifi_2015_rc
