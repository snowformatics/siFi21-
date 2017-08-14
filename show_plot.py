from PyQt4 import QtCore, QtGui
from matplotlib.backends.backend_qt4agg \
    import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg \
    import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg \
    import NavigationToolbar2QT
from matplotlib.patches import Rectangle
import matplotlib.lines as mlines

import tempfile
import general_helpers
import numpy as np
import seaborn as sns
from collections import Counter
import shutil
import os


class NavigationToolbar(NavigationToolbar2QT):
    # We want our own toolbar
    toolitems = [t for t in NavigationToolbar2QT.toolitems if t[0] in ('Home', 'Zoom', 'Save')]


class DrawPlot(QtGui.QMainWindow):
    def __init__(self, sirna_size, query_sequence, sifi_data, off_target_pos_list, region_plot_lst, temp_location,
                 scoret_lst, main_targets, f_in, mode, table_data, parent = None):
        """Class for drawing the plots inside a window."""
        super(DrawPlot, self).__init__(parent)

        # Some plotting settings
        sns.set(style="white", palette="muted", color_codes=True)


        self.sirna_size = sirna_size                    # Size of siRNA chosen by user.
        self.query_sequence = query_sequence            # Query sequence.
        self.query_length = len(self.query_sequence)    # Size of the query sequence.
        self.sifi_data = sifi_data                      # siFi dara in json format file
        self.off_target_pos_list = off_target_pos_list  # Off target position list for design plot
        self.region_plot_lst = region_plot_lst          # Plot for suggested RNAi design region
        self.score_lst = scoret_lst
        self.temp_location = temp_location
        self.main_targets = main_targets
        self.f_in = f_in
        self.mode = mode
        self.home_location = str(QtGui.QDesktopServices.storageLocation(QtGui.QDesktopServices.HomeLocation))
        self.table_data = table_data

        if self.mode == 0:
            self.setWindowTitle('RNAi design plot')
        else:
            self.setWindowTitle('Off-target plot')

        # Image tempfile
        self.temp_img_file = tempfile.mkstemp()

        self.printer = QtGui.QPrinter()
        # Single query mode, get all data
        query = list(general_helpers.iterparse(self.sifi_data))[0]
        # Get number of hits per query
        hit_counter = Counter(player['hit_name'] for player in query)
        hit_overview = hit_counter.most_common()

        dpi = 80
        if self.mode == 0:
            figsize = (13, 4)
        else:
            if len(hit_overview) == 1:
                figsize = (13,5)
            elif len(hit_overview) == 2:
                figsize = (13,6)
            elif len(hit_overview) == 3:
                figsize = (13,8)
            elif len(hit_overview) == 4:
                figsize = (13,10)
            elif len(hit_overview) > 5:
                figsize = (13,20)
                dpi = 65

        self.figure = Figure(figsize =figsize, dpi = dpi, facecolor = '#FFFFFF')
        self.canvas = FigureCanvas(self.figure)
        # use addToolbar to add toolbars to the main window directly!
        #print 'ok'
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.addToolBar(self.toolbar)

        self.main_widget = QtGui.QWidget(self)
        self.setCentralWidget(self.main_widget)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.canvas)

        self.createActions()
        self.createMenus()

        self.main_widget.setLayout(layout)

    def plot_design(self):
        """ Plotting the design plot."""

        # Create main figure
        ax1 = self.figure.add_subplot(111)
        def format_coord(x, y):
            return 'Sequence position = %d , Nr. of siRNAs = %d '%(x, y)
        ax1.format_coord = format_coord

        off_target_dict, main_target_dict, efficient_dict, main_hits_histo = general_helpers.get_target_data(self.f_in, self.sirna_size)

        # If there are no hits, we show only the efficient sirnas for designing a construct
        if self.main_targets == []:
            off_target_dict = {}
            main_target_dict = {}
            main_hits_histo = []

        # Off-target position list
        off_targets_pos = set()
        for i in off_target_dict.values():
            off_targets_pos = off_targets_pos | i

        # Main-target position list
        main_targets_plot = set()
        for i in main_target_dict.values():
            main_targets_plot = main_targets_plot | i

        # Efficient position list
        eff_sirna_plot = []
        for i in efficient_dict.values():
            eff_sirna_plot.extend(i)
        eff_sirna_plot.sort()

        # Draw efficiency plot
        eff_sirna_histo = np.bincount(eff_sirna_plot, minlength=self.query_length)
        ax1.plot(eff_sirna_histo, 'r-', label='Efficient siRNA hits')

        # Draw main target histogram
        main_histo = np.bincount(main_hits_histo, minlength=self.query_length)
        ax1.plot(main_histo, 'b-', label='Main target siRNA hits')

        # Draw main targets as green rectangles inside plot figure
        for region in main_targets_plot:
            someX, someY = region, 0.
            currentAxis = self.figure.gca()
            tri1 = Rectangle((someX, someY), 1, self.sirna_size + 3, color='g', alpha=0.2)#, label='green')
            currentAxis.add_patch(tri1)
            tri1.set_clip_on(False)

        # Draw off-targets as red rectangles inside plot figure
        for off_target in off_targets_pos:
            someX, someY = off_target, 0.
            currentAxis = self.figure.gca()
            tri2 = Rectangle((someX, someY), 1, self.sirna_size + 3, color='r', alpha=0.2)#, label='red')
            currentAxis.add_patch(tri2)
            tri2.set_clip_on(False)

        if max(eff_sirna_histo) > self.sirna_size:
            max_y = max(eff_sirna_histo) + 3
        else:
            max_y = self.sirna_size + 3


        ax1.set_ylim([0, max_y])
        x_text = 'mRNA sequence position\n\n'
        ax1.set_xlabel(x_text, fontsize=12)
        ax1.set_ylabel('Nr of siRNAs', color='b', fontsize=12)
        ax1.set_xlim([0, self.query_length])
        ax1.set_title("RNAi design plot\n\n", fontsize=14)
        self.figure.tight_layout()
        p1 = Rectangle((0, 0), 1, 1, fc="g", alpha=0.2)
        p2 = Rectangle((0, 0), 1, 1, fc="r", alpha=0.2)
        p3 = mlines.Line2D([], [], color='r')
        p4 = mlines.Line2D([], [], color='b')
        ax1.legend([p1,p2, p3, p4], ["Main target", "Off target", 'Efficient siRNAs', 'All siRNAs'],
                    bbox_to_anchor=(0., 1.02, 1., .102), loc=4, ncol=4, mode="", borderaxespad=0.5, frameon=True)

        self.figure.savefig(self.temp_img_file[1] + '.png')
        self.canvas.draw()
        return self.temp_img_file[1] + '.png', off_target_dict, main_target_dict

    def plot_offtarget(self):
        """Plotting the off-target plot."""
        max_sirna_count = self.sirna_size

        # Single query mode, get all data
        query = list(general_helpers.iterparse(self.sifi_data))[0]
        # Get number of hits per query
        hit_counter = Counter(player['hit_name'] for player in query)
        efficicent_counter = Counter(player['hit_name'] for player in query if player['is_efficient'])
        hit_overview = hit_counter.most_common()

        table_data = []
        for x in hit_counter.most_common():
            for y in efficicent_counter.most_common():
                if x[0] == y[0]:
                    table_data.append([x[0], x[1], y[1]])
            if x[0] not in list(efficicent_counter):
                table_data.append([x[0], x[1], 0])

        if len(hit_overview) > 5:
            data = table_data[:5]
            data.append([str(len(table_data)-5) + ' More targets', '...', '...'])
        else:
            data = table_data

        nrows = len(hit_overview[:5]) + 1
        ncols = 1

        ### Create main figure
        axes = self.figure.add_subplot(nrows, 1, 1)


        # First plot is table
        axes.axis('off')

        collabel = ("Targets", "Total siRNA hits", "Efficient siRNA hits")
        the_table = axes.table(cellText=data, colLabels=collabel, loc='center', cellLoc='left')
        cellDict = the_table.get_celld()
        for row, col in cellDict:
            cellDict[(row, col)].set_width(0.08)
        the_table.set_fontsize(12)
        the_table.scale(2.1, 2.1)

        p3 = mlines.Line2D([], [], color='r')
        p4 = mlines.Line2D([], [], color='b')

        # Extract all siRNA position, separate for all and efficient siRNAs
        # Show only up to five targets in plot, otherwise it become to much
        counter = 2

        for hit, nr_hits in hit_overview[:5]:
            all_sirna_plot = []
            eff_sirna_plot = []
            for data in query:
                if hit == data['hit_name']:
                    all_sirna_plot.extend(range(int(data['sirna_position']) + 1, int(data['sirna_position']) + self.sirna_size + 1))
                    if data['is_efficient']:
                        eff_sirna_plot.extend(range(int(data['sirna_position']) + 1, int(data['sirna_position']) + self.sirna_size + 1))
            axes = self.figure.add_subplot(nrows, 1, counter)
            # Create a new plot for each hit
            axes.set_title(str(hit), loc='left')#, color='k', fontsize=14, fontweight='bold')

            # Create a histogram of siRNA positions
            all_sirna_histo = np.bincount(all_sirna_plot, minlength=self.query_length)
            eff_sirna_histo = np.bincount(eff_sirna_plot, minlength=self.query_length)

            # Plot both histograms into one plot
            axes.plot(all_sirna_histo, color="b", label='Total siRNA hits')
            axes.plot(eff_sirna_histo, color="r", label='Strand selected siRNAs')

            # Adjust axis
            if np.amax(all_sirna_histo) > max_sirna_count:
                max_sirna_count = np.amax(all_sirna_histo)

            axes.set_ylim([0, max_sirna_count + 3])
            axes.set_xlim([0, self.query_length + 5])

            # Next plot
            counter += 1

        # Some figure settings
        bottom_offset = (nrows-1)
        if bottom_offset > 5:
            bottom_offset = 5

        axes.set_xlabel('\nRNAi trigger sequence position', fontsize=12)
        axes.set_ylabel('siRNA counts per position', fontsize=12)
        def format_coord(x, y):
            return 'Sequence position = %d , Nr. of siRNAs = %d '%(x, y)
        axes.format_coord = format_coord

        self.figure.subplots_adjust(left=0.05, bottom=0.15/bottom_offset, right=0.98, top=0.99, wspace=None, hspace=0.35)
        self.figure.legend([p4, p3], ['All siRNAs', 'Efficient siRNAs'], loc = "lower right", ncol=2, frameon=True,)
        self.figure.tight_layout()
        self.figure.savefig(self.temp_img_file[1] + '.png')
        self.canvas.draw()
        return self.temp_img_file[1] + '.png', table_data


    def print_(self):
        fileName = self.temp_img_file[1] + '.png'

        if fileName:
            image = QtGui.QImage(fileName)

        self.imageLabel = QtGui.QLabel()
        self.imageLabel.setBackgroundRole(QtGui.QPalette.Base)
        self.imageLabel.setSizePolicy(QtGui.QSizePolicy.Ignored,
                QtGui.QSizePolicy.Ignored)
        self.imageLabel.setScaledContents(True)
        self.imageLabel.setPixmap(QtGui.QPixmap.fromImage(image))
        self.scaleFactor = 1.0

        dialog = QtGui.QPrintDialog(self.printer, self)
        if dialog.exec_():
            painter = QtGui.QPainter(self.printer)
            rect = painter.viewport()
            size = self.imageLabel.pixmap().size()
            size.scale(rect.size(), QtCore.Qt.KeepAspectRatio)
            painter.setViewport(rect.x(), rect.y(), size.width(), size.height())
            painter.setWindow(self.imageLabel.pixmap().rect())
            painter.drawPixmap(0, 0, self.imageLabel.pixmap())

    def export_image(self):
        """Export the images as png."""
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save File', self.home_location, '*.png')
        if filename:
            if not str(filename).endswith('.png'):
                filename += filename + '.png'
            shutil.copy(self.temp_img_file[1] + '.png', filename)
            if os.path.exists(filename):
                message = 'Image successfully saved!'
            else:
                message = 'Could not save image!'
            self.show_info_message(message)

    def export_table(self):
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Export table', self.home_location, '*.csv')
        if filename:
            if not str(filename).endswith('.csv'):
                filename += filename + '.csv'
            try:
                f_in = open(filename, 'w')
                f_in.write("Targets" + ';' + "Total siRNA hits" + ';' + "Efficient siRNA hits" + '\n')
                for data in self.table_data:
                    f_in.write(str(data[0]) + ';' + str(data[1]) + ';' + str(data[2]) + '\n')
                f_in.close()
                if os.path.exists(filename):
                    message = 'Table successfully saved!'
                else:
                    message = 'Could not save table!'
            except IOError:
                message = 'Permission denied! Please close the file.'

            self.show_info_message(message)

    def show_info_message(self, message):
        """Pop up an info message."""
        QtGui.QMessageBox.information(self,
                                      u"Information",
                                      message)


    def createActions(self):
        self.printAct = QtGui.QAction("&Print...", self, enabled=True, triggered=self.print_)
        self.exitAct = QtGui.QAction("E&xit", self, triggered=self.close)

        self.imgAct = QtGui.QAction("Image file", self, triggered=self.export_image)
        self.tableAct = QtGui.QAction("Table (CSV)", self, triggered=self.export_table)

    def createMenus(self):
        self.fileMenu = QtGui.QMenu("&File", self)
        self.fileMenu.addAction(self.printAct)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAct)

        # Export menu
        self.exportMenu = QtGui.QMenu("&Save as", self)
        self.exportMenu.addAction(self.imgAct)

        self.exportMenu.addAction(self.tableAct)

        self.menuBar().addMenu(self.fileMenu)
        self.menuBar().addMenu(self.exportMenu)
