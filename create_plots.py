import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import tempfile
from matplotlib.patches import Rectangle
import matplotlib.lines as mlines

import general_helpers


class CreatePlots(object):
    def __init__(self, sirna_size, query_sequence, sifi_data, off_target_pos_list, region_plot_lst, temp_location,
                 scoret_lst, main_targets, f_in):
        """Class which hold all plotting functions.
           Two types of plots, one for off-target prediction and one for RNAi desgin."""

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
        # Some plotting settings
        sns.set(style="white", palette="muted", color_codes=True)

    def create_design_plot(self):
        """Creates a RNAi design plot."""

        # Temp file for storing
        temp_img_file = tempfile.mkstemp()

        # Create main figure
        fig, ax1 = plt.subplots(figsize=(15, 4))

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
            currentAxis = fig.gca()
            tri1 = Rectangle((someX, someY), 1, self.sirna_size + 3, color='g', alpha=0.2)#, label='green')
            currentAxis.add_patch(tri1)
            tri1.set_clip_on(False)
        
        # Draw off-targets as red rectangles inside plot figure
        for off_target in off_targets_pos:
            someX, someY = off_target, 0.
            currentAxis = fig.gca()
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
        ax1.xaxis.set_ticks(np.arange(0, self.query_length, 50))
        ax1.set_ylabel('Nr of siRNAs', color='b', fontsize=12)
        ax1.set_xlim([0, self.query_length])
        plt.title("RNAi design plot\n\n", fontsize=14)
        plt.tight_layout()
        p1 = Rectangle((0, 0), 1, 1, fc="g", alpha=0.2)
        p2 = Rectangle((0, 0), 1, 1, fc="r", alpha=0.2)
        p3 = mlines.Line2D([], [], color='r')
        p4 = mlines.Line2D([], [], color='b')
        plt.legend([p1,p2, p3, p4], ["Main target", "Off target", 'Efficient siRNAs', 'All siRNAs'],
                    bbox_to_anchor=(0., 1.02, 1., .102), loc=4, ncol=4, mode="", borderaxespad=0.5, frameon=True)
        plt.savefig(temp_img_file[1] + '.png')

        #plt.show()
        plt.close(fig)
        return temp_img_file[1] + '.png', off_target_dict, main_target_dict

    def create_off_target_plot(self):
        """Creates an off-target plot.
           The plot shows two lines, one line show all siRNA counts, the other line show only
           the efficient siRNA counts based on strand selection.
           One plot per hit found in database as subplot per figure.
           Works also in batch mode for queries > 1 (one figure per query)."""

        # Temp file for storing
        temp_img_file = tempfile.mkstemp()

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
        fig, axes = plt.subplots(nrows,ncols,figsize=(12,2.9*nrows))

        # First plot is table
        axes[0].axis('off')

        collabel = ("Targets", "Total siRNA hits", "Efficient siRNA hits")
        the_table = axes[0].table(cellText=data, colLabels=collabel, loc='center', cellLoc='left')
        cellDict = the_table.get_celld()
        for row, col in cellDict:
            cellDict[(row, col)].set_width(0.08)
        the_table.set_fontsize(12)
        the_table.scale(2.1, 2.1)

        # Extract all siRNA position, separate for all and efficient siRNAs
        # Show only up to five targets in plot, otherwise it become to much
        counter = 1
        for hit, nr_hits in hit_overview[:5]:
            all_sirna_plot = []
            eff_sirna_plot = []
            for data in query:
                if hit == data['hit_name']:
                    all_sirna_plot.extend(range(int(data['sirna_position']) + 1, int(data['sirna_position']) + self.sirna_size + 1))
                    if data['is_efficient']:

                    #if data['strand_selection'] or data['end_stability'] or data['target_site_accessibility']:
                        eff_sirna_plot.extend(range(int(data['sirna_position']) + 1, int(data['sirna_position']) + self.sirna_size + 1))

            # Create a new plot for each hit
            axes[counter].set_title(str(hit), loc='left')#, color='k', fontsize=14, fontweight='bold')

            # Create a histogram of siRNA positions
            all_sirna_histo = np.bincount(all_sirna_plot, minlength=self.query_length)
            eff_sirna_histo = np.bincount(eff_sirna_plot, minlength=self.query_length)

            # Plot both histograms into one plot
            axes[counter].plot(all_sirna_histo, color="b", label='Total siRNA hits')
            axes[counter].plot(eff_sirna_histo, color="r", label='Strand selected siRNAs')

            # Adjust axis
            if np.amax(all_sirna_histo) > max_sirna_count:
                max_sirna_count = np.amax(all_sirna_histo)

            axes[counter].set_ylim([0, max_sirna_count + 3])
            axes[counter].set_xlim([0, self.query_length + 5])

            # Next plot
            counter += 1

        ### Some figure settings
        bottom_offset = (nrows-1)
        if bottom_offset > 5:
            bottom_offset = 5

        p3 = mlines.Line2D([], [], color='r')
        p4 = mlines.Line2D([], [], color='b')

        axes[-1].set_xlabel('\nRNAi trigger sequence position', fontsize=12)
        axes[-1].xaxis.set_ticks(np.arange(0, self.query_length, 50))
        axes[-1].set_ylabel('siRNA counts per position', fontsize=12)
        axes[1].legend([p4, p3], ['All siRNAs', 'Efficient siRNAs'], loc=4, ncol=4, mode="", borderaxespad=.5, frameon=True,
                        bbox_to_anchor=(0., 1.02, 1., .102))
        plt.subplots_adjust(left=0.05, bottom=0.25/bottom_offset, right=0.98, top=0.96, wspace=None, hspace=0.35)

        plt.savefig(temp_img_file[1] + '.png')
        plt.close()
        return temp_img_file[1] + '.png', table_data
