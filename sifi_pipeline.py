import subprocess
import tempfile
from Bio import SeqIO
import os
import numpy as np
import json
from collections import Counter
from Bio.Seq import Seq

import free_energy
#import create_plots
import popup
import general_helpers
import show_plot


class SifiPipeline(object):
    def __init__(self, bowtie_db, db_location, query_sequences, sirna_size, mismatches, accessibility_check,
                 accessibility_window, rnaplfold_location, bowtie_location, mode, strand_check, end_check,
                 end_stability_treshold, target_site_accessibility_treshold, temp_location, terminal_check,
                 no_efficience):

        """Class for si-Fi pipeline:

            si-Fi pipeline for RNAi design.
           1. Save query sequence as fasta file.
           2. Split query sequence into xmers and store as fasta file.
           3. Run BOWTIE against DB and extract all positions of off-targets and main targets.
           4. For each siRNA, do strand selection, end stability and target site accessibility for efficiency.
           5. Start RNAplfold to get pair probabilities.
           6. Store data into json file for plotting.

           si-Fi pipeline for off-target prediction.
           1. Save query sequence as fasta file.
           2. Split query sequence into xmers and store as fasta file.
           3. Run BOWTIE against DB and extract all information.
           4. For each hit, do strand selection for efficiency.
           5. Store data into json file for plotting."""

        self.bowtie_db = bowtie_db                                                                                      # Bowtie DB
        self.query_sequences = query_sequences                                                                          # List of all query sequences in multi fasta format
        self.sirna_size = sirna_size                                                                                    # siRNA size
        self.mismatches = mismatches                                                                                    # Allowed mismatches
        self.db_location = db_location                                                                                  # DB path
        self.rnaplfold_location = rnaplfold_location                                                                    # Rnaplfold path
        self.bowtie_location = bowtie_location                                                                          # Bowtie path
        self.temp_location = temp_location                                                                              # Temp location
        self.mode = mode                                                                                                # Mode, either RNAi design or off-target prediction

        self.strand_check = strand_check                                                                                # Strand selection is enabled or disabled
        self.end_check = end_check                                                                                      # End stability selection is enabled or disabled
        self.accessibility_check = accessibility_check                                                                  # Target site accessibility is enabled or disabled
        self.accessibility_window = accessibility_window                                                                # Accessibility window
        self.end_stability_treshold = end_stability_treshold                                                            # End stability treshold
        self.ts_accessibility_treshold = target_site_accessibility_treshold                                             # Target site accessibility threshold
        self.terminal_check = terminal_check
        self.no_efficience = no_efficience

        # Some constants
        self.winsize = 80                                                                                               # Average the pair probabil. over windows of given size
        self.span = 40                                                                                                  # Set the maximum allowed separation of a base pair to span
        self.temperature = 22                                                                                           # Temperature for calculation free energy
        self.sirna_start_position = 0
        self.overhang = 2                                                                                               # siRNA overhang
        self.end_nucleotides = 3                                                                                        # siRNA end nucleotides

    @property
    def run_pipeline(self):
        """Start the si-Fi pipeline either in off-target or design mode."""

        # Iterate over all query sequences
        for seq_record in SeqIO.parse(self.query_sequences, "fasta"):
            # Store ID and sequence
            query_name = seq_record.id
            query_sequence = str(seq_record.seq)
            self.len_seq = len(query_sequence)

            # Create a temp single fasta file of the query sequence
            query_seq_temp_file = self.create_single_fasta_file(query_name, query_sequence)

            # Create all siRNAs of size "sirna_size" and save them into multi fasta file
            sirna_file, tab_file_name = self.create_sirnas(query_sequence, self.sirna_size)

            # Run BOWTIE against DB
            bowtie_data = self.run_bowtie(sirna_file, self.bowtie_db, self.mismatches)
            self.bowtie_data_l = self.bowtie_to_lst(bowtie_data)

            # Run RNAplfold
            lunp_data = self.run_rnaplfold(query_name, query_seq_temp_file)

            # Design mode
            if self.mode == 0:
                # Get main targets
                main_targets = self.get_main_target()
                if main_targets == None:
                    # We still want to see an efficiency plot, even without hits
                    no_target = True
                else:
                    no_target = False
            # Off target mode, we don't need main targets
            else:
                main_targets = None


            self.sirna_l = []
            f_in = open(tab_file_name, 'r')
            sirna_data = f_in.readlines()
            for x in sirna_data:
                x = x.split('\t')
                self.sirna_l.append(x)

            sirna_sequence_n2 = self.sirna_l

            # In case we get hits to DB
            if self.bowtie_data_l:
                # Both modes with hits
                input_data = self.bowtie_data_l
                no_target = False


            # No hits to DB, off-target plot will be not shown, design plot only show efficient sirnas
            else:
                # Design mode with no hits, we show only efficient siRNAs
                if self.mode == 0:
                    f_in = open(tab_file_name, 'r')
                    sirna_data = f_in.readlines()
                    l = []
                    for x in sirna_data:
                        x = x.split('\t')
                        l.append(x)
                    input_data = l
                    f_in.close()
                    no_target = True
                else:
                    # Off target mode without hits, no plot
                    input_data = None

            if main_targets == "Canceled":
                return None, None, None, None, None, 'Canceled by user.'
            else:
                if input_data:
                    # Json temp
                    temp_json_file = tempfile.mkstemp()
                    out_file = open(temp_json_file[1] + '.json', "w")
                    json_lst = self.data_to_json(query_name, input_data, no_target, lunp_data, main_targets)
                    json.dump(json_lst, out_file, indent=4)
                    out_file.close()

                    if json_lst:
                        sifi_data = open(temp_json_file[1] + '.json', "r").read()
                        table_data = general_helpers.get_table_data(temp_json_file[1] + '.json')
                        for x in table_data:
                            print(x[0], x[1], x[2])
                        plot = show_plot.DrawPlot(self.sirna_size, query_sequence, sifi_data, None, None,
                                                self.temp_location, None, main_targets,
                                                temp_json_file[1] + '.json', self.mode, table_data)
                        plot.show()

                        if self.mode == 0:
                            temp_img_file, off_target_dict, main_target_dict = plot.plot_design()
                        else:
                            temp_img_file, table_data = plot.plot_offtarget()
                            # Only for design
                            main_target_dict = None
                            off_target_dict = None
                    else:
                        return None, None, None, None, None, 'No targets found. Please make sure that the query and/or database sequences are in correct orientation.'
                    return temp_img_file, temp_json_file[1] + '.json', table_data, main_target_dict, off_target_dict, None
                else:
                    return None, None, None, None, None, 'No targets found. Please make sure that the query and/or database sequences are in correct orientation'

    def bowtie_to_lst(self, bowtie_data):
        """Converts Bowtie data into lists of list. Just for convenience."""
        bowtie_lst = []
        for sirna_nr in range(0, len(bowtie_data)):
            bowtie_data_split = bowtie_data[sirna_nr].strip()
            bowtie_data_split = bowtie_data_split.split('\t')
            # If we have missmatches, we will append a value for each entry.
            if len(bowtie_data_split) == 7:
                bowtie_data_split.append(0)
            bowtie_lst.append(bowtie_data_split)
        return bowtie_lst

    def data_to_json(self, query_name, input_data, no_target, lunp_data, main_targets):
        """Extracts the data from bowtie results and put everything into json format.
           Efficiency is calculated fro each siRNA.
           If no target is found, for design mode the siRNA fasta file is used instead of Bowtie data."""

        json_lst = []
        for data_split in input_data:
            if not no_target:
                sirna_name = data_split[0]
                strand = data_split[1]
                hit_name = data_split[2]
                # Target sequence position coming from Bowtie (0-based offset)
                reference_strand_pos = int(data_split[3])
                # Position on query sequence starting from 1
                query_position = int(sirna_name.split('sirna')[1])
                sirna_sequence = data_split[4]
                missmatches = data_split[7]

                if self.mode == 0:
                    if hit_name in main_targets:
                        off_target = False
                    else:
                        off_target = True
                else:
                    off_target = None
            else:
                # We use the siRNA fasta file to get the efficiency information for each siRNA
                sirna_name = data_split[0]
                query_position = int(sirna_name.split('sirna')[1])
                sirna_sequence = data_split[1]
                # We don't have this information because we got no Bowtie hits
                off_target = False
                strand = None
                hit_name = False
                reference_strand_pos = None
                missmatches = None

            if strand == '+':
                print('ok')
                # We need the antisense siRNA (c_seq) for the energy
                # Antisense sequence = sequence position - 2
                # First two siRNas are ignored
                # Example, siRNA3 -> Antisense siRNA1 as input for c_seq
                #antisense_sequence = self.sirna_l[query_position-3][1]

                # Calculate strand selection for each siRNA
                #if self.strand_check:
                    # We must ignore the first two siRNAs because we can not calculate free energy
                if query_position == 1 or query_position == 2:
                    sirna_sequence_n2 = None
                else:
                    sirna_sequence_n2 = self.sirna_l[query_position-3][1].strip()
            #else:
                #sirna_sequence_n2 = None

                # query_position = self.len_seq - query_position - self.sirna_size + 2
                # print query_position, self.len_seq - self.sirna_size
                #
                # if query_position > self.len_seq - self.sirna_size - 3:
                #     sirna_sequence_n2 = None
                # else:
                #     s = sirna_sequence
                #     sirna_sequence_n2 = Seq(s)
                #     sirna_sequence_n2 = str(sirna_sequence_n2.reverse_complement())
                #     sirna_sequence = self.sirna_l[query_position+3][1].strip()

            # print sirna_name, sirna_sequence_n2, sirna_sequence, query_position

                lunp_data_xmer = lunp_data[int(sirna_name.split('sirna')[1])-1, :].astype(np.float).tolist()[self.accessibility_window]

                is_efficient, strand_selection, end_stability, target_site_accessibility, thermo_effcicient = self.calculate_efficiency(sirna_sequence, sirna_sequence_n2, lunp_data_xmer)

                json_dict = {"query_name": query_name, "sirna_name":sirna_name,
                             "sirna_position": query_position, "sirna_sequence": sirna_sequence,
                             "is_efficient": is_efficient,
                             "strand_selection": strand_selection, "end_stability": end_stability,
                             "target_site_accessibility": target_site_accessibility,
                             "accessibility_value": lunp_data_xmer, "is_off_target": off_target,
                             "hit_name": hit_name, "reference_strand_pos":reference_strand_pos,
                             "strand": strand, "mismatches": missmatches}
                json_lst.append(json_dict)

        return json_lst

    def get_main_target(self):
        """ The user must choose the main targets from all bowtie hit names.
           All remaining hits will be marked as off-targets."""

        # Get all hit names and number of hits
        all_targets = Counter(player[2] for player in self.bowtie_data_l)
        # Open popup to choose main target
        lb = popup.ListSelection(sorted(all_targets.items(), key=lambda x: x[1], reverse=True))
        if lb.exec_():
            main_targets = lb.get_seletced(sorted(all_targets.items(), key=lambda x: x[1], reverse=True))
            return main_targets
        else:
            return "Canceled"

    def create_sirnas(self, query_sequence, sirna_size):
        """Create siRNA's of size "sirna_size" of a sequence.
        Return a path with the siRNA multi fasta temp file"""
        # Fasta format
        temp_query_file = tempfile.mkstemp()
        # Tab format
        temp_tab_file = tempfile.mkstemp()
        # Slice over sequence and split into xmers.
        start = self.sirna_start_position
        end = sirna_size
        seq_list = []
        for dummy_x in range(len(query_sequence)):
            if len(query_sequence[start:end]) == sirna_size:
                # Position of siRNA
                sirna_position = start
                # Original sequence
                sirna = query_sequence[start:end]
                sirna.upper()
                seq_list.append((('sirna'+ str(sirna_position+1)), sirna))
                start += 1
                end += 1
        # Store siRNAs in multiple fasta format
        sirna_file_name, tab_file_name = self.create_multi_fasta_file(seq_list, temp_query_file[1], temp_tab_file[1])
        return sirna_file_name, tab_file_name

    def create_single_fasta_file(self, query_name, query_sequence):
        """Create a temp file of the query sequence in fasta format."""
        temp_query_file = tempfile.mkstemp()
        f_temp = open(temp_query_file[1], 'w')
        f_temp.write('>' + str(query_name) + '\n')
        f_temp.write(str(query_sequence) + '\n')
        f_temp.close()
        return temp_query_file[1]

    def create_multi_fasta_file(self, sequence_list, file_name, tab_file_name):
        """Create a multiple fasta temp file."""
        f_fasta = open(file_name, 'w')
        f_tab = open(tab_file_name, 'w')
        for name, sequence in sequence_list:
            f_fasta.write('> ' + str(name) + '\n')
            f_fasta.write(str(sequence) + '\n')
            # Tab format
            f_tab.write(str(name) + '\t' + str(sequence) + '\n')
        f_fasta.close()
        return file_name, tab_file_name

    def run_bowtie(self, sequence, database_name, mismatches):
        """Run BOWTIE alignment."""
        temp_bowtie_file = tempfile.mkstemp()
        os.chdir(self.bowtie_location)
        process = subprocess.Popen(["bowtie", "-a", "-v", str(mismatches),  "-y",
                                    self.db_location / database_name, "-f",
                                    sequence, temp_bowtie_file[1]])
        process.wait()

        if os.path.exists(temp_bowtie_file[1]):
            bowtie_data = open(temp_bowtie_file[1], 'r').readlines()
        else:
            bowtie_data = ''
        return bowtie_data

    def run_rnaplfold(self, query_name, sequence_file):
        """Run RNAplfold."""
        os.chdir(self.rnaplfold_location)
        #sequence_file = "c:\\users\\lueck~1.ta-\\appdata\local\\temp\\tmpl90aux"
        seq = open(sequence_file, 'r').read()

        cwd = tempfile.mkdtemp()
        prc_stdout = subprocess.PIPE
        prc = subprocess.Popen(['RNAplfold', '-W', '%d'%self.winsize,'-L', '%d'% self.span, '-u', '%d'%self.sirna_size, '-T', '%.2f'%self.temperature], stdin=subprocess.PIPE, stdout=prc_stdout, cwd=cwd)
        prc.stdin.write(seq)
        prc.stdin.write('\n')
        prc.communicate()

        if os.path.exists(cwd + '/' + query_name + '_lunp'):
            lunp_file = cwd + '/' + query_name + '_lunp'
            lunp_data = np.loadtxt(lunp_file, dtype='str')
            # Delete first lines self.sirna_size-1 because they are not complete
            lunp_data = np.delete(lunp_data, np.r_[:self.sirna_size-1], 0)
        else:
            lunp_file = ''
        return lunp_data

    def free_energy3(self, sirna_sequence):
        """Calculate the free energy of a sequence.

           Code was taken from Biopython.
           http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-pysrc.html

           Example for 21mer
           siRNA GGGATGGCTCAAAGGCGTAGT
           Sense5prime_MFE siRNA position [0,1,2] -> GGG
           Antisense5prime_MFE siRNA position [17,18,29] -> GTA"""

        # Sense5_MFE
        sense_five_seq = sirna_sequence[self.sirna_start_position:self.end_nucleotides]
        # Anitsense5_MFE
        antisense_five_seq = sirna_sequence[self.sirna_size-self.overhang-self.end_nucleotides:self.sirna_size-self.overhang]

        #sense5_MFE_enegery = free_energy.calculate_free_energy(sense_five_seq)
        sense5_MFE_enegery = free_energy.calculate_free_energy(sense_five_seq)
        anti_sense5_MFE_enegery = free_energy.calculate_free_energy(antisense_five_seq)

        #print sense_five_seq, sense5_MFE_enegery, antisense_five_seq, anti_sense5_MFE_enegery
        return sense5_MFE_enegery, anti_sense5_MFE_enegery

    def free_energy_dangling_ends(self, sirna_sequence, sirna_sequence_n2):
        """Calculate the free energy of a sequence.
           c_seq is for dangling ends"""
        #print 'rc', Seq(sirna_sequence_n2).reverse_complement()

        # Sense5_MFE
        #sense_five_seq = sirna_sequence[self.sirna_start_position:self.end_nucleotides+1]
        sense_five_seq = sirna_sequence[self.sirna_start_position:self.end_nucleotides]
        #sense_c_seq = Seq(sirna_sequence_n2).reverse_complement().strip()[self.sirna_size-6:self.sirna_size-1]
        sense_c_seq = Seq(sirna_sequence_n2).reverse_complement().strip()[self.sirna_size-5:self.sirna_size-1]

        # Anitsense5_MFE
        antisense_five_seq = Seq(sirna_sequence_n2).reverse_complement().strip()[self.sirna_start_position:self.end_nucleotides]
        antisense_c_seq = sirna_sequence[self.sirna_size-5:self.sirna_size-1]

        #print 'sense ',  sirna_sequence, sense_five_seq, sense_c_seq[::-1]
        sense5_MFE_enegery = free_energy.calculate_free_energy(sense_five_seq, check=True, strict=True, c_seq=sense_c_seq[::-1], shift=1)
       # print 'G ', sense5_MFE_enegery

        #print 'antisense ', sirna_sequence_n2, antisense_five_seq, antisense_c_seq[::-1]
        anti_sense5_MFE_enegery = free_energy.calculate_free_energy(antisense_five_seq, check=True, strict=True, c_seq=antisense_c_seq[::-1], shift=1)
        #print 'G ', anti_sense5_MFE_enegery

        return sense5_MFE_enegery, anti_sense5_MFE_enegery

    def strand_selection(self, sirna_sequence, sirna_sequence_n2):
        """Returns whether the strand will be selected (True) or not (False) based on energy rules."""
        # For siRNA n>3 we calculate with dangling ends
        if sirna_sequence_n2 != None:
            sense5_MFE_enegery, anti_sense5_MFE_enegery = self.free_energy_dangling_ends(sirna_sequence, sirna_sequence_n2)
        else:
            #pass
            # For the first two siRNAs no dangling ends
            sense5_MFE_enegery, anti_sense5_MFE_enegery = self.free_energy3(sirna_sequence)

        if anti_sense5_MFE_enegery >= sense5_MFE_enegery:
            strand_selection = True
        else:
            strand_selection = False
        #print "G_S/G_AS ", sense5_MFE_enegery, anti_sense5_MFE_enegery
        return strand_selection

    def end_stability(self, sirna_sequence, sirna_sequence_n2):
        """Calculate whether the end stability is higher or equal threshold (default=1).
           Return True if yes and False if it is lower."""
        # For siRNA n>3 we calculate with dangling ends
        if sirna_sequence_n2 != None:
            sense5_MFE_enegery, anti_sense5_MFE_enegery = self.free_energy_dangling_ends(sirna_sequence, sirna_sequence_n2)
        else:
            #pass
            # For the first two siRNAs no dangling ends
            sense5_MFE_enegery, anti_sense5_MFE_enegery = self.free_energy3(sirna_sequence)

        # End stability
        if (anti_sense5_MFE_enegery - sense5_MFE_enegery) >= self.end_stability_treshold:
            end_stability = True
        else:
            end_stability = False
        return end_stability

    def pair_probability(self, lunp_data_xmer):
        """Calculates whether the pair probability the siRNA at a certain window (default 8) is higher or equal
           the threshold (default=0.1)
           Return True if yes and False if it is lower."""
        # Get pair probability of the accessibility window (chosen by user) of siRNA
        #print "accessibility ", lunp_data_xmer
        if lunp_data_xmer >= self.ts_accessibility_treshold:
            return True
        else:
            return False

    def check_efficient(self, strand_selection, end_stability, target_site_accessibility):
        """Calculates whether a siRNA is efficient or not. Default priority rules (if checked by user):
           1. Strand selection must be True.
           2. End stability must be higher or equal threshold (default=1).
           3. Pair probability (lunp data) of xmer (chosen by user, default 8) must be higher or equal threshold (default=0.1).
           If all rules apply, a siRNA is efficient."""

        is_efficient = None
        if self.strand_check and self.end_check and self.accessibility_check:
            if strand_selection and end_stability and target_site_accessibility:
                is_efficient = True
            else:
                is_efficient = False
        if self.strand_check and self.end_check and not self.accessibility_check:
            if strand_selection and end_stability:
                is_efficient = True
            else:
                is_efficient = False
        if self.strand_check and self.accessibility_check and not self.end_check:
            if strand_selection and target_site_accessibility:
                is_efficient = True
            else:
                is_efficient = False
        if self.end_check and self.accessibility_check and not self.strand_check:
            if end_stability and target_site_accessibility:
                is_efficient = True
            else:
                is_efficient = False
        if self.strand_check and not self.end_check and not self.accessibility_check:
            if strand_selection:
                is_efficient = True
            else:
                is_efficient = False
        if self.end_check and not self.strand_check and not self.accessibility_check:
            if end_stability:
                is_efficient = True
            else:
                is_efficient = False
        if self.accessibility_check and not self.strand_check and not self.end_check:
            if target_site_accessibility:
                is_efficient = True
            else:
                is_efficient = False
        if not self.accessibility_check and not self.strand_check and not self.end_check:
            is_efficient = True

        return is_efficient



    def calculate_efficiency(self, sirna_sequence, sirna_sequence_n2, lunp_data_xmer):
        """"""
        is_efficient = None

        if self.no_efficience:
            is_efficient = False
            strand_selection = None
            end_stability = None
            target_site_accessibility = None
            thermo_effcicient = None
        else:
            strand_selection = self.strand_selection(sirna_sequence, sirna_sequence_n2)
            end_stability = self.end_stability(sirna_sequence, sirna_sequence_n2)
            target_site_accessibility = self.pair_probability(lunp_data_xmer)
            thermo_effcicient = self.check_efficient(strand_selection, end_stability, target_site_accessibility)

            if self.terminal_check:
                if sirna_sequence[self.sirna_size-3] == 'A' or sirna_sequence[self.sirna_size-3] == 'T':
                    #print "A/T at S19"
                    if sirna_sequence[1] == 'A' or sirna_sequence[1] == 'T':
                        #print "A/T at S1"
                        if thermo_effcicient:
                            is_efficient = True
                        else:
                            is_efficient = False
                    else:
                        if self.accessibility_check:
                            if target_site_accessibility:
                                is_efficient = True
                            else:
                                is_efficient = False
                        else:
                            is_efficient = True
                else:
                    if sirna_sequence[1] == 'G' or sirna_sequence[1] == 'C':
                        #print "C/G at S1"
                        if thermo_effcicient:
                            is_efficient = True
                        else:
                            is_efficient = False
                    else:
                        is_efficient = False
            else:
                if thermo_effcicient:
                    is_efficient = True
                else:
                    is_efficient = False
        return is_efficient, strand_selection, end_stability, target_site_accessibility, thermo_effcicient
        #print sirna_sequence, strand_selection, end_stability, target_site_accessibility, thermo_effcicient, is_efficient
