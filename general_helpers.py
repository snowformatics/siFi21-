import logging
import platform
import os
import time
import shutil
import re
import tempfile
import json
from operator import itemgetter
from itertools import groupby
from collections import Counter

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature


def start_logging(method, log_file):
    """Writes into a log file the current status and OS."""
    logging.basicConfig(filename=log_file, level=logging.DEBUG)
    logging.exception('Got exception on main handler')
    logging.debug(time.strftime("%d.%m.%Y um %H:%M:%S Uhr"))
    logging.debug(str(platform.system()+platform.release()))
    logging.debug(str(method))


def create_folders(folder_lst, file_lst):
    """Create folders of si-Fi at the application data location."""
    for folder in folder_lst:
        if not os.path.exists(folder):
            os.mkdir(folder)

    for fi in file_lst:
        if not os.path.exists(fi):
            f_in = open(fi, 'w')
            f_in.close()

def copying_files(app_location):
    """Copy important files for si-Fi into the application data location."""
    to_copy_folders = os.listdir(os.getcwd() + '/ToCopy/')
    shutil.copyfile(os.getcwd() + '/ToCopy/Images/about.tif', app_location + '/Images/about.tif')
    for folder in to_copy_folders:
        to_copy_files = os.listdir(os.getcwd() + '/ToCopy/' + folder)
        for files in to_copy_files:
            if not folder.startswith('.'):
                if not files.startswith('.'):
                    if not os.path.exists(app_location + '/' + folder + '/' + files):
                        shutil.copyfile(os.getcwd() + '/ToCopy/' + folder + '/' + files, app_location + '/' + folder + '/' + files)

def validate_seq(sequence):
    """Validate plain text DNA sequence without header."""
    sequence = sequence.strip()
    sequence = sequence.replace(" ", "")
    sequence = sequence.replace("\n", "")
    sequence.upper()
    regex = re.compile('^[ACTGNRYSWKMBDHVEFILPQSXZ]*$', re.I)
    if regex.search(sequence) is not None:
        return True
    else:
        return False


def validate_fasta_seq(sequence):
    """Validate sequence in single Fasta format."""
    sequence = sequence.replace(" ", "")
    sequence.upper()
    regex = re.compile('>\S*\n[ACTGNRYSWKMBDHVEFILPQSXZ]*', re.MULTILINE)
    if len(regex.findall(sequence)) > 0:
        return len(regex.findall(sequence))
    else:
        return False


def save_seq_file(sequence):
    """Saves a sequence as fasta format into a file."""
    temp_out = tempfile.mkstemp()
    f_seq = open(temp_out[1], 'w')
    f_seq.write(str(sequence))
    f_seq.close()
    return temp_out[1]


def iterparse(j):
    """Work around because json can not load multiple objects.
       Taken from http://stackoverflow.com/questions/22112439/valueerror-extra-data-while-loading-json."""
    nonspace = re.compile(r'\S')
    decoder = json.JSONDecoder()
    pos = 0
    while True:
        matched = nonspace.search(j, pos)
        if not matched:
            break
        pos = matched.start()
        decoded, pos = decoder.raw_decode(j, pos)
        yield decoded


def prepare_json_data(f_in):
    """Prepares the json file for easier parsing."""
    sifi_data = open(f_in, "r").read()
    data = list(iterparse(sifi_data))[0]
    return data

def get_table_data(f_in):
    """Extracts a summary of all and efficient hits."""
    query = prepare_json_data(f_in)

    # Get number of hits per query
    hit_counter = Counter(player['hit_name'] for player in query)
    efficicent_counter = Counter(player['hit_name'] for player in query if player['is_efficient'])
    #hit_overview = hit_counter.most_common()
    table_data = []
    for x in hit_counter.most_common():
        for y in efficicent_counter.most_common():
            if x[0] == y[0]:
                #print y[1]
                table_data.append([x[0], x[1], y[1]])
        if x[0] not in list(efficicent_counter):
            table_data.append([x[0], x[1], 0])
    return table_data

def get_target_data(f_in, sirna_size):
    """Extracts the target hits and positions."""
    query = prepare_json_data(f_in)
    off_target_positions = set()
    off_target_dict = {}
    main_target_positions = set()
    main_target_dict = {}
    efficient_dict = {}
    main_hits_histo = []
    ready_sirnas = []
    main_target_ready = []
    for data in query:
        #print data['sirna_name'], data['sirna_position'], data['sirna_position']+21, data['is_efficient'], data['sirna_sequence']
        plot_range = set(range(int(data['sirna_position']), int(data['sirna_position']) + sirna_size))
        plot_range_l = range(int(data['sirna_position']), int(data['sirna_position']) + sirna_size)
        if data['is_efficient']:
            if data['sirna_name'] not in ready_sirnas:
                if data['hit_name'] in efficient_dict.keys():
                    efficient_dict[data['hit_name']].extend(plot_range_l)
                else:
                    efficient_dict[data['hit_name']] = plot_range_l
                ready_sirnas.append(data['sirna_name'])
        if data['is_off_target']:
            off_target_positions = off_target_positions | plot_range
            off_target_dict[data['hit_name']] = off_target_positions
        else:
            if (data['sirna_name'], data['sirna_position']) not in main_target_ready:
                main_target_ready.append((data['sirna_name'], data['sirna_position']))
                main_target_positions = main_target_positions | plot_range
                main_target_dict[data['hit_name']] = main_target_positions - off_target_positions
                main_hits_histo.extend(plot_range)

    return off_target_dict, main_target_dict, efficient_dict, main_hits_histo


def group_ranges(data):
    ranges = []
    for k, g in groupby(enumerate(data), lambda (i,x):i-x):
        group = map(itemgetter(1), g)
        ranges.append((group[0], group[-1]))
    return ranges


def create_gbk(main_target_dict, off_target_dict, seq_file, out_file):
    """Create a genbank file."""

    shutil.copy(seq_file, seq_file + '.fasta')
    data = SeqIO.read(open(seq_file + '.fasta'), "fasta")
    my_sequence = Seq(str(data.seq))
    my_sequence_record = SeqRecord(my_sequence)
    my_sequence_record.seq.alphabet = generic_dna

    for target, positions in main_target_dict.iteritems():
        group = group_ranges(positions)
        for x in group:
            my_feature_location = FeatureLocation(x[0], x[1])
            my_feature_type = "MT " + str(target)
            my_feature = SeqFeature(my_feature_location, type=my_feature_type)
            my_sequence_record.features.append(my_feature)

    for target, positions in off_target_dict.iteritems():
        group = group_ranges(positions)
        for x in group:
            my_feature_location = FeatureLocation(x[0], x[1])
            my_feature_type = "OT " + str(target)
            my_feature = SeqFeature(my_feature_location, type=my_feature_type)
            my_sequence_record.features.append(my_feature)

    today = time.strftime("%d") + '-' + time.strftime("%b").upper() + '-' + time.strftime("%Y")
    my_sequence_record.annotations["date"] = today
    gbk_file = open(out_file, 'w')
    SeqIO.write(my_sequence_record, gbk_file, "genbank")
