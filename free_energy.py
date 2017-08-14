

import math
import warnings

from Bio import SeqUtils, Seq
from Bio import BiopythonWarning


### Taken from Biopython source code:
# http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-pysrc.html

# Thermodynamic lookup tables (dictionaries):
# Enthalpy (dH) and entropy (dS) values for nearest neighbors and initiation
# process. Calculation of duplex initiation is quite different in several
# papers; to allow for a general calculation, all different initiation
# parameters are included in all tables and non-applicable parameters are set
# to zero.
# The key is either an initiation type (e.g., 'init_A/T') or a nearest neighbor
# duplex sequence (e.g., GT/CA, to read 5'GT3'-3'CA5'). The values are tuples
# of dH (kcal/mol), dS (cal/mol K).
# RNA/RNA
 # Freier et al. (1986), Proc Natl Acad Sci USA 83: 9373-9377
RNA_NN1 = {
       'init': (0, -10.8), 'init_A/T': (0, 0), 'init_G/C': (0, 0),
      'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0),
       'sym': (0, -1.4),
      'AA/TT': (-6.6, -18.4), 'AT/TA': (-5.7, -15.5), 'TA/AT': (-8.1, -22.6),
      'CA/GT': (-10.5, -27.8), 'GT/CA': (-10.2, -26.2), 'CT/GA': (-7.6, -19.2),
     'GA/CT': (-13.3, -35.5), 'CG/GC': (-8.0, -19.4), 'GC/CG': (-14.2, -34.9),
      'GG/CC': (-12.2, -29.7)}

 # Xia et al (1998), Biochemistry 37: 14719-14735
RNA_NN2 = {
      'init': (3.61, -1.5), 'init_A/T': (3.72, 10.5), 'init_G/C': (0, 0),
      'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0),
      'sym': (0, -1.4),
      'AA/TT': (-6.82, -19.0), 'AT/TA': (-9.38, -26.7), 'TA/AT': (-7.69, -20.5),
     'CA/GT': (-10.44, -26.9), 'GT/CA': (-11.40, -29.5),
      'CT/GA': (-10.48, -27.1), 'GA/CT': (-12.44, -32.5),
     'CG/GC': (-10.64, -26.7), 'GC/CG': (-14.88, -36.9),
     'GG/CC': (-13.39, -32.7)}

# Chen et al. (2012), Biochemistry 51: 3508-3522
RNA_NN3 = {
    'init': (6.40, 6.99), 'init_A/T': (3.85, 11.04), 'init_G/C': (0, 0),
    'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0),
    'sym': (0, -1.4),
    'AA/TT': (-7.09, -19.8), 'AT/TA': (-9.11, -25.8), 'TA/AT': (-8.50, -22.9),
    'CA/GT': (-11.03, -28.8), 'GT/CA': (-11.98, -31.3),
    'CT/GA': (-10.90, -28.5), 'GA/CT': (-13.21, -34.9),
    'CG/GC': (-10.88, -27.4), 'GC/CG': (-16.04, -40.6),
    'GG/CC': (-14.18, -35.0), 'GT/TG': (-13.83, -46.9),
    'GG/TT': (-17.82, -56.7), 'AG/TT': (-3.96, -11.6),
    'TG/AT': (-0.96, -1.8), 'TT/AG': (-10.38, -31.8), 'TG/GT': (-12.64, -38.9),
    'AT/TG': (-7.39, -21.0), 'CG/GT': (-5.56, -13.9), 'CT/GG': (-9.44, -24.7),
    'GG/CT': (-7.03, -16.8), 'GT/CG': (-11.09, -28.8)}


# Internal mismatch and inosine table (DNA)
# Allawi & SantaLucia (1997), Biochemistry 36: 10581-10594
# Allawi & SantaLucia (1998), Biochemistry 37: 9435-9444
# Allawi & SantaLucia (1998), Biochemistry 37: 2170-2179
# Allawi & SantaLucia (1998), Nucl Acids Res 26: 2694-2701
# Peyret et al. (1999), Biochemistry 38: 3468-3477
# Watkins & SantaLucia (2005), Nucl Acids Res 33: 6258-6267
DNA_IMM1 = {
    'AG/TT': (1.0, 0.9), 'AT/TG': (-2.5, -8.3), 'CG/GT': (-4.1, -11.7),
    'CT/GG': (-2.8, -8.0), 'GG/CT': (3.3, 10.4), 'GG/TT': (5.8, 16.3),
    'GT/CG': (-4.4, -12.3), 'GT/TG': (4.1, 9.5), 'TG/AT': (-0.1, -1.7),
    'TG/GT': (-1.4, -6.2), 'TT/AG': (-1.3, -5.3), 'AA/TG': (-0.6, -2.3),
    'AG/TA': (-0.7, -2.3), 'CA/GG': (-0.7, -2.3), 'CG/GA': (-4.0, -13.2),
    'GA/CG': (-0.6, -1.0), 'GG/CA': (0.5, 3.2), 'TA/AG': (0.7, 0.7),
    'TG/AA': (3.0, 7.4),
    'AC/TT': (0.7, 0.2), 'AT/TC': (-1.2, -6.2), 'CC/GT': (-0.8, -4.5),
    'CT/GC': (-1.5, -6.1), 'GC/CT': (2.3, 5.4), 'GT/CC': (5.2, 13.5),
    'TC/AT': (1.2, 0.7), 'TT/AC': (1.0, 0.7),
    'AA/TC': (2.3, 4.6), 'AC/TA': (5.3, 14.6), 'CA/GC': (1.9, 3.7),
    'CC/GA': (0.6, -0.6), 'GA/CC': (5.2, 14.2), 'GC/CA': (-0.7, -3.8),
    'TA/AC': (3.4, 8.0), 'TC/AA': (7.6, 20.2),
    'AA/TA': (1.2, 1.7), 'CA/GA': (-0.9, -4.2), 'GA/CA': (-2.9, -9.8),
    'TA/AA': (4.7, 12.9), 'AC/TC': (0.0, -4.4), 'CC/GC': (-1.5, -7.2),
    'GC/CC': (3.6, 8.9), 'TC/AC': (6.1, 16.4), 'AG/TG': (-3.1, -9.5),
    'CG/GG': (-4.9, -15.3), 'GG/CG': (-6.0, -15.8), 'TG/AG': (1.6, 3.6),
    'AT/TT': (-2.7, -10.8), 'CT/GT': (-5.0, -15.8), 'GT/CT': (-2.2, -8.4),
    'TT/AT': (0.2, -1.5),
    'AI/TC': (-8.9, -25.5), 'TI/AC': (-5.9, -17.4), 'AC/TI': (-8.8, -25.4),
    'TC/AI': (-4.9, -13.9), 'CI/GC': (-5.4, -13.7), 'GI/CC': (-6.8, -19.1),
    'CC/GI': (-8.3, -23.8), 'GC/CI': (-5.0, -12.6),
    'AI/TA': (-8.3, -25.0), 'TI/AA': (-3.4, -11.2), 'AA/TI': (-0.7, -2.6),
    'TA/AI': (-1.3, -4.6), 'CI/GA': (2.6, 8.9), 'GI/CA': (-7.8, -21.1),
    'CA/GI': (-7.0, -20.0), 'GA/CI': (-7.6, -20.2),
    'AI/TT': (0.49, -0.7), 'TI/AT': (-6.5, -22.0), 'AT/TI': (-5.6, -18.7),
    'TT/AI': (-0.8, -4.3), 'CI/GT': (-1.0, -2.4), 'GI/CT': (-3.5, -10.6),
    'CT/GI': (0.1, -1.0), 'GT/CI': (-4.3, -12.1),
    'AI/TG': (-4.9, -15.8), 'TI/AG': (-1.9, -8.5), 'AG/TI': (0.1, -1.8),
    'TG/AI': (1.0, 1.0), 'CI/GG': (7.1, 21.3), 'GI/CG': (-1.1, -3.2),
    'CG/GI': (5.8, 16.9), 'GG/CI': (-7.6, -22.0),
    'AI/TI': (-3.3, -11.9), 'TI/AI': (0.1, -2.3), 'CI/GI': (1.3, 3.0),
    'GI/CI': (-0.5, -1.3)}

# Terminal mismatch table (DNA)
# SantaLucia & Peyret (2001) Patent Application WO 01/94611
DNA_TMM1 = {
    'AA/TA': (-3.1, -7.8), 'TA/AA': (-2.5, -6.3), 'CA/GA': (-4.3, -10.7),
    'GA/CA': (-8.0, -22.5),
    'AC/TC': (-0.1, 0.5), 'TC/AC': (-0.7, -1.3), ' CC/GC': (-2.1, -5.1),
    'GC/CC': (-3.9, -10.6),
    'AG/TG': (-1.1, -2.1), 'TG/AG': (-1.1, -2.7), 'CG/GG': (-3.8, -9.5),
    'GG/CG': (-0.7, -19.2),
    'AT/TT': (-2.4, -6.5), 'TT/AT': (-3.2, -8.9), 'CT/GT': (-6.1, -16.9),
    'GT/CT': (-7.4, -21.2),
    'AA/TC': (-1.6, -4.0), 'AC/TA': (-1.8, -3.8), 'CA/GC': (-2.6, -5.9),
    'CC/GA': (-2.7, -6.0), 'GA/CC': (-5.0, -13.8), 'GC/CA': (-3.2, -7.1),
    'TA/AC': (-2.3, -5.9), 'TC/AA': (-2.7, -7.0),
    'AC/TT': (-0.9, -1.7), 'AT/TC': (-2.3, -6.3), 'CC/GT': (-3.2, -8.0),
    'CT/GC': (-3.9, -10.6), 'GC/CT': (-4.9, -13.5), 'GT/CC': (-3.0, -7.8),
    'TC/AT': (-2.5, -6.3), 'TT/AC': (-0.7, -1.2),
    'AA/TG': (-1.9, -4.4), 'AG/TA': (-2.5, -5.9), 'CA/GG': (-3.9, -9.6),
    'CG/GA': (-6.0, -15.5), 'GA/CG': (-4.3, -11.1), ' GG/CA': (-4.6, -11.4),
    'TA/AG': (-2.0, -4.7), 'TG/AA': (-2.4, -5.8),
    'AG/TT': (-3.2, -8.7), 'AT/TG': (-3.5, -9.4), 'CG/GT': (-3.8, -9.0),
    'CT/GG': (-6.6, -18.7), 'GG/CT': (-5.7, -15.9), 'GT/CG': (-5.9, -16.1),
    'TG/AT': (-3.9, -10.5), 'TT/AG': (-3.6, -9.8)}

# Dangling ends table (DNA)
# Bommarito et al. (2000), Nucl Acids Res 28: 1929-1934
DNA_DE1 = {
    'AA/.T': (0.2, 2.3), 'AC/.G': (-6.3, -17.1), 'AG/.C': (-3.7, -10.0),
    'AT/.A': (-2.9, -7.6), 'CA/.T': (0.6, 3.3), 'CC/.G': (-4.4, -12.6),
    'CG/.C': (-4.0, -11.9), 'CT/.A': (-4.1, -13.0), 'GA/.T': (-1.1, -1.6),
    'GC/.G': (-5.1, -14.0), 'GG/.C': (-3.9, -10.9), 'GT/.A': (-4.2, -15.0),
    'TA/.T': (-6.9, -20.0), 'TC/.G': (-4.0, -10.9), 'TG/.C': (-4.9, -13.8),
    'TT/.A': (-0.2, -0.5),
    '.A/AT': (-0.7, -0.8), '.C/AG': (-2.1, -3.9), '.G/AC': (-5.9, -16.5),
    '.T/AA': (-0.5, -1.1), '.A/CT': (4.4, 14.9), '.C/CG': (-0.2, -0.1),
    '.G/CC': (-2.6, -7.4), '.T/CA': (4.7, 14.2), '.A/GT': (-1.6, -3.6),
    '.C/GG': (-3.9, -11.2), '.G/GC': (-3.2, -10.4), '.T/GA': (-4.1, -13.1),
    '.A/TT': (2.9, 10.4), '.C/TG': (-4.4, -13.1), '.G/TC': (-5.2, -15.0),
    '.T/TA': (-3.8, -12.6)}

# Dangling ends table (RNA)
# Turner & Mathews (2010), Nucl Acids Res 38: D280-D282
RNA_DE1 = {
    'AA/T.': (-4.9, -13.2), 'AC/T.': (-0.9, -1.3), 'AG/T.': (-5.5, -15.1),
    'AT/T.': (-2.3, -5.5),
    'CA/G.': (-9.0, -23.5), 'CC/G.': (-4.1, -10.6), 'CG/G.': (-8.6, -22.2),
    'CT/G.': (-7.5, -20.31),
    'GA/C.': (-7.4, -20.3), 'GC/C.': (-2.8, -7.7), 'GG/C.': (-6.4, -16.4),
    'GT/C.': (-3.6, -9.7),
    'GA/T.': (-4.9, -13.2), 'GC/T.': (-0.9, -1.3), 'GG/T.': (-5.5, -15.1),
    'GT/T.': (-2.3, -5.5),
    'TA/A.': (-5.7, -16.1), 'TC/A.': (-0.7, -1.9), 'TG/A.': (-5.8, -16.4),
    'TT/A.': (-2.2, -6.8),
    'TA/G.': (-5.7, -16.1), 'TC/G.': (-0.7, -1.9), 'TG/G.': (-5.8, -16.4),
    'TT/G.': (-2.2, -6.8),
    'A./TA': (-0.5, -0.6), 'A./TC': (6.9, 22.6), 'A./TG': (0.6, 2.6),
    'A./TT': (0.6, 2.6),
    'C./GA': (-1.6, -4.5), 'C./GC': (0.7, 3.2), 'C./GG': (-4.6, -14.8),
    'C./GT': (-0.4, -1.3),
    'G./CA': (-2.4, -6.1), 'G./CC': (3.3, 11.6), 'G./CG': (0.8, 3.2),
    'G./CT': (-1.4, -4.2),
    'G./TA': (-0.5, -0.6), 'G./TC': (6.9, 22.6), 'G./TG': (0.6, 2.6),
    'G./TT': (0.6, 2.6),
    'T./AA': (1.6, 6.1), 'T./AC': (2.2, 8.1), 'T./AG': (0.7, 3.5),
    'T./AT': (3.1, 10.6),
    'T./GA': (1.6, 6.1), 'T./GC': (2.2, 8.1), 'T./GG': (0.7, 3.5),
    'T./GT': (3.1, 10.6)}

RNA_DE2 = {
    '.T/AA': (-4.9, -13.2), '.T/CA': (-0.9, -1.3), '.T/GA': (-5.5, -15.1),
    '.T/TA': (-2.3, -5.5),
    '.G/AC': (-9.0, -23.5), '.G/CC': (-4.1, -10.6), '.G/GC': (-8.6, -22.2),
    '.G/TC': (-7.5, -20.31),
    '.C/AG': (-7.4, -20.3), '.C/CG': (-2.8, -7.7), '.C/GG': (-6.4, -16.4),
    '.C/TG': (-3.6, -9.7),
    '.T/AG': (-4.9, -13.2), '.T/CG': (-0.9, -1.3), '.T/GG': (-5.5, -15.1),
    '.T/TG': (-2.3, -5.5),
    '.A/AT': (-5.7, -16.1), '.A/CT': (-0.7, -1.9), '.A/GT': (-5.8, -16.4),
    '.A/TT': (-2.2, -6.8),
    '.G/AT': (-5.7, -16.1), '.G/CT': (-0.7, -1.9), '.G/GT': (-5.8, -16.4),
    '.G/TT': (-2.2, -6.8),
    'AT/.A': (-0.5, -0.6), 'CT/.A': (6.9, 22.6), 'GT/.A': (0.6, 2.6),
    'TT/.A': (0.6, 2.6),
    'AG/.C': (-1.6, -4.5), 'CG/.C': (0.7, 3.2), 'GG/.C': (-4.6, -14.8),
    'TG/.C': (-0.4, -1.3),
    'AC/.G': (-2.4, -6.1), 'CC/.G': (3.3, 11.6), 'GC/.G': (0.8, 3.2),
    'TC/.G': (-1.4, -4.2),
    'AT/.G': (-0.5, -0.6), 'CT/.G': (6.9, 22.6), 'GT/.G': (0.6, 2.6),
    'TT/.G': (0.6, 2.6),
    'AA/.T': (1.6, 6.1), 'CA/.T': (2.2, 8.1), 'GA/.T': (0.7, 3.5),
    'TA/.T': (3.1, 10.6),
    'AG/.T': (1.6, 6.1), 'CG/.T': (2.2, 8.1), 'GG/.T': (0.7, 3.5),
    'TG/.T': (3.1, 10.6)}


def _check(seq, method):
    """Return a sequence which fullfils the requirements of the given method (PRIVATE)."""
    seq = ''.join(seq.split()).upper()
    seq = str(Seq.Seq(seq).back_transcribe())
    if method == 'Tm_Wallace':
        return seq
    if method == 'Tm_GC':
        baseset = ('A', 'B', 'C', 'D', 'G', 'H', 'I', 'K', 'M', 'N', 'R', 'S',
                   'T', 'V', 'W', 'X', 'Y')
    if method == 'Tm_NN':
        baseset = ('A', 'C', 'G', 'T', 'I')
    seq = ''.join([base for base in seq if base in baseset])
    return seq


def salt_correction(Na=0, K=0, Tris=0, Mg=0, dNTPs=0, method=1, seq=None):
    """Calculate a term to correct Tm for salt ions."""

    if method in (5, 6, 7) and not seq:
        raise ValueError('sequence is missing (is needed to calculate ' +
                         'GC content or sequence length).')
    if seq:
        seq = str(seq)
    corr = 0
    if not method:
        return corr
    Mon = Na + K + Tris / 2.0  # Note: all these values are millimolar
    mg = Mg * 1e-3             # Lowercase ions (mg, mon, dntps) are molar
    # Na equivalent according to von Ahsen et al. (2001):
    if sum((K, Mg, Tris, dNTPs)) > 0 and not method == 7:
        if dNTPs < Mg:
            # dNTPs bind Mg2+ strongly. If [dNTPs] is larger or equal than
            # [Mg2+], free Mg2+ is considered not to be relevant.
            Mon += 120 * math.sqrt(Mg - dNTPs)
    mon = Mon * 1e-3
    # Note: math.log = ln(), math.log10 = log()
    if method in range(1, 7) and not mon:
        raise ValueError('Total ion concentration of zero is not allowed in ' +
                         'this method.')
    if method == 1:
        corr = 16.6 * math.log10(mon)
    if method == 2:
        corr = 16.6 * math.log10((mon) / (1.0 + 0.7 * (mon)))
    if method == 3:
        corr = 12.5 * math.log10(mon)
    if method == 4:
        corr = 11.7 * math.log10(mon)
    if method == 5:
        corr = 0.368 * (len(seq) - 1) * math.log(mon)
    if method == 6:
        corr = (4.29 * SeqUtils.GC(seq) / 100 - 3.95) * 1e-5 * math.log(mon) + \
            9.40e-6 * math.log(mon) ** 2
    if method == 7:
        a, b, c, d = 3.92, -0.911, 6.26, 1.42
        e, f, g = -48.2, 52.5, 8.31
        if dNTPs > 0:
            dntps = dNTPs * 1e-3
            ka = 3e4  # Dissociation constant for Mg:dNTP
            # Free Mg2+ calculation:
            mg = (-(ka * dntps - ka * mg + 1.0) +
                  math.sqrt((ka * dntps - ka * mg + 1.0) ** 2 + 4.0 * ka * mg)) / (2.0 * ka)
        if Mon > 0:
            R = math.sqrt(mg) / mon
            if R < 0.22:
                corr = (4.29 * SeqUtils.GC(seq) / 100 - 3.95) * \
                    1e-5 * math.log(mon) + 9.40e-6 * math.log(mon) ** 2
                return corr
            elif R < 6.0:
                a = 3.92 * (0.843 - 0.352 * math.sqrt(mon) * math.log(mon))
                d = 1.42 * (1.279 - 4.03e-3 * math.log(mon) -
                            8.03e-3 * math.log(mon) ** 2)
                g = 8.31 * (0.486 - 0.258 * math.log(mon) +
                            5.25e-3 * math.log(mon) ** 3)
        corr = (a + b * math.log(mg) + (SeqUtils.GC(seq) / 100) *
                (c + d * math.log(mg)) + (1 / (2.0 * (len(seq) - 1))) *
                (e + f * math.log(mg) + g * math.log(mg) ** 2)) * 1e-5
    if method > 7:
        raise ValueError('Allowed values for parameter \'method\' are 1-7.')
    return corr


def calculate_free_energy(seq, check=True, strict=True, c_seq=None, shift=0, nn_table=RNA_NN3,
          tmm_table=DNA_TMM1, imm_table=DNA_IMM1, de_table=RNA_DE2,
          dnac1=25, dnac2=0, selfcomp=False, Na=20, K=50, Tris=0, Mg=0,
          dNTPs=0, saltcorr=5):
    """Return the delatG using nearest neighbor thermodynamics."""
    #print shift
    print seq
    seq = str(seq)
    if not c_seq:
        # c_seq must be provided by user if dangling ends or mismatches should
        # be taken into account. Otherwise take perfect complement.
        c_seq = Seq.Seq(seq).complement()
    c_seq = str(c_seq)
    if check:
        seq = _check(seq, 'Tm_NN')
        c_seq = _check(c_seq, 'Tm_NN')
    tmpseq = seq
    tmp_cseq = c_seq
    deltaH = 0
    deltaS = 0
    dH = 0  # Names for indexes
    dS = 1  # 0 and 1
    #print tmpseq, tmp_cseq
    # Dangling ends?
    if shift or len(seq) != len(c_seq):
        # Align both sequences using the shift parameter
        if shift > 0:
            tmpseq = '.' * shift + seq
        if shift < 0:
            tmp_cseq = '.' * abs(shift) + c_seq

        if len(tmp_cseq) > len(tmpseq):
            tmpseq += (len(tmp_cseq) - len(tmpseq)) * '.'
        if len(tmp_cseq) < len(tmpseq):
            tmp_cseq += (len(tmpseq) - len(tmp_cseq)) * '.'
        # Remove 'over-dangling' ends
        while tmpseq.startswith('..') or tmp_cseq.startswith('..'):
            tmpseq = tmpseq[1:]
            tmp_cseq = tmp_cseq[1:]
        while tmpseq.endswith('..') or tmp_cseq.endswith('..'):
            tmpseq = tmpseq[:-1]
            tmp_cseq = tmp_cseq[:-1]
        #print tmpseq, tmp_cseq
        # Now for the dangling ends
        if tmpseq.startswith('.') or tmp_cseq.startswith('.'):
            left_de = tmpseq[:2] + '/' + tmp_cseq[:2]
            #print 'left ', left_de
            deltaH += de_table[left_de][dH]
            deltaS += de_table[left_de][dS]
            tmpseq = tmpseq[1:]
            tmp_cseq = tmp_cseq[1:]
        if tmpseq.endswith('.') or tmp_cseq.endswith('.'):
            right_de = tmp_cseq[-2:][::-1] + '/' + tmpseq[-2:][::-1]
            deltaH += de_table[right_de][dH]
            deltaS += de_table[right_de][dS]
            tmpseq = tmpseq[:-1]
            tmp_cseq = tmp_cseq[:-1]
    # Now for terminal mismatches
    left_tmm = tmp_cseq[:2][::-1] + '/' + tmpseq[:2][::-1]
    if left_tmm in tmm_table:
        deltaH += tmm_table[left_tmm][dH]
        deltaS += tmm_table[left_tmm][dS]
        tmpseq = tmpseq[1:]
        tmp_cseq = tmp_cseq[1:]
    right_tmm = tmpseq[-2:] + '/' + tmp_cseq[-2:]
    if right_tmm in tmm_table:
        deltaH += tmm_table[right_tmm][dH]
        deltaS += tmm_table[right_tmm][dS]
        tmpseq = tmpseq[:-1]
        tmp_cseq = tmp_cseq[:-1]

    # Now everything 'unusual' at the ends is handled and removed and we can
    # look at the initiation.
    # One or several of the following initiation types may apply:

    # Type: General initiation value
    deltaH += nn_table['init'][dH]
    deltaS += nn_table['init'][dS]

    # Type: Duplex with no (allA/T) or at least one (oneG/C) GC pair
    if SeqUtils.GC(seq) == 0:
        deltaH += nn_table['init_allA/T'][dH]
        deltaS += nn_table['init_allA/T'][dS]
    else:
        deltaH += nn_table['init_oneG/C'][dH]
        deltaS += nn_table['init_oneG/C'][dS]

    # Type: Penalty if 5' end is T
    if seq.startswith('T'):
        deltaH += nn_table['init_5T/A'][dH]
        deltaS += nn_table['init_5T/A'][dS]
    if seq.endswith('A'):
        deltaH += nn_table['init_5T/A'][dH]
        deltaS += nn_table['init_5T/A'][dS]

    # Type: Different values for G/C or A/T terminal basepairs
    ends = seq[0] + seq[-1]
    AT = ends.count('A') + ends.count('T')
    GC = ends.count('G') + ends.count('C')
    deltaH += nn_table['init_A/T'][dH] * AT
    deltaS += nn_table['init_A/T'][dS] * AT
    deltaH += nn_table['init_G/C'][dH] * GC
    deltaS += nn_table['init_G/C'][dS] * GC

    # Finally, the 'zipping'
    for basenumber in range(len(tmpseq) - 1):
        neighbors = tmpseq[basenumber:basenumber + 2] + '/' + \
            tmp_cseq[basenumber:basenumber + 2]
        #print neighbors
        if neighbors in imm_table:
            deltaH += imm_table[neighbors][dH]
            deltaS += imm_table[neighbors][dS]
        elif neighbors[::-1] in imm_table:
            deltaH += imm_table[neighbors[::-1]][dH]
            deltaS += imm_table[neighbors[::-1]][dS]
        elif neighbors in nn_table:
            deltaH += nn_table[neighbors][dH]
            deltaS += nn_table[neighbors][dS]
        elif neighbors[::-1] in nn_table:
            deltaH += nn_table[neighbors[::-1]][dH]
            deltaS += nn_table[neighbors[::-1]][dS]
        else:
            # We haven't found the key...
            if strict:
                raise ValueError('no data for neighbors \'' + neighbors + '\'')
            else:
                warnings.warn('no data for neighbors \'' + neighbors +
                              '\'. Calculation will be wrong',
                              BiopythonWarning)

    k = (dnac1 - (dnac2 / 2.0)) * 1e-9
    if selfcomp:
        k = dnac1 * 1e-9
        deltaH += nn_table['sym'][dH]
        deltaS += nn_table['sym'][dS]
    R = 1.987  # universal gas constant in Cal/degrees C*Mol
    if saltcorr:
        corr = salt_correction(Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs,
                               method=saltcorr, seq=seq)
    if saltcorr == 5:
        deltaS += corr

    tao = 273.15 + 22 # Constant temperature tao in Kelvin
    deltaG = (deltaH * 1000 - tao * deltaS) / 1000
    #print deltaG

    Tm = (1000 * deltaH) / (deltaS + (R * (math.log(k)))) - 273.15

    if saltcorr in (1, 2, 3, 4):
        Tm += corr
    if saltcorr in (6, 7):
        Tm = (1 / (1 / (Tm + 273.15) + corr) - 273.15)

    return deltaG