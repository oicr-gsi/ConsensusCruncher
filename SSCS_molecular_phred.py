#!/usr/bin/env python3

###############################################################
#
#                        Median SSCS Maker
#
# Author: Nina Wang
# Date Created: Jun 23, 2016
###############################################################
# Function: 
# Written for Python 3.5.1
#
# Inputs: 
# 1. A position-sorted paired-end BAM file containing reads with a duplex tag 
#    in the header.  
#
# Outputs:
# 1. A paired-end BAM file containing single stranded consensus sequences.
# 2. A singleton BAM file containing all read families with single reads.
# 3. A tag family size distribution plot (x-axis: family size, y-axis: number of reads).
# 4. A text file containing summary statistics (Total reads, SSCS reads, 
#    singletons, rescued reads)
#    - Read family: reads that share the same molecular barcode, chr, and start
#                   coordinates for Read1 and Read2
#    - Singleton: a read family containing only one member (a single read)
#
# usage: SSCS_maker.py [--cutoff CUTOFF] [--Ncutoff NCUTOFF] [--infile INFILE] 
#                      [--outfile OUTFILE]
# optional arguments:
# --cutoff CUTOFF     Percentage of nucleotides at a given position in a 
#                     sequence required to be identical for a consensus [0.7]
#                        Example (--cutoff = 0.7):
#                           Four reads (readlength = 10) are as follows:
#                              Read 1: ACTGATACTT
#                              Read 2: ACTGAAACCT
#                              Read 3: ACTGATACCT
#                              Read 4: ACTGATACTT
#                           The resulting SSCS is: ACTGATACNT    
#
# --Ncutoff NCUTOFF   Percentage of Ns allowed in a consensus sequence [0.3]
#                        Example (--ncutoff = 0.3):
#                           SSCS 1: ACGTGANCTAGTNCTNTACC
#                           SSCS 2: GATCTAGTNCATGACCGATA
#                        SSCS 2 passes the n filter (10%) with 1/20 = 5% Ns, 
#                        while SSCS 1 does not with 3/20 = 15% Ns.
#
# --infile INFILE     input BAM file
# --outfile OUTFILE   output BAM file
#
# python3 SSCS_maker.py --Ncutoff 0.3 --cutoff 0.7 --infile MEM-001_KRAS.bam --outfile MEM-001_KRAS.sscs.bam
#
###############################################################

import pysam # Need to install
import collections
import re
import array
from random import *
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
import statistics
import math
import time
import inspect
import os

from consensus_helper import *

###############################
##         Functions         ##
###############################

def mismatch_pos(cigar, mismatch_tag):
    '''(list, str) -> lst
    Return 0-based index of mismatch positions in sequence (including insertions,
    ignoring deletions). List of positions to check phred quality score.
    
    cigar: soft clips, hard clips, insertions, deletions
    mismatch_tag: position of mismatch
    
    E.g. mismatch_pos([(0, 19), (2, 1), (0, 79)], '19^A8G70')
    [27]
    First mismatch "G' at position 19 + 8
    
    Test cases:
    >>> mismatch_pos([(0, 98)], '31A0T22T41G') 
    [31, 32, 55, 97]
    >>> mismatch_pos([(0, 98)], '31A0T22T42')
    [31, 32, 55]
    >>> mismatch_pos([(0, 98)], 'G30A0T65')
    [0, 31, 32]
    
    In/del examples (1 = insertion, 2 = deletion in cigar):
    >>> mismatch_pos('19M1D26M1D28M25S', '19^A8G5T11^T8C0T12T5') # two deletions
    [27, 33, 53, 54, 67]
    >>> mismatch_pos([(0, 19), (2, 1), (0, 54), (4, 25)], '19^A8G5T10T8C0T12T5')
    [27, 33, 44, 53, 54, 67]
    >>> mismatch_pos([(0, 28), (1, 1), (0, 69)], '19T74G2') # one insertion 
    [19, 28, 94]
    >>> mismatch_pos([(0, 26), (1, 1), (0, 5), (2, 1), (0, 66)], '21C4G0G3^C14C6G5T33T4')
    [21, 26, 27, 45, 52, 58, 92]
    >>> mismatch_pos([(0, 74), (2, 2), (0, 3), (1, 2), (0, 19)], '70T1A1^GC22') # cases when insertions and deletions are in the same read
    [70, 72, 77]    
    >>> mismatch_pos([(0, 74), (2, 2), (0, 3), (2, 1), (1, 2), (0, 19)], '70T1A1^GC22')
    [70, 72, 77]
    >>> mismatch_pos('3M1I48M15D8M1D38M', '8G4C7G1G23C3^GAATTAAGAGAAGCA8^G38')
    >>> mismatch_pos([(0, 3), (1, 1), (0, 48), (2, 15), (0, 8), (2, 1), (0, 38)], '8G4C7G1G23C3^GAATTAAGAGAAGCA8^G38')
    [3, 8, 13, 21, 23, 47]
    
    Hard clip examples:
    >>> mismatch_pos([(5, 65), (0, 33)], '31A1')
    [31]
    >>> mismatch_pos([(0, 37), (5, 61)], '37')
    []
    '''
    mismatches = re.split('[^0-9, \^]+', mismatch_tag) #split by letters
    mis_pos = []
    index = 0    
    
    del_pos = 0 # keep track of pos # before deletion
    prev_del = False    
    
    for i in range(len(mismatches)-1):  
        # SNP in the first position         
        if mismatches[i] == '':
            mis_pos.append(0)
            continue
        # Ignore deletions, add pos num to subsequent mismatches
        if '^' in mismatches[i]:
            del_pos += int(mismatches[i][:-1])
            prev_del = True
            continue 
        
        mis_pos.append(int(mismatches[i]))
        
        # If prev pos contains deletion, add prev pos to current as we're ignoring deletions
        if prev_del:
            mis_pos[-1] += del_pos
            prev_del = False
            del_pos = 0

        if len(mismatches) >= 2:
            if i != 0 and len(mis_pos) > 1:
                # Need to add 1 to positions for correct indexing (except for first position)
                mis_pos[-1] += mis_pos[-2] + 1
    
    # Incorporate insertions 
    insert = 0
    for i in cigar:
        if i[0] == 1:
            # If there's multiple insertions
            for j in range(i[1]):
                mis_pos.append(insert)
                insert += 1
        elif i[0] == 0:
            # Keep track of positional num
            insert += i[1]
        else:
            pass
        
    mis_pos = list(set(mis_pos)) # Incase insertion and SNP share same position and there's repeat pos
    mis_pos.sort()
    
    return mis_pos


def query_seq_pos(cigar, readLength):
    '''(list of tuples, int) -> tuples
    Return tuple of seq position excluding soft clips and hard clips (0-based).
    
    0 = Match/mismatch
    1 = Insertion
    2 = Deletion
    4 = Soft clip
    5 = Hard clip
    
    >>> query_seq_pos([(4, 73), (0, 20), (4, 5)], 98)
    (73, 93)
    >>> query_seq_pos([(4, 6), (0, 92)], 98)
    (6, 98)
    >>> query_seq_pos([(0, 23), (4, 75)], 98)
    (0, 23)
    >>> query_seq_pos([(0, 37), (5, 61)], 98)
    (0, 37)
    '''
    start = 0 
    end = readLength
    
    if cigar[0][0] == 4 or cigar[0][0] == 5:
        start += cigar[0][1]
        
    if cigar[-1][0] == 4 or cigar[-1][0] == 5:
        end -= cigar[-1][1]

    return start, end


def consensus_maker(readList, readLength, cutoff):
    '''(list, int, int) -> str
    Return consensus sequence (without soft-clips) and quality score consensus.
    
    Majority rules concept where if no majority is reached above the cutoff, an 'N' is assigned to the position. 
    
    - Add N's for soft clipped regions so it aligns with full length sequences
    
    - At each position, add quality score to list corresponding to nucleotide. 
      Take max quality score of nucleotide with highest frequency
    '''
    nuc_lst = ['A', 'C', 'G', 'T', 'N']
    consensus_read = ''
    quality_consensus = array.array('B')
    
    mismatch_pos_lst = []
    
    for read in readList:
        mismatch_pos_lst.append(mismatch_pos(read.cigar, read.get_tag('MD')))

    for i in range(readLength):
        position_score = [0, 0 ,0, 0, 0] # A, C, G, T, N 
        #quality_score = [[], [], [], [], []] 
        #quality_score = [0, 0, 0, 0, 0]
        phred_fail = 0
        
        ### HOW MANY Ns ARE IN THE FINAL CONSENSUS AND HOW MANY TIE BREAKING EVENTS? ###
        
        for j in range(len(readList)):
            # === Find position of sequence without soft clips ===
            query_pos = query_seq_pos(readList[j].cigar, readLength)
            # if seq length < or > region of query seq, add 1 to N and set qual score 0
            if i < query_pos[0] or i >= query_pos[1]:
                position_score[4] += 1
                #quality_score[4].append(0)
                continue

            # === Phred filter mismatch positions ===     
            if i in mismatch_pos_lst[j]:
                if readList[j].query_alignment_qualities[i] < 30: # Phred cutoff of 30
                    phred_fail += 1
                    continue
            
            i = i - query_pos[0] # account for seq with soft/hard clips
            # index subtract clipped bps to iterate through sequence
            # (e.g. 2S96M -> indexes 0 and 1 are N,
            # but at index 2 actual position in query seq is 0)
            # IF you use query_sequence (which includes soft clips), then you don't need to offset the query pos => however, it will throw off consensus making between those with and without soft clips (e.g. readA: ACGTT, readB: TGACGTT (TG soft clips) but pos based comparison will show they're off) 

            # If pass filter, add 1 to nuc
            nuc = readList[j].query_alignment_sequence[i]
            nuc_index = nuc_lst.index(nuc)
        
            position_score[nuc_index] += 1  
            #quality_score[nuc_index] += 1
            #quality_score[nuc_index].append(readList[j].query_alignment_qualities[i])
            
            i = i + query_pos[0]                        


        try:
            # Find most common nuc #
            max_nuc_pos = [f for f, k in enumerate(position_score) if k == max(position_score)]
            # If there's more than one max, randomly select nuc
            max_nuc_pos = max_nuc_pos[randint(0, len(max_nuc_pos)-1)]
            
            # === Molecular Phred Quality Score ===
            # error = num variant bases (not most freq base)
            error_bases = position_score[:max_nuc_pos] + position_score[(max_nuc_pos +1):]
            # probability of observed bases
            P = sum(error_bases)/sum(position_score)
            if P == 0:
                Q = 62
            else:
                Q = round(-10 * math.log10(P))
                
            consensus_read += nuc_lst[max_nuc_pos]
            quality_consensus.append(Q)
        
            # frequency of nuc at position > cutoff 
            #if max(position_score)/(len(readList) - phred_fail) > cutoff:
                #consensus_read += nuc_lst[max_nuc_pos]
                #quality_consensus.append(max_qual)
            #else:
                #raise ValueError
                            
        except:
            # For cases when # matches fail phred > # reads
            consensus_read += 'N'
            quality_consensus.append(0)            

    return consensus_read, quality_consensus


def chr_arm_pos(chr_lst, chr_len):
    '''(list, list) -> list
    Return list of int indicating chromosomal arm positions given a list of chromosomes and their lengths.
    
    ChrM not divided by arms.
    
    Chrm arm positions are used to separate bam file reads into more manageable chunks, so dictionaries don't take up too much memory.
    
    Input: 
    - chr_lst 
    ['chrM', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    
    - chr_len
    [16571, 249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566]
        
    '''
    chr_arm_coor = collections.OrderedDict()
    
    if 'chrM' in chr_lst:
        chr_arm_coor['chrM'] = (0, chr_len[chr_lst.index('chrM')])
    
    filepath = os.path.abspath(inspect.getfile(inspect.currentframe())).rsplit('/', 1)[0]
    with open(filepath + '/cytoBand.txt') as f:
        next(f) # Skip header
        for line in f:
            chr_arm = line.split('\t')
            chr_key = '{}_{}'.format(chr_arm[0], chr_arm[3])
            start = int(chr_arm[1]) # python is 0-based (start is usually 1)
            end = int(chr_arm[2])
            chr_val = (start, end)

            chr_arm_coor[chr_key] = chr_val
    
    return chr_arm_coor


def reverse_seq(seq):
    '''(str) -> str
    Return reverse compliment of sequence (used for writing rev comp sequences to fastq files).
    
    >>> reverse_seq('TCAGCATAATT')
    'AATTATGCTGA'
    >>> reverse_seq('ACTGNN')
    'NNCAGT'
    '''
    rev_comp = ''
    nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    for base in seq:
        rev_comp = nuc[base] + rev_comp
        
    return rev_comp


def sscs_qname(tag, flag):
    '''(str, int) -> str
    Return new tag/queryname for consensus sequences (barcode_chr_start_chr_end).
    
    * Since multiple reads go into making a consensus, a new unique identifier is required to match up the read with its pair *
    
    Note: chr of start and end are used in case of translocations. In addition, coordinates are orientated as lower_chr -> higher_chr number for translocations.
    (e.g. TTCA_10_135461271_0_2364_fwd_R2 -> TTCA_0_2364_10_135461271_trans)
    
    Examples:
    (+)                  [Flag]
    ACGT_1_1_1_230_fwd_R1 [99] --->     ACGT_chr1_1_chr1_230_pos
    ACGT_1_230_1_1_rev_R2 [147]
    
    (-)
    ACGT_1_230_1_1_rev_R1 [83] --->     ACGT_chr1_1_chr1_230_neg
    ACGT_1_1_1_230_fwd_R2 [163]

    
    ATGT_1_249239818_1_10060_fwd_R1 [65] --->     ATGT_1_10060_1_249239818_pos
    ATGT_1_10060_1_249239818_fwd_R2 [129]
    
    
    Special cases (duplex and pair reads all in the same orientation):
    ['AGAG_3_178919046_8_75462483_rev_R1', 'AGAG_3_178919046_8_75462483_rev_R2', 'AGAG_8_75462483_3_178919046_rev_R2', 'AGAG_8_75462483_3_178919046_rev_R1']
    - Use pos/neg to differentiate between strand

    Test cases:
    >>> sscs_qname('CCCC_12_25398064_12_25398156_fwd_R1', 99)
    'CCCC_12_25398064_12_25398156_99_147
    >>> sscs_qname('CCCC_12_25398156_12_25398064_rev_R2', 147)
    'CCCC_12_25398064_12_25398156_99_147'
    >>> sscs_qname('CCCC_12_25398156_12_25398064_rev_R1', 83)
    'CCCC_12_25398064_12_25398156_83_163'
    >>> sscs_qname('CCCC_12_25398064_12_25398156_fwd_R2', 163)
    'CCCC_12_25398064_12_25398156_83_163'

    
    Translocation:
    >>> sscs_qname('TGGT_1_21842527_13_72956752_rev_R1', 113)
    'TGGT_1_21842527_13_72956752_113_177'
    >>> sscs_qname('TGGT_13_72956752_1_21842527_rev_R2', 177)
    'TGGT_1_21842527_13_72956752_113_177'
    >>> sscs_qname('TTCA_0_2364_10_135461271_fwd_R1', 65)
    'TTCA_0_2364_10_135461271_65_129'
    >>> sscs_qname('TTCA_10_135461271_0_2364_fwd_R2', 129)
    'TTCA_0_2364_10_135461271_65_129'
    '''
    flag_pairings = {99:147, 147:99, 83:163, 163:83, \
                     # mapped within insert size, but wrong orientation (++, --)
                     67:131, 131:67, 115:179, 179:115, \
                     ## === translocations ===
                     # mapped uniquely, but wrong insert size
                     81:161, 161:81, 97:145, 145:97, \
                     # wrong insert size and wrong orientation
                     65:129, 129:65, 113:177, 177:113
                     }

    start_chr = tag.split("_")[1]
    stop_chr = tag.split("_")[3]
    start_coor = tag.split("_")[2]
    stop_coor = tag.split("_")[4]

    # Order by chr coordinate
    if (int(start_coor) > int(stop_coor) and start_chr == stop_chr) or \
       (start_chr != stop_coor and int(start_chr) > int(stop_chr)):
        new_tag = tag.split("_")
        new_tag[1] = stop_chr
        new_tag[3] = start_chr
        new_tag[2] = stop_coor
        new_tag[4] = start_coor
        new_tag = "_".join(new_tag)[:-7]
        # Determine strand 
        if 'R1' in tag: # rev_R1
            new_tag = new_tag + '_neg'
        else: # rev_R2
            new_tag = new_tag + '_pos'
    else:
        # Use flags to determine strand for reads with the same coordinate as mate 
        if start_chr == stop_chr and start_coor == stop_coor:
            if flag in [99, 147]:
                new_tag = tag[:-7] + '_pos'
            elif flag in [83, 163]:
                new_tag = tag[:-7] + '_neg'
        elif 'R1' in tag: # fwd_R1
            new_tag = tag[:-7] + '_pos'
        else: # fwd_R2
            new_tag = tag[:-7] + '_neg'
            
    # Group pairs by strand unless translocation
    if flag < flag_pairings[flag]:
        new_tag = '{}_{}_{}'.format(new_tag, flag, flag_pairings[flag])
    else:
        new_tag = '{}_{}_{}'.format(new_tag, flag_pairings[flag], flag)
    
    return new_tag
    

def main():
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--cutoff", action = "store", dest="cutoff", help="nucleotide base % cutoff", required = True)
    parser.add_argument("--Ncutoff", action = "store", dest="Ncutoff", help="N % cutoff", required = True)
    parser.add_argument("--infile", action = "store", dest="infile", help="input BAM file", required = True)
    parser.add_argument("--outfile", action = "store", dest="outfile", help="output SSCS BAM file", required = True)
    args = parser.parse_args()
    
    
    start_time = time.time()
    
    # ===== Initialize input and output bam files =====
    bamfile = pysam.AlignmentFile(args.infile, "rb")
    SSCS_bam = pysam.AlignmentFile(args.outfile, "wb", template = bamfile)
    stats = open('{}.stats.txt'.format(args.outfile.split('.sscs')[0]), 'w')
    singleton_bam = pysam.AlignmentFile('{}.singleton.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)    
    doubleton_bam = pysam.AlignmentFile('{}.doubleton.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)
    
    # setup fastq files
    fastqFile1 = open('{}.sscs_R1.fastq.gz'.format(args.outfile.split('.sscs')[0]), 'w')
    fastqFile2 = open('{}.sscs_R2.fastq.gz'.format(args.outfile.split('.sscs')[0]), 'w')
    
    doubleton_fastqFile1 = open('{}.doubleton.sscs_R1.fastq.gz'.format(args.outfile.split('.sscs')[0]), 'w')
    doubleton_fastqFile2 = open('{}.doubleton.sscs_R2.fastq.gz'.format(args.outfile.split('.sscs')[0]), 'w')    
    
    time_tracker = open('{}.time_tracker.txt'.format(args.outfile.split('.sscs')[0]), 'w')
    
    bam_dict = collections.OrderedDict() # dict subclass that remembers order entries were added
    tag_dict = collections.defaultdict(int)
    paired_dict = collections.OrderedDict()
        
    tag_quality_dict = collections.defaultdict(list)
    quality_dict = collections.defaultdict(list)
    
    unmapped = 0
    unmapped_flag = 0
    bad_reads = 0
    counter = 0  
    doubletons = 0
    singletons = 0
    SSCS_reads = 0      
        
    chrm = [x['SN'] for x in bamfile.header['SQ']]
    chr_len = [x['LN'] for x in bamfile.header['SQ']]
    
    chr_arm_coor = chr_arm_pos(chrm, chr_len)
    #print(chr_arm_coor)
    
    #unmapped_key = []
    
    for x in chr_arm_coor.keys():
        bamLines = bamfile.fetch(reference = x.split('_')[0], start = chr_arm_coor[x][0], end =  chr_arm_coor[x][1]) # genomic start and end 0-based       
        # Create dictionary for each chrm
        for line in bamLines:
            counter += 1
            strand = 'fwd'
            if line.is_reverse:
                strand = 'rev'
                
            read = 'R1'
            if line.is_read2:
                read = 'R2'              
            
            # KEEP UNMAPPED READS THAT HAVE PAIR THATS MAPPED
            # if only map one side of it, chuck both
            ## chuck unmapped reads and their pairs!!! TODAY throw them for now
            # Later, if unmapped, use absolute sequence to build families 
            # if coordinates unavailable, compare sequence fix 
            
            #good_flags = [99, 147, 83, 163]
            
            # do we want to keep reads that have the pair mapped to reverse?????
            
            # Read 1 vs Read 2
            # Unique mol: read1 = 7975, read2 = 7934
            # original bam (test file): read1 = 229324, read2 = 229517
            # original bam (full): read1 = 35929483, read 2 = 35934622
            
            ## JUST FILTER OUT UNMAPPED as I won't have a R2 
            
            bad_flags = [73, 133, 89, 121, 165, 181, 101, 117, 153, 185, 69, 137, 77, 141]
            
            if line.is_unmapped:
                unmapped += 1
                continue
            elif line.is_secondary:
                bad_reads += 1
                continue
            elif line.is_supplementary:
                bad_reads += 1
                continue
            elif line.flag in bad_flags :
                unmapped_flag += 1
                continue        
        
            tag = '{}_{}_{}_{}_{}_{}_{}'.format(line.qname.split("|")[1], # mol barcode
                                          line.reference_id, # chr num
                                          line.reference_start, # start R1 (0-based)
                                          line.next_reference_id,
                                          line.next_reference_start, # start R2
                                          strand, # strand direction
                                          read # read num
                                          ) 
                        
            try:
                line.get_tag('MD')
            except:
                print(line)
                continue
                
            tag_dict[tag] += 1     
            #paired_dict[] 
            consensus_tag = sscs_qname(tag, int(line.flag))
            
            if tag not in bam_dict:
                bam_dict[tag] =[line]
                
                # Only add tag to paired_dict once (aka first time creating tag) 
                if consensus_tag not in paired_dict:
                    paired_dict[consensus_tag] = [tag]
                else:
                    paired_dict[consensus_tag].append(tag)                
        
            else:
                bam_dict[tag].append(line)              
        
        
        # ===== Create consenus seq for reads in each chrm arm and reset =====
        if bool(bam_dict):
            rand_key = choice(list(bam_dict.keys()))
            readLength = bam_dict[rand_key][0].infer_query_length()        
        
        ## WHY ARE MATES OF READS WITH TRANSLOCATIONS NOT GROUPED TOGETHER? 
        
        for readPair in paired_dict.keys():
            # Check pairing
            if len(paired_dict[readPair]) != 2:
                print('uh oh pairing problem!!!')
                print(readPair)
                print(paired_dict[readPair])
            else:
                for tag in paired_dict[readPair]:
                    # Check for singletons
                    if tag_dict[tag] < 2:
                        singletons += 1
                        singleton_bam.write(bam_dict[tag][0])
                    else:
                        SSCS = consensus_maker(bam_dict[tag], readLength, float(args.cutoff))
                        
                        quality_dict[tag] += [SSCS[1]]
                        tag_quality_dict[tag_dict[tag]] += [round(np.mean(SSCS[1]))]
                        
                        SSCS_read = create_aligned_segment(bam_dict[tag], SSCS[0], SSCS[1])
                        
                        #new_tag = sscs_qname(tag, SSCS_read.flag)
                        query_name = readPair + ':' + str(tag_dict[tag])   
                        SSCS_read.query_name = query_name
                        
                        # Use aligned sequence in case of soft clips
                        aligned_seq = SSCS_read.query_alignment_sequence
                        if aligned_seq.count('N')/len(aligned_seq) > float(args.Ncutoff):
                            print('uh oh too many Ns')
                            print(aligned_seq)
                            print(SSCS_read)
                            continue                        

                        if tag_dict[tag] == 2:
                            doubletons += 1
                            doubleton_bam.write(SSCS_read)
                        else:
                            SSCS_bam.write(SSCS_read)    
                            SSCS_reads += 1
                
                        # ===== write as fastq file =====                        
                        #if SSCS_read.is_reverse and SSCS_read.is_read1:
                        #fastq_seq = SSCS_read.query_sequence.decode("utf-8") 
                        fastq_seq = SSCS_read.query_sequence 
                        
                        
                        if 'rev' in tag:
                            fastq_seq = reverse_seq(fastq_seq)
                            fastq_qual = pysam.qualities_to_qualitystring(reversed(SSCS_read.query_qualities))
                        else:
                            fastq_qual = pysam.qualities_to_qualitystring(SSCS_read.query_qualities)
                        
                        if tag_dict[tag] == 2:
                            if 'R1' in tag:
                                doubleton_fastqFile1.write('@{}\n{}\n+\n{}\n'.format(SSCS_read.query_name, fastq_seq, fastq_qual))
                            else:
                                doubleton_fastqFile2.write('@{}\n{}\n+\n{}\n'.format(SSCS_read.query_name, fastq_seq, fastq_qual))                            
                        else:
                            if 'R1' in tag:
                                fastqFile1.write('@{}\n{}\n+\n{}\n'.format(SSCS_read.query_name, fastq_seq, fastq_qual))
                            else:
                                fastqFile2.write('@{}\n{}\n+\n{}\n'.format(SSCS_read.query_name, fastq_seq, fastq_qual))
                        
        import pickle
        read_dict_file = open(args.outfile.split('.sscs')[0] + '.read_dict.txt', 'ab+')
        pickle.dump(bam_dict, read_dict_file)
        read_dict_file.close()
        
        # reset dictionary            
        bam_dict = collections.OrderedDict() # dict subclass that remembers order entries were added        
        paired_dict = collections.OrderedDict()
        try:
            time_tracker.write(x + ': ')
            time_tracker.write(str((time.time() - start_time)/60) + '\n')
            
        except:
            continue
    
    # ===== write tag family size dictionary to file ===== 
    # (key = tags, value = int [number of reads in that family]) 
    tag_file = open(args.outfile.split('.sscs')[0] + '.read_families.txt', 'ab+')
    pickle.dump(tag_dict, tag_file)
    tag_file.close()
    
    summary_stats = '''Total reads: {} \n
Unmapped reads: {} \n
Unmapped flag reads: {} \n
Secondary/Supplementary reads: {} \n
SSCS reads: {} \n
singletons: {} \n
doubletons: {} \n
'''.format(counter, unmapped, unmapped_flag, bad_reads, SSCS_reads, singletons, doubletons)

    stats.write(summary_stats)
    
    time_tracker.close()
    stats.close()
    bamfile.close()
    SSCS_bam.close()
    doubleton_bam.close()
    fastqFile1.close()
    fastqFile2.close()
    doubleton_fastqFile1.close()
    doubleton_fastqFile2.close()    
    
    
    # ===== Create tag family size plot =====
    # Count number of families containing each number of read (e.g. Counter({1: 3737, 32: 660... -> 3737 families are singletons)
    fam_per_read_group = collections.Counter([i for i in tag_dict.values()])
    lst_fam_per_read = list(fam_per_read_group.items()) # conver to list
    
    total_reads = sum(tag_dict.values())
    # Multiply number of families by read num to get total number of reads in that read group, divide it by total reads to obtain read fraction
    read_fraction = [(i*j)/total_reads for i,j in lst_fam_per_read] 
    
    plt.bar(list(fam_per_read_group), read_fraction)
    #plt.locator_params(axis = 'x', nbins=lst_fam_per_read[-1][0]//5)
    plt.xlabel('Tag family size (# of reads per family)')
    plt.ylabel('Fraction of total reads')
    
    plt.savefig(args.outfile.split('.sscs')[0]+'_tag_fam_size.png')      
    
    
    ## READ IN BAM FILE LATER AND CREATE THESE PLOTS
    # ===== Create quality score plot =====
    
    # Should we take average quality score of every read??????
    import itertools
    
    qual_dist_lst = [list(i[0]) for i in quality_dict.values()]
    qual_dist = list(itertools.chain.from_iterable(qual_dist_lst))
    count_qual = collections.Counter(qual_dist).most_common()
    count_qual_sorted = sorted(count_qual, key=lambda tup: tup[0])
    x = [x for x,y in count_qual_sorted]
    y = [y for x,y in count_qual_sorted]    
    
    fig = plt.figure(figsize=(7.195, 3.841), dpi=100)  
    ax = fig.add_subplot(111)
    
    plt.plot(x, y, "-o", markersize=np.sqrt(5), linewidth = 0.5)
    plt.yscale('log')
    plt.xticks(np.arange(0, max(x)+5, 5.0), fontsize = 6)
    plt.yticks(fontsize = 6)    
    plt.xlabel('Phred Quality Score', fontsize = 8)
    plt.ylabel('Number of bases', fontsize = 8)
    
    for xy in zip(x, y):                                     
        ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data', fontsize = 4)     
    plt.grid()     
    
    plt.savefig(args.outfile.split('.sscs')[0]+'.phred_density.png', dpi = 1000) 
    
    
    # ===== Create quality score plot for each tag family =====
    
    # Figure out how to plot tag_quality_dict....
    
    # PLOT BASES -> density plot phred vs bases
    
    x = [[i]*len(tag_quality_dict[i]) for i in tag_quality_dict]
    y = list(tag_quality_dict.values())
    
    x = list(itertools.chain.from_iterable(x))
    y = list(itertools.chain.from_iterable(y))
    
    plt.figure(figsize=(7.195, 3.841), dpi=100)
    
    plt.plot(x, y, "o", markersize=np.sqrt(5))
    
    plt.xticks(np.arange(0, max(x)+1, 5.0), fontsize = 6)
    plt.yticks(np.arange(0, max(y)+5, 5.0), fontsize = 6)
    
    plt.ylabel('Average Molecular Q / Read', fontsize = 8)
    plt.xlabel('Tag family size', fontsize = 8)    
    
    #tick.label.set_fontsize(8)
    
    plt.grid() 
    
    plt.savefig(args.outfile.split('.sscs')[0]+'.tag_fam_quality.png', dpi = 1000)
    

###############################
##           Main            ##
###############################
if __name__ == "__main__": 
    import time
    start_time = time.time()
    main()  
    print((time.time() - start_time)/60)  


# use python, write to temp files, then merge together using samtools or bamtools and then delete temp files




#if line.is_unmapped:
    #unmapped += 1
    #mate_read = 'R1'
    #if read == 'R1':
        #mate_read = 'R2'
    #unmapped_key += [tag[:-2] + mate_read]
    
    ##try:
        ##unmapped_reads += [bamfile.mate(line)]
    ##except:
        ##print('uh oh, this read has no mate! {}'.format(line))
    ##print(line)
    ##print(tag)

    ##print(bamfile.mate(line))
    ##break 
    #continue
#elif tag in unmapped_key:
    #pair_unmapped += 1
    #continue
#elif line.reference_start == line.next_reference_start:
    #try:
        #if bamfile.mate(line).is_unmapped:
            #pair_unmapped += 1
            #continue
    #except:
        #print(line)
        #pair_unmapped += 1
        #continue
        
    ##print(line)
    ###print(bamfile.mate(line))
    ##if bamfile.mate(line).is_unmapped:
        ##continue
    ##else:
        ##print('uh oh')
