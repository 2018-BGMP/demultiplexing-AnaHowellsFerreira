#!/usr/bin/env python

# constants
#quality threshold cutoff
QUAL_THR = 30

##############
# functions  #
#############

#create function to convert Quality scores:
def convert_phred(letter):
    """Converts a single character into a phred score"""
    phred_scores = ord(letter) - 33
    return phred_scores

#apply convert_phred function to calculate mean of Quality score

def mean_qual(string):
    qual = 0
    length = len(string)
    for i in string:
        qual += convert_phred(i)
    return qual/length

# tests for convert_phred function
assert convert_phred("A") == 32, "Phred score incorrect"
assert convert_phred("@") == 31, "Phred score incorrect"
assert convert_phred("#") == 2, "Phred score incorrect"

#create function to reverse complement barcode:
def reverse_complement(bases):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return ''.join([complement[base] for base in bases[::-1]])

#tests for reverse complement function:
assert reverse_complement("ATTGC") == "GCAAT", "Incorrect complement base"
assert reverse_complement("A") != "A", "Incorrect complement base"


def write_fastq_record(fh, header, index, seq, qual):
    fh.write("{}{}\n{}\n+\n{}\n".format(header,index, seq, qual))


##############
# list
#############
# create a list with expected known barcodes

index_list = ["GTAGCGTA", "TATGGCAC",
              "CGATCGAT", "TGTTCCGT",
              "GATCAAGG", "GTCCTAAG",
              "AACAGCGA", "TCGACAAG",
              "TAGCCATG", "TCTTCGAC",
              "CGGTAATC", "ATCATGCG",
              "CTCTGGAT", "ATCGTGGT",
              "TACCGGAT", "TCGAGAGT",
              "CTAGCTCA", "TCGGATTC",
              "CACTTCAC", "GATCTTGC",
              "GCTACTCT", "AGAGTCCA",
              "ACGATCAG","AGGATAGC"]




################
# python script  #########################
################

#define forward files:
read1 = "/projects/bgmp/ahowells/demultiplexing/demultiplexing_part2/1294_S1_L008_R1_001.fastq.gz"
index1 ="/projects/bgmp/ahowells/demultiplexing/demultiplexing_part2/1294_S1_L008_R2_001.fastq.gz"

#define reverse files:
read2 = "/projects/bgmp/ahowells/demultiplexing/demultiplexing_part2/1294_S1_L008_R4_001.fastq.gz"
index2 = "/projects/bgmp/ahowells/demultiplexing/demultiplexing_part2/1294_S1_L008_R3_001.fastq.gz"


#bad index reads that will be discarded (didn't pass quality tests)
bad_fw = "bad_fw.fastq"
bad_rv = "bad_rv.fastq"

#key = string represeting a barcode
#values = file handle for that barcode
bar_R1_fh_dict= {}
bar_R2_fh_dict= {}

import gzip

#Open all files at once:
with    gzip.open(index1, "rt") as i1_fh, \
        gzip.open(index2, "rt") as i2_fh, \
        gzip.open(read1, "rt") as r1_fh, \
        gzip.open(read2, "rt") as r2_fh, \
        open(bad_fw, "w") as bad_fw_fh, \
        open(bad_rv, "w") as bad_rv_fh:

    #open 48 fastq files to be written to
    for index in index_list:
        index_R1_fh = open("{}_R1.fastq".format(index), "w")
        index_R2_fh = open("{}_R2.fastq".format(index), "w")
        bar_R1_fh_dict[index] = index_R1_fh
        bar_R2_fh_dict[index] = index_R2_fh


    #start counters for desire output tracking of index hoping
    contain_N = 0 #containing N values
    low_qual_index= 0  #low phred Score < 30 for any bp (not sure how to do this!)
    mismatched_list = 0
    index_hoping = 0  #reads that do not match each other (read2 and reverse_complement read3)
    correct_index = 0   #reads that pass quality score and match each other to library

    #print(index_list)
    while True:
        #read index1 record
        header_i1 = i1_fh.readline().strip()
        if not header_i1:
            break   #tells to end when reach last line in file
        seq_i1 = i1_fh.readline().strip()
        plus_i1 = i1_fh.readline().strip()
        qual_i1 = i1_fh.readline().strip()
        #read index2 record...
        header_i2 = i2_fh.readline().strip()
        if not header_i2:
            break
        seq_i2 = i2_fh.readline().strip()
        plus_i2 = i2_fh.readline().strip()
        qual_i2 = i2_fh.readline().strip()
        #read read1 record
        header_r1 = r1_fh.readline().strip()
        if not header_r1:
            break
        seq_r1 = r1_fh.readline().strip()
        plus_r1 = r1_fh.readline().strip()
        qual_r1 = r1_fh.readline().strip()
        #read read2 record...
        header_r2 = r2_fh.readline().strip()
        if not header_r2:
            break
        seq_r2 = r2_fh.readline().strip()
        plus_r2 = r2_fh.readline().strip()
        qual_r2 = r2_fh.readline().strip()

        #apply reverse complement function to read3 file (referred as i2_fh)
        rev_seq_i2 = reverse_complement(seq_i2)

        #Remove indexes with N on the sequence and direct them to output bad files
        #print(header_i1)
        #print(seq_i1)
        if "N" in seq_i1 or "N" in seq_i2:
            #print("DEBUG AAA {} and {}".format(seq_i1, rev_seq_i2))
            write_fastq_record(bad_fw_fh, header_r1, seq_i1, seq_r1, qual_r1)
            write_fastq_record(bad_rv_fh, header_r2, rev_seq_i2, seq_r2, qual_r2)
            contain_N += 1
            continue

        #set a qual_threshhold to evaluate my poor quality indexes

        if mean_qual(qual_i1) <= QUAL_THR or mean_qual(qual_i2) <= QUAL_THR:
            #print("DEBUG BBB {} and {}".format(seq_i1, rev_seq_i2))
            #print("made it to low qual index 1")
            #print("made it to low qual index 2")
            write_fastq_record(bad_fw_fh, header_r1, seq_i1, seq_r1, qual_r1)
            write_fastq_record(bad_rv_fh, header_r2, rev_seq_i2, seq_r2, qual_r2)
            low_qual_index += 1
            continue


        #write to bad files those indexes that do not match expected/index in dictionary
        if seq_i1 not in bar_R1_fh_dict or rev_seq_i2 not in bar_R2_fh_dict:
            #print("DEBUG CCC {} and {}".format(seq_i1, rev_seq_i2))
            write_fastq_record(bad_fw_fh, header_r1, seq_i1, seq_r1, qual_r1)
            write_fastq_record(bad_rv_fh, header_r2, rev_seq_i2, seq_r2, qual_r2)
            mismatched_list += 1
            continue

        #write to bad files those indexes that do not match each other(hoping)
        if seq_i1 != rev_seq_i2:
            #print("DEBUG DDD {} and {}".format(seq_i1, rev_seq_i2))
            write_fastq_record(bad_fw_fh, header_r1, seq_i1, seq_r1, qual_r1)
            write_fastq_record(bad_rv_fh, header_r2, rev_seq_i2, seq_r2, qual_r2)
            index_hoping += 1
            continue

        #write to 48 files created
        if seq_i1 == rev_seq_i2:
            #print("DEBUG EEE {} and {}".format(seq_i1, rev_seq_i2))
            write_fastq_record(bar_R1_fh_dict[seq_i1], header_r1, seq_i1, seq_r1, qual_r1)
            write_fastq_record(bar_R2_fh_dict[rev_seq_i2], header_r2, rev_seq_i2, seq_r2, qual_r2)
            correct_index += 1

print("Total Index with N:", contain_N, "\t" "Total Low Quality Index:", low_qual_index, "\t" \
"Mismatched to List:", mismatched_list, "\t" "Total Index Hoping:", index_hoping, "\t" \
"Total Correct Index:",correct_index)
