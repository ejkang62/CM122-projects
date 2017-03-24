import sys
import os

from collections import defaultdict
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
import numpy as np
from os.path import join
import time
from BIOINFO_M260B.helpers import read_reads, read_reference, pretty_print_aligned_reads_with_ref

def hash_ref(ref_reads, partial_len):
    reads_table = defaultdict(list)
    for i in range(0,len(ref_reads)-partial_len):
        seq = ref_reads[i:i+partial_len]
        reads_table[seq].append(i)
    return reads_table

def faster_algorithm(all_reads, ref):

    all_read_alignment_locations = []
    output_read_pairs = []
    count = 0
    part_len = 10
    # start = time.clock()
    ref_table = hash_ref(ref, part_len)
    for read in all_reads:
        count += 1
        read_alignment_locations = []
        output_read_pair = []
        # if count % 10 == 0:
        #     time_passed = (time.clock()-start)/60
        #     print '{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed)
        #     remaining_time = time_passed/count*(len(all_reads)-count)
        #     print 'Approximately {:.3} minutes remaining'.format(remaining_time)
        #for read in read_pair:
        max_mismatches = 2
        min_mismatch_location = -1

        # split read into parts
        for part in range(0,part_len*5,part_len):
            section = read[part:part+part_len]
            if section in ref_table:
                for i in ref_table[section]:
                    mismatches = [1 if read[j] != ref[i - part + j] else 0 for j in range(len(read))]
                    n_mismatches = sum(mismatches)
                    if n_mismatches < max_mismatches:
                        min_mismatch_location = i - part

        rev_read = read[::-1]
        # print 'read is'
        # print read
        # print min_mismatch_location

        for part in range(0, part_len * 5, part_len):
            section = rev_read[part:part + part_len]
            if section in ref_table:
                for i in ref_table[section]:
                    mismatches = [1 if rev_read[j] != ref[i - part + j] else 0 for j in range(len(rev_read))]
                    n_mismatches = sum(mismatches)
                    if n_mismatches < max_mismatches:
                        min_mismatch_location = i - part
                        read = rev_read

        read_alignment_locations.append(min_mismatch_location)
        output_read_pair.append(read)
        all_read_alignment_locations.append(read_alignment_locations)
        output_read_pairs.append(output_read_pair)

    return all_read_alignment_locations, output_read_pairs

# read in strains into matrix
def read_strains(read_fn):
    f = open(read_fn, 'r')
    first_line = True
    all_reads = []
    for line in f:
        if first_line:
            first_line = False
            continue
        line = line.strip()
        for read in line.split('\n'):
            all_reads.append(read)
    return all_reads

# compare other 3 strains to a reference strain
def compare_strains(strains):
    comparisons = []
    for track in range(0,len(strains[0])):
        col = []
        col.append(1)
        if strains[0][track] == strains[1][track] and strains[1][track] == strains[2][track] and strains[2][track] == strains[3][track]:
            continue
        for i in range(1,4):
            if strains[0][track] != strains[i][track]:
                col.append(0)
            else:
                col.append(1)
        comparisons.append(col)
    return comparisons

def snp_frequency(read_align, reads, strain):
    snps_file = open('./hw4_W_0/hw4_W_0_snps.txt')
    snps = []
    for line in snps_file:
        line = line.strip()
        snps.append(int(line))

    nearby_snps = defaultdict(list)

    for index, read in zip(read_align, reads):
        for ind, snp in enumerate(snps):
            if snp >= index[0] and snp < index[0]+50:
                snpPosition = snp - index[0]
                nearby_snps[ind].append(1 if read[0][snpPosition] == strain[snp] else 0)

    # print reads
    return nearby_snps

# from python notebook
def get_freqs_from_strains_and_counts(input_snp_counts, input_strains):
    #     print input_snp_count
    snp_freqs = []
    for k in range(len(input_snp_counts)):
        raw_snps = input_snp_counts[k]
        snp_freq = float(sum(raw_snps)) / len(raw_snps)
        snp_freqs.append(snp_freq)

    snp_freqs = np.array(snp_freqs)

    strain_matrix = np.array(input_strains)
    print strain_matrix

    strain_freqs = np.linalg.lstsq(strain_matrix.T, snp_freqs)
    for freq, strain in zip(strain_freqs[0], input_strains):
        print freq, strain
    return

if __name__ == "__main__":
    # data_folder = 'hw4_W_0'
    # input_folder = join('./', data_folder)
    # reads_fn = join(input_folder, '{}_strains.txt'.format(input_folder))
    #
    # all_strains = read_strains(reads_fn)
    # snps_matrix = compare_strains(all_strains)
    # snps_matrix = np.transpose(snps_matrix)

    # [['1' '1' '1' '1' '1' '1' '1']
    #  ['0' '0' '0' '1' '1' '0' '0']
    #  ['1' '1' '0' '0' '1' '0' '1']
    #  ['0' '0' '1' '0' '0' '1' '1']

    #  [[1 1 1 1 1 1 1]
    #   [0 0 0 1 1 0 0]
    #   [1 1 0 0 1 0 1]
    #   [0 0 1 0 0 1 1]]

    #  ['1320' '1451' '4575' '4646' '7358' '8288' '8347']]


    data_folder = 'hw4_W_0'
    input_folder = join('./', data_folder)
    reads_fn =  join(input_folder,'{}_reads.txt'.format(input_folder))
    start = time.clock()
    input_reads = read_reads(reads_fn)

    reference_fn = join(input_folder, '{}_strains.txt'.format(input_folder))
    reference = read_reference(reference_fn)
    alignments, reads = faster_algorithm(input_reads, reference[0])

    strain_ref = read_strains(reference_fn)
    snps = snp_frequency(alignments, reads, strain_ref[0])

    # snp_freq = defaultdict(list)
    # for index, values in snp_freqs.iteritems():
    #     snp_freq[index].append(float(sum(values))/len(values))

    # defaultdict( < type 'list' >,
    #           {0: [0.6601941747572816],
    #            1: [0.6235294117647059],
    #            2: [0.4074074074074074],
    #            3: [0.8155339805825242],
    #            4: [0.9393939393939394],
    #            5: [0.47619047619047616],
    #            6: [0.7471264367816092]})

    # for getting snps from strains
    all_strains = read_strains(reference_fn)
    snps_matrix = compare_strains(all_strains)
    snps_matrix = np.transpose(snps_matrix)

    get_freqs_from_strains_and_counts(snps, snps_matrix)