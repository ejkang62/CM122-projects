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

def faster_algorithm(paired_end_reads, ref):
    """
    :param paired_end_reads: Paired-end reads generated from read_reads
    :param ref: A reference genome generated from read_reference
    :return: 2 lists:
                1) a list of alignment locations for each read (all_alignment_locations).
                    The list gives the starting position of the minimum-mismatch alignment of both reads.
                2) a list of the paired-end reads set so that both reads are in their optimal orientation
                   with respect to the reference genome.
    """
    all_read_alignment_locations = []
    output_read_pairs = []
    count = 0
    part_len = 16
    start = time.clock()
    ref_table = hash_ref(reference, part_len)
    for read_pair in paired_end_reads:
        count += 1
        read_alignment_locations = []
        output_read_pair = []
        if count % 10 == 0:
            time_passed = (time.clock()-start)/60
            print '{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed)
            remaining_time = time_passed/count*(len(paired_end_reads)-count)
            print 'Approximately {:.3} minutes remaining'.format(remaining_time)
        for read in read_pair:
            max_mismatches = 2
            min_mismatch_location = -1

            # include another if-else part that only checks reverse if original didn't match?
            # split read into 3 parts
            for part in range(0,part_len*3,part_len):
                section = read[part:part+part_len]
                if section in ref_table:
                    for i in ref_table[section]:
                        mismatches = [1 if read[j] != ref[i - part + j] else 0 for j in range(len(read))]
                        n_mismatches = sum(mismatches)
                        if n_mismatches < max_mismatches:
                            min_mismatch_location = i - part

            rev_read = read[::-1]

            for part in range(0,part_len*3,part_len):
                section = rev_read[part:part+part_len]
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


if __name__ == "__main__":
    data_folder = 'practice_W_1'
    input_folder = join('./', data_folder)
    f_base = '{}_chr_1'.format(data_folder)
    reads_fn = join(input_folder, 'reads_{}.txt'.format(f_base))
    start = time.clock()
    input_reads = read_reads(reads_fn)

    reference_fn = join(input_folder, 'ref_{}.txt'.format(f_base))
    reference = read_reference(reference_fn)
    alignments, reads = faster_algorithm(input_reads, reference)
    print alignments
    print reads
    output_str = pretty_print_aligned_reads_with_ref(reads, alignments, reference)
    output_fn = join(input_folder, 'aligned_{}.txt'.format(f_base))
    with(open(output_fn, 'w')) as output_file:
        output_file.write(output_str)
