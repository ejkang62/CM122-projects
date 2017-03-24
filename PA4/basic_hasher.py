from collections import defaultdict, Counter
import cPickle as pickle
from os.path import join, exists, splitext
import time
import os
import sys
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
from BIOINFO_M260B.helpers import read_reads, read_reference, pretty_print_aligned_reads_with_ref


def hash_end(end, genome_ht):
    """
    Uses hashing to identify the set of locations spanned by
    a read.

    :param end: A single end of a read
    :param genome_ht: A hash of the genome with uniform key length
    :return:
    """
    key_length = len(genome_ht.keys()[0])
    end_pieces = [end[i * key_length: (i + 1) * key_length]
                  for i in range(len(end) / key_length)]

    hashed_read_locations = [genome_ht[read_piece]
                             for read_piece in end_pieces]
    start_positions = [[x - i * key_length for x in hashed_read_locations[i]]
                       for i in range(len(hashed_read_locations))]
    start_counter = Counter()

    for position_list in start_positions:
        start_counter.update(position_list)

    if not start_counter:
        return -1, 0
    else:
        best_alignment_location, best_alignment_count = \
            start_counter.most_common(1)[0]

    if best_alignment_count < 2:
        return -1, best_alignment_count
    else:
        return best_alignment_location, best_alignment_count


def hash_read(read, genome_ht):
    """
    Uses hashing to identify the set of locations spanned by
    a read.

    :param read: A single read
    :param genome_ht: A hash of the genome with uniform key length
    :return:
    """

    oriented_reads = [(read[0][::i], read[1][::j]) for i, j in ((1, -1), (-1, 1))]
    ## Either one end is forward and the other end is reversed, or vice versa.

    best_score = -1
    best_alignment_locations = (-1, -1)
    best_oriented_read = ('', '')
    for oriented_read in oriented_reads:
        hash_results = [hash_end(_, genome_ht) for _ in oriented_read]
        hash_locations = [_[0] for _ in hash_results]
        hash_score = sum([_[1] for _ in hash_results])
        if hash_score > best_score:
            best_alignment_locations = hash_locations
            best_oriented_read = oriented_read
    return best_oriented_read, best_alignment_locations


def make_genome_hash(reference, key_length):
    """

    :param reference: The reference as a string stored
    :param key_length: The length of keys to use.
    :return:
    """
    genome_hash = defaultdict(list)
    for i in range(len(reference) - key_length):
        ref_piece = reference[i: i + key_length]
        genome_hash[ref_piece].append(i)
    return genome_hash


def build_hash_and_pickle(ref_fn, key_length, force_rebuild=False):
    reference_hash_pkl_fn = '{}_hash_keylength_{}.pkl'.format(splitext(ref_fn)[0], key_length)
    if exists(reference_hash_pkl_fn) and not force_rebuild:
        ref_genome_hash = pickle.load(open(reference_hash_pkl_fn, 'rb'))
        if len(ref_genome_hash.keys()[0]) == key_length:
            return ref_genome_hash
        else:
            pass
    else:
        pass
    reference = read_reference(ref_fn)
    ref_genome_hash = make_genome_hash(reference, key_length)
    pickle.dump(ref_genome_hash, open(reference_hash_pkl_fn, 'wb'))
    return ref_genome_hash


def hashing_algorithm(paired_end_reads, genome_ht):
    """

    :param paired_end_reads:
    :param genome_ht:
    :return:
    """
    alignments = []
    genome_aligned_reads = []
    count = 0
    start = time.clock()

    for read in paired_end_reads:
        alignment, genome_aligned_read = hash_read(read, genome_ht)
        alignments.append(alignment)
        genome_aligned_reads.append(genome_aligned_read)
        count += 1
        if count % 100 == 0:
            time_passed = (time.clock()-start)/60
            print '{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed)
            remaining_time = time_passed/count*(len(paired_end_reads)-count)
            print 'Approximately {:.3} minutes remaining'.format(remaining_time)
    return alignments, genome_aligned_reads



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
                col.append('0')
            else:
                col.append('1')
        col.append(track)
        comparisons.append(col)
    return comparisons

if __name__ == "__main__":
    genome_name = 'hw4_W_0'
    input_folder = join('./', genome_name)
    reads_fn = join(input_folder, '{}_reads.txt'.format(input_folder))
    reads = read_reads(reads_fn)
    key_length = 7
    start = time.clock()

    ref_fn = join(input_folder, '{}_strains_1.txt'.format(input_folder))

    genome_hash_table = build_hash_and_pickle(reads_fn, key_length)
    ref = read_reference(ref_fn)
    genome_aligned_reads, alignments = hashing_algorithm(reads, genome_hash_table)

    strain_ref = read_strains(ref_fn)

    output_str = pretty_print_aligned_reads_with_ref(genome_aligned_reads, alignments, ref)
    print output_str[:5000]
