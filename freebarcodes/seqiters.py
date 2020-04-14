import itertools
import random
import numpy as np
from seqtools import bases, dna_rev_comp, dna2num

def possible_barcode_iterator(k, AT_max, GC_max):
    """
    Iterates all k-mers such that:
        AT count <= AT_max
        GC count <= GC_max
        No triplet homopolymers

        k :int:         barcode k-mer length
        AT_max :int:
        GC_max :int:
    """
    bases_and_incs = [('A', np.array([1, 0, 0, 0])),
                      ('C', np.array([0, 1, 0, 0])),
                      ('G', np.array([0, 0, 1, 0])),
                      ('T', np.array([0, 0, 0, 1]))]

    def recursive_extension(prev_seq, prev_cnt):  # No defaults
        bad_bases = ''
        if prev_seq[-2] == prev_seq[-1]:  # Don't allow triplets
            bad_bases += prev_seq[-1]
            if prev_seq[-1] == 'G':  # Illumina has higher errors with GGC
                bad_bases += 'C'
        if prev_cnt[0] + prev_cnt[3] == AT_max:  # Enforce AT/GC content within bounds
            bad_bases += 'AT'
        elif prev_cnt[1] + prev_cnt[2] == GC_max:
            bad_bases += 'CG'
        for i in range(len(prev_seq)-4):  # Don't allow rev-comp seqs of 3+ bp
            if dna_rev_comp(prev_seq[i+1:i+3]) == prev_seq[-2:]:
                bad_bases += dna_rev_comp(prev_seq[i])


        if len(prev_seq) + 1 == k:
            for base in bases:
                if base in bad_bases:
                    continue
                else:
                    yield prev_seq + base
        else:
            for base, inc in bases_and_incs:
                if base in bad_bases:
                    continue
                for seq in recursive_extension(prev_seq + base, prev_cnt + inc):
                    yield seq

    def iterate_seqs():
        for b1, c1 in bases_and_incs:
            for b2, c2 in bases_and_incs:
                for seq in recursive_extension(b1 + b2, c1 + c2):
                    yield seq

    return iterate_seqs


def idx_possible_barcode_iterator(k, AT_max, GC_max):
    def iterate_seqs():
        for seq in possible_barcode_iterator(k, AT_max, GC_max)():
            yield dna2num(seq)
    return iterate_seqs


def possible_barcode_foursomes_iterator(k, AT_max, GC_max):
    """
    Iterates all sets of four k-mers such that:
        ACGT each appear in one of the four k-mers at each position
        AT count <= AT_max
        GC count <= GC_max
        No quadruplet homopolymers

        k :int:         barcode k-mer length
        AT_max :int:
        GC_max :int:
    """
    bases_and_incs = [('A', np.array([1, 0, 0, 0])),
                      ('C', np.array([0, 1, 0, 0])),
                      ('G', np.array([0, 0, 1, 0])),
                      ('T', np.array([0, 0, 0, 1]))]

    all_bases_and_incs = list(itertools.permutations(bases_and_incs))

    def recursive_extension(prev_seqs, prev_cnts):  # No defaults
        all_bad_bases = ['' for _ in range(4)]
        for i, (prev_seq, prev_cnt) in enumerate(zip(prev_seqs, prev_cnts)):
            if prev_seq[-2] == prev_seq[-1]:  # Don't allow triplets
                all_bad_bases[i] += prev_seq[-1]
                if prev_seq[-1] == 'G':  # Illumina has higher errors with GGC
                    all_bad_bases[i] += 'C'
            if prev_cnt[0] + prev_cnt[3] == AT_max:  # Enforce AT/GC content within bounds
                all_bad_bases[i] += 'AT'
            elif prev_cnt[1] + prev_cnt[2] == GC_max:
                all_bad_bases[i] += 'CG'
            for i in range(len(prev_seq)-4):  # Don't allow rev-comp seqs of 3+ bp
                if dna_rev_comp(prev_seq[i+1:i+3]) == prev_seq[-2:]:
                    all_bad_bases[i] += dna_rev_comp(prev_seq[i])

        good_bases_and_incs = [
            b_and_i for b_and_i in all_bases_and_incs
            if all(base not in bad_bases
                   for (base, _), bad_bases in zip(b_and_i, all_bad_bases))
        ]
        #random.shuffle(good_bases_and_incs)


        if len(prev_seqs[0]) + 1 == k:
            for b_and_i in good_bases_and_incs:
                yield [prev_seq + base for prev_seq, (base, inc) in zip(prev_seqs, b_and_i)]
        else:
            for b_and_i in good_bases_and_incs:
                new_prev_seqs = [prev_seq + base for prev_seq, (base, inc) in zip(prev_seqs, b_and_i)]
                new_prev_incs = [prev_cnt + inc for prev_cnt, (base, inc) in zip(prev_cnts, b_and_i)]
                for seqs in recursive_extension(new_prev_seqs, new_prev_incs):
                    yield seqs

    def iterate_seq_sets():
        # Note the first bases are always ACGT in that order to remove redundancy
        for b_and_i_2 in all_bases_and_incs:
            for b_and_i_3 in all_bases_and_incs:
                prev_seqs = [b1 + b2 + b3 for (b1, c1), (b2, c2), (b3, c3)
                             in zip(bases_and_incs, b_and_i_2, b_and_i_3)]
                prev_cnts = [c1 + c2 + c3 for (b1, c1), (b2, c2), (b3, c3)
                             in zip(bases_and_incs, b_and_i_2, b_and_i_3)]
                for seqs in recursive_extension(prev_seqs, prev_cnts):
                    yield seqs

    return iterate_seq_sets


def idx_possible_barcode_foursomes_iterator(k, AT_max, GC_max):
    def iterate_seq_sets():
        for seqs in possible_barcode_foursomes_iterator(k, AT_max, GC_max)():
            yield [dna2num(seq) for seq in seqs]
    return iterate_seq_sets


