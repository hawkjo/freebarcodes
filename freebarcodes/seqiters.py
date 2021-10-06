import itertools
import random
import numpy as np
from .seqtools import bases, dna_rev_comp, dna2num

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


def seq_iterator_avoiding_prev_bases(k, AT_max, GC_max, first_base, seqs_so_far=[]):
    """
    Iterates all k-mers such that:
        The first base is as specified
        k-mer[i] != seq[i] for any seq in seqs_so_far
        AT count <= AT_max
        GC count <= GC_max
        No triplet homopolymers

        k :int:         barcode k-mer length
        AT_max :int:
        GC_max :int:
        first_base :str:
        seqs_so_far :list:
    """
    bases_and_incs = [('A', np.array([1, 0, 0, 0])),
                      ('C', np.array([0, 1, 0, 0])),
                      ('G', np.array([0, 0, 1, 0])),
                      ('T', np.array([0, 0, 0, 1]))]

    def recursive_extension(prev_seq, prev_cnt):  # No defaults
        bad_bases = ''.join([other_seq[len(prev_seq)] for other_seq in seqs_so_far])
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
        if len(prev_seq) >= 4:
            for base in [b for b in bases if b not in bad_bases]:
                rc_seq = dna_rev_comp(prev_seq[-4:] + base)
                if any(rc_seq in ssf for ssf in seqs_so_far):
                    bad_bases += base


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
        b1 = first_base
        c1 = next(inc for base, inc in bases_and_incs if base == b1)
        bad_bases = ''.join([other_seq[1] for other_seq in seqs_so_far])
        for b2, c2 in bases_and_incs:
            if b2 in bad_bases:
                continue
            for seq in recursive_extension(b1 + b2, c1 + c2):
                yield seq

    return iterate_seqs


def idx_iterator(seq_iterator_generator, *args, **kw_args):
    def iterate_idxs():
        for seq in seq_iterator_generator(*args, **kw_args)():
            yield dna2num(seq)
    return iterate_idxs


def idx_possible_barcode_iterator(k, AT_max, GC_max):
    return idx_iterator(possible_barcode_iterator, k, AT_max, GC_max)


def idx_seq_iterator_avoiding_prev_bases(*args, **kw_args):
    return idx_iterator(seq_iterator_avoiding_prev_bases, *args, **kw_args)
