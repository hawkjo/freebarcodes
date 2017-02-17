import numpy as np


bases = 'ACGT'


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
        if prev_seq[-2] == prev_seq[-1]:
            bad_bases += prev_seq[-1]
        if prev_cnt[0] + prev_cnt[3] == AT_max:
            bad_bases += 'AT'
        elif prev_cnt[1] + prev_cnt[2] == GC_max:
            bad_bases += 'CG'

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

    return iterate_seqs()


def dna2num(s):
    """
    Convert dna to number where dna is considered base 4 with '0123' = 'ACGT'.

        s :str:     Given dna string
    """
    return sum(bases.index(c) << 2*i for i, c in enumerate(s[::-1]))


def num2dna(n, dnalen):
    """
    Convert number to dna of given length where dna is considered base 4 with '0123' = 'ACGT'

        n :int:         Numerical representation of dna string
        dnalen :int:    Length of dna string
    """
    return ''.join(bases[(n & (3 << i)) >> i] for i in xrange(2*dnalen-2, -1, -2))

