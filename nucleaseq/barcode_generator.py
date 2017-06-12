import sys
import os
import numpy as np
import SeqlevSphere
import seqlev_dist
import string


bases = 'ACGT'

dna_complements = string.maketrans('acgtnACGTN', 'tgcanTGCAN')
def dna_rev_comp(dna_string):
    return dna_string.translate(dna_complements)[::-1]


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


class SeqlevBarcodeGenerator(object):
    def __init__(self, bc_len, min_dist, seq_idx_iter_func=None):
        self.bc_len = bc_len
        self.min_dist = min_dist
        self._codewords = set()
        self.barcodes = set()
        self.manual_codewords = set()
        self.reserved_words = np.zeros((4**self.bc_len, ), dtype=np.uint8)
        if seq_idx_iter_func is not None:
            self.seq_idx_iter_func = seq_idx_iter_func
        else:
            self.seq_idx_iter_func = lambda : xrange(4**self.bc_len)

    def _add_codeword(self, idx):
        assert isinstance(idx, int), '{} is not a valid codeword. Must be int'.format(idx)
        self._codewords.add(idx)
        word = num2dna(idx, self.bc_len)
        for seq_idx in SeqlevSphere.SeqlevSphere(word, self.min_dist-1).parallel_num_iterator():
            self.reserved_words[seq_idx] = 1

    def _add_barcode(self, seq_idx):
        assert not self._idx_is_reserved(seq_idx), seq_idx
        self.barcodes.add(seq_idx)
        self._codewords.add(seq_idx)
        self._add_codeword(seq_idx)

    def add_dnastr_nonbarcode_codeword(self, dnastring):
        seq_idx = dna2num(dnastring)
        self.add_idx_nonbarcode_codeword(seq_idx)
        
    def add_idx_nonbarcode_codeword(self, seq_idx):
        self.manual_codewords.add(seq_idx)
        self._add_codeword(seq_idx)

    def _idx_is_reserved(self, idx):
        return self.reserved_words[idx]

    def Conway_closure(self):
        for seq_idx in self.seq_idx_iter_func():
            if not self._idx_is_reserved(seq_idx):
                self._add_barcode(seq_idx)
                print len(self.barcodes),

    @property
    def dna_barcodes(self):
        return (num2dna(seq_idx, self.bc_len) for seq_idx in self.barcodes)

    def manual_barcodes_test(self):
        bc_list = list(self.barcodes)
        for i in range(len(self.barcodes)):
            bc1 = num2dna(bc_list[i], self.bc_len)
            for j in range(i+1, len(self.barcodes)):
                bc2 = num2dna(bc_list[j], self.bc_len)
                dist = seqlev_dist.seqlev_dist(bc1, bc2)
                if dist < self.min_dist:
                    print '!'*10 + ' FAIL ' + '!'*10
                    print 'Distance {} between {} and {}.'.format(dist, bc1, bc2)
                    return
        print 'Barcodes Pass Manual Check'

def write_barcodes(bc_len, min_dist, dpath):
    import time
    start_time = time.time()
    fpath = os.path.join(dpath, 'barcodes{}-{}.txt'.format(bc_len, min_dist))
    GC_max = min(range(bc_len), key=lambda x: abs(float(x)/bc_len-0.6))
    print 'AT/GC max:', GC_max
    bc_iter = idx_possible_barcode_iterator(bc_len, GC_max, GC_max)
    sbg = SeqlevBarcodeGenerator(bc_len, min_dist, bc_iter)
    sbg.Conway_closure()
    with open(fpath, 'w') as out:
        out.write('\n'.join(sorted(sbg.dna_barcodes)))
    comp_time = time.time() - start_time
    print
    print 'Barcode generation time:', comp_time
    stats_fpath = os.path.join(dpath, 'barcodes{}-{}_stats.txt'.format(bc_len, min_dist))
    with open(stats_fpath, 'w') as out:
        out.write('Barcode length:\t{:d}\n'.format(bc_len))
        out.write('Min-dist:\t{:d}\n'.format(min_dist))
        out.write('Barcode generation time:\t{:.2f} seconds\n'.format(comp_time))
