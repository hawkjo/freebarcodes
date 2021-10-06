import os
import numpy as np
import psutil
import logging
from . import FreeDivSphere
from . import editmeasures
from .seqtools import bases, dna2num, num2dna, frac_to_int_max_GC
from .seqiters import idx_possible_barcode_iterator, idx_seq_iterator_avoiding_prev_bases

log = logging.getLogger(__name__)


class FreeDivBarcodeGenerator(object):
    """
    A class to generate a set of barcodes via the FreeDiv semimetric.
    """
    # Algorithm explanation:
    #   
    # In a metric space, one can reserve spheres around each new codeword of radius 2*max_err, and
    # any remaining word not reserved by any codeword is a valid new codeword. The FreeDiv
    # semimetric does not give that guarantee.  However, the above algorithm is still useful in
    # that any word so reseved is guaranteed to not be a valid codeword. Hence, we perform the
    # following steps:
    #
    #   1. Alphabetically iterate (for the Conway closure) through the words of acceptable
    #       composition until finding a word still marked with int(0). For that word, check its
    #       decode sphere. If all words in the decode sphere are marked < int(2), add codeword.
    #   2. After adding a codeword, reseve two surrounding spheres of different radii:
    #       -A sphere of radius max_err, the decode sphere, marked with int(2)
    #       -A sphere of radius 2*max_err, the approx encode sphere, marked with int(1), careful to
    #       not overwrite any 2's
    #   3. Repeat until space exhausted.

    def __init__(self, bc_len, max_err, seq_idx_iter_func=None):
        self.bc_len = bc_len
        self.max_err = max_err
        self._codewords = set()
        self.barcodes = set()
        self.manual_codewords = set()
        self.dna_barcode_4sets = []
        needed_bytes = 4**self.bc_len
        available_bytes = psutil.virtual_memory().available
        if needed_bytes > available_bytes:
            raise RuntimeError('Not enough memory. {:,d} bytes needed, {:,d} bytes available'.format(
                needed_bytes,
                available_bytes
            ))

        self.reserved_words = np.zeros((needed_bytes, ), dtype=np.uint8)
        if seq_idx_iter_func is not None:
            self.seq_idx_iter_func = seq_idx_iter_func
        else:
            self.seq_idx_iter_func = lambda : range(4**self.bc_len)

    def _add_codeword(self, cw_idx):
        assert isinstance(cw_idx, int), '{} is not a valid codeword. Must be int'.format(cw_idx)
        self._codewords.add(cw_idx)
        for seq_idx in self.iterate_approx_encode_sphere(cw_idx):
            if self.reserved_words[seq_idx] == 0:
                # Important to not overwrite previous 2's
                self.reserved_words[seq_idx] = 1
        for seq_idx in self.iterate_decode_sphere(cw_idx):
            self.reserved_words[seq_idx] = 2

    def iterate_decode_sphere(self, center_idx):
        word = num2dna(center_idx, self.bc_len)
        for seq_idx in FreeDivSphere.FreeDivSphere(word, self.max_err).parallel_num_iterator():
            yield seq_idx

    def iterate_approx_encode_sphere(self, center_idx):
        word = num2dna(center_idx, self.bc_len)
        for seq_idx in FreeDivSphere.FreeDivSphere(word, 2*self.max_err, min_r=self.max_err+1).parallel_num_iterator():
            yield seq_idx

    def _add_barcode(self, seq_idx):
        assert self._idx_is_available(seq_idx), seq_idx
        self.barcodes.add(seq_idx)
        self._codewords.add(seq_idx)
        self._add_codeword(seq_idx)

    def add_dnastr_nonbarcode_codeword(self, dnastring):
        seq_idx = dna2num(dnastring)
        self.add_idx_nonbarcode_codeword(seq_idx)
        
    def add_idx_nonbarcode_codeword(self, seq_idx):
        self.manual_codewords.add(seq_idx)
        self._add_codeword(seq_idx)

    def dnastr_codeword_is_available(self, dnastring):
        seq_idx = dna2num(dnastring)
        self._idx_is_available(seq_idx)
        
    def _idx_is_available(self, test_idx):
        if self.reserved_words[test_idx] != 0:
            return False
        else:
            for seq_idx in self.iterate_decode_sphere(test_idx):
                if self.reserved_words[seq_idx] == 2:
                    return False
            return True

    def Conway_closure(self, tmp_fpath=None):
        for seq_idx in self.seq_idx_iter_func():
            if self._idx_is_available(seq_idx):
                self._add_barcode(seq_idx)
                log.info('Found barcode {}: {}'.format(len(self.barcodes),
                                                       num2dna(seq_idx, self.bc_len)))
                if tmp_fpath:
                    with open(tmp_fpath, 'a') as out:
                        out.write('{}\n'.format(num2dna(seq_idx, self.bc_len)))

    def Conway_closure_until_satisfied(self, n_desired_barcodes):
        for seq_idx in self.seq_idx_iter_func():
            if self._idx_is_available(seq_idx):
                self._add_barcode(seq_idx)
                log.info('Found barcode {}: {}'.format(len(self.barcodes),
                                                       num2dna(seq_idx, self.bc_len)))
                if len(self.barcodes) >= n_desired_barcodes:
                    return

    def restart_Conway_closure(self, prev_fpath, tmp_fpath=None):
        log.info('Restarting {}...'.format(prev_fpath))
        prev_bc_idxs = [dna2num(line.strip()) for line in open(prev_fpath)]
        for bc_idx in prev_bc_idxs:
            self._add_barcode(bc_idx)
            log.info('Adding previous {}: {}'.format(len(self.barcodes),
                                                     num2dna(bc_idx, self.bc_len)))
        with open(tmp_fpath, 'w') as out:
            for bc_idx in prev_bc_idxs:
                out.write('{}\n'.format(num2dna(bc_idx, self.bc_len)))

        seq_iter = self.seq_idx_iter_func()
        max_prev = max(prev_bc_idxs)
        bc_idx = next(seq_iter)
        while bc_idx < max_prev:
            bc_idx = next(seq_iter)
        log.info('Reached last previous barcode: {}'.format(num2dna(max_prev, self.bc_len)))
        log.info('Restarting after {}'.format(num2dna(bc_idx, self.bc_len)))

        for seq_idx in seq_iter:
            if self._idx_is_available(seq_idx):
                self._add_barcode(seq_idx)
                log.info('Found barcode {}: {}'.format(len(self.barcodes),
                                                       num2dna(seq_idx, self.bc_len)))
                if tmp_fpath:
                    with open(tmp_fpath, 'a') as out:
                        out.write('{}\n'.format(num2dna(seq_idx, self.bc_len)))


    def exclude_barcodes(self, exclude_fpath):
        log.info('Excluding barcodes in {}...'.format(exclude_fpath))
        exclude_bc_idxs = [dna2num(line.strip()) for line in open(exclude_fpath)]
        for seq_idx in exclude_bc_idxs:
            assert self.reserved_words[seq_idx] == 0, num2dna(seq_idx, self.bc_len)
            self.reserved_words[seq_idx] = 1


    def find_barcode_4sets(
            self,
            AT_max,
            GC_max,
            seqs_so_far=[],
            prev_spheres=[],
            tmp_fpath=None,
            last_prev_Aidx=None,
            ):
        """
        A barcode 4-set is here defined as a set of four barcodes such that no two barcodes have
        the same base in the same position. I.e., all four bases are in each position in exactly
        one barcode.
        """
        first_base = bases[len(seqs_so_far)]

        seq_idx_iter_func = idx_seq_iterator_avoiding_prev_bases(
            self.bc_len,
            AT_max,
            GC_max,
            first_base,
            seqs_so_far
        )

        # Restart previous run (For 'A' seqs only)
        seq_idx_iter = seq_idx_iter_func()
        if last_prev_Aidx is not None:
            for seq_idx in seq_idx_iter:
                if seq_idx >= last_prev_Aidx:
                    break

        for seq_idx in seq_idx_iter:
            if self._idx_is_available(seq_idx):
                seq_sphere = set(self.iterate_decode_sphere(seq_idx))
                for prev_sphere in prev_spheres:
                    if prev_sphere & seq_sphere:
                        break
                else:
                    new_seqs = seqs_so_far + [num2dna(seq_idx, self.bc_len)]
                    new_spheres = prev_spheres + [seq_sphere]
                    if len(new_seqs) == 4:
                        assert first_base == bases[-1], new_seqs
                        for seq in new_seqs:
                            self._add_barcode(dna2num(seq))
                        if tmp_fpath:
                            with open(tmp_fpath, 'a') as out:
                                out.write('\n'.join(new_seqs) + '\n\n')
                        self.dna_barcode_4sets.append(new_seqs)
                        log.info(
                            'Found barcode set {}: {}'.format(len(self.dna_barcode_4sets), new_seqs)
                        )
                        return
                    else:
                        self.find_barcode_4sets(
                            AT_max,
                            GC_max,
                            new_seqs,
                            new_spheres,
                            tmp_fpath
                        )
                        if len(new_seqs) > 1:
                            return
            elif self.reserved_words[seq_idx] == 0:
                self.reserved_words[seq_idx] = 1


    def restart_barcode_4sets(self, AT_max, GC_max, prev_bc4sets_fpath, tmp_fpath=None):
        with open(prev_bc4sets_fpath) as f:
            while True:
                try:
                    seq_4set = [next(f).strip() for _ in range(4)]
                except StopIteration:
                    break
                assert all(len(seq) == self.bc_len for seq in seq_4set), 'Prev bcs not specified length'
                self.dna_barcode_4sets.append(seq_4set)
                assert not next(f).strip(), next(f)
        
        log.info('Read previous barcodes file')
        
        for i, seq_4set in enumerate(sorted(self.dna_barcode_4sets)):
            for seq in seq_4set:
                self._add_barcode(dna2num(seq))
            log.info('Added prev set {}: {}'.format(i+1, seq_4set))

        if tmp_fpath:
            with open(tmp_fpath, 'w') as out:
                for seq_4set in sorted(self.dna_barcode_4sets):
                    out.write('\n'.join(seq_4set) + '\n\n')

        Aseqs = [seq for seq in self.dna_barcodes if seq.startswith('A')]
        Aseqs.sort()
        last_prev_Aidx = dna2num(Aseqs[-1])

        self.find_barcode_4sets(
            AT_max,
            GC_max,
            tmp_fpath=tmp_fpath,
            last_prev_Aidx=last_prev_Aidx,
        )


    @property
    def dna_barcodes(self):
        return (num2dna(seq_idx, self.bc_len) for seq_idx in self.barcodes)

    def manual_barcodes_test(self):
        bc_list = list(self.barcodes)
        for i in range(len(self.barcodes)):
            bc1 = num2dna(bc_list[i], self.bc_len)
            for j in range(i+1, len(self.barcodes)):
                bc2 = num2dna(bc_list[j], self.bc_len)
                dist = editmeasures.free_divergence(bc1, bc2)
                if dist < self.max_err:
                    log.error('!'*10 + ' FAIL ' + '!'*10)
                    log.error('Distance {} between {} and {}.'.format(dist, bc1, bc2))
                    return
        log.info('Barcodes Pass Manual Check')


def generate_barcodes(arguments):
    if arguments.generate_as_4sets:
        generate_barcode_4sets(arguments)
    else:
        generate_normal_barcodes(arguments)


def generate_normal_barcodes(arguments):
    import time
    start_time = time.time()
    fpath = os.path.join(arguments.output_dir,
                         'barcodes{}-{}.txt'.format(arguments.barcode_length,
                                                    arguments.num_errors))
    tmp_fpath = os.path.join(arguments.output_dir,
                             'barcodes{}-{}.txt.tmp'.format(arguments.barcode_length,
                                                            arguments.num_errors))
    GC_max = frac_to_int_max_GC(arguments.barcode_length, 0.6)
    log.info('Barcode length: {}'.format(arguments.barcode_length))
    log.info('AT/GC max: {}'.format(GC_max))
    bc_iter = idx_possible_barcode_iterator(arguments.barcode_length, GC_max, GC_max)
    sbg = FreeDivBarcodeGenerator(arguments.barcode_length,
                                  arguments.num_errors,
                                  bc_iter)
    if arguments.exclude_bc_fpath:
        sbg.exclude_barcodes(arguments.exclude_bc_fpath)
    if arguments.prev_bc_fpath:
        sbg.restart_Conway_closure(arguments.prev_bc_fpath, tmp_fpath=tmp_fpath)
    else:
        sbg.Conway_closure(tmp_fpath=tmp_fpath)
    with open(fpath, 'w') as out:
        out.write('\n'.join(sorted(sbg.dna_barcodes)))
    os.remove(tmp_fpath)
    comp_time = time.time() - start_time
    log.info('Barcode generation time: {}'.format(comp_time))


def generate_barcode_4sets(arguments):
    import time
    start_time = time.time()
    fpath = os.path.join(arguments.output_dir,
                         'barcode_4sets_{}-{}.txt'.format(arguments.barcode_length,
                                                          arguments.num_errors))
    alph_fpath = os.path.join(
        arguments.output_dir,
        'barcode_4sets_alphabetical_{}-{}.txt'.format(arguments.barcode_length,
                                                      arguments.num_errors)
    )
    tmp_fpath = os.path.join(arguments.output_dir,
                             'barcode_4sets{}-{}.txt.tmp'.format(arguments.barcode_length,
                                                                 arguments.num_errors))
    GC_max = frac_to_int_max_GC(arguments.barcode_length, 0.6)
    log.info('Barcode length: {}'.format(arguments.barcode_length))
    log.info('AT/GC max: {}'.format(GC_max))
    sbg = FreeDivBarcodeGenerator(arguments.barcode_length,
                                  arguments.num_errors)
    if arguments.exclude_bc_fpath:
        sbg.exclude_barcodes(arguments.exclude_bc_fpath)
    if arguments.prev_bc_fpath:
        sbg.restart_barcode_4sets(GC_max, GC_max, arguments.prev_bc_fpath, tmp_fpath)
    else:
        sbg.find_barcode_4sets(GC_max, GC_max, tmp_fpath=tmp_fpath)

    with open(fpath, 'w') as out:
        for bc_4set in sbg.dna_barcode_4sets:
            out.write('\n'.join(bc_4set) + '\n\n')
    with open(alph_fpath, 'w') as out:
        out.write('\n'.join(sorted([bc for bc_4set in sbg.dna_barcode_4sets for bc in bc_4set])))
    os.remove(tmp_fpath)
    comp_time = time.time() - start_time
    log.info('Barcode generation time: {}'.format(comp_time))


