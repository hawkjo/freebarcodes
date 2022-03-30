import os
import numpy as np
import h5py
import time
import re
import psutil
import logging
from Bio import SeqIO
from . import FreeDivSphere
from . import seqtools
from . import editmeasures

log = logging.getLogger(__name__)


class FreeDivBarcodeDecoder(object):
    """
    A class for decoding freediv barcodes after possible introduction of errors.
    """
    # This decoder is a simple codebook lookup. For short enough barcodes (as determined by
    # available memory), this is easily done.
    #
    # The basic strategy is to encode each possible barcode as an integer. We then use that integer
    # as an index in our codebook array. At the indexed location is a reference to the barcode of
    # interest. Specifically, each location contains an integer which is then the (index plus one)
    # of the decoded valid codeword position in the codework list. The reason for (index plus one)
    # is that we reserve zero to indicate undecodable words.

    def __init__(self):
        pass

    def build_codebook_from_cw_fpath(self, cw_fpath, max_err_decode):
        """
        Builds codebook given path to file with list of one barcode per line.

            cw_fpath :str: barcode file path
            max_err_decode :int: max correctible error
        """
        codewords = [line.strip() for line in open(cw_fpath)]
        self.build_codebook_from_codewords(codewords, max_err_decode)

    def build_codebook_from_codewords(self, codewords, max_err_decode):
        """
        Builds codebook given list or set of codewords

            codewords :iterable: list or set of codewords
            max_err_decode :int: max correctible error
        """
        self.max_err_decode = max_err_decode
        self._codewords = list(codewords)
        self._codewords.sort()
        self._set_cw_len()
        self._initialize_codebook()

        for i, cw in enumerate(self._codewords):
            cw_idx = i + 1
            for seq in FreeDivSphere.FreeDivSphere(cw, self.max_err_decode):
                seq_idx = seqtools.dna2num(seq)
                self._codebook[seq_idx] = cw_idx

    def build_codebook_from_random_codewords(
            self,
            codewords,
            max_err_decode,
            reject_delta=0):
        """
        Builds codebook given list or set of random, undesigned codewords

            codewords :iterable: list or set of codewords
            max_err_decode :int: max correctible error
            reject_delta :int: additional radius to find and reject conflicts
        """
        assert isinstance(max_err_decode, int) and max_err_decode >= 0
        assert isinstance(reject_delta, int) and reject_delta >= 0
        self.max_err_decode = max_err_decode
        self.reject_delta = reject_delta
        self._codewords = list(codewords)
        self._codewords.sort()
        self._set_cw_len()
        self._initialize_codebook(self.max_err_decode)

        for decode_radius in range(self.max_err_decode, 0, -1): # not including zero
            # Assign decode spheres proper index, rejecting conflicts
            reject_idx = len(self._codewords) + decode_radius
            for i, cw in enumerate(self._codewords):
                cw_idx = i + 1
                for seq in FreeDivSphere.FreeDivSphere(cw, decode_radius):
                    seq_idx = seqtools.dna2num(seq)
                    val = self._codebook[seq_idx] 
                    if val == 0 or val > reject_idx: # unassigned or rejected in a previous round (try again)
                        self._codebook[seq_idx] = cw_idx
                    elif 0 < val <= len(self._codewords) and val != cw_idx: # assigned to other bc (reject)
                        self._codebook[seq_idx] = reject_idx
                    # else val matches cw_idx or was rejected already this round

            if self.reject_delta != 0:
                # Iterate reject spheres, rejecting conflicts
                for i, cw in enumerate(self._codewords):
                    cw_idx = i + 1
                    for seq in FreeDivSphere.FreeDivSphere(
                            cw, 
                            min_r=decode_radius+1,
                            r=decode_radius+self.reject_delta,
                            ):
                        seq_idx = seqtools.dna2num(seq)
                        val = self._codebook[seq_idx] 
                        if 0 < val <= len(self._codewords) and val != cw_idx:
                            self._codebook[seq_idx] = reject_idx
                        # else val unassigned, matches current, or already rejected

        for i, cw in enumerate(self._codewords):
            # Set all barcodes with zero errors to themselves
            cw_idx = i + 1
            seq_idx = seqtools.dna2num(cw)
            self._codebook[seq_idx] = cw_idx


    def analyze_random_codeword_codebook(self):
        """
        Return stats_given_cw with entries 'good', 'bad', 'total', and 'self'
        """
        stats_given_cw = {}
        for i, cw in enumerate(self._codewords):
            cw_idx = i + 1
            stats = {'good': 0, 'bad': 0}
            if self._codebook[seqtools.dna2num(cw)] == cw_idx:
                stats['self'] = 'good'
            else:
                stats['self'] = 'bad'
            for seq in FreeDivSphere.FreeDivSphere(cw, self.max_err_decode):
                seq_idx = seqtools.dna2num(seq)
                if self._codebook[seq_idx] == cw_idx: 
                    stats['good'] += 1
                else:
                    stats['bad'] += 1
            stats['total'] = stats['good'] + stats['bad']
            stats_given_cw[cw] = stats
        return stats_given_cw


    def _initialize_codebook(self, n_reject_levels=0):
        # need space for codewords + 2 (unassigned, conflict => ignore)
        needed_vals = len(self._codewords) + 1 + n_reject_levels
        if needed_vals <= 2**8:
            dtype = np.uint8
            cw_bytes = 1
        elif needed_vals <= 2**16:
            dtype = np.uint16
            cw_bytes = 2
        elif needed_vals <= 2**32:
            dtype = np.uint32
            cw_bytes = 4
        elif needed_vals <= 2**64:
            dtype = np.uint64
            cw_bytes = 8
        else:
            raise ValueError('More than 2^64 barcodes currently not supported')

        space_size = 4**self.cw_len + 1
        needed_bytes = space_size * cw_bytes
        available_bytes = psutil.virtual_memory().available
        if needed_bytes > available_bytes:
            raise RuntimeError('Not enough memory. {:,d} bytes needed, {:,d} bytes available'.format(
                needed_bytes,
                available_bytes
            ))
                
        self._codebook = np.zeros((space_size,), dtype=dtype)


    def _set_cw_len(self):
        assert len(set(map(len, self._codewords))) == 1, \
                'Barcodes not all constant length in {}'.format(cw_fpath)
        self.cw_len = len(self._codewords[0])


    def save_codebook(self, out_fpath):
        dtype = 'S{:d}'.format(self.cw_len)
        with h5py.File(out_fpath, 'w') as f:
            f.create_dataset('codewords', (len(self._codewords),), dtype, self._codewords)
            f.create_dataset('codebook', (len(self._codebook),), self._codebook.dtype, self._codebook)
            f.attrs['max_err_decode'] = self.max_err_decode


    def load_codebook(self, codebook_fpath):
        with h5py.File(codebook_fpath, 'r') as f:
            self.max_err_decode = int(f.attrs['max_err_decode'])
            self._codewords = [bc.decode('ascii') for bc in f['codewords']]
            self._codebook = np.array(f['codebook'])
        self._set_cw_len()


    def decode(self, seq):
        """
        Return barcode if decoded, None if outside all spheres, or negative number by conflict level.
        """
        seq_idx = seqtools.dna2num(seq)
        cw_idx = self._codebook[seq_idx]
        if cw_idx == 0:
            return
        cw_idx -= 1
        if cw_idx >= len(self._codewords):
            return -1 - cw_idx + len(self._codewords)
        return self._codewords[cw_idx]


    def time_decoder(self, n_decodes=1000):
        ground_truth = np.random.choice(self._codewords, n_decodes)
        bcs = [seqtools.add_random_freediv_errors(cw, self.max_err_decode) for cw in ground_truth]

        start_time = time.time()
        decoded = list(map(self.decode, bcs))
        end_time = time.time()
        decode_time = end_time - start_time

        for i, (gt, bc, dw) in enumerate(zip(ground_truth, bcs, decoded)):
            if not gt == dw:
                raise RuntimeError('Decoding Errors in Test: {} -> {} -> {}'.format(gt, bc, dw))

        log.info('Decoding time: {:.2f} seconds'.format(decode_time))
        return decode_time


def decode_fastqs(arguments):
    """
    A wrapper function for decoding barcodes in fastq files. 
    """
    # For speed, the actual method for decoding each file depends on the number of:
    #   prefixes (zero, one, many) 
    #   barcodes (one, many)

    if len(arguments.prefixes) != len(arguments.max_prefix_err):
        raise ValueError('Number of prefixes does not match number of max prefix errors.')

    count_tuple = (
        len(arguments.prefixes),
        min(2, len(arguments.barcode_files))
    )

    decode_fastq_func = {
        (0, 1): decode_no_prefix_one_barcode,
        (0, 2): decode_no_prefix_mult_barcode,
        (1, 1): decode_one_prefix_one_barcode,
        (1, 2): decode_one_prefix_mult_barcode,
        (2, 1): decode_mult_prefix_one_barcode,
        (2, 2): decode_mult_prefix_mult_barcode,
    }
    log.debug('Decode function: {}'.format(decode_fastq_func[count_tuple].__name__))

    decoders = list(map(load_or_build_and_save_decoder, arguments.barcode_files))
    for fastq_fpath in arguments.fastq_files:
        fastq_fname = os.path.split(fastq_fpath)[1]
        fastq_bname = os.path.splitext(fastq_fname)[0]
        out_fpath = os.path.join(arguments.output_dir, fastq_bname + '_decoded.txt')

        log.info('Decoding {}'.format(fastq_fname))
        start_time = time.time()
        decode_fastq_func[count_tuple](
            arguments,
            decoders,
            fastq_fpath,
            out_fpath
        )
        log.debug('Decoded {} in {:.1f} seconds'.format(fastq_fname, time.time() - start_time))


def load_or_build_and_save_decoder(bc_fpath):
    bd = FreeDivBarcodeDecoder()
    bc_fpath_re = re.compile('barcodes(\d+)-(\d+).txt$')
    codebook_fpath = re.sub('barcodes(\d+)-(\d+).txt$', r'codebook\1-\2.txt', bc_fpath)
    log.debug('Codebook path: {}'.format(codebook_fpath))
    start_time = time.time()
    if os.path.exists(codebook_fpath):
        log.info('Loading barcode codebook')
        bd.load_codebook(codebook_fpath)
    else:
        log.info('First use of {}. Building and saving codebook.'.format(bc_fpath))
        num_errors = int(bc_fpath_re.search(bc_fpath).group(2))
        bd.build_codebook_from_cw_fpath(bc_fpath, num_errors)
        bd.save_codebook(codebook_fpath)
    log.debug('Codebook prep time: {:.1f} s'.format(time.time() - start_time))
    return bd


def process_multiple_prefixes(arguments, seq):
    # Multiple options for prefixes, but only one per seq
    for prefix, max_err in zip(arguments.prefixes, arguments.max_prefix_err):
        res = editmeasures.prefix_identification(prefix, seq, max_err)
        if res:
            return prefix, res[0]
    return None, None


def process_multiple_barcodes(decoders, seq, start=0):
    # Process concatenated barcodes in a single sequence
    bcs = []
    for i, decoder in enumerate(decoders):
        bc = decoder.decode(seq[start:start + decoder.cw_len])
        if not bc:
            return None
        bcs.append(bc)
        if i + 1 == len(decoders):
            return bcs
        res = editmeasures.prefix_identification(bc, seq[start:], decoder.max_err_decode)
        if not res:
            return None
        start += res[0]


def decode_no_prefix_one_barcode(arguments, decoders, fastq_fpath, out_fpath):
    decoder = decoders[0]
    with open(out_fpath, 'w') as out:
        for rec in SeqIO.parse(open(fastq_fpath), 'fastq'):
            seq = str(rec.seq)
            bc = decoder.decode(seq[:decoder.cw_len])
            if bc:
                out.write('\t'.join([rec.id, bc, seq]) + '\n')


def decode_no_prefix_mult_barcode(arguments, decoders, fastq_fpath, out_fpath):
    with open(out_fpath, 'w') as out:
        for rec in SeqIO.parse(open(fastq_fpath), 'fastq'):
            seq = str(rec.seq)
            bcs = process_multiple_barcodes(decoders, seq)
            if bcs:
                out.write('\t'.join([rec.id] + bcs + [seq]) + '\n')


def decode_one_prefix_one_barcode(arguments, decoders, fastq_fpath, out_fpath):
    prefix = arguments.prefixes[0]
    max_prefix_err = arguments.max_prefix_err[0]
    decoder = decoders[0]
    with open(out_fpath, 'w') as out:
        for rec in SeqIO.parse(open(fastq_fpath), 'fastq'):
            seq = str(rec.seq)
            res = editmeasures.prefix_identification(prefix, seq, max_prefix_err)
            if not res:
                continue
            prefix_len, _ = res
            bc = decoder.decode(seq[prefix_len:prefix_len + decoder.cw_len])
            if bc:
                out.write('\t'.join([rec.id, prefix, bc, seq]) + '\n')


def decode_one_prefix_mult_barcode(arguments, decoders, fastq_fpath, out_fpath):
    prefix = arguments.prefixes[0]
    max_prefix_err = arguments.max_prefix_err[0]
    with open(out_fpath, 'w') as out:
        for rec in SeqIO.parse(open(fastq_fpath), 'fastq'):
            seq = str(rec.seq)
            res = editmeasures.prefix_identification(prefix, seq, max_prefix_err)
            if not res:
                continue
            prefix_len, _ = res
            bcs = process_multiple_barcodes(decoders, seq, prefix_len)
            if bcs:
                out.write('\t'.join([rec.id, prefix] + bcs + [seq]) + '\n')


def decode_mult_prefix_one_barcode(arguments, decoders, fastq_fpath, out_fpath):
    decoder = decoders[0]
    with open(out_fpath, 'w') as out:
        for rec in SeqIO.parse(open(fastq_fpath), 'fastq'):
            seq = str(rec.seq)
            prefix, prefix_len = process_multiple_prefixes(arguments, seq)
            if prefix is None:
                continue
            bc = decoder.decode(seq[prefix_len:prefix_len + decoder.cw_len])
            if bc:
                out.write('\t'.join([rec.id, prefix, bc, seq]) + '\n')


def decode_mult_prefix_mult_barcode(arguments, decoders, fastq_fpath, out_fpath):
    with open(out_fpath, 'w') as out:
        for rec in SeqIO.parse(open(fastq_fpath), 'fastq'):
            seq = str(rec.seq)
            prefix, prefix_len = process_multiple_prefixes(arguments, seq)
            if prefix is None:
                continue
            bcs = process_multiple_barcodes(decoders, seq, prefix_len)
            if bcs:
                out.write('\t'.join([rec.id, prefix] + bcs + [seq]) + '\n')
