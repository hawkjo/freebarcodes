import numpy as np
import FreeDivSphere
import seqtools
import h5py
import time
import psutil


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

    def build_codebook_from_cw_fpath(self, cw_fpath, max_err):
        """
        Builds codebook given path to file with list of one barcode per line.

            cw_fpath :str: barcode file path
            max_err :int: max correctible error
        """
        codewords = [line.strip() for line in open(cw_fpath)]
        self.build_codebook_from_codewords(codewords, max_err)

    def build_codebook_from_codewords(self, codewords, max_err):
        """
        Builds codebook given list or set of codewords

            codewords :iterable: list or set of codewords
            max_err :int: max correctible error
        """
        self.max_err = max_err
        self._codewords = list(codewords)
        self._codewords.sort()
        self._set_cw_len()

        if len(self._codewords) < 2**8:
            dtype = np.uint8
            cw_bytes = 1
        elif len(self._codewords) < 2**16:
            dtype = np.uint16
            cw_bytes = 2
        elif len(self._codewords) < 2**32:
            dtype = np.uint32
            cw_bytes = 4
        elif len(self._codewords) < 2**64:
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

        for i, cw in enumerate(self._codewords):
            cw_idx = i + 1
            for seq in FreeDivSphere.FreeDivSphere(cw, self.max_err):
                seq_idx = seqtools.dna2num(seq)
                self._codebook[seq_idx] = cw_idx


    def _set_cw_len(self):
        assert len(set(map(len, self._codewords))) == 1, \
                'Barcodes not all constant length in {}'.format(cw_fpath)
        self.cw_len = len(self._codewords[0])


    def save_codebook(self, out_fpath):
        dtype = 'S{:d}'.format(self.cw_len)
        with h5py.File(out_fpath) as f:
            f.create_dataset('codewords', (len(self._codewords),), dtype, self._codewords)
            f.create_dataset('codebook', (len(self._codebook),), self._codebook.dtype, self._codebook)
            f.attrs['max_err'] = self.max_err


    def load_codebook(self, codebook_fpath):
        with h5py.File(codebook_fpath) as f:
            self.max_err = int(f.attrs['max_err'])
            self._codewords = list(f['codewords'])
            self._codebook = np.array(f['codebook'])
        self._set_cw_len()


    def decode(self, seq):
        seq_idx = seqtools.dna2num(seq)
        cw_idx = self._codebook[seq_idx]
        if cw_idx == 0:
            return
        else:
            cw_idx -= 1
            return self._codewords[cw_idx]


    def time_decoder(self, n_decodes=1000, verbose=False):
        ground_truth = np.random.choice(self._codewords, n_decodes)
        bcs = [seqtools.add_random_freediv_errors(cw, self.max_err) for cw in ground_truth]

        start_time = time.time()
        decoded = map(self.decode, bcs)
        end_time = time.time()
        decode_time = end_time - start_time

        for i, (gt, bc, dw) in enumerate(zip(ground_truth, bcs, decoded)):
            if not gt == dw:
                raise RuntimeError('Decoding Errors in Test: {} -> {} -> {}'.format(gt, bc, dw))

        if verbose:
            print 'Decoding time: {:.2f} seconds'.format(decode_time)
        return decode_time
