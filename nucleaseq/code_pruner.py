import sys
import numpy as np
import SeqlevSphere
import seqtools
import cPickle
import time
import psutil


class SeqlevCodePruner(object):
    """
    A class for pruning codes which were made assuming SeqLev is a metric.
    """
    # We use the codebook idea to stake a claim to space for codes and prune new ones with
    # overlapping decode spheres.

    def __init__(self):
        pass

    def prune_and_build_codebook_from_cw_fpath(self, cw_fpath, max_err):
        """
        Builds and prunes codebook given path to file with list of one barcode per line.

            :str cw_fpath: barcode file path
            :int max_err: max correctible error
        """
        self._codewords = [line.strip() for line in open(cw_fpath)]
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
        free_bytes = psutil.virtual_memory().free
        if needed_bytes > free_bytes:
            raise RuntimeError('Not enough memory. {:,d} bytes needed, {:,d} bytes free'.format(
                needed_bytes,
                free_bytes
            ))
                
        self._codebook = np.zeros((space_size,), dtype=dtype)

        def decode_sphere_unclaimed(cw):
            for seq in SeqlevSphere.SeqlevSphere(cw, max_err):
                seq_idx = seqtools.dna2num(seq)
                if self._codebook[seq_idx] != 0:
                    return False
            return True


        def claim_decode_sphere(cw, cw_idx):
            for seq in SeqlevSphere.SeqlevSphere(cw, max_err):
                seq_idx = seqtools.dna2num(seq)
                self._codebook[seq_idx] = cw_idx

        i = 0
        while i < len(self._codewords):
            cw = self._codewords[i]
            if decode_sphere_unclaimed(cw):
                cw_idx = i + 1
                claim_decode_sphere(cw, cw_idx)
                i += 1
            else:
                self._codewords.pop(i)


    def _set_cw_len(self):
        assert len(set(map(len, self._codewords))) == 1, \
                'Barcodes not all constant length in {}'.format(cw_fpath)
        self.cw_len = len(self._codewords[0])

    def save_codebook(self, out_fpath):
        with open(out_fpath, 'w') as out:
            cPickle.dump((self._codewords, self._codebook), out)

    def load_codebook(self, codebook_fpath):
        with open(codebook_fpath) as f:
            self._codewords, self._codebook = cPickle.load(f)
        self._set_cw_len()

    def save_codewords(self, out_fpath):
        with open(out_fpath, 'w') as out:
            out.write('\n'.join(self._codewords))


if __name__ == '__main__':
    usg = '{} <input_fpath> <output_fpath>'.format(sys.argv[0])
    if len(sys.argv) != len(usg.split()):
        sys.exit('Usage: ' + usg)

    in_fpath, out_fpath = sys.argv[1:]
    print in_fpath
    print '{} barcodes before pruning'.format(sum(1 for line in open(in_fpath)))

    pruner = SeqlevCodePruner()
    pruner.prune_and_build_codebook_from_cw_fpath(in_fpath, 2)
    pruner.save_codewords(out_fpath)
    print '{} barcodes after pruning'.format(len(pruner._codewords))
