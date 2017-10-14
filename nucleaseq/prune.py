import sys
import os
import generate
import seqtools
import time

bases = 'ACGT'

def make_is_good_seq(GC_min, GC_max):
    def is_good_seq(seq):
        # Don't allow triplets
        for b in bases:
            if b*3 in seq:
                return False
        # Don't allow GGC
        if 'GGC' in seq:
            return False
        # Enforce GC content within bounds
        if not (GC_min <= seq.count('G') + seq.count('C') <= GC_max):
            return False
        # Don't allow rev-comps of 3+ bp
        seq_rc = seqtools.dna_rev_comp(seq)
        for i in range(len(seq)-2):
            if seq[i:i+2] in seq_rc[:-i-2]:
                return False
        return True
    return is_good_seq
    

def make_iterator(raw_fpath):
    bc_list = [line.strip() for line in open(raw_fpath)]
    bc_len = len(bc_list[0])
    assert all(len(bc) == bc_len for bc in bc_list), set(map(len, bc_list))
    GC_max = min(range(bc_len), key=lambda x: abs(float(x)/bc_len-0.6))
    is_good_seq = make_is_good_seq(bc_len - GC_max, GC_max)

    print 'Barcode length:', bc_len
    print 'AT/GC max:', GC_max
    print 'Starting list size:', len(bc_list)
    bc_list = [bc for bc in bc_list if is_good_seq(bc)]
    print 'Valid sequences:', len(bc_list)
    def iterate_good_barcodes():
        for seq in bc_list:
            yield seqtools.dna2num(seq)
    return iterate_good_barcodes, bc_len


def prune_barcode_file(raw_fpath, max_err, dpath):
    start_time = time.time()
    bc_iter, bc_len = make_iterator(raw_fpath)
    iter_made_time = time.time()
    print 'First filter time:', iter_made_time - start_time

    fpath = os.path.join(dpath, 'barcodes{}-{}.txt'.format(bc_len, max_err))

    sbg = barcode_generator.FreeDivBarcodeGenerator(bc_len, max_err, bc_iter)
    sbg.Conway_closure()
    with open(fpath, 'w') as out:
        out.write('\n'.join(sorted(sbg.dna_barcodes)))
    comp_time = time.time() - start_time
    print
    print 'Barcode pruning time:', time.time() - iter_made_time
    print 'Total time:', comp_time
    stats_fpath = os.path.join(dpath, 'barcodes{}-{}_stats.txt'.format(bc_len, max_err))
    with open(stats_fpath, 'w') as out:
        out.write('Barcode length:\t{:d}\n'.format(bc_len))
        out.write('Max-error:\t{:d}\n'.format(max_err))
        out.write('Barcode pruning time:\t{:.2f} seconds\n'.format(comp_time))


if __name__ == '__main__':
    usg = '{} <bc_fpath> <max_err> <out_dir>'.format(sys.argv[0])
    if len(sys.argv) != len(usg.split()):
        sys.exit('Usage: {}'.format(usg))
    bc_fpath, max_err, dpath = sys.argv[1:]
    max_err = int(max_err)
    prune_barcode_file(bc_fpath, max_err, dpath)
