import os
import re
import time
import logging
from .seqtools import bases, dna_rev_comp

log = logging.getLogger(__name__)


bad_triplets = [b*3 for b in bases] + ['GGC']
def go_together(bc, other_bcs, bc_rev_comp_triplets):
    # Check for bad triplets in only two positions with new triplets at the boundary
    rest_with_overhang = bc[-2:] + ''.join(other_bcs)
    if rest_with_overhang[:3] in bad_triplets or rest_with_overhang[1:4] in bad_triplets:
        return False

    # Check for reverse complementary triplets, careful to not overfilter at the boundary
    for bc_rc_trip in bc_rev_comp_triplets[:-2]:
        if bc_rc_trip in rest_with_overhang:
            return False
    if (bc_rev_comp_triplets[-2] in rest_with_overhang[1:]
        or bc_rev_comp_triplets[-1] in rest_with_overhang[2:]):
        return False
    return True


def multiple_barcodes_generator(bc_lists, r):
    if r == 1:
        for bc in bc_lists[0]:
            yield [bc]
    elif r > 1:
        for bc in bc_lists[r-1]:
            bc_rev_comp_triplets = [dna_rev_comp(bc[i:i+3]) for i in range(len(bc)-3)]
            for other_bcs in multiple_barcodes_generator(bc_lists, r=r-1):
                if go_together(bc, other_bcs, bc_rev_comp_triplets):
                    yield other_bcs + [bc] 
    else:
        raise ValueError('r < 1 encountered: {}'.format(r))


def concatenate_barcodes(arguments):
    if len(arguments.barcode_files) <= 1:
        raise ValueError('Concatenate requires more than one barcode file.')
    bc_len_re = re.compile('barcodes(\d+)-\d+.txt')
    for bc_fpath in arguments.barcode_files:
        if int(bc_len_re.search(bc_fpath).group(1)) < 5:
            raise ValueError('Concatenated sub-barcodes must be at least 5 bp long.')

    log.info('Loading sub-barcodes')
    bc_lists = [[line.strip() for line in open(bc_fpath)]
                for bc_fpath in arguments.barcode_files]
    # The output name simply lists all the length and error-correction tuples separated with _'s
    fpath_re = re.compile('barcodes(\d+-\d+).txt$')
    out_fname = 'barcodes{}.txt'.format(
        '_'.join(fpath_re.search(bc_fpath).group(1) for bc_fpath in arguments.barcode_files)
    )
    out_fpath = os.path.join(arguments.output_dir, out_fname)

    log.info('Writing to {}'.format(out_fpath))
    lim = arguments.max_bc or 1e9
    start_time = time.time()
    with open(out_fpath, 'w') as out:
        for i, bcs in enumerate(multiple_barcodes_generator(bc_lists, len(bc_lists))):
            if i >= lim:
                break
            out.write('\t'.join(bcs) + '\n')
    log.info('Concatenation time: {:.1f} seconds'.format(time.time() - start_time))
