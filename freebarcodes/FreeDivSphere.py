import itertools
import numpy as np
import logging
from pathos.multiprocessing import ProcessPool
from . import editmeasures
from . import seqtools

log = logging.getLogger(__name__)


bases = 'ACGT'

class FreeDivSphere(object):
    """
    Iterator over the sphere centered at dna c with radius r using Sequence-Levenstein semimetric.
        All iterated seqs have same length as c.

        *** Warning: Some repetition will occur, the amount dependent on center string ***
    
        c :str:     Dna string at center of sphere
        r :int:     Radius of sphere
    """
    def __init__(self, c, r, min_r=0):
        c = c.upper()
        assert set(c) <= set(bases), 'Center object must be DNA.'
        assert r >= 0, 'Radius must be >= 0.'
        assert min_r <= r, 'Min radius must be <= sphere radius.'
        self.c = c
        self.r = r
        self.min_r = min_r

    def __iter__(self):
        for nsub, ndel, nins in self._nsub_ndel_nins_iterator():
            for seq in self._freediv_subsphere_given_counts(nsub, ndel, nins):
                yield seq
    
    def _nsub_ndel_nins_iterator(self):
        for nsub in range(self.r+1):
            for ndel in range(self.r+1 - nsub):
                ins_start = max(0, self.min_r - nsub - ndel)
                ins_end = self.r+1 - nsub - ndel
                for nins in range(ins_start, ins_end):
                    yield nsub, ndel, nins

    def _freediv_subsphere_given_counts(self, nsub, ndel, nins):
        """
        Iterates through seqs with given number of errors from c.
    
            c :str:     Dna string at center
            nsub :int:  Number substitutions
            ndel :int:  Number deletions
            nins :int:  Number insertions
        """
        #--------------------------------------------------------------------------------
        # Notes:
        #
        # All changes in last max(0, nins - ndel) positions will be cleaved off,
        #   so cleave them off ahead of time to save computation.
        # Due to this cleavage, insertions after last cleaved position must be included.
        # Deletions in last position are equivalent to substitutions, so skip them.
        # Deletions next to insertions are equivalent to substitutions, so skip them.
        # Ins-Sub = Sub-Ins, so skip the former.
        # If nins < ndel, all possible additional bases must be added.
        #
        # delidxs are given in original coordinates
        # subidxs are given in post-deletion coords
        # insidxs take work to generate, but at the end are given in post-deletion coords
        #
        # Hence, order of errors must be del-sub-ins
        #--------------------------------------------------------------------------------

        cleaved_k = len(self.c) - max(0, nins - ndel)
        cleaved_seq = self.c[:cleaved_k]

        # We first create the fastest function possible to iterate all full-len seqs 
        if nins >= ndel:
            # Longer seqs already cleaved appropriately. 
            def filled_end_seqs(seq):
                yield seq
        else:
            # Must fill in end
            def filled_end_seqs(seq):
                for end_bases in itertools.product(bases, repeat=ndel-nins):
                    yield seq + ''.join(end_bases)

        # Now find the appropriate coordinate sets and iterate seqs
        all_delidxs = list(range(cleaved_k-1))
        for delidxs in itertools.combinations(all_delidxs, r=ndel):
            # We need only create the deletion seq once per set of indices
            delseq = self._deletion_seq(cleaved_seq, delidxs)

            # Find postdel idxs
            postdel_k = cleaved_k - ndel
            postdel_delidxs = set(idx - sum(1 for didx in delidxs if didx < idx)
                                  for idx in delidxs)

            all_postdel_subidxs = list(range(postdel_k))

            # Don't insert at deletion sites
            all_postdel_insidxs = (set(range(postdel_k + 1))
                                   - postdel_delidxs)

            # Iterate remaining indices
            for subidxs in itertools.combinations(all_postdel_subidxs, r=nsub):
                # Don't insert before substition because Ins-Sub = Sub-Ins
                all_postdel_postsub_insidx_combs = set(
                    tuple(sorted(tup))  # all unique w/multiple insertions per site
                    for tup in itertools.product(all_postdel_insidxs - set(subidxs), repeat=nins)
                )
                for insidxs in all_postdel_postsub_insidx_combs:

                    # All indices established, iterate seqs
                    for subdelseq in self._substitution_seqs(delseq, subidxs):
                        for isd_seq in self._insertion_seqs(subdelseq, insidxs):
                            for seq in filled_end_seqs(isd_seq):
                                yield seq

    def _deletion_seq(self, seq, idxtup):
        """Returns sequence with given deletions from given seq."""
        if not idxtup:
            return seq
        newseq = seq[:idxtup[0]]
        for i, j in zip(idxtup, idxtup[1:]):
            newseq += seq[i+1:j]
        newseq += seq[idxtup[-1]+1:]
        return newseq
    
    def _substitution_seqs(self, seq, idxtup):
        """Iterates all sequences with mutations at given idxs from given seq."""
        if not idxtup:
            yield seq
        else:
            all_mm_bases = [bases.replace(seq[i], '') for i in idxtup]
            for mm_bases in itertools.product(*all_mm_bases):
                newseq = seq[:idxtup[0]]
                for i, c in enumerate(mm_bases[:-1]):
                    newseq += c + seq[idxtup[i] + 1:idxtup[i+1]]
                newseq += mm_bases[-1] + seq[idxtup[-1]+1:]
                yield newseq
    
    def _insertion_seqs(self, seq, idxtup):
        """Iterates all sequences with insertions at given idxs from given seq."""
        if not idxtup:
            yield seq
        else:
            for ins_bases in itertools.product(bases, repeat=len(idxtup)):
                newseq = seq[:idxtup[0]]
                for base_idx, (i, j) in enumerate(zip(idxtup, idxtup[1:])):
                    newseq += ins_bases[base_idx] + seq[i:j]
                newseq += ins_bases[-1] + seq[idxtup[-1]:]
                yield newseq


    def iterator_test(self, iterator='self'):
        log.info('Generating self set...')
        if iterator == 'self':
            self_set = set(self)
        elif iterator == 'parallel_num':
            self_set = set(seqtools.num2dna(seq, len(self.c))
                           for seq in self.parallel_num_iterator())
        else:
            raise ValueError('Invalid iterator to test: {}'.format(iterator))

        log.info('Generating brute force set...')
        bf_set = set(''.join(tup) for tup in itertools.product(bases, repeat=len(self.c))
                     if self.min_r <= editmeasures.free_divergence(self.c, ''.join(tup)) <= self.r)
        log.info('Comparing...')
        if self_set == bf_set:
            log.info('PASS')
        else:
            log.error('#### FAIL ####')
            log.error('{} missing seqs, {} extra seqs'.format(
                len(bf_set - self_set), len(self_set - bf_set)
            ))

    def parallel_num_iterator(self, num_proc=None):
        nerr_tups = list(self._nsub_ndel_nins_iterator())
        def dna_nums_given_nerr_tup(nerr_tup):
            nsub, ndel, nins = nerr_tup
            return [seqtools.dna2num(seq)
                    for seq in self._freediv_subsphere_given_counts(nsub, ndel, nins)]

        pl = ProcessPool(num_proc)
        results = pl.map(dna_nums_given_nerr_tup, nerr_tups)
        for num in itertools.chain(*results):
            yield num
        pl._clear()
        del pl
