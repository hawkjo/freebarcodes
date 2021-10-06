import sys
import itertools
import numpy as np
from collections import defaultdict, Counter
import random
import string


bases='ACGT'


dna_complements = str.maketrans('acgtnACGTN', 'tgcanTGCAN')
def dna_rev_comp(dna_string):
    """
    Reverse complement.
    """
    return dna_string.translate(dna_complements)[::-1]


dna2num_tr = str.maketrans('ACGT', '0123')
def dna2num(dna_string):
    """
    Convert dna to number where dna is considered base 4 with '0123' = 'ACGT'.

        s :str:     Given dna string
    """
    return int(dna_string.translate(dna2num_tr), 4)


digs = string.digits + string.ascii_lowercase
def to_base(n, b):
    """
    Convert integer n to a string in base b.

        n :int:     Integer to convert
        b :int:     New base
    """
    res = ""
    while n:
        res += digs[n%b]
        n //= b
    return res[::-1] or "0"


num2dna_tr = str.maketrans('0123', 'ACGT')
def num2dna(n, dnalen):
    """
    Convert number to dna of given length where dna is considered base 4 with '0123' = 'ACGT'

        n :int:         Numerical representation of dna string
        dnalen :int:    Length of dna string
    """
    s = to_base(n, 4).translate(num2dna_tr)
    padlen = dnalen - len(s)
    return 'A'*padlen + s


def frac_to_int_max_GC(k, frac_max_GC):
    """
    Takes a fractional max GC value, like 0.6, and returns the integer maximum number of GC bases
    in a k-mer, length k.
    """
    return min(list(range(k)), key=lambda x: abs(float(x)/k-frac_max_GC))


def mm_names(ref, seq):
    mms = []
    for i, (c1, c2) in enumerate(zip(ref, seq)):
        if c1 != c2:
            mms.append('{}{}{}'.format(c1, i, c2))
    return ','.join(mms)


def get_deletion_seqs(seq, ndel):
    """Returns set of all sequences with ndel deletions from given seq."""
    outset = set()
    for tup in itertools.combinations(list(range(len(seq))), r=ndel):
        newseq = seq[:tup[0]]
        for i, j in zip(tup, tup[1:]):
            newseq += seq[i+1:j]
        newseq += seq[tup[-1]+1:]
        assert len(newseq) == len(seq) - ndel, (tup, newseq)
        outset.add(newseq)
    return outset


def get_contiguous_insertion_seqs(seq, len_ins):
    """Returns set of all sequences with single insertions of length len_ins from given seq."""
    outset = set()
    all_insertions = [''.join(tup) for tup in itertools.product(bases, repeat=len_ins)]
    for i in range(1, len(seq) + 1):
        outset.update([seq[:i] + insertion + seq[i:] for insertion in all_insertions])
    assert all([len(outseq) == len(seq) + len_ins for outseq in outset])
    return outset


def get_insertion_seqs(seq, nins):
    """Returns set of all sequences with nins insertions from given seq."""
    outset = set()
    for tup in itertools.combinations(list(range(1, len(seq) + 1)), r=nins):
        for ins_bases in itertools.product(bases, repeat=nins):
            assert len(ins_bases) == len(tup), (tup, ins_bases)
            newseq = seq[:tup[0]]
            for base_idx, (i, j) in enumerate(zip(tup, tup[1:])):
                newseq += ins_bases[base_idx] + seq[i:j]
            newseq += ins_bases[-1] + seq[tup[-1]:]
            assert len(newseq) == len(seq) + nins, (tup, newseq)
            outset.add(newseq)
    return outset


def get_mismatch_seqs(seq, num_mm):
    """Returns set of all sequences with num_mm mutations from given seq."""
    outset = set()
    for tup in itertools.combinations(list(range(len(seq))), r=num_mm):
        all_mm_bases = [bases.replace(seq[i], '') for i in tup]
        for mm_bases in itertools.product(*all_mm_bases):
            newseq = seq[:tup[0]]
            for i, c in enumerate(mm_bases[:-1]):
                newseq += c + seq[tup[i] + 1:tup[i+1]]
            newseq += mm_bases[-1] + seq[tup[-1]+1:]
            assert len(newseq) == len(seq), '{}\n{}'.format(seq, newseq)
            outset.add(newseq)
    return outset


complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def forward_complement(seq):
    return ''.join([complements[c] for c in seq])


def switch_end_to_complement(seq, num_bases):
    """Replaces num_bases bases of tail with complement"""
    if num_bases <= 0:
        return seq
    return seq[:-num_bases] + forward_complement(seq[-num_bases:])


def get_stretch_of_complement_seqs(seq, num_bases):
    """Returns all seqs with num_bases-long stretches of sequence replaced with complements"""
    outset = set()
    for i in range(len(seq)-num_bases+1):
        outset.add(seq[:i] + forward_complement(seq[i:i+num_bases]) + seq[i+num_bases:])
    return outset


def get_randomized_stretch_seqs(seq, num_bases):
    """Returns all seqs with num_bases-long stretches of randomized nt."""
    outset = set()
    all_randomized = [''.join(tup) for tup in itertools.product(bases, repeat=num_bases)]
    for i in range(len(seq) - num_bases + 1):
        outset.update([seq[:i] + rand_seq + seq[i+num_bases:] for rand_seq in all_randomized])
    return outset


def get_randomized_pam_seqs(seq, num_pam_bases, num_randomized_bases, end='5p'):
    """Returns set of sequences with randomized pam and leading bases, at preferred end."""
    assert num_randomized_bases >= num_pam_bases
    all_randomized = (''.join(tup) for tup in itertools.product(bases, repeat=num_randomized_bases))
    if end == '5p':
        return set([rand_seq + seq[num_pam_bases:] for rand_seq in all_randomized])
    else:
        assert end == '3p', end
        return set([seq[:-num_pam_bases] + rand_seq for rand_seq in all_randomized])


def get_randomized_region_seqs(seq, start, end):
    """Returns set of sequences where seq[start:end] is randomized."""
    assert start < end, (start, end)
    all_randomized = (''.join(tup) for tup in itertools.product(bases, repeat=end-start))
    return set([seq[:start] + rand_seq + seq[end:] for rand_seq in all_randomized])


def get_mismatches_in_region(seq, start, end, num_mm):
    """Return all seqs with given number of mismatches in given region."""
    return set([seq[:start] + mm_seq + seq[end:] for mm_seq in get_mismatch_seqs(seq[start:end], num_mm)])


def get_complementary_bundle_sets(seq):
    """
    Return all sequences with combinations of stretches set to complementary sequence.
    
    For instance, take a sequence of length 13. Considering bundles of length 3 will 
    produce the following bundles:
    
        ... ... ... ....
        
    Then forward-complimenting 2 bundles at a time will produce the following set of sequences:
        
        *** *** ... ....
        *** ... *** ....
        *** ... ... ****
        ... *** *** ....
        ... *** ... ****
        ... ... *** ****
        
    Note the last bundle includes left-over nucleotides. 
    """
    outset = set()
    for bundle_len in range(3, 11, 2):  # only consider bundles up to length 10
        if bundle_len * 2 > len(seq):
            bundles = [(0, len(seq))]
        else:
            bundles = []
            for start in range(0, len(seq) - len(seq) % bundle_len - bundle_len, bundle_len):
                bundles.append((start, start + bundle_len))
            # Extend last bundle to end of sequence
            bundles.append((len(seq) - len(seq) % bundle_len - bundle_len, len(seq)))
            
        for num_bundles in range(2, 4):  
            for selected_bundles in itertools.combinations(bundles, r=num_bundles):
                if len(seq) - sum(end - start for start, end in selected_bundles) <= len(seq)/3.0:
                    # skip if there will be less than a third of the original sequence
                    continue
                newseq = seq[:]
                for start, end in selected_bundles:
                    newseq = newseq[:start] + forward_complement(newseq[start:end]) + newseq[end:]
                assert len(newseq) == len(seq)
                outset.add(newseq)
    return outset


def simple_hamming_distance(s1, s2): 
    return sum(1 for c1, c2 in zip(s1, s2) if c1 != c2)


def build_read_names_given_seq(target,
                               read_names_by_seq_fpath, 
                               allowed_read_names_set, 
                               is_interesting_seq,
                               max_ham,
                               verbose=True):
    interesting_reads = defaultdict(set)
    i = 0
    for i, line in enumerate(open(read_names_by_seq_fpath)):
        if verbose and i % 10000 == 0:
            sys.stdout.write('.')
            sys.stdout.flush()

        words = line.strip().split()
        seq = words[0]
        if is_interesting_seq(seq):
            read_names = set(words[1:]) & allowed_read_names_set
            interesting_reads[seq].update(read_names)
            last_start = len(seq) - len(target)
            if last_start < 0:
                continue
            min_ham_idx = min(list(range(0, last_start+1)),
                              key=lambda i: simple_hamming_distance(target, seq[i:i+len(target)]))
            min_ham = simple_hamming_distance(target, seq[min_ham_idx:min_ham_idx+len(target)])
            if min_ham <= max_ham:
                min_ham_seq = seq[min_ham_idx:min_ham_idx+len(target)]
                interesting_reads[min_ham_seq].update(read_names)
    return interesting_reads


def fill_or_truncate(seq, slen):
    if len(seq) >= slen:
        return seq[:slen]
    else:
        fill = ''.join([random.choice(bases) for _ in range(slen - len(seq))])
        return seq + fill


def add_random_mismatch(seq):
    i = random.randrange(len(seq))
    b = random.choice(bases.replace(seq[i], ''))
    return seq[:i] + b + seq[i+1:]


def add_random_deletion(seq):
    i = random.randrange(len(seq))
    return seq[:i] + seq[i+1:]


def add_random_filled_deletion(seq):
    b = random.choice(bases)
    return add_random_deletion(seq) + b


def add_random_insertion(seq):
    i = random.randrange(len(seq))
    b = random.choice(bases)
    return seq[:i] + b + seq[i:]


def add_random_truncated_insertion(seq):
    return add_random_insertion(seq)[:len(seq)]


def add_random_mismatches(seq, nerr):
    for i in random.sample(list(range(len(seq))), nerr):
        b = random.choice(bases.replace(seq[i], ''))
        seq = seq[:i] + b + seq[i+1:]
    return seq


def add_random_deletions(seq, nerr):
    idxs = random.sample(list(range(len(seq))), nerr)
    idxs.sort(reverse=True)
    for i in idxs:
        seq = seq[:i] + seq[i+1:]
    return seq


def add_random_filled_deletions(seq, nerr):
    return fill_or_truncate(add_random_deletions(seq, nerr), len(seq))


def add_random_insertions(seq, nerr):
    idxs = random.sample(nerr*list(range(len(seq))), nerr)
    idxs.sort(reverse=True)
    for i in idxs:
        b = random.choice(bases)
        seq = seq[:i] + b + seq[i:]
    return seq


def add_random_truncated_insertions(seq, nerr):
    return add_random_insertions(seq, nerr)[:len(seq)]


def add_random_freediv_errors(seq, nerr):
    n_mm = random.randint(0, nerr)
    n_del = random.randint(0, nerr - n_mm)
    n_ins = nerr - n_mm - n_del
    return fill_or_truncate(
        add_random_deletions(
            add_random_insertions(
                add_random_mismatches(seq, n_mm),
                n_ins),
            n_del),
        len(seq)
    )


def add_correlated_errors(seq, psub, pdel, pins, corr_factor):
    # Adds errors according to given probability, where if the previous position contained an error
    # the current position has all probabilities multiplied by corr_factor.

    pi = [pins, pins * corr_factor]
    pd = [pdel, pdel * corr_factor]
    psd = [psub + pdel, (psub + pdel) * corr_factor] # Del and sub chosen

    prev_err = 0
    out = ''
    for c in seq:
        while random.random() < pi[prev_err]:
            out += random.choice(bases)
            prev_err = 1

        v = random.random()
        if v < pd[prev_err]:
            # out += ''
            prev_err = 1
        elif v < psd[prev_err]:
            out += random.choice(bases.replace(c, ''))
            prev_err = 1
        else:
            out += c
            prev_err = 0
    return out


def sample_and_add_errors(seqs, lbda, psub, pdel, pins, corr_factor):
    out_seqs = []
    for seq in seqs:
        ncopies = np.random.poisson(lbda)
        for _ in range(ncopies):
            corrupt_seq = add_correlated_errors(seq, psub, pdel, pins, corr_factor)
            if random.random() < 0.5:
                corrupt_seq = dna_rev_comp(corrupt_seq)
            out_seqs.append(corrupt_seq)
    random.shuffle(out_seqs)
    return out_seqs
