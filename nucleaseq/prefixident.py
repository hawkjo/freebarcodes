import numpy as np


def prefix_identification(prefix_seq, observed_seq, max_err):
    """
    Identifies whether observed seq starts with given prefix seq.

    Returns (optimal observed length, number of edits), or None if edits > max.
    
    Args:
        prefix_seq :str: The expected prefix sequence
        observed_seq :str: The observed full length sequence
        max_err :str: Maximum allowed error
    """
    # Rows associated with prefix, cols with observed sequence
    ni = len(prefix_seq) + 1                        # num rows
    nj = min(len(observed_seq) + 1, ni + max_err)   # num cols
    max_ins = max(0, nj - ni)
    min_del = max(0, ni - nj)

    M = np.zeros((ni, nj), dtype=np.uint8)

    # First row and column
    M[:, 0] = np.arange(ni)
    M[0, :] = np.arange(nj)

    # Rest of matrix
    for i in xrange(1, ni):
        for j in xrange(1, nj):
            if prefix_seq[i-1] == observed_seq[j-1]:
                diag_penalty = M[i-1, j-1]
            else:
                diag_penalty = M[i-1, j-1] + 1
            M[i, j] = min(diag_penalty,
                          M[i-1, j] + 1,
                          M[i, j-1] + 1)

    # Penalize lengths further from original length
    last_row = np.array(M[ni-1, :], dtype=np.float)
    penalty = 1.0/(2*nj) # Choose penalty that will never be >= 1
    for ndel in range(min_del, ni):
        last_row[ni - 1 - ndel] += penalty * ndel
    for nins in range(1, max_ins + 1):
        last_row[ni - 1 + nins] += 1.1 * penalty * nins # Prefer deletions to insertions
    
    # Test for score within tolerance
    best_len = last_row.argmin()
    n_err = M[ni-1, best_len]
    if n_err <= max_err:
        return best_len, n_err
    else:
        return None
