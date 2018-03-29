import numpy as np
cimport numpy as np
DTYPE = np.uint8
ctypedef np.uint8_t DTYPE_t


def free_divergence_matrix(s1, s2):
    assert len(s1) == len(s2), 'Free divergence requires strings of equal length.'
    cdef int n = len(s1) + 1
    cdef np.ndarray[DTYPE_t, ndim=2] M = np.zeros((n, n), dtype=DTYPE)
    cdef int i, j

    # First row and column
    M[0, :] = np.arange(n)
    M[:, 0] = np.arange(n)

    # Inner square
    for i in range(1, n-1):
        for j in range(1, n-1):
            if s1[i-1] == s2[j-1]:
                diag_penalty = M[i-1, j-1]
            else:
                diag_penalty = M[i-1, j-1] + 1
            M[i, j] = min(diag_penalty,
                          M[i-1, j] + 1,
                          M[i, j-1] + 1)

    # Last row
    for i in range(1, n-1):
        if s1[n-2] == s2[i-1]:
            diag_penalty = M[n-2, i-1]
        else:
            diag_penalty = M[n-2, i-1] + 1
        M[n-1, i] = min(diag_penalty,
                        M[n-2, i] + 1,
                        M[n-1, i-1])  # No penalty along last row

    # Last column
    for i in range(1, n-1):
        if s1[i-1] == s2[n-2]:
            diag_penalty = M[i-1, n-2]
        else:
            diag_penalty = M[i-1, n-2] + 1
        M[i, n-1] = min(diag_penalty,
                        M[i-1, n-1],  # No penalty along last col
                        M[i, n-2] + 1)

    # Last elt
    if s1[n-2] == s2[n-2]:
        diag_penalty = M[n-2, n-2]
    else:
        diag_penalty = M[n-2, n-2] + 1
    M[n-1, n-1] = min(diag_penalty,
                      M[n-2, n-1],  # No penalty along last col
                      M[n-1, n-2])  # No penalty along last row
    return M

def free_divergence(s1, s2):
    M = free_divergence_matrix(s1, s2)
    return M[-1, -1]



def prefix_identification(prefix_seq, observed_seq, int max_err):
    """
    Identifies whether observed seq starts with given prefix seq.

    Returns (optimal observed length, number of edits), or None if edits > max.
    
    Args:
        prefix_seq :str: The expected prefix sequence
        observed_seq :str: The observed full length sequence
        max_err :str: Maximum allowed error
    """
    # Rows associated with prefix, cols with observed sequence
    cdef int ni = len(prefix_seq) + 1                        # num rows
    cdef int nj = min(len(observed_seq) + 1, ni + max_err)   # num cols
    cdef int max_ins = max(0, nj - ni)
    cdef int min_del = max(0, ni - nj)
    cdef float penalty = 1.0/(2*nj) # Choose penalty that will never be >= 1
    cdef int i, j, ndel, nins, best_len, n_err

    cdef np.ndarray[DTYPE_t, ndim=2] M = np.zeros((ni, nj), dtype=DTYPE)
    cdef np.ndarray[np.float_t, ndim=1] last_row = np.zeros((nj,), dtype=np.float)

    # First row and column
    M[:, 0] = np.arange(ni)
    M[0, :] = np.arange(nj)

    # Rest of matrix
    for i in range(1, ni):
        for j in range(1, nj):
            if prefix_seq[i-1] == observed_seq[j-1]:
                diag_penalty = M[i-1, j-1]
            else:
                diag_penalty = M[i-1, j-1] + 1
            M[i, j] = min(diag_penalty,
                          M[i-1, j] + 1,
                          M[i, j-1] + 1)

    # Penalize lengths further from original length
    last_row = np.array(M[ni-1, :], dtype=np.float)
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


cdef int hamming_distance(char *read,
                          char *adapter,
                          int read_length,
                          int adapter_length,
                          int start,
                         ):
    cdef int compare_length = min(adapter_length, read_length - start)
    cdef int mismatches = 0
    cdef int i

    for i in range(compare_length):
        if read[start + i] != adapter[i]:
            mismatches += 1

    return mismatches

def simple_hamming_distance(first, second):
    return hamming_distance(first, second, len(first), len(second), 0)

