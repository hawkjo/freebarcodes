import numpy as np


def seqlev_dist_matrix(s1, s2):
    n = min(len(s1), len(s2)) + 1
    M = np.zeros((n, n), dtype=np.uint8)

    # First row and column
    M[0, :] = np.arange(n)
    M[:, 0] = np.arange(n)

    # Inner square
    for i in xrange(1, n-1):
        for j in xrange(1, n-1):
            if s1[i-1] == s2[j-1]:
                diag_penalty = M[i-1, j-1]
            else:
                diag_penalty = M[i-1, j-1] + 1
            M[i, j] = min(diag_penalty,
                          M[i-1, j] + 1,
                          M[i, j-1] + 1)

    # Last row
    for i in xrange(1, n-1):
        if s1[n-2] == s2[i-1]:
            diag_penalty = M[n-2, i-1]
        else:
            diag_penalty = M[n-2, i-1] + 1
        M[n-1, i] = min(diag_penalty,
                      M[n-2, i] + 1,
                      M[n-1, i-1])  # No penalty along last row

    # Last column
    for i in xrange(1, n-1):
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

def seqlev_dist(s1, s2):
    M = seqlev_dist_matrix(s1, s2)
    return M[-1, -1]
