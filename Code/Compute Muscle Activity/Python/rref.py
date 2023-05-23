import numpy as np

def rref(A, tol = None):
    """
    RREF Reduced row echelon form.
    R = rref(A) produces the reduced row echelon form of A.
    Usage:
    [R, jb] = rref(A)
    [R, jb] = rref(A, tol) uses the given tolerance in the rank tests.

    INPUTS:
    - A: matrix to be treated and reduced to echelon form
    - tol: tolerance to be used to retain the significant pivots of R

    OUTPUTS:
    - R: reduced matrix
    - jb: indexes of the pivots of R (column where the non-null element of each row is found)
    - order_rows: order of the rows in R, wrt the initial matrix A (how the transformed rows of A are placed in R)
                  (only non-null rows in R are considered here)
    """
    # get the dimensions of the matrix
    m, n = np.shape(A)

    # set the tolerance if not specified by the user
    if tol is None:
        tol = max(m, n) * np.finfo(A.dtype).eps * np.linalg.norm(A, np.inf)

    # loop over the matrix
    i = 0
    j = 0
    jb = []
    rows_order = np.array(range(0, m, 1))

    while i < m and j < n:
        # find value and index of largest element in the remainder of column j
        p = max(abs(A[i:m, j]))
        k = np.where(abs(A[i:m, j]) == p)[0][0]
        k = k+i*1

        # check if p is smaller than the tolerance
        if p <= tol:
            A[i: m, j] = np.zeros(np.shape(A[i:m, j]))
            j = j + 1
        else:
            jb.append(j)
            # Swap i-th and k-th rows
            A[[i, k], j: n] = A[[k, i], j: n]
            rows_order[[i, k]] = rows_order[[k, i]]

            # divide the pivot row by the pivot element
            A[i, j: n] = A[i, j: n] / A[i, j]

            # subtract multiples of the pivot rows from all the other rows
            indexes = np.concatenate((range(0, i), range(i+1, m)))
            for ind in indexes:
                ind = int(ind)
                A[ind, j: n] = A[ind, j: n] - A[ind, j] * A[i, j: n]

            i = i+1
            j = j+1

    return np.array(A), jb, rows_order





