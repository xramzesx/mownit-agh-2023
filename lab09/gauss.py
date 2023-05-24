import numpy as np

def gauss_solver(A, B, float_type):
    AB = np.hstack((A,B))

    n = len(AB)
    m = len(AB[0])

    row = 0  # Initialization of the pivot row
    col = 0  # Initialization of the pivot column

    while row < n and col < m:
        # Find the k-th pivot:
        i_max = np.argmax(np.abs(AB[row:, col])) + row
        if AB[i_max, col] == 0:
            # No pivot in this column, pass to next column
            col += 1
        else:
            AB[[row, i_max], :] = AB[[i_max, row], :]  # Swap rows(h, i_max)
            # Do for all rows below pivot:
            for i in range(row + 1, n):
                f = AB[i, col] / AB[row, col]
                # Fill with zeros the lower part of pivot column:
                AB[i, col] = 0
                # Do for all remaining elements in current row:
                for j in range(col + 1, m):
                    AB[i, j] -= AB[row, j] * f
            # Increase pivot row and column
            row += 1
            col += 1
    
    # Backward substitution

    X = np.zeros((n, 1)).astype(float_type)

    def subsum_ax(i):
        nonlocal AB, X
        result = 0
        for j in range(i + 1, n):
            result += AB[i][j] * X[j]
        return result

    # x[n] = b'[n] / a'[n][n]
    X[-1] = AB[-1][-1] / AB[n - 1][n - 1]

    for i in range( n - 2, -1, -1):
        X[i] = (AB[i][-1] - subsum_ax(i)) / AB[i][i]

    return X


#### MATRIX GENERATIONS ####

def generate_x_vector(n, float_type):
    result = np.zeros((n, 1)).astype(float_type)
    for i in range(n):
        result[i] = (-1) ** i
    return result

def generate_first_matrix(n, float_type):
    A = np.zeros(shape = (n,n)).astype(float_type)

    # a[1][j]= 1
    for j in range(n):
        A[0][j] = 1
    
    # a[i][j] = 1 / (i + j - 1) for i != 1
    for i in range(1, n):
        _i = i + 1
        for j in range(n):
            _j = j + 1
            A[i][j] = 1 / (_i + _j - 1)
    return A

def generate_second_matrix(n, float_type):
    A = np.zeros(shape = (n,n)).astype(float_type)
    
    for i in range(n):
        _i = i + 1
        for j in range(i, n):
            _j = j + 1
            A[i][j] = 2 * _i / _j
            A[j][i] = A[i][j]
    return A

#### CONSTANTS ####

matrix_sizes = [
    3, 4, 5, 6, 7, 8, 9, 10, 
    11, 12, 13, 14, 15, 20, 
    30, 50, 70, 100, 200, 
    300, 500, 700 
]

float_types = [np.float16, np.float32, np.float64]

def count_condition(matrix):
    return np.linalg.norm(matrix) * np.linalg.norm(invert_matrix(matrix))

def test_matrix(generate_matrix, solver = gauss_solver):
    global matrix_sizes, float_types

    print("n", end="\t")
    for float_type in float_types:
        print(float_type.__name__, end="\t")
    print()

    for n in matrix_sizes:
        print(n, end="\t")
        for float_type in float_types:
            A = generate_matrix(n, float_type)
            X = generate_x_vector(n, float_type)
            B = A @ X
            solved_X = solver(A, B, float_type)
            norm = np.linalg.norm( X - solved_X )
            print(f"{norm:.6e}", end="\t")
        print()

def get_conditions(generate_matrix):
    global matrix_sizes, float_types

    print("n", end="\t")
    for float_type in float_types:
        print(float_type.__name__, end="\t")
    print()

    for n in matrix_sizes:
        print(n, end="\t")
        for float_type in float_types:

            print(f"{count_condition(generate_matrix(n, float_type)):.6e}", end="\t")
        print()

def invert_matrix(M):
    # store dimension
    n = M.shape[0]

    # A must be square with non-zero determinant
    # assert np.linalg.det(M) != 0

    # identity matrix with same shape as A
    I = np.identity(n=n)

    # form the augmented matrix by concatenating A and I
    M = np.concatenate((M, I), axis=1)

    # move all zeros to buttom of matrix
    M = np.concatenate((M[np.any(M != 0, axis=1)], M[np.all(M == 0, axis=1)]), axis=0)

    # iterate over matrix rows
    for i in range(0, n):

        # initialize row-swap iterator
        j = 1

        # select pivot value
        pivot = M[i][i]

        # find next non-zero leading coefficient
        while pivot == 0 and i + j < n:
            # perform row swap operation
            M[[i, i + j]] = M[[i + j, i]]

            # incrememnt row-swap iterator
            j += 1

            # get new pivot
            pivot = M[i][i]

        # if pivot is zero, remaining rows are all zeros
        if pivot == 0:
            # return inverse matrix
            return M[:, n:]

        # extract row
        row = M[i]

        # get 1 along the diagonal
        M[i] = row / pivot

        # iterate over all rows except pivot to get augmented matrix into reduced row echelon form
        for j in [k for k in range(0, n) if k != i]:
            # subtract current row from remaining rows
            M[j] = M[j] - M[i] * M[j][i]

    # return inverse matrix
    return M[:, n:]
