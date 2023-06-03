import numpy as np
import scipy as sp

import time

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


def thomas_solver(A, B, float_type):
    # AB = np.hstack((A,B)).astype(float_type)
    # count coefficients:
    n = len(B)
    
    C = np.zeros(n).astype(float_type)
    C[0] = A[0][1] / A[0][0]

    for i in range(1, n - 1):
        C[i] = A[i][i + 1] / (A[i][i] - A[i][i-1] * C[i-1])
    
    D = np.zeros(n).astype(float_type)
    D[0] = B[0] / A[0][0]

    for i in range(1, n):
        D[i] = (B[i] - A[i][i-1] * D[i - 1]) / (A[i][i] - A[i][i-1] * C[i-1])

    X = np.zeros((n, 1)).astype(float_type)
    X[-1] = D[-1]

    for i in range( n - 2, -1, -1 ):
        X[i] = D[i] - C[i] * X[i + 1]

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

def generate_third_matrix(n, float_type):

    k = m = 5.0
    
    main_diag = np.ones(n) * k
    upper_diag = np.zeros(n-1)
    bottom_diag = np.zeros(n-1)

    for i in range(n - 1):
        upper_diag[i] = 1 / (i + m)
        bottom_diag[i] = k / (i + m + 1)
    
    offsets = np.array([0,1,-1])

    return sp.sparse.diags([main_diag, upper_diag, bottom_diag], offsets, shape = (n,n)).toarray()
    

#### CONSTANTS ####

matrix_sizes = [
    3, 4, 5, 6, 7, 8, 9, 10, 
    11, 12, 13, 14, 15, 20, 
    30, 50, 70, 100, 200, 
    300, 500, 700 
]

float_types = [np.float16, np.float32, np.float64]

#### UTILS ####

def invert_matrix(matrix):
    # store dimension
    n = len(matrix)

    # identity matrix with same shape as A
    I = np.identity(n=n)

    # form the augmented matrix by concatenating A and I
    matrix = np.concatenate((matrix, I), axis=1)

    # move all zeros to buttom of matrix
    matrix = np.concatenate((matrix[np.any(matrix != 0, axis=1)], matrix[np.all(matrix == 0, axis=1)]), axis=0)

    for i in range(n):
        pivot = matrix[i][i]

        for j in range(1, n):
            if pivot != 0 or i + j >= n:
                break

            # swap rows
            matrix[[i, i + j]] = matrix[[i + j, i]]
            pivot = matrix[i][i]

        if pivot == 0:
            return matrix[:, n:]
        
        # extract row
        row = matrix[i]

        # get 1 along the diagonal
        matrix[i] = row / pivot

        for j in range(n):
            if j == i: continue
            matrix[j] = matrix[j] - matrix[i] * matrix[j][i]

    return matrix[:, n:]


def count_condition(matrix):
    return np.linalg.norm(matrix) * np.linalg.norm(invert_matrix(matrix))

#### TESTS ####

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

def time_test_matrix(generate_matrix, solver = gauss_solver):
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

            start_time = time.perf_counter()
            solved_X = solver(A, B, float_type)
            end_time = time.perf_counter()
            print(f"{end_time - start_time:.6e}", end="\t")
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


#### TASKS ####

n = 10

print("Task 1. cond")
get_conditions(generate_first_matrix)
print("Task 2. cond")
get_conditions(generate_second_matrix)

print("Task 1.")
test_matrix(generate_first_matrix)
print("Task 2.")
test_matrix(generate_second_matrix)

print("Task 3. cond")
get_conditions(generate_third_matrix)

print("Task 3. : Thomas")
test_matrix(generate_third_matrix, thomas_solver)

print("Task 3. : Gauss")
test_matrix(generate_third_matrix, gauss_solver)