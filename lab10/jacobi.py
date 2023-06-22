import numpy as np
import scipy as sp

import time

MAX_ITERATIONS = 1000

#### MATRIX GENERATIONS ####

def generate_x_vector(n, float_type):
    result = np.zeros((n, 1)).astype(float_type)
    for i in range(n):
        result[i] = (-1) ** i
    return result

def generate_first_matrix(n, float_type):
    A = np.zeros(shape = (n,n)).astype(float_type)
    k = 8
    m = 4

    # a[i][i]= k
    for i in range(n):
        A[i][i] = k

    # a[i][j] = 1 / (|i - j| + m), for i,j = 1, 2, ..., n
    for i in range(n):
        _i = i + 1
        for j in range(n):
            if i == j: continue
            
            _j = j + 1
            A[i][j] = 1 / (np.abs(_i - _j) + m)

    return A
    
#### CONSTANTS ####

matrix_sizes = [
    3, 4, 5, 6, 7, 8, 9, 10, 
    11, 12, 13, 14, 15, 20, 
    30, 50, 70, 100, 200, 
    300, 500, 700 
]

accuracies = [10 ** (-i) for i in range(8, 14)]

float_type = np.float64

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

#### STOP CONDITIONS ####

def stop_condition_first(prev_x, next_x, accuracy, **kwargs ):
    return np.linalg.norm(next_x - prev_x) < accuracy

def stop_condition_second(prev_x, A, B, accuracy, **kwargs):
    return np.linalg.norm( A @ prev_x - B) < accuracy

#### BASE VECTORS ####

def first_base_vector (n, float_type):
    return np.zeros((n, 1)).astype(float_type)

def second_base_vector (n, float_type):
    return generate_x_vector(n, float_type) / 2

#### SOLVER ####

def jacobi_solver(A, B, float_type, stop_condition, generate_base_vector, accuracy = 10 ** (-10)):
    global MAX_ITERATIONS
    n = len(A)

    prev_x = generate_base_vector(n, float_type)
    next_x = generate_base_vector(n, float_type)
    
    iterations = 0

    while iterations < MAX_ITERATIONS:
        for i in range(n):
            prev_x[i][0] = next_x[i][0]
        for i in range(n):
            sigma = 0
            for j in range(n):
                if j == i: continue
                sigma += A[i][j] * prev_x[j][0]
            next_x[i][0] = (B[i] - sigma) / A[i][i]
        iterations += 1
        if stop_condition(
            accuracy = accuracy,
            prev_x = prev_x, 
            next_x = next_x,
            A = A,
            B = B
        ):
            break

    return next_x, iterations

#### TESTS ####

def print_results(matrix):
    global matrix_sizes, float_type, accuracies

    print("n\\r", end="\t")
    for accuracy in accuracies:
        print(f"{accuracy:.0e}", end="\t")
    print()

    for i in range(len(matrix_sizes)):
        print(matrix_sizes[i], end="\t")
        for j in range(len(accuracies)):
            print(f"{matrix[i][j]:.6e}", end="\t")
        print()

def test_matrix(generate_matrix, solver, stop_condition = stop_condition_first, generate_base_vector = first_base_vector):
    global matrix_sizes, float_type, accuracies

    matrix_norms = np.zeros((len(matrix_sizes), len(accuracies)))
    matrix_iters = np.zeros((len(matrix_sizes), len(accuracies)))
    matrix_times = np.zeros((len(matrix_sizes), len(accuracies)))

    for i in range(len(matrix_sizes)):
        n = matrix_sizes[i]
        for j in range(len(accuracies)):
            accuracy = accuracies[j]
            
            A = generate_matrix(n, float_type)
            X = generate_x_vector(n, float_type)
            B = A @ X

            start_time = time.perf_counter()
            solved_X, iterations = solver(A, B, float_type, stop_condition, generate_base_vector, accuracy)
            end_time = time.perf_counter()

            matrix_norms[i][j] = np.linalg.norm( X - solved_X )
            matrix_iters[i][j] = iterations
            matrix_times[i][j] = end_time - start_time

    print("Solution errors:")
    print_results(matrix_norms)
    print("Solution iterations:")
    print_results(matrix_iters)
    print("Solution times:")
    print_results(matrix_times)

def time_test_matrix(generate_matrix, solver, stop_condition = stop_condition_first):
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
            solved_X = solver(A, B, float_type, stop_condition)
            end_time = time.perf_counter()
            print(f"{end_time - start_time:.6e}", end="\t")
        print()

def get_spectral_radius(matrix):
    diagonals = np.diagflat(np.diag(matrix))

    inverted_diagonals = invert_matrix(diagonals)

    iteration_matrix = inverted_diagonals @ (matrix - diagonals)

    return np.max(np.abs(np.linalg.eigvals(iteration_matrix)))

def get_spectral_radii(generate_matrix):
    print("n\tspectral radius")
    
    for n in matrix_sizes:
        print(f"{n}\t{get_spectral_radius(generate_matrix(n, float_type)):.6e}")

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

n = 5

print("Spectral radii:")
get_spectral_radii(generate_first_matrix)

print("First stop condition and first vector")
test_matrix(generate_first_matrix, jacobi_solver, stop_condition_first, first_base_vector)
print("First stop condition and second vector")
test_matrix(generate_first_matrix, jacobi_solver, stop_condition_first, second_base_vector)

print("Second stop condition and first vector")
test_matrix(generate_first_matrix, jacobi_solver, stop_condition_second, first_base_vector)
print("Second stop condition and second vector")
test_matrix(generate_first_matrix, jacobi_solver, stop_condition_second, second_base_vector)

