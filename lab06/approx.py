import os
import numpy as np
from collections import deque
from matplotlib import pyplot as plt 

#### TRANSLATION ####

translation = {
    "even" : "w. równoodległe",
    "chebyschev" : "Czebyszewa",
}

#### CONSTANTS ####

plots_dir_path = "plots"

if not os.path.exists(plots_dir_path):
    os.makedirs(plots_dir_path)

#### TEST FUNCTION ####

def fun(x):
    return np.exp(-np.sin(2*x)) + np.sin(2 * x) - 1

#### SPACES ####

def even_space(start, stop, n):
    return np.linspace(start, stop, n)

def chebyschev_space(start, stop, n):
    result = deque()

    for i in range(1, n+1):
        cosinus = np.cos( (2 * i - 1) / (2 * n) * np.pi )
        result.append( 1/2 * (start + stop) + 1/2 * (stop - start) * cosinus )
    
    return np.array(sorted(result))


#### VISUALISATION ####

def visualise(start, stop, n, m, function, title, type = "even", option = "save" ):
    global plots_dir_path

    plt.clf()

    #### GENERATE PROPER SPACE ####
    
    if type == "even":
        X = even_space(start, stop, n)
    elif type == "chebyschev":
        X = chebyschev_space(start, stop, n)
    else:
        print("Specify proper type")
        return
    
    Y = fun(X)
    points = np.column_stack((X, Y))

    #### GENERATE PLOT ####

    domain = even_space( start, stop, 10000 )

    plt.title(f'{title} - n = {n} - m={m}')
    plt.xlabel("x")
    plt.ylabel("y")

    #### PLOTS ####

    f = function(points, m)

    plt.plot(domain, fun(domain), label = "funkcja", color="red")
    plt.plot(domain, np.array([f(x) for x in domain]), label = "aproksymacja")
    
    #### INTERSECTIONS ####

    plt.scatter(X,Y)

    #### SHOW ####

    plt.legend()

    if option == 'save':
        plt.savefig(f'{plots_dir_path}/{title}.{type}.{n}.{m}.png')
    if option == 'show':
        plt.show()

#### ERRORS ####

def max_error(interpolation, start, stop, n = 501):

    error_function = np.vectorize(
        lambda x: np.abs(
            fun(x) - interpolation(x)
        ) 
    )
    return np.max(
        error_function( 
            even_space(start, stop, n) 
        )
    )

def sum_error(interpolation, start, stop, n = 501):
    error_function = np.vectorize(
        lambda x: (
            fun(x) - interpolation( x)
        ) ** 2
    )

    return np.sqrt( np.sum(
        error_function(
            even_space(start, stop, n)
        )
    )) / (n - 1)


def get_errors(
    interpolation, 
    start, 
    stop,
    type,
    base_n, 
    base_m,
    accuracy_n = 501
) -> tuple:

    if type == "even":
        X = even_space(start, stop, base_n)
    else:
        X = chebyschev_space(start, stop, base_n)

    Y = fun(X)
    points = np.column_stack((X, Y))

    f = interpolation(points, base_m)

    return (
        max_error(f, start, stop, accuracy_n),
        sum_error(f, start, stop, accuracy_n)
    )

def test_interpolation(interpolation, start, stop, point_counts, function_counts):

    #### CONSTANTS ####

    max_err = 0
    sum_err = 1
    
    w = len(point_counts)
    k = len(function_counts)

    max_errors = np.zeros((w, k))
    sum_errors = np.zeros((w, k))

    for i in range(w):
        n = point_counts[i]
        for j in range(k):
            m = function_counts[j]
            print(n, m)

            if n < 2 * m: continue

            max_errors[i][j], sum_errors[i][j] = get_errors(interpolation, start, stop, "even", n, m)

    print("n\m", end = "\t")

    for m in function_counts:
        print(m, end="\t")
    
    print()

    for i in range(w):
        print(point_counts[i], end="\t")
        for j in range(k):
            if point_counts[i] < 2 * function_counts[j]:
                print("------------", end="\t")
            else:
                print(f"{max_errors[i][j]:.6e}", end="\t")
        print()

    print("n\m", end = "\t")

    for m in function_counts:
        print(m, end="\t")
    
    print()

    for i in range(w):
        print(point_counts[i], end="\t")

        for j in range(k):
            if point_counts[i] < 2 * function_counts[j]:
                print("------------", end="\t")
            else:
                print(f"{sum_errors[i][j]:.6e}", end="\t")
        print()

#### APPROXIMATION ####

def solve_equation(A, B):
    return np.linalg.solve(A, B)

def algebraic_polynomial_approx(points : np.array, m):

    #### CONSTANTS ####
    _x, _y = (0, 1)
    n = len(points)

    X, Y = points.T

    ### WEIGHTS ###
    weights = np.empty(n)
    weights.fill(1)

    ### MATRIXES ####
    G = np.zeros((m, m))
    B = np.zeros(m)

    for i in range(m):
        for j in range(m):
            G[i][j] = np.sum(weights * X ** (i + j))

    for i in range(m):
        B[i] = np.sum(weights * Y * X ** (i))

    A = solve_equation(G, B)

    def approximation(x):
        nonlocal A
        result = 0
        for j in range(len(A)):
            result += A[j] * x ** j

        return result

    return approximation

def scale_points(x):
    global start, stop
    interval = stop - start
    return x * 2 * np.pi / interval - np.pi - (2 * np.pi * start / interval) 

def trigonometric_approx(points : np.array, m):

    #### CONSTANTS ####
    
    n = len(points)
    X, Y = points.T

    scaled = np.column_stack((scale_points(X), Y))

    A = np.array( [ 2 / n * sum([y * np.cos(x * j) for x, y in scaled]) for j in range(n)] )
    B = np.array( [ 2 / n * sum([y * np.sin(x * j) for x, y in scaled]) for j in range(n)] )

    def approximation(x):
        nonlocal A, B, m
        x = scale_points(x)
        return A[0] / 2 + sum([A[j] * np.cos(j * x) + B[j] * np.sin(j * x) for j in range(1, m + 1)])

    return approximation


#### TESTS ####

start = 0
stop = 3 * np.pi

test_count = 7
offset_n = 20

point_counts = [4,6,8,10,12,14,25,30,35,40,45,50,55,60,65,70,75,80,85]
function_counts = [3,4,5,8,9,10,11,12,15,20,25,30,35,40,45,50,55]

for n in point_counts:
    for m in function_counts:
        if m + 1 > n: continue
        visualise(start, stop, n, m, trigonometric_approx, "Aproksymacja trygonometryczna", "even")


new_points = sorted( list(set(point_counts + [ offset_n + 5 * n for n in range(test_count, 2 * test_count)] )))
new_functions = sorted( list(set(function_counts + [ 5 * n for n in range(test_count, 2 * test_count) ])))

test_interpolation(trigonometric_approx, start, stop, new_points, new_functions)