import os
import numpy as np
import random
import numdifftools as nd
from collections import deque
from matplotlib import pyplot as plt 

#### TRANSLATION ####

translation = {
    "even" : "równomierny",
    "chebyschev" : "Chebyschewa"
}

#### CONSTANTS ####

plots_dir_path = "plots"

if not os.path.exists(plots_dir_path):
    os.makedirs(plots_dir_path)

#### TEST FUNCTION ####

def fun(x):
    return np.exp(-np.sin(2*x)) + np.sin(2 * x) - 1

def d_fun(x):
    return 2 * (1 - np.exp( -np.sin( 2 * x ))) * np.cos(2 * x)

#### SPACES ####

def even_space(start, stop, n):
    return np.linspace(start, stop, n)

def chebyschev_space(start, stop, n):
    result = deque()

    for i in range(1, n+1):
        cosinus = np.cos( (2 * i - 1) / (2 * n) * np.pi )
        result.append( 1/2 * (start + stop) + 1/2 * (stop - start) * cosinus )
    
    return np.array(result)

#### DERIVATIVES ####

def get_derivatives(f, max_derivatives ):
    if max_derivatives == 2:
        print("uses analitical derivatives")
        return [fun, d_fun]

    return [nd.Derivative(f, n=i) for i in range(max_derivatives + 1)]



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
    
    nodes, derivatives = prepare_nodes(start, stop, n, m, type)

    interpolated_f = function(nodes, derivatives)

    #### GENERATE PLOT ####

    domain = even_space( start, stop, 10000 )
    plt.title(f'{title} - {translation[type]} - n = {n}')
    plt.xlabel("x")
    plt.ylabel("y")

    #### PLOTS ####

    plt.plot(domain, fun(domain), label = "Function",color="red")
    plt.plot(domain, interpolated_f(domain), label = "Interpolation")
    # plt.plot(domain, function(domain), label = "Interpolation")
    
    #### INTERSECTIONS ####

    Y = fun(X)

    plt.scatter(X,Y)

    #### SHOW ####

    plt.legend()

    if option == 'save':
        plt.savefig(f'{plots_dir_path}/{title}.{type}.{n}.png')
    if option == 'show':
        plt.show()

#### ERRORS ####

def max_error(interpolation, nodes, derivatives, start, stop, n = 501, base_m=2):
    error_function = np.vectorize(
        lambda x: np.abs(
            fun(x) - interpolation(nodes, derivatives)(x)
        ) 
    )
    return np.max(
        error_function( 
            even_space(start, stop, n) 
        )
    )

def sum_error(interpolation, nodes, derivatives, start, stop, n = 501, base_m = 2):
    error_function = np.vectorize(
        lambda x: (
            fun(x) - interpolation(nodes, derivatives)(x)
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
    base_m = 2,
    accuracy_n = 501
) -> tuple:

    nodes, derivatives = prepare_nodes(start, stop, base_n, base_m, type)

    return (
        max_error(interpolation, nodes, derivatives, start, stop, accuracy_n, base_m=base_m),
        sum_error(interpolation, nodes, derivatives, start, stop, accuracy_n, base_m=base_m)
    )

def test_interpolation(interpolation, start, stop, point_counts, m):

    #### CONSTANTS ####

    max_err = 0
    sum_err = 1
    
    print("n\teven max\tchebyschev max\teven square\tchebyschev square")

    for n in point_counts:
        interpolation_even = get_errors(interpolation, start, stop, "even", n, base_m=m)
        interpolation_chebyschev = get_errors(interpolation, start, stop, "chebyschev", n, base_m = m)

        print(
            f'{n}\t'
            f'{interpolation_even[max_err]:.6e}\t'
            f'{interpolation_chebyschev[max_err]:.6e}\t'
            f'{interpolation_even[sum_err]:.6e}\t'
            f'{interpolation_chebyschev[sum_err]:.6e}'
        )
    
#### TESTS ####

start = 0
stop = 3 * np.pi


n = 10
m = 2
def divided_differences(points, derivatives, size):
    diffs = [[0 for _ in range(size)] for _ in range(size)]

    _x, _k = (0 , 1)

    for i in range(size):
        diffs[i][0] = fun(points[i][_x])

    for j in range(1, size):
        for i in range(size - j):
            if points[i][_x] == points[i+j][_x]:
                diffs[i][j] = derivatives[points[j][_k]](points[j][_x]) \
                         / np.math.factorial(j)
            else:
                diffs[i][j] = (diffs[i + 1][j - 1] - diffs[i][j - 1]) / (
                    points[i + j][_x] - points[i][_x]
                )
    print(size, diffs[0][size - 1])
    return diffs[0][size - 1]

def hermite_interpolation(nodes, derivatives):
    points = []

    ### CONSTANTS ###
    _x, _k = (0,1)

    #### NODES TO POINTS ####
    # k_i -> derivative degree
    # x_i -> argument x
    for x_i, k_i in reversed(nodes):
        for m_i in range(int(k_i)):
            points.append((x_i, m_i))

    n = len(points)
    dy = [ [0 for _ in points] for _ in points ]
    
    # Rewrite y to differences array
    for i in range(n):
        dy[i][0] = fun( points[i][_x])
    
    # Count divided difference
    for j in range(1, n):
        for i in range(n - j):
            if points[i][_x] == points[i+j][_x]:
                dy[i][j] = derivatives[points[j][_k]](points[j][_x]) \
                         / np.math.factorial(j)
            else:
                dy[i][j] = (dy[i + 1][j - 1] - dy[i][j - 1]) / (
                    points[i + j][_x] - points[i][_x]
                )


    def result(x: float):
        nonlocal dy, _x, n

        multipler = 1
        res = 0
        for i in range(n):
            res += dy[0][i] * multipler
            multipler *= (x - points[i][_x])

        return res

    return result


def prepare_nodes( start, stop, n, m, type ):
    nodes = chebyschev_space(start, stop, n) if type == "chebyschev" else even_space(start, stop, n)
   
    initial_k = np.vectorize(lambda x: (x, m))
    nodes = np.column_stack(initial_k(nodes))

    derivatives = get_derivatives(fun, m)

    return (nodes, derivatives)

    remaining = n - m + 1
    maximum = 1

    ### RANDOMIZE ####

    while remaining > 0:
        i = random.randint(0, m - 1)
        nodes[i][1] += 1
        remaining -= 1

        maximum = max(maximum, nodes[i][1])
    derivatives = get_derivatives( fun, int(maximum) )
    print(nodes)
    return (nodes, derivatives)

point_counts = [3, 4, 5, 7,8, 9, 10, 15, 20]

for n in point_counts:
    visualise(start, stop, n, m, hermite_interpolation, "Hermite" , "even", "save")

    visualise(start, stop, n, m, hermite_interpolation, "Hermite" , "chebyschev", "save")
test_interpolation(hermite_interpolation, start, stop, point_counts, m)
