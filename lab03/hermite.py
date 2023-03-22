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

#### SPACES ####

def even_space(start, stop, n):
    return np.linspace(start, stop, n)

def chebyschev_space(start, stop, n):
    result = deque()

    for i in range(1, n+1):
        cosinus = np.cos( (2 * i - 1) / (2 * n) * np.pi )
        result.append( 1/2 * (start + stop) - 1/2 * (stop - start) * cosinus )
    
    return np.array(result)

#### DERIVATIVES ####

def get_derivatives(f, max_derivatives ):
    return [nd.Derivative(f, n=i) for i in range(max_derivatives + 1)]



#### VISUALISATION ####

def visualise(start, stop, n,m, function, title, type = "even", option = "save" ):
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
    
    nodes, derivatives = prepare_nodes(start, stop, n, m)

    interpolated_f = function(nodes, derivatives)
    # Y = fun(X)
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

    # plt.scatter(X,Y)

    #### SHOW ####

    plt.legend()

    if option == 'save':
        plt.savefig(f'{plots_dir_path}/{title}.{type}.{n}.png')
    if option == 'show':
        plt.show()

#### ERRORS ####

def max_error(interpolation, points, start, stop, n = 501):

    error_function = np.vectorize(
        lambda x: np.abs(
            fun(x) - interpolation(points, x)
        ) 
    )
    return np.max(
        error_function( 
            even_space(start, stop, n) 
        )
    )

def sum_error(interpolation, points, start, stop, n = 501):
    error_function = np.vectorize(
        lambda x: (
            fun(x) - interpolation(points, x)
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
    accuracy_n = 501
) -> tuple:

    if type == "even":
        X = even_space(start, stop, base_n)
    elif type == "chebyschev":
        X = chebyschev_space(start, stop, base_n)
    else:
        print("Specify proper type")
        return (None, None)

    Y = fun(X)
    points = np.column_stack((X, Y))

    return (
        max_error(interpolation, points, start, stop, accuracy_n),
        sum_error(interpolation, points, start, stop, accuracy_n)
    )

def test_interpolation(interpolation, start, stop, point_counts):

    #### CONSTANTS ####

    max_err = 0
    sum_err = 1
    
    print("n\teven max\tchebyschev max\teven square\tchebyschev square")

    for n in point_counts:
        interpolation_even = get_errors(interpolation, start, stop, "even", n)
        interpolation_chebyschev = get_errors(interpolation, start, stop, "chebyschev", n)

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


n = 20
m = 10
def hermite_interpolation(nodes, derivatives):
    points = []

    _x, _k = (0,1)

    for x_i, k_i in nodes:
        for m_i in range(int(k_i)):
            points.append((x_i, m_i))

    n = len(points)

    dy = [ [0 for _ in points] for _ in points ]

    # Rewrite y to differences array
    for i in range(n):
        dy[i][0] = points[i][_x]

    # Count divided difference
    for j in range(1, n):
        for i in range(n - j):
            if points[i][_x] == points[i+j][_x]:
                dy[i][j] = derivatives[points[j][_k]](points[j][_x]) / np.math.factorial(j)
            else:
                dy[i][j] = (dy[i + 1][j - 1] - dy[i][j - 1]) / (points[i + j][_x] - points[i][_x])

    def result(x: float):
        multipler = 1
        res = dy[0][0]

        for i in range( 1, n ):
            multipler *= (x - points[i - 1][_x])
            res += dy[0][i] * multipler
        return res

    return result
    multipler = 1

def prepare_nodes( start, stop, n, m ):
    nodes = chebyschev_space(start, stop, m)
    # nodes = even_space(start, stop, m)

    initial_k = np.vectorize(lambda x: (x, 1))
    nodes = np.column_stack(initial_k(nodes))

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

visualise(start, stop, n, m, hermite_interpolation ,"Hermit")
