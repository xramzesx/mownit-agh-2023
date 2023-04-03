import os
import numpy as np
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
        result.append( 1/2 * (start + stop) + 1/2 * (stop - start) * cosinus )
    
    return np.array(sorted(result))


#### VISUALISATION ####

def visualise(start, stop, n, function, title, type = "even", spline_type = "cubic", option = "save" ):
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

    plt.title(f'{title} - {translation[type]} - n = {n}')
    plt.xlabel("x")
    plt.ylabel("y")

    #### PLOTS ####

    f = function(points, spline_type)

    plt.plot(domain, fun(domain), label = "Function", color="red")
    plt.plot(domain, np.array([f(x) for x in domain]), label = "Interpolation")
    
    #### INTERSECTIONS ####

    plt.scatter(X,Y)

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
    
### GETTERS ###

def get_h( points ):
    return np.array( [ points[i + 1][0] - points[i][0] for i in range(len(points) - 1)] )

def get_delta( points, degree ):
    n = len(points)
    _x, _y = (0, 1)

    # Δi = f'(x)
    if degree == 1:
        return np.array([
            (points[i+1][_y] - points[i][_y]) /
            (points[i+1][_x] - points[i][_x]) 
            for i in range( n - 1 ) 
        ])

    # 2 * Δ²i = f''(x)    
    if degree == 2:
        delta = get_delta(points, 1)
        return np.array([
            (delta[i + 1] - delta[i]) /
            (points[i + 2][_x] - points[i][_x])
            for i in range( n - 2 )
        ])

    # 6 * Δ³i = f'''(x)
    if degree == 3:
        delta2 = get_delta(points, 2)
        return np.array([
            (delta2[i + 1] - delta2[i]) /
            (points[i + 3][_x] - points[i][_x])
            for i in range( n - 3 )
        ])

    return np.array([])

def get_z(points, type):
    cache = np.zeros(len(points))
    delta = get_delta(points, 1)

    if type == "natural":
        cache[0] = 0
    elif type == "cubic":
        cache[0] = delta[0]

    for i in range(1, len(points)):
        cache[i] = - cache[i - 1] + 2 * delta[i - 1]

    return cache

def generate_matrix(points, type = "cubic" ):
    n =     len(points)
    h =     get_h(points)
    delta  = get_delta(points, 1)
    delta3 = get_delta(points, 3)
    
    left_matrix = np.zeros((n,n))
    right_matrix = np.zeros(n)

    ### MAIN DIAGONAl ###

    for i in range(1, n - 1):
        left_matrix[i][i - 1] = h[i - 1]
        left_matrix[i][i + 1] = h[i]
        left_matrix[i][i] = 2 *( h[i - 1] + h[i] )

    for i in range( 1, n - 1 ):
        right_matrix[i] = delta[i] - delta[i - 1]
        
    ### BOUNDATIRIES ###

    if type == "cubic":

        left_matrix[0][0] =   h[0]
        left_matrix[0][1] = - h[1]

        left_matrix[-1][-2] =   h[-1]
        left_matrix[-1][-1] = - h[-1]

        right_matrix[0] = h[0] ** 2 * delta3[0]
        right_matrix[-1] = h[0] ** 2 * delta3[-1]



    return np.linalg.solve(left_matrix, right_matrix)
#### INTERPOLATION ####

def cubic_spline_interpolation(points, type):
    sigmas = generate_matrix(points, type)
    h = get_h(points)
    n = len(points)

    def interpolation(x):
        nonlocal sigmas, h, n
        _x, _y = (0, 1)
        for j in range( n - 1 ):
            if points[j][_x] <= x and x <= points[j+1][_x]:
                i = j
                break

        b = (points[i + 1][_y] - points[i][_y]) / h[i] - h[i] * ( sigmas[i + 1] + 2 * sigmas[i] )
        c = 3 * sigmas[i]
        d = (sigmas[i+1] - sigmas[i]) / h[i]

        return points[i][_y] + b * (x - points[i][_x]) + c * (x - points[i][_x]) ** 2 + d * (x - points[i][_x]) ** 3

    return interpolation


def quadratic_spline_interpolation( points, type ):
    z = get_z(points, type)
    n = len(points)

    def interpolation(x):
        nonlocal z, n
        _x, _y = (0, 1)

        for j in range( n - 1 ):
            if points[j][_x] <= x and x <= points[j+1][_x]:
                i = j
                break

        return (z[i + 1] - z[i]) / (2 * (points[i + 1][_x] - points[i][_x])) * (x - points[i][_x]) ** 2 + z[i] * (x - points[i][_x]) + points[i][_y]
    
    return interpolation
#### TESTS ####

start = 0
stop = 3 * np.pi
point_counts = [ 4, 5, 7,8, 9, 10, 11, 12, 15, 20, 100 ]

for n in point_counts:
    visualise(start, stop, n, cubic_spline_interpolation, "Spline sześcienny", "even", "cubic")
    visualise(start, stop, n, quadratic_spline_interpolation, "Spline kwadratowy", "even", "cubic")
