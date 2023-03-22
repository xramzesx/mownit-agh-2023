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
        result.append( 1/2 * (start + stop) - 1/2 * (stop - start) * cosinus )
    
    return np.array(result)

#### INTERPOLATIONS ####

def lagrange_interpolation(points, x):
    result = 0

    for x_i, y_i in points:
        tmp = y_i

        for x_j, y_j in points:
            if x_i == x_j:
                continue
            tmp *= ( x - x_j ) / (x_i - x_j)
        
        result += tmp

    return result

def newton_interpolation(points, x : float) -> float:
    dy = [ [None for _ in points] for _ in points ]
    n = len(points)
    _x, _y = (0,1)

    # Rewrite y to differences array
    for i in range(n):
        dy[i][0] = points[i][_y]

    # Count didived difference
    for j in range(1, n):
        for i in range(n - j):
            dy[i][j] = (dy[i + 1][j - 1] - dy[i][j - 1]) / (points[i + j][_x] - points[i][_x])
    
    multipler = 1
    result = dy[0][0]

    for i in range( 1, n ):
        multipler *= (x - points[i - 1][_x])
        result += dy[0][i] * multipler

    return result

#### VISUALISATION ####

def visualise(start, stop, n, function, title, type = "even", option = "save" ):
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

    plt.plot(domain, fun(domain), label = "Function",color="red")
    plt.plot(domain, function(points, domain), label = "Interpolation")
    
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
    
#### TESTS ####

start = 0
stop = 3 * np.pi

point_counts = [3, 4, 5, 7, 10, 15, 20, 100]

for n in point_counts:
    visualise(start, stop, n, lagrange_interpolation, "Lagrange" , "even")
    visualise(start, stop, n, newton_interpolation, "Newton" , "even")

    visualise(start, stop, n, lagrange_interpolation, "Lagrange" , "chebyschev")
    visualise(start, stop, n, newton_interpolation, "Newton" , "chebyschev")

new_points = sorted( point_counts + [20 + 10 * n for n in range(1, 8)] )

print(new_points)
test_interpolation(lagrange_interpolation, start, stop, new_points)
test_interpolation(newton_interpolation, start, stop, new_points)
