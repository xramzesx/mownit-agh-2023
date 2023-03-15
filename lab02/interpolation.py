import numpy as np
from collections import deque
from matplotlib import pyplot as plt 

#### TEST FUNCTION ####

def fun(x):
    return np.cos(x) * x

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

def visualise(start, stop, n, function, title, type = "even" ):
    
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

    #### GENERATE GRAPH ####

    domain = even_space( start, stop, 10000 )

    plt.title(f'Wykres {title} - {type} space')
    plt.xlabel("x")
    plt.ylabel("y")

    #### PLOTS ####
    
    plt.plot(domain, fun(domain), label = "Function",color="red")
    plt.plot(domain, function(points, domain), label = "Interpolation")
    
    #### INTERSECTIONS ####

    plt.scatter(X,Y, label="Points")

    #### SHOW ####

    plt.legend()
    plt.show()

#### TESTS ####

start = 0
stop = 20
n = 20

visualise(start, stop, n, lagrange_interpolation, "Lagrange'a" , "even")
visualise(start, stop, n, newton_interpolation, "Newton'a" , "even")

visualise(start, stop, n, lagrange_interpolation, "Lagrange'a" , "chebyschev")
visualise(start, stop, n, newton_interpolation, "Newton'a" , "chebyschev")
