import os
import numpy as np
from collections import deque
from matplotlib import pyplot as plt 

#### CONSTANTS ####

start = 0.1
stop  = 1.9

n = int((stop - start) / 0.1 + 2)

# epsilon == rho
epsilon = 1/1000
min_epsilon = 3
max_epsilon = 13

max_iterations = 1000

#### PLOTS ####

plots_dir_path = "plots"

if not os.path.exists(plots_dir_path):
    os.makedirs(plots_dir_path)


#### FUNCTIONS ####

def fun(x):
    return x**2 - 20 * (np.sin(x)) **12

def dfun(x):
    return 2 * (x - 120 * (np.sin (x)) ** 11 * np.cos(x))

#### NON-LINEAR SOLVERS ####

# x - start point
def newton_rhapsod(x, stop_condition, e = epsilon, max_iter = max_iterations):
    
    cached = [x]

    for _ in range(max_iter):
        tmp = x
        x -= fun(x) / dfun(x)        
        
        if stop_condition(tmp, x, e):
            break
        cached.append(x)
    
    return cached

def generate_secants(start=None, stop = None):

    def secants(x, stop_condition, e = epsilon, max_iter = max_iterations):
        if start is None:
            x1 = x
            x2 = stop
        else:
            x1 = start
            x2 = x

        cached = deque([x1,x2])
        for _ in range(max_iter - 2):
            x1, x2 = x2, x2 - (x2 - x1) / (fun(x2) - fun(x1)) * fun(x2)
            cached.append(x2)
            if stop_condition(x1, x2, e):
                break
            if fun(x1) == fun(x2):
                cached.pop()
                cached.append(np.nan)
                break
        return cached

    return secants


#### CONDITIONS ####

def stop_distance_condition(x1, x2, e = epsilon):
    return np.abs(x1 - x2) < e

def stop_residual_condition(x1, x2, e = epsilon):
    return np.abs(fun(x1)) < e

#### UTILIS ####

def even_space(start, stop, n):
    return np.round(np.linspace(start, stop, n), 1)

def print_row(x, result_arr):
    print(x, len(result_arr), f"{result_arr[-1]:.6e}", sep="\t")

def generate_result_table(solve_method, e = epsilon):
    global start, stop, n

    print("x\tn\tresult")
    for x in even_space(start, stop, n ):
        print_row(x, solve_method(x, stop_distance_condition, e=e))

def generate_result_matrix(solve_method, stop_condition, print_length = False):
    global start, stop, n, min_epsilon, max_epsilon

    print("x\\e", end="\t")

    ei_range = range(min_epsilon, max_epsilon + 1)
    x_range = even_space(start, stop, n)
    

    for ei in ei_range:
        print(f"10E{-ei}", end="\t")

    print()

    for x in x_range:
        print(x, end="\t")
        for ei in ei_range:
            if print_length:
                print(f"{len(solve_method(x, stop_condition, e = 10 ** (-ei)))}", end="\t")
            else:
                print(f"{solve_method(x, stop_condition, e = 10 ** (-ei))[-1]:.6e}", end="\t")
        print()

def visualise(start, stop, pack, option = "save" ):
    functions, title = pack

    global plots_dir_path

    plt.clf()
    num = 10000

    #### GENERATE PROPER SPACE ####
    X = np.linspace(start, stop, num)
    
    Y = fun(X)
    points = np.column_stack((X, Y))

    #### GENERATE PLOT ####

    domain = np.linspace(start, stop, num)

    plt.title(f'{title}')
    plt.xlabel("x")
    plt.ylabel("y")

    #### PLOTS ####

    for function, label, color in functions:
        plt.plot(domain, function(domain), label=label, color=color)

    plt.grid()

    #### SHOW ####

    plt.legend()

    if option == 'save':
        plt.savefig(f'{plots_dir_path}/{title}.png')
    if option == 'show':
        plt.show()



#### PRINT RESULTS ####

function_packs = [
    ([(fun, "funkcja", "red")], "Funkcja"),
    ([(dfun, "pochodna" ,"blue")], "Pochodna"),
    ([(fun, "funkcja", "red"),(dfun, "pochodna", "blue")], "Funkcja z pochodną"),
]

for pack in function_packs:
    visualise(start, stop, pack)

secants_start = generate_secants(start=start)
secants_stop  = generate_secants(stop=stop)

solve_functions = [newton_rhapsod, secants_start, secants_stop]

for solve_function in solve_functions:
    print(solve_function.__name__)
    generate_result_matrix(solve_function, stop_distance_condition)
    generate_result_matrix(solve_function, stop_distance_condition, True)
    generate_result_matrix(solve_function, stop_residual_condition)
    generate_result_matrix(solve_function, stop_residual_condition, True)
