from collections import deque
import numpy as np

#### CONSTANTS ####

start = 0.1
stop  = 1.9

n = int((stop - start) / 0.1 + 2)


# epsilon == rho
epsilon = 1/1000
min_epsilon = 3
max_epsilon = 13

max_iterations = 1000

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
        
        if stop_condition(tmp, x):
            break
        cached.append(x)
    
    return cached
    pass

def generate_secants(x0):
    def secants(x, stop_condition, e = epsilon, max_iter = max_iterations):
        x1 = x0
        x2 = x
        cached = deque([x1,x2])
        for _ in range(max_iter):
            x1, x2 = x2, (x2 - x1) / (fun(x2) - fun(x1)) * fun(x2)
            cached.append(x2)
            if stop_condition(x1, x2):
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
    global start, stop, n

    print("e\\x", end="\t")

    for x in even_space(start, stop, n):
        print(x,end="\t")
    print()
    
    for ei in range(min_epsilon, max_epsilon + 1):
        print(f"10E{-ei}", end="\t")
        for x in even_space(start, stop, n):
            if print_length:
                print(f"{len(solve_method(x, stop_condition, e = 10 ** (-ei)))}", end="\t")
            else:
                print(f"{solve_method(x, stop_condition, e = 10 ** (-ei))[-1]:.6e}", end="\t")
        print()


#### PRINT RESULTS ####

secants_start = generate_secants(start)
secants_stop  = generate_secants(stop)

solve_functions = [newton_rhapsod, secants_start, secants_stop]

for solve_function in solve_functions:
    print(solve_function.__name__)
    generate_result_matrix(solve_function, stop_distance_condition)
    generate_result_matrix(solve_function, stop_distance_condition, True)
    generate_result_matrix(solve_function, stop_residual_condition)
    generate_result_matrix(solve_function, stop_residual_condition, True)
