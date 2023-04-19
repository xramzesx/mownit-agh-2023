from matplotlib import pyplot as plt 

#### TRANSLATION ####

translation = {
    "even" : "równomierny",
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

    plt.title(f'{title} - {translation[type]} - n = {n}')
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
    
    print("n\tm\teven max\tchebyschev max\teven square\tchebyschev square")

    for i in range(len(point_counts)):
        n = point_counts[i]
        m = function_counts[i]
        interpolation_cubic = get_errors(interpolation, start, stop, "even", n, m)
        interpolation_natural = get_errors(interpolation, start, stop, "chebyschev", n, m)

        print(
            f'{n}\t'
            f'{m}\t'
            f'{interpolation_cubic[max_err]:.6e}\t'
            f'{interpolation_natural[max_err]:.6e}\t'
            f'{interpolation_cubic[sum_err]:.6e}\t'
            f'{interpolation_natural[sum_err]:.6e}'
        )

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


#### TESTS ####

start = 0
stop = 3 * np.pi

test_count = 7
offset_n = 20
point_counts = [offset_n + 5 * x for x in range(1, test_count + 1)] + [100, 120]
function_counts = [5 * x for x in range(1, test_count + 1)] + [50, 100]

for i in range(len(point_counts)):
    n = point_counts[i]
    m = function_counts[i]
    visualise(start, stop, n, m, algebraic_polynomial_approx, "Aproksymacja wielomianami algebraicznymi", "even")
    visualise(start, stop, n, m, algebraic_polynomial_approx, "Aproksymacja wielomianami algebraicznymi", "chebyschev")

new_points = sorted( point_counts + [ offset_n + 5 * n for n in range(test_count, 2 * test_count)] )
new_functions = sorted( function_counts + [ 5 * n for n in range(test_count, 2 * test_count) ])

test_interpolation(algebraic_polynomial_approx, start, stop, new_points, new_functions)