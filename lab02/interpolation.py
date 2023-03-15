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


points = [ (0,2), (1,3), (2,12), (5,147) ]
x = 12
print(f'Test {lagrange_interpolation(points, x)}')
print(f'Test {newton_interpolation(points, x)}')