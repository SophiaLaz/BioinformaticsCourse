def neighbor_joining(D, n):
    if n == 2:
        T = [[0, 1, D[0][1]]]
        return T

    D_prime = construct_neighbor_joining_matrix(D, n)
    i, j = find_min_non_diagonal_element(D_prime)
    delta = (total_distance(D, i) - total_distance(D, j)) / (n - 2)
    limb_length_i = (D[i][j] + delta) / 2
    limb_length_j = (D[i][j] - delta) / 2

    m = n
    D = add_row_column(D, m, i, j)
    D = remove_row_column(D, i, j)

    T = neighbor_joining(D, n - 1)

    # Add two new limbs to the tree
    T.append([m, i, limb_length_i])
    T.append([m, j, limb_length_j])

    return T


def construct_neighbor_joining_matrix(D, n):
    D_prime = [[0] * n for _ in range(n)]
    total_distances = [sum(row) for row in D]

    for i in range(n):
        for j in range(n):
            if i != j:
                D_prime[i][j] = (n - 2) * D[i][j] - total_distances[i] - total_distances[j]

    return D_prime


def find_min_non_diagonal_element(D_prime):
    min_value = float('inf')
    min_i, min_j = 0, 0
    n = len(D_prime)

    for i in range(n):
        for j in range(i + 1, n):
            if D_prime[i][j] < min_value:
                min_value = D_prime[i][j]
                min_i, min_j = i, j

    return min_i, min_j


def total_distance(D, i):
    return sum(D[i])


def add_row_column(D, m, i, j):
    n = len(D)
    new_row = [0] * (n + 1)
    new_col = [0] * (n + 1)

    for k in range(n):
        if k != i and k != j:
            new_row[k] = (D[k][i] + D[k][j] - D[i][j]) / 2
            new_col[k] = new_row[k]

    D.append(new_row)

    for row in D:
        row.append(new_col.pop(0))

    D[n][n] = 0

    return D


def remove_row_column(D, i, j):
    D.pop(j)

    for row in D:
        row.pop(j)

    D.pop(i)

    for row in D:
        row.pop(i)

    return D


# Чтение входных данных
n = int(input())
distance_matrix = []

for _ in range(n):
    row = list(map(float, input().split()))
    distance_matrix.append(row)

# Решение задачи и получение результата
tree = neighbor_joining(distance_matrix, n)

# Вывод результата
for edge in tree:
    print(f'{edge[0]}->{edge[1]}:{edge[2]:.3f}')
