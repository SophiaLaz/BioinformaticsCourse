def viterbi_algorithm(x, alphabet, states, transition_matrix, emission_matrix):
    emission_index = {symbol: index for index, symbol in enumerate(alphabet)}
    state_index = {state: index for index, state in enumerate(states)}

    n = len(x)
    m = len(states)

    # Инициализация таблицы вероятностей
    viterbi = [[0.0] * m for _ in range(n)]
    backtrack = [[0] * m for _ in range(n)]

    # Инициализация первого столбца
    for i in range(m):
        viterbi[0][i] = emission_matrix[i][emission_index[x[0]]]

    # Заполнение таблицы вероятностей и backtrack
    for j in range(1, n):
        for i in range(m):
            max_prob = -1
            max_state = -1
            for k in range(m):
                prob = viterbi[j - 1][k] * transition_matrix[k][i] * emission_matrix[i][emission_index[x[j]]]
                if prob > max_prob:
                    max_prob = prob
                    max_state = k
            viterbi[j][i] = max_prob
            backtrack[j][i] = max_state

    # Восстановление пути
    max_prob = -1
    last_state = -1
    for i in range(m):
        if viterbi[n - 1][i] > max_prob:
            max_prob = viterbi[n - 1][i]
            last_state = i

    path = [states[last_state]]
    for j in range(n - 1, 0, -1):
        last_state = backtrack[j][last_state]
        path.insert(0, states[last_state])

    return ''.join(path)


x = "xyxzzxyxyy"
alphabet = ['x', 'y', 'z']
states = ['A', 'B']
transition_matrix = [
    [0.641, 0.359],
    [0.729, 0.271]
]
emission_matrix = [
    [0.117, 0.691, 0.192],
    [0.097, 0.42, 0.483]
]

path = viterbi_algorithm(x, alphabet, states, transition_matrix, emission_matrix)
print(path)
