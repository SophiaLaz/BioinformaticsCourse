
def prob_HPATH(hidden_path, T_matrix, states):
    prob = 0.5
    prev = states.index(hidden_path[0])
    for state in hidden_path[1:]:
        prob *= T_matrix[prev][states.index(state)]
        prev = states.index(state)
    return prob


pi = 'AABBBAABABAAAABBBBAABBABABBBAABBAAAABABAABBABABBAB'
states = ['A', 'B']
T_matrix = [[0.194, 0.806], [0.273, 0.727]]

prob_pi = prob_HPATH(pi, T_matrix, states)
print(prob_pi)
