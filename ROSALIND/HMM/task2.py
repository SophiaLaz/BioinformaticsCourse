import numpy as np

def p_emissions(emissions, alphabet, hidden_path, states):
    hidden_path = [states.index(state) for state in hidden_path]
    emissions = [alphabet.index(em) for em in emissions]
    probability = np.prod(E_matrix[hidden_path, emissions])
    return probability

emissions = "xxyzyxzzxzxyxyyzxxzzxxyyxxyxyzzxxyzyzxzxxyxyyzxxzx"
alphabet = ['x', 'y', 'z']
pi = "BBBAAABABABBBBBBAAAAAABAAAABABABBBBBABAABABABABBBB"
states = ['A', 'B']
E_matrix = np.array([[0.612, 0.314, 0.074], [ 0.346, 0.317, 0.336]])

print(p_emissions(emissions, alphabet, pi, states))
