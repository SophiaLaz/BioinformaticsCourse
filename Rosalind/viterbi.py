import math


def viterbi(string: str, alphabet : list, states: list, 
            transition_matrix: dict, emission_matrix: dict) -> str:
    
    n = len(string)
    d = [[0.0 for _ in range(len(states))] for _ in range(n + 1)]
    best_point = [[0 for _ in range(len(states))] for _ in range(n + 1)]
    d[1] = [emission_matrix[(states[i], string[0])] for i in range(len(states))]
    for i, c in enumerate(string):
        if i == 0:
            continue
        for l in range(len(states)):
            max_prob, best_point[i + 1][l] = max(
                [(d[i][k] + math.log2(transition_matrix[(states[k], states[l])]), k) for k in range(len(states))]
            )
            d[i + 1][l] = max_prob + math.log2(emission_matrix[states[l], c])

    end_prob, end_state = max([(d[n][k], k) for k in range(len(states))])
    answer = [alphabet[0] for _ in range(n)]
    cur_state = end_state
    for i in range(n, 0, -1):
        answer[i - 1] = states[cur_state]
        cur_state = best_point[i][cur_state]
    return ''.join(answer)


if __name__ == '__main__':
    string = input()
    input()

    alphabet = input().split()
    input()

    states = input().split()
    input()
    input()

    transition_matrix = {}
    for i in range(len(states)):
        s = input().split()
        for j in range(1, len(states) + 1):
            transition_matrix[(s[0], states[j - 1])] = float(s[j])
    input()
    input()

    emission_matrix = {}
    for i in range(len(states)):
        s = input().split()
        for j in range(1, len(alphabet) + 1):
            emission_matrix[(s[0], alphabet[j - 1])] = float(s[j])
    
    print(viterbi(string, alphabet, states, transition_matrix, emission_matrix))
