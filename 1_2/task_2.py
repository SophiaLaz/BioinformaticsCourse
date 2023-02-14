import Bio.Align as Align


def Needleman_Wunsch(s1, s2, matrix, gap):

    alphabet = {matrix.alphabet[i]: i for i in range(len(matrix.alphabet))}

    n, m = len(s1), len(s2)
    D = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    D[0] = [i * gap for i in range(n + 1)]
    for i in range(1, m + 1):
        D[i][0] = i * gap
        for j in range(1, n + 1):
            replace = matrix[alphabet[s2[i-1]]][alphabet[s1[j-1]]]
            D[i][j] = max(D[i - 1][j - 1] + replace, D[i - 1][j] + gap, D[i][j - 1] + gap)

    i, j, answer1, answer2 = m, n, [], []
    while i > 0 and j > 0:
        if D[i][j] == D[i][j - 1]:
            answer1.append('_')
            answer2.append(s2[i-1])
            i -= 1
        elif D[i][j] == D[i - 1][j] + gap:
            answer1.append(s1[j-1])
            answer2.append('_')
            j -= 1
        else:
            answer1.append(s1[j-1])
            answer2.append(s2[i-1])
            i -= 1
            j -= 1

    while j > 0:
        answer1.append(s1[j-1])
        answer2.append('_')
        j -= 1
    while i > 0:
        answer1.append('_')
        answer2.append(s2[i-1])
        i -= 1

    print('\nОптимальное выравнивание 1 строки: ', *answer1[::-1], sep='')
    print('Оптимальное выравнивание 2 строки: ', *answer2[::-1], sep='')
    print('Вес выравнивания:', D[-1][-1])


def Gotoh(s1, s2, matrix, alpha, beta):

    alphabet = {matrix.alphabet[i]: i for i in range(len(matrix.alphabet))}

    n, m = len(s1), len(s2)
    D = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    A = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    B = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    D[0] = [alpha + (i - 1) * beta if i > 0 else 0 for i in range(n + 1)]
    B[0] = D[0].copy()
    for i in range(1, m + 1):
        D[i][0] = A[i][0] = alpha + (i - 1) * beta
        for j in range(1, n + 1):
            A[i][j] = max(A[i][j-1] + beta, D[i][j-1] + alpha)
            B[i][j] = max(B[i-1][j] + beta, D[i-1][j] + alpha)

            replace = matrix[alphabet[s2[i-1]]][alphabet[s1[j-1]]]
            D[i][j] = max(D[i - 1][j - 1] + replace, A[i][j], B[i][j])

    i, j, answer1, answer2 = m, n, [], []
    while i > 0 and j > 0:
        if D[i][j] == B[i][j]:
            answer1.append('_')
            answer2.append(s2[i-1])
            i -= 1
        elif D[i][j] == A[i][j]:
            answer1.append(s1[j-1])
            answer2.append('_')
            j -= 1
        else:
            answer1.append(s1[j-1])
            answer2.append(s2[i-1])
            i -= 1
            j -= 1

    while j > 0:
        answer1.append(s1[j-1])
        answer2.append('_')
        j -= 1
    while i > 0:
        answer1.append('_')
        answer2.append(s2[i-1])
        i -= 1

    print('\nОптимальное выравнивание 1 строки: ', *answer1[::-1], sep='')
    print('Оптимальное выравнивание 2 строки: ', *answer2[::-1], sep='')
    print('Вес выравнивания:', D[-1][-1])


if __name__ == "__main__":
    print('\ntask 1')
    s1, s2 = 'GATTACA', 'AAGAGTAC'
    matrix = Align.substitution_matrices.load("BLOSUM62")
    Needleman_Wunsch(s1, s2, matrix, matrix[-1][-2])

    s1, s2 = 'AGTA', 'ATA'
    Needleman_Wunsch(s1, s2, matrix, matrix[-1][-2])

    print('\ntask 2')
    s1, s2 = 'GATTACA', 'AAGAGTAC'
    Gotoh(s1, s2, matrix, matrix[-1][-2], matrix[-1][-2])

    s1, s2 = 'AGTA', 'ATA'
    Gotoh(s1, s2, matrix, matrix[-1][-2], matrix[-1][-2])
