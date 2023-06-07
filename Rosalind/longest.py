def needleman_wunsch(s1, s2, gap):
    n, m = len(s1), len(s2)
    D = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    D[0] = [i * gap for i in range(n + 1)]
    for i in range(1, m + 1):
        D[i][0] = i * gap
        for j in range(1, n + 1):
            if s2[i-1] == s1[j-1]:
                replace = 1
            else:
                replace = 0
            D[i][j] = max(D[i - 1][j - 1] + replace, D[i - 1][j] + gap, D[i][j - 1] + gap)

    res = []
    i, j, answer1, answer2 = m, n, [], []
    while i > 0 and j > 0:
        if D[i][j] == D[i][j - 1] + gap:
            j -= 1
        elif D[i][j] == D[i - 1][j] + gap:
            i -= 1
        else:
            res.append(s2[i-1])
            i -= 1
            j -= 1
    return "".join(reversed(res))


if __name__ == "__main__":
    s1 = input()
    s2 = input()
    print(Needleman_Wunsch(s1, s2, 0))
