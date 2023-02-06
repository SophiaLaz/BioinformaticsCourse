from Bio import SeqIO


def hamming(s1, s2):
    if len(s1) != len(s2):
        print("Different length of string!")
    else:
        return sum([1 if s1[i] != s2[i] else 0 for i in range(len(s1))])


def substring(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1
    m, n = len(s1), len(s2)
    count_min, index = m, 0
    count = hamming(s1, s2[:m])
    for i in range(n - m + 1):
        count = hamming(s1, s2[i:m + i])
        if count < count_min:
            count_min, index = count, i
    return index, s2[index:m + index], count


def levenshtein(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1
    m, n = len(s1), len(s2)
    dist = [[0 for _ in range(n+1)] for _ in range(2)]
    dist[0] = [i for i in range(n+1)]
    for i in range(1, m + 1):
        dist[1][0] = i
        for j in range(1, n + 1):
            mask = int(s1[i-1] != s2[j-1])
            dist[1][j] = min(dist[1][j-1] + 1, dist[0][j], dist[0][j-1] + mask)
        dist[0] = dist[1]
    return dist[-1][-1]


if __name__ == "__main__":
    print("Testing from data/f8.fasta:\n")
    data = list(SeqIO.parse('data/f8.fasta', 'fasta'))
    data = [str(data[0].seq), str(data[1].seq)]
    print("Nearest substring with Hamming's distance:")
    print("{0} - index, {1} - substring, \n{2} - Hamming's distance".format(*substring(*data)))
    print("Levenshtein's distance is ", levenshtein(*data))

    print("Testing from data/gattaca.fasta:\n")
    data = list(SeqIO.parse('data/gattaca.fasta', 'fasta'))
    data = [str(data[0].seq), str(data[1].seq)]
    print("Nearest substring with Hamming's distance:")
    print("{0} - index, {1} - substring, \n{2} - Hamming's distance".format(*substring(*data)))
    print("Levenshtein's distance is ", levenshtein(*data))

