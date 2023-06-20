from typing import List


class MultipleLocalAlignment:
    """
    input: two string: a, b
    ----------
    Class for a Hirshberg's algorithm
    ----------
    output: tree of local aligment
    """

    def __init__(self):
        self.a, self.b = None, None
        self.delete, self.ins, self.match, self.mismatch = [0] * 4

    def _parameters(self, delete: float = -1, insert: float = -1,
                    match: float = 1, mismatch: float = -1):
        self.delete = delete
        self.ins = insert
        self.match = match
        self.mismatch = mismatch

    def _needleman_wunsch(self, a: str, b: str):
        n, m = len(b), len(a)
        D = [[i * self.ins for i in range(n + 1)], [0 for _ in range(n + 1)]]
        for i in range(1, m + 1):
            D[1][0] = i * self.delete
            for j in range(1, n + 1):
                replace = self.match if b[j - 1] == a[i - 1] else self.mismatch
                D[1][j] = max(D[0][j - 1] + replace, D[0][j] + self.delete, D[1][j - 1] + self.ins)
            D[0] = D[1].copy()
        return D[-1]

    def _needleman_wunsch_distance(self, a: str, b: str):
        return self._needleman_wunsch(a, b)[-1]

    def _one_iteration_of_hirshberg(self, a: str, b: str):
        """
        Returns
        -------
        i, j - indexes for a, b
        """
        i = round(len(a) / 2)
        D_1 = self._needleman_wunsch(a[:i], b)
        D_2 = self._needleman_wunsch(a[:i-1:-1], b[::-1])
        D = [D_1[i] + D_2[-1-i] for i in range(len(D_1))]

        return [i, D.index(max(D))]

    def _hirshberg_consensus(self, a: str, b: str):
        i, j = self._one_iteration_of_hirshberg(a, b)

        if len(b[:j]) > 1 and len(a[:i]) > 1:
            left = self._hirshberg_consensus(a[:i], b[:j])
        else:
            left = b[:j] if len(b[:j]) > len(a[:i]) else a[:i]

        if len(b[j:]) > 1 and len(a[i:]) > 1:
            right = self._hirshberg_consensus(a[i:], b[j:])
        else:
            right = b[j:] if len(b[j:]) > len(a[i:]) else a[i:]

        return left + right

    def start(self, s: List[str]):
        self._parameters()
        while len(s) > 1:
            best_indexes, best_distance = None, None
            for i in range(len(s)):
                for j in range(i):
                    distance = self._needleman_wunsch_distance(s[i], s[j])
                    if best_distance is None or distance > best_distance:
                        best_distance = distance
                        best_indexes = i, j
            consensus = self._hirshberg_consensus(s[best_indexes[0]], s[best_indexes[1]])
            s.pop(max(*best_indexes))
            s.pop(min(*best_indexes))
            s.append(consensus)
        return s[0]


if __name__ == "__main__":
    s = ['GATTACA', 'AAGAGTAC', 'GATT', 'AAGGG', 'AAGGGTGCTGTGA', 'TGCTGAGA']
    print("Input:", *s)

    loc_al = MultipleLocalAlignment()
    result = loc_al.start(s)
    print("Result:", result)
