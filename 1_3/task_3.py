from binarytree import Node as leave_tree


class node:
    def __init__(self, value: list):
        self.left = None
        self.right = None
        self.value = value

    def add_left(self, val):
        self.left = val

    def add_right(self, val):
        self.right = val


class local_alignment:
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

    def _parameters(self, delete: float = -2, insert: float = -2,
                    match: float = 2, mismatch: float = -1):
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

    def _hirshberg(self, a: str, b: str, leave: node):
        i, j = self._one_iteration_of_hirshberg(a, b)
        left = node([[leave.value[0][0], leave.value[0][0] + i], [leave.value[1][0], leave.value[1][0] + j]])
        leave.add_left(left)
        right = node([[leave.value[0][0] + i, leave.value[0][1]], [leave.value[1][0] + j, leave.value[1][1]]])
        leave.add_right(right)

        if len(b[:j]) > 1 and len(a[:i]) > 1:
            self._hirshberg(a[:i], b[:j], left)
        if len(b[j:]) > 1 and len(a[i:]) > 1:
            self._hirshberg(a[i:], b[j:], right)

    def start(self, a: str, b: str):
        self.a, self.b = a, b
        result = node([[0, len(a)], [0, len(b)]])
        self._parameters()
        self._hirshberg(self.a, self.b, result)

        def disp(edge, leave):
            if edge.right is not None:
                right = edge.right.value
                leave.right = leave_tree('({0}, {1})'.format(self.a[right[0][0]:right[0][1]], self.b[right[1][0]:right[1][1]]))
                disp(edge.right, leave.right)

            if edge.left is not None:
                left = edge.left.value
                leave.left = leave_tree('({0}, {1})'.format(self.a[left[0][0]:left[0][1]], self.b[left[1][0]:left[1][1]]))
                disp(edge.left, leave.left)

        root = result.value
        tree = leave_tree('({0}, {1})'.format(self.a[root[0][0]:root[0][1]], self.b[root[1][0]:root[1][1]]))
        disp(result, tree)
        print(tree)


if __name__ == "__main__":
    s1, s2 = 'GATTACA', 'AAGAGTAC'
    a = 'AGTACGCA'
    b = 'TATGC'
    loc_al = local_alignment()
    loc_al.start(a, b)
    loc_al.start(s1, s2)
