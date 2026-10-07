"""Exact root certificates for exceptional cosets carrying Gern's D6 interval.

Run with Python 3, no dependencies. Matrices are stored by columns in the basis
of simple roots. E_n numbering: 1--3--4--5--...--n, with 2 attached to 4.
This checks Coxeter length/descent statements, not KL polynomial computations.
"""


def certificate(n, prescribed_strip):
    edges = [(0, 2), (2, 3), (1, 3)] + [(i, i + 1) for i in range(3, n - 1)]
    cartan = [[2 * (i == j) for j in range(n)] for i in range(n)]
    for i, j in edges:
        cartan[i][j] = cartan[j][i] = -1
    identity = tuple(tuple(int(i == j) for i in range(n)) for j in range(n))

    def right_multiply(matrix, s):
        image = matrix[s]
        return tuple(
            tuple(column[i] - cartan[s][j] * image[i] for i in range(n))
            for j, column in enumerate(matrix)
        )

    def is_positive(root):
        return all(coefficient >= 0 for coefficient in root)

    def is_negative(root):
        return all(coefficient <= 0 for coefficient in root)

    def longest(nodes):
        matrix, word = identity, []
        while True:
            choices = [j for j in nodes if is_positive(matrix[j])]
            if not choices:
                return matrix, word
            s = min(choices)
            matrix = right_multiply(matrix, s)
            word.append(s)

    w0, full_word = longest(range(n))
    nodes = list(range(1, 7))
    _, parabolic_word = longest(nodes)
    assert len(parabolic_word) == 30
    a = w0
    for s in parabolic_word:
        assert is_negative(a[s])
        a = right_multiply(a, s)
    length_a = len(full_word) - 30
    assert all(is_positive(a[j]) for j in nodes)
    remainder = a
    for j in prescribed_strip:
        assert is_negative(remainder[j - 1]), (n, j, remainder[j - 1])
        remainder = right_multiply(remainder, j - 1)
    adjacent_descents = [
        (s + 1, t + 1)
        for s, t in edges
        if is_negative(remainder[s]) and is_negative(remainder[t])
    ]
    assert adjacent_descents
    print(f"E{n}: length(a)={length_a}; lengths(a*x6,a*w6)="
          f"({length_a + 4},{length_a + 15}); interval rank=11")
    print("  stripped right descents:", prescribed_strip)
    print("  adjacent right descents remaining:", adjacent_descents)
    for s, t in adjacent_descents:
        print(f"  image(alpha_{s}) = {remainder[s - 1]}")
        print(f"  image(alpha_{t}) = {remainder[t - 1]}")


if __name__ == "__main__":
    certificate(7, [1, 3, 4, 2, 5, 4, 3, 1, 6, 5, 4, 2, 3, 4, 5])
    certificate(8, [1, 3, 4, 2, 5, 4, 3, 1, 6])
