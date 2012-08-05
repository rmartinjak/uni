#!/usr/bin/env python3

MIN_DIST = 3

def fasta_read(filename):
    """Generator function to read all IDs and sequences from a fasta file."""
    with open(filename) as f:
        line = f.readline()
        while 1:
            if not line or line[0] != '>':
                raise StopIteration

            f_id = line[1:].rstrip()

            line = f.readline()
            while line[0] == ';':
                line = f.readline()

            f_data = ""
            while line and line[0] != '>':
                f_data += line.strip()
                line = f.readline()
                
            yield f_id, f_data


def pairing_score(x, y):
    """Returns the 'pairing score' of two nucleotides."""
    x = x.upper()
    y = y.upper()

    if x not in 'ACGU' or y not in 'ACGU':
        raise ValueError("invalid sequence data")

    if x not in 'AG':
        x, y = y, x

    if x == 'A' and y == 'U':
        return 1
    if x == 'G' and y in 'CU':
        return 1

    return 0


def nussinov(seq):
    """Returns a list of tuples with the indexes of paired nucleotides."""

    seq_length = len(seq)

    # create empty matrix
    mat = []
    for i in range(seq_length):
        mat.append([0] * seq_length)

    # fill matrix
    for n in range(MIN_DIST-1, seq_length):
        for j in range(n, seq_length):
            i = j-n
            mat[i][j] = max(
                mat[i+1][j-1] + pairing_score(seq[i], seq[j]),
                mat[i+1][j],
                mat[i][j-1],
                max(mat[i][k] + mat[k+1][j] for k in range(i+1, j))
                )

    # traceback
    result = []
    stack = [ (0, seq_length-1) ]
    while len(stack) > 0:
        i, j = stack.pop()

        if  j - i <= MIN_DIST:
            continue

        x = mat[i][j]

        if x == mat[i+1][j-1] + pairing_score(seq[i], seq[j]):
            if pairing_score(seq[i], seq[j]) > 0:
                result.append((i, j))
            stack.append((i+1, j-1))

        elif x == mat[i+1][j]:
            stack.append((i+1, j))

        elif x == mat[i][j-1]:
            stack.append((i, j-1))

        else:
            # find k in (i, j) such that mat[i][k] + mat[k+1][j] is maximal
            _, k = max((mat[i][k] + mat[k+1][j], k) for k in range(i+1, j))

            stack.append((i, k))
            stack.append((k+1, j))

    return result



if __name__ == '__main__':
    from sys import argv, exit

    if len(argv) < 2:
        print("usage: %s <file> [<file2>...]" % argv[0])
        exit(1)

    for f in argv[1:]:
        for seq_id, seq_data in fasta_read(f):
            print(seq_id)

            p = [ '-' ] * len(seq_data)
            for opn, cls in nussinov(seq_data):
                p[opn] = '('
                p[cls] = ')'

            print(''.join(p))
            print(seq_data)
