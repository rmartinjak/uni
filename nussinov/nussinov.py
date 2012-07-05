#!/usr/bin/env python3

MIN_DIST = 3

class Stack(object):
    def __init__(self):
        self.items = []

    def push(self, data):
        self.items.append(data)

    def pop(self):
        return self.items.pop()

    def __iter__(self):
        while len(self.items) > 0:
            yield self.items.pop()


def fasta_read(filename):
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

def pair2(x, y):
    if x == 'A' and y == 'U':
        return 1
    if x == 'G':
        if y == 'C' or y == 'U':
            return 1
    return 0

def pair(x, y):
    x = x.upper()
    y = y.upper()

    if x not in 'ACGU' or y not in 'ACGU':
        raise ValueError("invalid sequence data")

    return max(pair2(x,y), pair2(y, x))


def nussinov(seq):
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
                mat[i+1][j-1] + pair(seq[i], seq[j]),
                mat[i+1][j],
                mat[i][j-1],
                max([mat[i][k] + mat[k+1][j] for k in range(i+1, j)])
                )

    # traceback
    s = Stack()
    result = []
    s.push((0, seq_length-1))
    for p in s:
        i, j = p
        if  j - i <= MIN_DIST:
            continue

        x = mat[i][j]

        if x == mat[i+1][j-1] + pair(seq[i], seq[j]):
            if pair(seq[i], seq[j]) > 0:
                result.append(p)
            s.push((i+1, j-1))

        elif x == mat[i+1][j]:
            s.push((i+1, j))

        elif x == mat[i][j-1]:
            s.push((i, j-1))

        else:
            mx = 0
            mx_k = -1
            for k in range(i+1, j):
                tmp = mat[i][k] + mat[k+1][j] 
                if tmp > mx:
                    mx = tmp
                    mx_k = k

            if mx_k == -1:
                raise ValueError
            s.push((i, k))
            s.push((k+1, j))

    return result

if __name__ == '__main__':
    from sys import argv, exit

    if len(argv) < 2:
        print("usage: %s <file> [<file2>...]" % argv[0])
        exit(1)

    for f in argv[1:]:
        for fastadata in fasta_read(f):
            p = ""
            (seq_id, seq) = fastadata
            print(seq_id)
            pairs = nussinov(seq)

            for i in range(0, len(seq)):
                opn = [True for t in pairs if t[0] == i]
                cls = [True for t in pairs if t[1] == i]

                if opn != []:
                    p += '('
                elif cls != []:
                    p += ')'
                else:
                    p += '-'

            print(p)
            print(seq)
