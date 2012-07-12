def pcp(pairs):
    queue = []
    for v1, v2 in pairs:
        queue.insert(0, (1, v1, v2))

    while (len(queue) > 0):
        n, w1, w2 = queue.pop()
        for v1, v2 in pairs:
            s1 = w1 + v1
            s2 = w2 + v2
            if (s1 == s2):
                return (n+1, s1)
            elif s1.startswith(s2) or s2.startswith(s1):
                queue.insert(0, (n+1, s1, s2))

if __name__ == '__main__':
    pairs = [ ("001", "0"), ("01", "011"), ("01", "101"), ("10", "001") ]
    n, string = pcp(pairs)
    print("{}: {}".format(n, string))
