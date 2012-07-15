from collections import deque

def pcp(pairs):
    queue = deque()
    for v1, v2 in pairs:
        queue.append((1, v1, v2))

    while (len(queue) > 0):
        n, w1, w2 = queue.popleft()
        for v1, v2 in pairs:
            s1 = w1 + v1
            s2 = w2 + v2
            if (s1 == s2):
                return (n+1, s1)
            elif s1.startswith(s2) or s2.startswith(s1):
                queue.append((n+1, s1, s2))

if __name__ == '__main__':
    pairs = [ ("001", "0"), ("01", "011"), ("01", "101"), ("10", "001") ]
    n, string = pcp(pairs)
    print("{}: {}".format(n, string))
