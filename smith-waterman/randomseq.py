#!/usr/bin/python

import sys
from random import choice, randint

letters='CSTPAGNDEQHRKMILVFYW'

if (len(sys.argv) < 4):
    print("usage: %s <count> <length_min> <length_max>" % sys.argv[0])
    sys.exit(1)

c = int(sys.argv[1])
mi = int(sys.argv[2])
ma = int(sys.argv[3])

for i in range(c):
    out = ''
    print(">random%d" % (i+1))
    for k in range(randint(mi, ma)):
        out += choice(letters)
    print(out)
