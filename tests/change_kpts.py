from __future__ import division
import sys

N = int(sys.argv[1])
print_w = False
if len(sys.argv) >= 3:
    print_w = True

for i in range(N):
    for j in range(N):
        for k in range(N):
            if print_w:
                print i/N, '  ', j/N, '  ', k/N, '  ', 1
            else:
                print i/N, '  ', j/N, '  ', k/N, '  '
