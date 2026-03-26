from dataclasses import dataclass
import numpy as np
import pandas as pd

@dataclass
class Cursor:
    d: list         # Location of last match
    y: list         # Current value in sort order
    a: list         # Index back to original order
    M: int          # Number of samples * 2 (biallel)
    N: int          # Number of sites

class PBWT:
    def __init__(self, records, order):
        self.records = np.array(records)
        self.M = len(records[0])
        self.N = len(records)
        self.cursor = Cursor(
            M = self.M,
            N = self.N,
            d = [0] * (self.M + 1),
            y = [0] * self.M,
            a = order
        )

    def pbwtCursorForwardsAD(self, x, k):
        """
        algorithm 2 in the manuscript
        src/pbwtCore.c:487
        """

        p = k + 1
        q = k + 1

        b = [0] * x.M
        e = [0] * x.M
        
        u = 0
        v = 0

        for i in range(x.M):

            if x.d[i] > p:
                p = x.d[i]
            if x.d[i] > q:
                q = x.d[i]

            if x.y[i] == 0:         # NB x[a[i]] = y[i] in manuscript
                x.a[u] = x.a[i]
                x.d[u] = p
                u += 1
                p = 0
            else:                   # y[i] == 1, since bi-allelic
                b[v] = x.a[i]       
                e[v] = q
                v += 1
                q = 0

        # copy 1-alleles to tail (memcpy equivalent)
        for i in range(v):
            x.a[u + i] = b[i]
            x.d[u + i] = e[i]

        # sentinels
        x.d[0] = k + 2
        x.d[x.M] = k + 2

    def matchLongWithin2(self, T):
        # src/pbwtMatch.c:85
        i0 = 0
        na = 0
        nb = 0

        ia = 0

        u = self.cursor
        report = []

        for k in range(self.N):

            u.y = self.records[k][u.a]

            for i in range(self.M):
                if (u.d[i] > k-T):
                    if na and nb:   # then there is something to report
                        for _ in range(i0, i):
                            ia += 1
                    
                    dmin = 0
                    for ib in range(ia+1, i):
                        if (u.d[ib] > dmin): 
                            dmin = u.d[ib]
                        if (u.y[ib] != u.y[ia]):
                            report.append([u.a[ia], u.a[ib], dmin, k])
                    
                    na = 0
                    nb = 0
                    i0 = i
                
                if (u.y[i] == 0):
                    na += 1
                else:
                    nb += 1
            
            self.pbwtCursorForwardsAD(u, k)


        return report

if __name__ == "__main__":
    records = [                     # 6 haps x 3 sites
        [0,0,1,0,1,1],
        [0,1,1,0,0,0],
        [1,0,0,1,1,0],
        [0,1,1,0,1,1]
    ]

    order = [5, 1, 2, 0, 3, 4]         # Read from the pbwt file (applied the algorithm 1)

    pbwt = PBWT(records, order)

    T = 2
    report = pbwt.matchLongWithin2(T)

    print(report)