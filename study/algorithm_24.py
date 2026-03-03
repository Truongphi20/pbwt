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


    def matchMaximalWithin(self):
        """
        algorithm 4 in the paper
        src/pbwtMatch.c:115
        """
        x = self.cursor
        report = []

        for k in range(self.N):

            # IMPORTANT: y must reflect current site
            x.y = self.records[k][x.a]

            for i in range(x.M):

                m = i - 1
                n = i + 1

                # ---- left scan ----
                if x.d[i] <= x.d[i + 1]:
                    while m >= 0 and x.d[m + 1] <= x.d[i]:

                        if k < x.N - 1 and x.y[m] == x.y[i]:
                            break
                        m -= 1

                # ---- right scan ----
                if x.d[i] >= x.d[i + 1]:
                    while n < x.M and x.d[n] <= x.d[i + 1]:

                        if k < x.N - 1 and x.y[n] == x.y[i]:
                            break
                        n += 1

                # (Here is where matches would be reported)
                for j in range(m+1, i):
                    report.append([x.a[i], x.a[j], x.d[i], k])
                
                for j in range(i+1, n):
                    report.append([x.a[i], x.a[j], x.d[i+1], k])


            # move PBWT forward
            self.pbwtCursorForwardsAD(x, k)
        
        return pd.DataFrame(
                report, 
                columns=["cor_hap", "match_hap", "start", "length"]
            ).sort_values("cor_hap")\
            .reset_index(drop=True)



if __name__ == "__main__":
    records = [                     # 5 haps x 3 sites
        [0,0,1,0,1],
        [0,1,1,0,0],
        [1,0,0,1,1]
    ]

    order = [1, 2, 0, 3, 4]         # Read from the pbwt file (applied the algorithm 1)

    pbwt = PBWT(records, order)
    report = pbwt.matchMaximalWithin()

    print("Final permutation a:", pbwt.cursor.a)
    print("Final divergence d:", pbwt.cursor.d)
    
    print(report.loc[report["length"]!=0])

    