from dataclasses import dataclass
import numpy as np

@dataclass
class Cursor:
    d: list         # Location of last match
    y: list         # Current value in sort order
    a: list         # Index back to original order
    M: int          # Number of samples * 2 (biallel)
    N: int          # Number of sites
    c: int          # number of 0s in y
    u: list         # number of 0s up to and including this position

class PBWT:
    def __init__(self, base_records, base_order, query_records):
        self.base_records = np.array(base_records)
        self.query_records = np.array(query_records)
        self.M = len(base_records[0])
        self.N = len(base_records)
        self.cursor = Cursor(
            M = self.M,
            N = self.N,
            d = [0] * (self.M + 1),
            y = [0] * self.M,
            a = base_order,
            c = 0,
            u = [0] * (self.M + 1)
        )

    def pbwtCursorCalculateU(self, x):
        """
        src/pbwtCore.c:512
        Count number of 0s to now and redefine number of 0s in y
        """

        u = 0
        for i in range(self.M):
            x.u[i] = u 
            if (x.y[i] == 0): 
                u += 1
        
        # need one off the end of update intervals
        x.c = u
        x.u[-1] = u

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

    def matchSequencesIndexed(self):
        """
        algorithm 5 in the article
        src/pbwtMatch.c:255
        """

        up = self.cursor
        report = []
        totLen = 0
        nTot = 0

        ## build indexes
        a = np.zeros((self.N+1, self.M), dtype=int)
        d = np.zeros((self.N+1, self.M+1), dtype=int)
        u = np.zeros((self.N, self.M+1), dtype=int)
        cc = [0] * self.N
        
        for k in range(0, self.N):
            up.y = self.base_records[k][up.a]

            a[k,:] = up.a
            d[k,:] = up.d
            cc[k] = (self.base_records[k] == 0).sum()       # Count number of 0s
            self.pbwtCursorCalculateU(up)
            u[k,:] = up.u
            self.pbwtCursorForwardsAD(up, k)
        
        a[self.N] = up.a
        d[self.N] = up.d

        ## match each query in turn
        for j in range(len(self.query_records[0])):
            x = self.query_records[:,j]

            ## start of match, and pbwt interval as in algorithm 5
            e = 0
            f = 0
            g = self.M - 1

            ## next versions of the above, e' etc in algorithm 5
            e1 = 0
            f1 = 0 
            g1 = 0

            for k in range(self.N):     # use classic FM updates to extend [f,g) interval to next position
                num_one_f = (f + 1 - u[k][f+1])
                num_one_f_start = num_one_f if num_one_f > 0 else (num_one_f + 1)

                num_one_g = (g + 1 - u[k][g+1])
                num_one_g_start = num_one_g if num_one_g > 0 else (num_one_g + 1)

                f1 = (cc[k] + num_one_f_start if x[k] else u[k][f+1]) - 1
                g1 = (cc[k] + num_one_g_start if x[k] else u[k][g+1]) - 1

                # if the interval is non-zero we can just proceed
                if g1 > f1:
                    f = f1
                    g = g1
                else:   ## we have reached a maximum | src/pbwtMatch.c:304 
                    for i in range(f, g):
                        report.append([j, a[k][i], e, k])
                    
                    nTot += 1
                    totLen += k-e 

                    e1 = d[k+1][f1] - 1   # y[f1] and y[f1-1] diverge here, so upper bound for e

                    if (f1 == self.M) or ((f1 > 0) and (x[e1] == 0)):
                        f1 = g1 - 1
                        y = self.base_records[:,a[k+1][f1]]
                        while e1 > 0 and x[e1-1] == y[e1-1]:
                            e1 -= 1
                        
                        while f1 > 0 and d[k+1][f1] <= e1:
                            f1 -= 1
                    elif (f1 < self.M):
                        g1 = f1 + 1
                        y = self.base_records[:,a[k+1][f1]]
                        while e1 > 0 and x[e1-1] == y[e1-1]:
                            e1 -= 1
                        
                        while ((g1 < self.M) and (d[k+1][g1] <= e1)):
                            g1 += 1

                    e = e1
                    f = f1
                    g = g1

            ## report the maximal matches to the end
            last_k = self.N - 1
            for i in range(f, g):
                report.append([j, a[last_k][i], e, last_k])
            
            nTot += 1
            totLen += k-e

        return report, nTot/len(self.query_records[0]), totLen/nTot


if __name__ == "__main__":
    base_records = [                     # 6 haps x 3 sites
        [0,0,1,0,1,1],
        [0,1,1,0,0,0],
        [1,0,0,1,1,0]
    ]
    base_order = [0, 1, 2, 3, 4, 5]      # Read from the pbwt file (applied the algorithm 1)

    query_records = [                    # 3 haps x 3 sites
        [1,0,1],
        [1,1,0],
        [0,1,0]
    ]

    report, avg_best_matches, avg_length = PBWT(base_records, base_order, query_records).matchSequencesIndexed()

    print(report)
    print(f"Average number of best matches: {avg_best_matches}")
    print(f"Average length: {avg_length}")