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
            a = order
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

    def matchSequencesIndexed(self):
        """
        algorithm in the article
        src/pbwtMatch.c:255
        """

        up = self.cursor

        ## stored indexes
        a = np.empty((self.N+1, self.M), dtype=int)
        d = np.empty((self.N+1, self.M+1), dtype=int)
        u = np.empty((self.N, self.M+1), dtype=int)
        cc = [0] * self.N
        
        for k in range(0, self.N):
            up.y = self.records[k][up.a]

            a[k,:] = up.a
            d[k,:] = up.d
            cc[k] = (self.base_records[k] == 0).sum()       # Count number of 0s
            self.pbwtCursorCalculateU(up)

        pass


if __name__ == "__main__":
    base_records = [                     # 5 haps x 3 sites
        [0,0,1,0,1],
        [0,1,1,0,0],
        [1,0,0,1,1]
    ]
    base_order = [1, 2, 0, 3, 4]         # Read from the pbwt file (applied the algorithm 1)

    query_records = [                    # 3 haps x 3 sites
        [1,0,1],
        [1,1,0],
        [0,1,0]
    ]
    query_order = [2, 0, 1]              # Read from the pbwt file (applied the algorithm 1)