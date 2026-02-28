from dataclasses import dataclass
import numpy as np

@dataclass
class Cursor:
    d: list         # Location of last match
    y: list         # Current value in sort order
    M: int          # Number of samples * 2 (biallel)
    N: int          # Number of sites

class PBWT:
    def __init__(self, records, order):
        self.records = np.array(records)
        self.M = len(records[0])
        self.N = len(records)
        self.order = order
        self.cursor = Cursor(
            M = self.M,
            N = self.N,
            d = [0] * self.M,
            y = self.records[:,order]
        )

    def pbwtCursorForwardsAD(self, k):
        """
        algorithm 2 in the manuscript
        src/pbwtCore.c:487
        """

        pass


    def matchMaximalWithin(self):
        """
        algorithm 4 in the paper
        src/pbwtMatch.c:115
        """

        m = 0
        n = 0

        for k in range(self.N + 1):
            for i in range(self.M):
                m = i-1
                n = i+1

                skip_i = False

                if (self.cursor.d[i] <= self.cursor.d[i+1]):
                    while (self.cursor.d[m+1] <= uself.cursor.d[i]):
                        if ((self.cursor.y[m] == self.cursor.y[i]) and (k < self.N)):
                            m -= 1
                            skip_i = True
                            break
                        m -= 1
    
                    if skip_i: continue
                
                if (self.cursor.d[i] >= self.cursor.d[i+1]):
                    while (self.cursor.d[n] <= self.cursor.d[i+1]):
                        if ((self.cursor.y[n] == self.cursor.y[i]) and (k < self.N)):
                            n += 1
                            skip_i = True
                            break
                        n += 1
                    
                    if skip_i: continue

            self.pbwtCursorForwardsAD(k)



if __name__ == "__main__":
    records = [
        [0,0,1,0,1],
        [0,1,1,0,0],
        [1,0,0,1,1]
    ]

    a = [1, 2, 0, 3, 4]

    pbwt = PBWT(records, a)
    print(pbwt.cursor)