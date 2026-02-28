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


    def matchMaximalWithin(self):
        """
        algorithm 4 in the paper
        """



        pass



if __name__ == "__main__":
    records = [
        [0,0,1,0,1],
        [0,1,1,0,0],
        [1,0,0,1,1]
    ]

    a = [1, 2, 0, 3, 4]

    pbwt = PBWT(records, a)
    print(pbwt.cursor)