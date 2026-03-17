


def pbwtCursorForwardsA(record, a):
    # Algorithm 1 in manuscript
    # src/pbwtCore.c:460
    # Index to create *.pbwt file
    M = len(a)
    zeros = []
    ones = []
    
    for i in range(M):
        if record[a[i]] == 0:   # IMPORTANT
            zeros.append(a[i])
        else:
            ones.append(a[i])
    
    return zeros + ones



if __name__ == "__main__":
    records = [
        [0,0,1,0,1,1],
        [0,1,1,0,0,0],
        [1,0,0,1,1,0],
        [0,1,1,0,1,1]
    ]

    a = [0,1,2,3,4,5]

    result = []
    for record in records:
        a = pbwtCursorForwardsA(record, a)
        result.append(a.copy())

    print(result)
