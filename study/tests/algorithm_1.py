


def pbwtCursorForwardsA(record, a):
    # src/pbwtCore.c:460
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
        [0,0,1,0,1],
        [0,1,1,0,0],
        [1,0,0,1,1]
    ]

    a = [0,1,2,3,4]

    result = []
    for record in records:
        a = pbwtCursorForwardsA(record, a)
        result.append(a.copy())

    print(result)
