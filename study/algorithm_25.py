
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