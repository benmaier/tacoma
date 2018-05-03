import tacoma as tc

# This is a test whether or not the contact persistence distribution is rightfully
# adjusted compared to the pair persistence distribution

L = tc.edge_lists()
L.N = 5
L.t = [0., 1., 2., 3. ]
L.tmax = 4.
L.edges = [ 
            [],
            [
                (0, 1), (2,3), (3, 4), (2, 4),
            ],
            [
                (2,3), (3, 4), (2, 4),
            ],
            [],
        ]

C = tc.convert(L)

L_result = tc.measure_group_sizes_and_durations(L)
C_result = tc.measure_group_sizes_and_durations(C)

for res in [L_result, C_result]:
    contact_durations = res.contact_durations
    pair_durations = res.group_durations[2]

    print("contact durations =", contact_durations)
    print("(should be [ 1.0, 2.0, 2.0, 2.0])")
    print("pair durations =", pair_durations)
    print("(should be [ 1.0 ])")
    print()


