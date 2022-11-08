def get_n_string(atom_list, i, m):
    output = []
    for j in range(m):
        output.append(atom_list[i+j])
    return tuple(output)

def get_n_strings(atom_list, n):
    output = []
    # number of times we can fit an n_string into our atom list
    m = len(atom_list) + 1 - n
    for i in range(m):
        output.append(get_n_string(atom_list, i, n))
# naive approach
def get_ic_idcs(atom_list):
    n = len(atom_list)
    bonds = get_n_strings(atom_list, 2)
    angles = get_n_strings(atom_list, 3)
    dihedrals = get_n_strings(atom_list, 4)
    return bonds, angles, dihedrals