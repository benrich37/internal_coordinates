import numpy as np
import copy

radii_dict = {
    1: 0.31,
    2: 0.25,
    3: 1.28,
    4: 0.96,
    5: 0.84,
    6: 0.76,
    7: 0.71,
    8: 0.66
}


def is_bonded(posn1, posn2, num1, num2, margin=0.1):
    dist = abs(np.linalg.norm(posn2 - posn1))
    maxdist = radii_dict[num1] + radii_dict[num2]
    maxdist = maxdist + margin * maxdist
    return maxdist >= dist and dist != 0.


def atoms_is_bonded(atoms_obj, id1, id2):
    nums = atoms_obj.get_atomic_numbers()
    posns = atoms_obj.get_positions()
    return is_bonded(posns[id1], posns[id2], nums[id1], nums[id2])


def get_bonds(idx, posns, nums):
    n = len(nums)
    output = []
    main_posn = posns[idx]
    main_num = nums[idx]
    for i in range(n):
        if is_bonded(main_posn, posns[i], main_num, nums[i]):
            output.append(i)
    return output


def bonds_dict(atoms_obj):
    return_dict = {}
    nums = atoms_obj.get_atomic_numbers()
    posns = atoms_obj.get_positions()
    n = len(nums)
    for i in range(n):
        bonds = get_bonds(i, posns, nums)
        return_dict[i] = bonds
    return return_dict


def coordination_numbers(atoms_obj, sort=True):
    bond_dict = bonds_dict(atoms_obj)
    return_list_of_tuples = []
    for i in range(len(bond_dict)):
        return_list_of_tuples.append(tuple([i, len(bond_dict[i])]))
    if sort:
        return_list_of_tuples = sorted(return_list_of_tuples, key=lambda x: -x[1])
    return return_list_of_tuples


def get_bond_disp_vec(posns, idcs):
    posn1 = posns[idcs[0]]
    posn2 = posns[idcs[1]]
    bvec = posn2 - posn1
    posn2_disp = posn2 - (bvec / 2)
    posn1_disp = posn1 + (bvec / 2)
    posns_new = copy.copy(posns)
    posns_new[idcs[0]] = posn1_disp
    posns_new[idcs[1]] = posn2_disp
    posns_new = posns - posns_new
    return posns_new / np.linalg.norm(posns_new)


def stretch_posn_vec(posn_vec):
    n = len(posn_vec)
    output = np.empty(n * 3)
    for i in range(n):
        for j in range(3):
            output[3 * i + j] = posn_vec[i][j]
    return output


def bond_vec(posns, idcs):
    vec = get_bond_disp_vec(posns, idcs)
    return stretch_posn_vec(vec)


def get_unique_bonds(atoms_obj):
    bond_dict = bonds_dict(atoms_obj)
    return_dict = {}
    tmp_list = []
    for i in range(len(bond_dict)):
        return_dict[i] = []
        bonds = bond_dict[i]
        for a in bonds:
            if a < i:
                if not i in return_dict[a]:
                    return_dict[i].append(a)
            else:
                return_dict[i].append(a)
    return return_dict


def unique_bonds(atoms_obj):
    u_bond_dict = get_unique_bonds(atoms_obj)
    output = []
    for i in range(len(u_bond_dict)):
        for j in u_bond_dict[i]:
            output.append(tuple([i, j]))
    return output


def get_empty_vec(n):
    return np.zeros(n)


def bond_disp_vecs_sq(atoms_obj):
    posns = atoms_obj.get_positions()
    u_bonds = unique_bonds(atoms_obj)
    n_atoms = len(atoms_obj.get_atomic_numbers())
    m_bonds = len(u_bonds)
    output = np.empty(tuple([3 * n_atoms, 3 * n_atoms]))
    for m in range(m_bonds):
        output[m] = bond_vec(posns, [u_bonds[m][0], u_bonds[m][1]])
    for s in range(3 * n_atoms - m_bonds):
        output[m_bonds + s] = get_empty_vec(3 * n_atoms)
    return output


def bond_disp_vecs(atoms_obj):
    posns = atoms_obj.get_positions()
    u_bonds = unique_bonds(atoms_obj)
    n_atoms = len(atoms_obj.get_atomic_numbers())
    m_bonds = len(u_bonds)
    output = np.empty(tuple([m_bonds, 3 * n_atoms]))
    for m in range(m_bonds):
        output[m] = bond_vec(posns, [u_bonds[m][0], u_bonds[m][1]])
    return output