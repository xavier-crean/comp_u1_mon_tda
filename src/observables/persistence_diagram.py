import numpy as np
import itertools
from gtda.homology import CubicalPersistence
import h5py
from multiprocessing import Pool
from sys import argv
from configuration import *

"""
On command line: 
argv[1] = Ls, argv[2] = beta, argv[3] = N, argv[4] = sweep_step, argv[5] = NP
e.g. 'py persistence_diagram.py 6 0.9000 200 350000 5'
"""

PLANE_LABELS = {(0, 1): 0, (0, 2): 1, (0, 3): 2, (1, 2): 3, (1, 3): 4, (2, 3): 5}


def diracStringNum(plaq_ang):
    """Take a plaquette angle and return the number of Dirac strings passing through as int"""
    a1 = [plaq_ang > -4 * np.pi, plaq_ang <= -3 * np.pi]
    a2 = [plaq_ang > -3 * np.pi, plaq_ang <= -1 * np.pi]
    a3 = [plaq_ang > -1 * np.pi, plaq_ang <= np.pi]
    a4 = [plaq_ang > np.pi, plaq_ang <= 3 * np.pi]
    a5 = [plaq_ang > 3 * np.pi, plaq_ang <= 4 * np.pi]
    if all(a1):
        return -2
    elif all(a2):
        return -1
    elif all(a3):
        return 0
    elif all(a4):
        return 1
    elif all(a5):
        return 2
    else:
        raise Exception("Plaq value too large/small.")


def levi_cevita(dim):
    """Generates Levi-Cevita symbol for dimensions dim"""
    arr = np.zeros(tuple([dim for _ in range(dim)]))
    for x in itertools.permutations(tuple(range(dim))):
        mat = np.zeros((dim, dim), dtype=np.int32)
        for i, j in zip(range(dim), x):
            mat[i, j] = 1
        arr[x] = int(np.linalg.det(mat))
    return arr


def mon(plang):
    """
    Input: array (dtype=float64) of plaquette angles
    Return: array (dtype=int32) of monopole currents on dual lattice
    """

    # Dirac string array
    n = np.zeros(plang.shape[:4] + (4, 4), dtype=np.int32)
    for t, x, y, z in itertools.product(*[range(s) for s in plang.shape[:4]]):
        for k, v in PLANE_LABELS.items():
            n[t, x, y, z, k[0], k[1]] = diracStringNum(plang[t, x, y, z, v])

    # monopole 3-cube array
    levi = levi_cevita(4)
    M = np.zeros(plang.shape[:4] + (4,), dtype=np.int32)
    for t, x, y, z in itertools.product(*[range(s) for s in plang.shape[:4]]):
        # free index
        for r in range(4):
            total = 0
            for s in range(4):
                if s == 0:
                    ds = [1, 0, 0, 0]
                elif s == 1:
                    ds = [0, 1, 0, 0]
                elif s == 2:
                    ds = [0, 0, 1, 0]
                elif s == 3:
                    ds = [0, 0, 0, 1]
                n_diff = (
                    n[
                        (t + ds[0]) % (plang.shape[0]),
                        (x + ds[1]) % (plang.shape[0]),
                        (y + ds[2]) % (plang.shape[0]),
                        (z + ds[3]) % (plang.shape[0]),
                    ]
                    - n[t, x, y, z]
                )
                total += np.einsum("mn,mn->", levi[r, s], n_diff)
            M[t, x, y, z][r] = total

    # monopole dual 1-cube array
    j = np.zeros(plang.shape[:4] + (4,), dtype=np.int32)
    for t, x, y, z in itertools.product(*[range(s) for s in plang.shape[:4]]):
        for r in range(4):
            if r == 0:
                dr = [1, 0, 0, 0]
            elif r == 1:
                dr = [0, 1, 0, 0]
            elif r == 2:
                dr = [0, 0, 1, 0]
            elif r == 3:
                dr = [0, 0, 0, 1]
            j[t, x, y, z][r] = M[
                (t + dr[0]) % (plang.shape[0]),
                (x + dr[1]) % (plang.shape[0]),
                (y + dr[2]) % (plang.shape[0]),
                (z + dr[3]) % (plang.shape[0]),
            ][r]
    return j


# +- unit vector in each direction
pm = np.array(
    [
        [[1, 0, 0, 0], [-1, 0, 0, 0]],
        [[0, 1, 0, 0], [0, -1, 0, 0]],
        [[0, 0, 1, 0], [0, 0, -1, 0]],
        [[0, 0, 0, 1], [0, 0, 0, -1]],
    ]
)


def cubicalFiltration(m):
    """
    Trivial filtration:
        s=-1: non-zero current lines and their boundary vertices
        s=0: all other d-cubes enter completing the 4-torus
    """
    filt = np.zeros(tuple([2 * s for s in m.shape[:4]]))

    for t, x, y, z in itertools.product(*[range(2 * s) for s in m.shape[:4]]):
        """1-cubes"""
        dim = t % 2 + x % 2 + y % 2 + z % 2
        if dim == 1:
            parDirs = tuple([i for i in range(4) if [t, x, y, z][i] % 2 == 1])
            j = m[
                (t - (t % 2)) // 2,
                (x - (x % 2)) // 2,
                (y - (y % 2)) // 2,
                (z - (z % 2)) // 2,
                parDirs[0],
            ]
            if j != 0:
                filt[t, x, y, z] = -1

    for t, x, y, z in itertools.product(*[range(2 * s) for s in m.shape[:4]]):
        """
        place 0-cubes in the correct location in filtered cubical complex array
        """
        dim = t % 2 + x % 2 + y % 2 + z % 2
        if dim == 0:
            dirs = np.concatenate(pm[[i for i in range(4) if [t, x, y, z][i] % 2 == 0]])
            filt[t, x, y, z] = min(
                [
                    filt[
                        (t + dt) % (2 * m.shape[0]),
                        (x + dx) % (2 * m.shape[1]),
                        (y + dy) % (2 * m.shape[2]),
                        (z + dz) % (2 * m.shape[3]),
                    ]
                    for [dt, dx, dy, dz] in dirs
                ]
            )
    return filt


def main():
    # Set lattice size
    Ls = argv[1]

    # Set beta
    b = argv[2]

    # Number of sampled configurations
    N = argv[3]

    # Set number of sweeps between measurement in the Markov chain
    sweep_step = argv[4]

    # Markov chain recording length
    rec_len = N * sweep_step

    # Define number of parallel computations
    NP = argv[5]

    # Homology is computed via a Cubical Persistence method from giotto-tda
    # We set
    # - periodic boundary conditions in t,x,y,z
    # - Z/2Z coefficient field
    # - homology dimensions to [0,1] only
    cp = CubicalPersistence(
        homology_dimensions=[0, 1],
        coeff=2,
        periodic_dimensions=np.array([True, True, True, True]),
        reduced_homology=False,
        infinity_values=np.inf,
        n_jobs=NP,
    )

    confs_b = []
    for r in range(sweep_step, (rec_len + 1), sweep_step):
        print(r)
        fn = f"data/configurations/{Ls}.{Ls}.{Ls}.{Ls}/{b:.4f}/conf.dat{r}"
        with open(fn, "rb") as f:
            HL = len(f.readline())
            f.seek(0, 0)
            confs_b.append(parseConfig(f.read(), HL, Ls))
    with Pool(NP) as p:
        plaqs = p.map(plaquettes, confs_b)
        mons = p.map(mon, plaqs)
        filts = p.map(cubicalFiltration, mons)
        pers = cp.fit_transform(filts)

    filename = f"data/observables/pd/{Ls}.{Ls}.{Ls}.{Ls}/pers_mon_Ns={Ls}{Ls}{Ls}{Ls}_b={b:.4f}.h5"
    with h5py.File(filename, "w") as hf:
        hf.create_dataset("persistence", data=pers)


if __name__ == "__main__":
    main()
