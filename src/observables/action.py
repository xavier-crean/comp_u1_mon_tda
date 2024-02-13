import numpy as np
import h5py
from sys import argv
from configuration import *

"""
On command line: 
argv[1] = Ls, argv[2] = beta, argv[3] = N, argv[4] = sweep_step
e.g. 'py action.py 6 0.9000 200 350000'
"""


def action(plaq):
    """Compute total action for a configuration, given an array of plaquettes shape (4,4,4,4,6)"""
    return np.sum(np.cos(plaq))


def main():
    # Set lattice size
    Ls = argv[1]

    # Set beta
    b = float(argv[2])

    # Number of sampled configurations
    N = argv[3]

    # Set number of sweeps between measurement in the Markov chain
    sweep_step = argv[4]

    # Markov chain recording length
    rec_len = N * sweep_step

    confs_b = []
    for r in range(sweep_step, (rec_len + 1), sweep_step):
        fn = f"data/configurations/{Ls}.{Ls}.{Ls}.{Ls}/{b:.4f}/conf.dat{r}"
        with open(fn, "rb") as f:
            HL = len(f.readline())
            f.seek(0, 0)
            confs_b.append(parseConfig(f.read(), HL, Ls))
    plaqs = map(plaquettes, confs_b)
    act = np.array(list(map(action, plaqs)))

    filename = f"data/observables/action/{Ls}.{Ls}.{Ls}.{Ls}/action_mon_Ns={Ls}{Ls}{Ls}{Ls}_b={b:.4f}.h5"
    with h5py.File(filename, "w") as hf:
        hf.create_dataset("action", data=act)


if __name__ == "__main__":
    main()
