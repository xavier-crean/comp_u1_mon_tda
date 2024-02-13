from struct import unpack
import itertools
import numpy as np

PLANE_LABELS = {(0, 1): 0, (0, 2): 1, (0, 3): 2, (1, 2): 3, (1, 3): 4, (2, 3): 5}


def parseConfig(config_buffer, HL, L):
    def cNumFromBytes(buffer):
        cnum = np.array(unpack(">dd", buffer))
        return cnum[0] + 1j * cnum[1]

    HEAD_LEN, DTYPE_LEN = HL, 16
    STDIM, DSIZE_T, DSIZE_X, DSIZE_Y, DSIZE_Z = 4, L, L, L, L
    config_array = np.array(
        [
            [
                [
                    [
                        [
                            cNumFromBytes(
                                config_buffer[
                                    HEAD_LEN
                                    + DTYPE_LEN
                                    * (
                                        STDIM
                                        * (
                                            t
                                            + (x * DSIZE_T)
                                            + (y * DSIZE_T * DSIZE_X)
                                            + (z * DSIZE_T * DSIZE_X * DSIZE_Y)
                                        )
                                        + link
                                    ) : HEAD_LEN
                                    + DTYPE_LEN
                                    * (
                                        STDIM
                                        * (
                                            t
                                            + (x * DSIZE_T)
                                            + (y * DSIZE_T * DSIZE_X)
                                            + (z * DSIZE_T * DSIZE_X * DSIZE_Y)
                                        )
                                        + link
                                        + 1
                                    )
                                ]
                            )
                            for link in range(STDIM)
                        ]
                        for z in range(DSIZE_Z)
                    ]
                    for y in range(DSIZE_Y)
                ]
                for x in range(DSIZE_X)
            ]
            for t in range(DSIZE_T)
        ]
    )
    return config_array


def wilsonLoop(conf, t, x, y, z, directions):
    d_mu = [[1 if i == d else 0 for i in range(4)] for d in directions]
    confang = np.angle(conf)
    ang = (
        confang[t, x, y, z, directions[0]]
        + confang[
            (t + d_mu[0][0]) % conf.shape[0],
            (x + d_mu[0][1]) % conf.shape[1],
            (y + d_mu[0][2]) % conf.shape[2],
            (z + d_mu[0][3]) % conf.shape[3],
            directions[1],
        ]
        - confang[
            (t + d_mu[1][0]) % conf.shape[0],
            (x + d_mu[1][1]) % conf.shape[1],
            (y + d_mu[1][2]) % conf.shape[2],
            (z + d_mu[1][3]) % conf.shape[3],
            directions[0],
        ]
        - confang[t, x, y, z, directions[1]]
    )
    return np.array(ang)


def plaquettes(conf):
    plaq = np.zeros(conf.shape[:4] + (6,), dtype=np.double)
    for t, x, y, z in itertools.product(*[range(s) for s in conf.shape[:4]]):
        for k, v in PLANE_LABELS.items():
            plaq[t, x, y, z, v] = wilsonLoop(conf, t, x, y, z, k)
    return plaq
