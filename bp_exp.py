import numpy as np
from bp import BeliefPropagation, TannerGraph
from typing import Optional
import functools

from tqdm import tqdm
from multiprocessing import Pool


def llr_mean(snr: float):
    return 4 * 0.5 * 10 ** (snr / 10)


def llr_var(snr: float):
    return 8 * 0.5 * 10 ** (snr / 10)


def log_density(x: float, snr: float):
    return -((x - llr_mean(snr)) ** 2) / (2 * llr_var(snr))


def awgn_llr(snr: float):
    return lambda y: log_density(-y, snr) - log_density(y, snr)


def BP(H, codeword, snr, max_iter=10) -> Optional[np.ndarray]:
    tg = TannerGraph.from_biadjacency_matrix(H, channel_model=awgn_llr(snr))
    bp = BeliefPropagation(tg, H, max_iter=max_iter)
    estimate, _, decode_success = bp.decode(codeword)
    if not decode_success:
        return None
    return estimate


def transmit(c, snr):
    var = llr_var(snr)
    mean = llr_mean(snr)
    std = np.sqrt(var)
    return (c * 2 - 1) * np.random.normal(mean, std, (len(c),))


def read_H(filepath):
    with open(filepath) as f:
        lines = f.readlines()
    n, m = map(int, lines[0].split())
    H = np.zeros((n, m))
    for i, l in enumerate(lines[4:4 + n]):
        for j in map(int, l.split()):
            if j == 0:
                continue
            H[i, j - 1] = 1
    return H.T


def read_codewords(filepath):
    with open(filepath) as f:
        lines = f.readlines()
    n = int(lines[0])
    codewords = []
    for l in lines[1:]:
        codeword = []
        for c in l.strip():
            codeword.append(int(c))
        codewords.append(np.array(codeword))
    return codewords


def test(c, H, snr):
    tr = transmit(c, snr)
    # print(c)
    # print((tr > 0).astype(int) - c)
    # print()
    res = BP(H, tr, snr)
    # print('  Input:', ''.join([str(f) for f in c]))
    # if res is not None:
    #     print('Decoded:', ''.join([str(f) for f in res]))
    #     print()
    # else:
    #     print('Not found')
    return res is not None


if __name__ == '__main__':
    np.random.seed(239)

    snr = 1.0
    H = read_H('H.txt')
    codewords = read_codewords('codewords.txt')

    with Pool(26) as p:
        res = list(tqdm(p.imap(functools.partial(test, H=H, snr=snr), codewords), total=len(codewords)))
    print(sum(res) / len(res))
