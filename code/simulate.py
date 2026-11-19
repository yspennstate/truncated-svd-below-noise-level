"""Finite-N check of the truncation limits.

W = S + R, where R is N x M with independent N(0, 1/N) entries, M/N = c, and S has rank 3 with singular values
3e, 1.15e, 0.9e (e = 1 + sqrt(c)) and random orthonormal singular vectors drawn afresh in each draw. W is truncated
at a cutoff T (singular triples with singular value > T are kept). For each draw the script records the operator
norm of W_T - S, the norm of (W_T - S) v_2, and the limits predicted by the paper: y(t) = sqrt((t^2+1)(t^2+c))/t,
a^2 = (1-c/t^4)/(1+1/t^2), b^2 = (1-c/t^4)/(1+c/t^2), and the largest singular value of the 2x2 block
y (a, sqrt(1-a^2))^T (b, sqrt(1-b^2)) - diag(t, 0) over the kept spikes, together with t over the dropped ones.

usage: python simulate.py N seeds out.tsv
"""
import sys
import numpy as np


def pred(t, c):
    if t ** 4 <= c:
        return None
    y = np.sqrt((t * t + 1) * (t * t + c)) / t
    a = np.sqrt((1 - c / t ** 4) / (1 + 1 / t ** 2))
    b = np.sqrt((1 - c / t ** 4) / (1 + c / t ** 2))
    A = y * np.outer([a, np.sqrt(1 - a * a)], [b, np.sqrt(1 - b * b)]) - np.diag([t, 0.0])
    return y, a * a, b * b, np.linalg.norm(A, 2), np.sqrt(t * t + y * y * b * b - 2 * t * y * a * b)


N, SEEDS, OUT = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
cols = ["c", "T_name", "T", "seed", "t1", "t2", "t3", "kept", "ret1", "ret2", "ret3", "loss", "e", "max_discarded_t",
        "claim_lhs", "claim_rhs", "agree", "pred_op", "meas_v1", "pred_v1", "meas_v2", "pred_v2", "meas_v3", "pred_v3",
        "meas_y", "pred_y", "meas_a2", "pred_a2", "meas_b2", "pred_b2"]
rows = ["\t".join(cols)]
nan = float("nan")
for c in (0.25, 1.0, 4.0):
    M = int(round(c * N))
    e = 1 + np.sqrt(c)
    t = np.array([0.9 * e, 1.15 * e, 3.0 * e])
    pr = [pred(t[i], c) for i in range(3)]
    for seed in range(SEEDS):
        rng = np.random.default_rng(1000 * seed + int(100 * c))
        U, _ = np.linalg.qr(rng.standard_normal((N, 3)))
        V, _ = np.linalg.qr(rng.standard_normal((M, 3)))
        S = (U * t) @ V.T
        W = S + rng.standard_normal((N, M)) / np.sqrt(N)
        P, sig, Qt = np.linalg.svd(W, full_matrices=False)
        for name, T in (("edge+0.1", e + 0.1), ("above_low", sig[2] + 1e-6), ("above_mid", sig[1] + 1e-6)):
            k = int(np.sum(sig > T))
            Wk = (P[:, :k] * sig[:k]) @ Qt[:k]
            D = Wk - S
            loss = np.linalg.norm(D, 2)
            ret = [float(np.sum((P[:, :k].T @ U[:, i]) ** 2)) > 0.5 for i in range(3)]
            disc = [float(t[i]) for i in range(3) if not ret[i]]
            lhs = loss < e
            rhs = all(ret[i] for i in range(3) if t[i] >= e)
            pop = max([pr[i][3] for i in range(3) if ret[i] and pr[i]] + disc + [0.0])
            ex = []
            for i in range(3):
                ex += [np.linalg.norm(D @ V[:, i]) if ret[i] else nan, pr[i][4] if (ret[i] and pr[i]) else nan]
            ex += [sig[0], pr[2][0], float((P[:, 0] @ U[:, 2]) ** 2), pr[2][1], float((Qt[0] @ V[:, 2]) ** 2), pr[2][2]]
            rows.append("%g\t%s\t%.6f\t%d\t%.4f\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%d\t%d\t%d\t%.6f" % (
                c, name, T, seed, t[0], t[1], t[2], k, ret[0], ret[1], ret[2], loss, e, max(disc) if disc else 0.0,
                lhs, rhs, lhs == rhs, pop) + "".join("\t%.6f" % x for x in ex))
        print("c", c, "seed", seed, "done", flush=True)
open(OUT, "w").write("\n".join(rows) + "\n")
body = [r.split("\t") for r in rows[1:]]
print("rows", len(body), "claim agrees on", sum(int(b[16]) for b in body), "of", len(body))
