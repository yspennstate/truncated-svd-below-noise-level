# Truncating Singular Values Below the Noise Level

Source, PDF and simulation code for the paper "Truncating Singular Values Below the Noise Level" by Yitzchak Shmalo.

Let W = S + R, where S is an N x M matrix of fixed rank and R has independent N(0, 1/N) entries, with M/N tending to
c in (0, infinity), so that ||R|| tends to 1 + sqrt(c). Let W_T keep the singular values of W above a cutoff
T > 1 + sqrt(c) and set the others to zero. The almost sure limit of ||W_T - S|| follows from results of
Benaych-Georges and Nadakuditi and of Gavish and Donoho. We show that a singular value of S that produces a singular
value of W above T contributes less than 1 + sqrt(c) to this limit, whatever its size. Hence the limit, and the limits
of ||(W_T - S)v|| and ||u^T (W_T - S)|| for the singular vectors u, v of S, are below 1 + sqrt(c) if and only if every
singular value of S that is at least 1 + sqrt(c) produces a singular value of W above T. For every cutoff in
(1 + sqrt(c), y(1 + sqrt(c))), where y(t) = sqrt((t^2 + 1)(t^2 + c))/t, the limit of ||W_T - S||, where it exists, is
below 1 + sqrt(c) for every S, with a bound that does not depend on S; for every cutoff above y(1 + sqrt(c)) some S of
rank one gives a larger limit. The interval, (2, 5/2) for square matrices, contains the cutoff that Gavish and Donoho
found optimal for the mean squared error. For Gaussian noise this answers a question of Berlyand, Sandier, Shmalo and
Zhang on pruning random weights.

## Contents

- `main.tex`, `main.pdf`: the paper.
- `code/simulate.py`: draws W for c = 0.25, 1, 4 and three cutoffs, and records the measured errors beside the limits.
  Usage: `python code/simulate.py N seeds out.tsv` (numpy).
- `code/table.py`: summarizes a draws file as in Table 1. Usage: `python code/table.py data/draws_N1000.tsv`.
- `data/draws_N1000.tsv` (40 draws per value of c) and `data/draws_N2000.tsv` (20 draws per value of c): the draws
  behind Table 1 and the N = 2000 comparison in Section 3.

## License

Text, figures and data CC BY 4.0; code MIT.
