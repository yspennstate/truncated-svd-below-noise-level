# Truncating Singular Values Below the Noise Level

Source, PDF and simulation code for the paper "Truncating Singular Values Below the Noise Level" by Yitzchak Shmalo.

Let W = S + R, where S is an N x M matrix of fixed rank and R has independent N(0, 1/N) entries, with M/N tending to
c in (0, infinity), so that ||R|| tends to 1 + sqrt(c). Let W_T keep the singular values of W above a cutoff
T > 1 + sqrt(c) and set the others to zero. Each singular value theta > c^(1/4) of S produces a singular value of W near
y(theta) = sqrt((theta^2 + 1)(theta^2 + c))/theta. When the singular values of S are distinct and no theta > c^(1/4)
satisfies y(theta) = T, the almost sure limit of ||W_T - S||, which the argument of Gavish and Donoho gives from results
of Benaych-Georges and Nadakuditi, is the largest of several numbers, one for each singular value of S. We show that the
number belonging to a theta > c^(1/4) with y(theta) > T is less than 1 + sqrt(c), whatever its size. Hence the limit,
and the limits of ||(W_T - S)v|| and ||u^T (W_T - S)|| for the singular vectors u, v of S, are below 1 + sqrt(c) if and
only if every singular value theta >= 1 + sqrt(c) of S has y(theta) > T. For every cutoff in
(1 + sqrt(c), y(1 + sqrt(c))) the limit of ||W_T - S|| is below 1 + sqrt(c) for every such S, with a bound that does
not depend on S; for every cutoff above y(1 + sqrt(c)) some S of rank one gives a limit above 1 + sqrt(c). The
interval, (2, 5/2) for square matrices, contains the cutoff that Gavish and Donoho found optimal for the mean squared
error. For Gaussian noise this answers a question of Berlyand, Sandier, Shmalo and Zhang on truncating random weight
matrices.

## Contents

- `main.tex`, `main.pdf`: the paper.
- `code/simulate.py`: draws W for c = 0.25, 1, 4 and three cutoffs, and records the measured errors beside the limits.
  Usage: `python code/simulate.py N seeds out.tsv` (numpy).
- `code/table.py`: summarizes a draws file as in Table 1. Usage: `python code/table.py data/draws_N1000.tsv`.
- `data/draws_N1000.tsv` (40 draws per value of c) and `data/draws_N2000.tsv` (20 draws per value of c): the draws
  behind Table 1 and the N = 2000 comparison in Section 3.

## License

Text, figures and data CC BY 4.0; code MIT.
