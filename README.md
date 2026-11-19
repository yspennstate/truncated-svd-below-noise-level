# Truncating singular values below the noise level

Source, PDF and simulation code for the paper "Truncating singular values below the noise level" by Yitzchak Shmalo.

Let W = S + R, where S is an N x M matrix of fixed rank and R has independent N(0, 1/N) entries, with M/N tending to c.
Keep the singular values of W above a cutoff T > 1 + sqrt(c) and set the others to zero. The paper computes the almost
sure limit of the operator-norm error and of the error along each singular vector of S, and shows that every cutoff in
the interval (1 + sqrt(c), y(1 + sqrt(c))), where y(t) = sqrt((t^2+1)(t^2+c))/t, gives an error below 1 + sqrt(c) for
every S. For square matrices the interval is (2, 5/2).

## Contents

- `main.tex`, `main.pdf`: the paper.
- `code/simulate.py`: draws W for c = 0.25, 1, 4 and three cutoffs, and records the measured errors beside the limits.
  Usage: `python code/simulate.py N seeds out.tsv` (numpy).
- `code/table.py`: summarizes a draws file as in Table 1. Usage: `python code/table.py data/draws_N1000.tsv`.
- `data/draws_N1000.tsv` (40 draws per value of c) and `data/draws_N2000.tsv` (20 draws per value of c): the draws
  behind Table 1 and the N = 2000 comparison in Section 3.

## License

Text, figures and data CC BY 4.0; code MIT.
