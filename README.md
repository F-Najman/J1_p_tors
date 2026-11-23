# J1_p_tors

Magma code for the computations in Section 6 of

> F. Najman, *The rational torsion subgroup of J_1(p)*.

The proofs in the paper do not use any computer calculations. The computations here are consistency
checks of Theorem 1.3 (injectivity of reduction modulo 2), Proposition 4.2 (the index of the Eisenstein
ideal at 2) and Corollary 5.3, for the primes 11 <= p <= 139. The conjecture of Conrad, Edixhoven and
Stein was already known for p <= 157, so these computations do not prove anything new.

## Files

| file | what it computes | log |
|---|---|---|
| `orders.m` | the order `c(p)` of the rational cuspidal subgroup and `#J_1(p)(F_2)` for all primes `11 <= p <= 139` | `logs/orders.log` |
| `two_rank.m` | `dim_{F_2} J_1(p)(F_2)[2]` for the 12 primes of Table 1 | `logs/two_rank.log` |
| `congruence_module.m` | the 2-primary part of Ohta's congruence module, i.e. of the cokernel of the matrix `H` of Section 6 | `logs/congruence_module.log` |

## Claims

Every claim of Section 6 is labelled in the code by a comment `// Claim n: ...`, checked by `assert`,
and confirmed by the line `Claim n successfully verified` in the log.

| claim | statement in the paper | file |
|---|---|---|
| 1 | Section 6: `c(p)` divides `#J_1(p)(F_2)` for all primes `11 <= p <= 139` | `orders.m` |
| 2 | Table 1: the primes `11 <= p <= 139` with `c(p)` even, and the columns `v_2(c(p))` and `v_2(N_2)` | `orders.m` |
| 3 | Table 1: the column `d_2 = dim_{F_2} J_1(p)(F_2)[2]` | `two_rank.m` |
| 4 | Table 1: the column `coker(H) (x) Z_2`, and `v_2(det(pH)) = v_2(c(p))` | `congruence_module.m` |
| 5 | Section 6: for `p <= 113`, `coker(H) (x) Z_2` agrees with De Leo and Stoll, Table 4.5 | `congruence_module.m` |

The other statements of Section 6 (for example, which primes have `v_2(c(p)) = v_2(N_2)`, and the
comparison of 2-ranks) are read off from Table 1.

For Claim 5, the 2-primary parts of `J_1(p)(Q)_tors` for `p <= 113` are typed into
`congruence_module.m` from Table 4.5 of D. De Leo and M. Stoll, *On some open cases of a conjecture of
Conrad, Edixhoven and Stein*, arXiv:2505.10777v2.

## Running

From the repository root:

```sh
./run_all.sh
```

or run the scripts separately, e.g. `magma -b orders.m > logs/orders.log`. A single prime can be
done with `magma -b P:=73 two_rank.m`, and a range of primes with
`magma -b Pmin:=11 Pmax:=50 orders.m` (Claim 2 is only checked for the full range).

The logs in `logs/` were produced with Magma V2.30-1 on the *Mordell* server at the Department of
Mathematics, University of Zagreb (AMD EPYC 9175F, 384 GB of RAM). `congruence_module.m` and `orders.m`
take a few seconds of CPU time each, and `two_rank.m` about 13 minutes, of which about 6 minutes are
for `p = 137` (the CPU times per prime are in the logs).
