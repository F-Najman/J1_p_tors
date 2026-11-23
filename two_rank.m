// Section 6 of "The rational torsion subgroup of J_1(p)": the 2-rank of J_1(p)(F_2).
//
// Run from the repository root:   magma -b two_rank.m > logs/two_rank.log
// (or for a single prime:          magma -b P:=73 two_rank.m)
//
// By the isomorphism of Cartier and Serre, J_1(p)(F_2)[2] is isomorphic to the space of regular
// differentials on X_mu(p) over F_2 fixed by the Cartier operator. On q-expansions at the cusp
// infinity, the Cartier operator on S_2(Gamma_1(p); F_2) is U_2: sum a_k q^k -> sum a_{2k} q^k, so
//     dim_{F_2} J_1(p)(F_2)[2] = dim ker(U_2 - 1 | S_2(Gamma_1(p); F_2)).
// S_2(Gamma_1(p); Z) is the lattice of cusp forms with integral q-expansions. We compute it as the
// saturation in Z^(prec-1) of the lattice spanned by the q-expansions of Magma's basis (which is not
// the integral lattice for p >= 73), with prec = 2*(Sturm bound) + 4, and reduce it modulo 2. The
// matrix of U_2 is obtained by solving for the coefficients a_{2k}, 1 <= k <= (prec-1)/2, in terms of
// the first (prec-1)/2 coefficients, on which the reduced basis is checked to be linearly independent.

SetColumns(0);
SetMemoryLimit(40*10^9);

function TwoRankF2(p)
  S := CuspForms(Gamma1(p));
  B := Basis(S);
  g := #B;
  prec := 2*PrecisionBound(S) + 4;
  Q := [qExpansion(f, prec) : f in B];
  Mq := Matrix(Rationals(), [[Coefficient(Q[i], k) : k in [1..prec-1]] : i in [1..g]]);
  d := LCM([Denominator(c) : c in Eltseq(Mq)]);
  Mz := ChangeRing(d*Mq, Integers());
  assert Rank(Mz) eq g;
  H := EchelonForm(Saturation(Mz));            // a Z-basis of S_2(Gamma_1(p); Z), truncated
  assert Nrows(H) eq g;
  F2 := GF(2);
  H2 := ChangeRing(H, F2);
  half := (prec-1) div 2;
  V := Submatrix(H2, 1, 1, g, half);
  assert Rank(V) eq g;                          // reductions determined by the first half
  rowsU := [];
  for i in [1..g] do
    ok, c := IsConsistent(V, Vector(F2, [H2[i, 2*k] : k in [1..half]]));
    assert ok;                                  // U_2 preserves S_2(Gamma_1(p); F_2)
    Append(~rowsU, Eltseq(c));
  end for;
  U2 := Matrix(F2, rowsU);
  return Dimension(Kernel(U2 - 1)), g;
end function;

// Table 1, column d_2 = dim_{F_2} J_1(p)(F_2)[2].
Table1 := [ <17,1>, <29,5>, <31,2>, <41,2>, <43,2>, <73,3>, <89,5>, <97,1>,
            <109,6>, <113,9>, <127,8>, <137,2> ];

if assigned P then
  ps := [StringToInteger(P)];
else
  ps := [t[1] : t in Table1];
end if;

for p in ps do
  t0 := Cputime();
  d2, g := TwoRankF2(p);
  printf "p=%o genus=%o dim_F2 J1(p)(F2)[2]=%o (%.1os)\n", p, g, d2, Cputime(t0);
  exp := [t[2] : t in Table1 | t[1] eq p];
  if #exp eq 1 then
    // Claim 3: Table 1, column d_2.
    assert d2 eq exp[1];
  end if;
end for;
if not assigned P then print "Claim 3 successfully verified"; end if;
quit;
