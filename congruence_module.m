// Section 6 of "The rational torsion subgroup of J_1(p)": the 2-primary part of Ohta's congruence
// module, which is isomorphic to T/I and hence (Corollary 5.3) to J_1(p)(Q)[2^oo].
//
// Run from the repository root:   magma -b congruence_module.m > logs/congruence_module.log
//
// Let n = (p-1)/2, g a generator of (Z/pZ)^*, and beta(+-a) = a(a-p)/(2p) for 1 <= a <= n
// (Lemma 4.1). In the integral coordinates v(g^j), 1 <= j <= n-1, of the lattice L_1 of
// [Ohta, Lemma 2.4.6], the residue map of [Ohta, (2.2.5)] has matrix
//     H[i,j] = beta(g^-i) - beta(g^(j-i)),        1 <= i, j <= n-1,
// in the basis [g^-i]_oo - [1]_oo, so the congruence module over Z_2 is coker(H) tensor Z_2. As p is
// odd, this is the 2-primary part of coker(p*H), and p*H has integer entries.

SetColumns(0);

function Beta(a, p)
  a := a mod p; n := (p-1) div 2;
  if a gt n then a := p - a; end if;
  return (a*(a-p))/(2*p);
end function;

function TwoPartCokernel(p)
  n := (p-1) div 2; g := PrimitiveRoot(p);
  H := Matrix(Rationals(), n-1, n-1,
         [Beta(Modexp(g, -i, p), p) - Beta(Modexp(g, j-i, p), p) : j in [1..n-1], i in [1..n-1]]);
  pH := ChangeRing(p*H, Integers());
  D := ElementaryDivisors(pH);
  return [2^Valuation(d, 2) : d in D | Valuation(d, 2) gt 0], Valuation(Abs(Determinant(pH)), 2);
end function;

// Table 1, column coker(H) tensor Z_2, as invariant factors.
Table1 := [ <17,[8]>, <29,[4,4,4]>, <31,[2,2]>, <41,[16]>, <43,[2,2]>, <73,[2,2,2]>, <89,[2,2,2]>,
            <97,[16]>, <109,[4,4]>, <113,[4,4,16,16,16,16]>, <127,[2,2,2,2,2,2,2,2]>, <137,[32]> ];

// The 2-primary parts of J_1(p)(Q)_tors for p <= 113 with c(p) even, read off from
// De Leo and Stoll, Table 4.5 (arXiv:2505.10777v2).
DeLeoStoll := [ <17,[8]>, <29,[4,4,4]>, <31,[2,2]>, <41,[16]>, <43,[2,2]>, <73,[2,2,2]>,
                <89,[2,2,2]>, <97,[16]>, <109,[4,4]>, <113,[4,4,16,16,16,16]> ];

// v_2(c(p)) from Table 1 (see orders.m).
V2c := [ <17,3>, <29,6>, <31,2>, <41,4>, <43,2>, <73,3>, <89,3>, <97,4>, <109,4>, <113,20>,
         <127,8>, <137,5> ];

for t in Table1 do
  p := t[1];
  two, v2det := TwoPartCokernel(p);
  printf "p=%o 2-part of coker(pH)=%o v2(det(pH))=%o\n", p, two, v2det;
  // Claim 4: Table 1, column coker(H) tensor Z_2, and v_2(det(pH)) = v_2(c(p)).
  assert two eq t[2];
  assert v2det eq [s[2] : s in V2c | s[1] eq p][1];
end for;
print "Claim 4 successfully verified";

// Claim 5: for p <= 113 the cokernel agrees with De Leo-Stoll, Table 4.5.
for s in DeLeoStoll do
  assert [t[2] : t in Table1 | t[1] eq s[1]][1] eq s[2];
end for;
print "Claim 5 successfully verified";
quit;
