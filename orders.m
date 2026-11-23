// Section 6 of "The rational torsion subgroup of J_1(p)": the orders c(p) and #J_1(p)(F_2).
//
// Run from the repository root:   magma -b orders.m > logs/orders.log
//
// c(p) = p * prod_{chi even, chi != 1} B_{2,chi}/4 is, up to sign, the order of the rational cuspidal
// subgroup C_1(p) (Kubert-Lang). By Lemma 4.1, B_{2,chi}/4 = sum_{a=1}^{n} beta(a) chi(a) with
// n = (p-1)/2 and beta(a) = a(a-p)/(2p), so |c(p)| = p * |Res(Phi, F)|, where Phi = (x^n-1)/(x-1),
// F = sum_{i=0}^{n-1} beta(g^i) x^i, g a generator of (Z/pZ)^*, and beta(+-a) := beta(a) for 1 <= a <= n.
//
// #J_1(p)(F_2) is computed with the Eichler-Shimura relation: it is the product, over the Galois
// orbits of even Dirichlet characters e modulo p, of the norms of det(1 + 2 e(2) - T_2) on the
// cuspidal modular symbols of weight 2 and character e.

SetColumns(0);
SetMemoryLimit(40*10^9);
if not assigned Pmin then Pmin := 11; elif Type(Pmin) eq MonStgElt then Pmin := StringToInteger(Pmin); end if;
if not assigned Pmax then Pmax := 139; elif Type(Pmax) eq MonStgElt then Pmax := StringToInteger(Pmax); end if;

R<x> := PolynomialRing(Rationals());

function Beta(a, p)
  a := a mod p; n := (p-1) div 2;
  if a gt n then a := p - a; end if;
  return (a*(a-p))/(2*p);
end function;

function CuspOrder(p)
  // |c(p)|
  n := (p-1) div 2; g := PrimitiveRoot(p);
  F := &+[ Beta(Modexp(g, i, p), p) * x^i : i in [0..n-1] ];
  Phi := (x^n - 1) div (x - 1);
  return Integers()!Abs(p*Resultant(Phi, F));
end function;

function NumPointsF2(p)
  // #J_1(p)(F_2)
  D := DirichletGroup(p, CyclotomicField(p-1));
  tot := 1;
  for e0 in GaloisConjugacyRepresentatives(D) do
    if not IsEven(e0) then continue; end if;
    e := MinimalBaseRingCharacter(e0);
    S := CuspidalSubspace(ModularSymbols(e, 2, 1));
    d := Dimension(S);
    if d eq 0 then continue; end if;
    K := BaseRing(S);
    v := Determinant(ScalarMatrix(K, d, K!(1 + 2*Evaluate(e, 2))) - HeckeOperator(S, 2));
    tot *:= (K eq Rationals()) select v else Norm(v);
  end for;
  return Integers()!Abs(tot);
end function;

// The rows of Table 1: <p, v_2(c(p)), v_2(#J_1(p)(F_2))> for the primes 11 <= p <= 139 with c(p) even.
Table1 := [ <17,3,3>, <29,6,8>, <31,2,2>, <41,4,8>, <43,2,2>, <73,3,3>, <89,3,11>, <97,4,4>,
            <109,4,10>, <113,20,32>, <127,8,8>, <137,5,10> ];

rows := [];
for p in PrimesInInterval(Pmin, Pmax) do
  t0 := Cputime();
  C := CuspOrder(p);
  N2 := NumPointsF2(p);
  printf "p=%o |c(p)|=%o #J1(p)(F2)=%o v2(c(p))=%o v2(#J1(p)(F2))=%o (%.1os)\n",
    p, C, N2, Valuation(C, 2), Valuation(N2, 2), Cputime(t0);
  // Claim 1: c(p) divides #J_1(p)(F_2) for all primes 11 <= p <= 139.
  assert N2 mod C eq 0;
  if IsEven(C) then Append(~rows, <p, Valuation(C, 2), Valuation(N2, 2)>); end if;
end for;
print "Claim 1 successfully verified for all primes in", [Pmin, Pmax];

if Pmin le 11 and Pmax ge 139 then
  // Claim 2: Table 1, columns p, v_2(c(p)) and v_2(N_2).
  assert rows eq Table1;
  print "Claim 2 successfully verified";
end if;
quit;
