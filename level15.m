// Section 6 of "The rational torsion subgroup of J_1(p)": a test of Theorem 1.3 at the composite level N = 15.
//
// Run from the repository root:   magma -b level15.m > logs/level15.log
//
// Theorem 1.3 says that the reduction map J_1(N)(Q)_tors -> J_1(N)(F_2) is injective for every odd N >= 5.
// The curve X_1(15) has genus 1 and rational points (its rational cusps), so J_1(15) is isomorphic to X_1(15).
//
// X_1(15) is computed from the Tate normal form E(b,c): y^2 + (1-c)xy - by = x^3 - bx^2, P = (0,0). Every pair
// (E,P) of an elliptic curve and a point of order at least 4 is isomorphic to a unique (E(b,c),(0,0)), so the curve
// of pairs (b,c) for which P has order 15 is a plane model of X_1(15). The condition 15P = 0, P != 0, is
// x(7P) = x(8P). The polynomial x(7P) - x(8P) has two linear factors, b (where E(b,c) is singular and P has
// order 3) and b - c (where P has order 5), and the remaining factor defines X_1(15).
//
// On a minimal Weierstrass model, a nonzero point reduces to 0 modulo 2 if and only if its x-coordinate is not
// 2-integral. The Cremona labels are those of J. E. Cremona, Algorithms for modular elliptic curves.

SetColumns(0);

function ReductionMod2IsInjective(E)
  // E: a minimal model of an elliptic curve over Q with good reduction at 2.
  T, m := TorsionSubgroup(E);
  return forall{ t : t in T | t eq T!0 or Valuation(m(t)[1], 2) ge 0 };
end function;

// --- X_1(15) from the Tate normal form.
F<b,c> := FunctionField(Rationals(), 2);
E := EllipticCurve([1-c, -b, -b, 0, 0]);
P := E![0,0];
R<B,C> := PolynomialRing(Rationals(), 2);
f := R!Numerator((7*P)[1] - (8*P)[1]);
fac := Factorization(f);
lin := [ g[1] : g in fac | TotalDegree(g[1]) eq 1 ];
big := [ g[1] : g in fac | TotalDegree(g[1]) gt 1 ];
// the factorization drops the leading unit, so compare up to scalars
assert #fac eq 3 and #lin eq 2 and #big eq 1;
assert { Normalize(g) : g in lin } eq { Normalize(B), Normalize(B - C) };
X := Curve(AffineSpace(R), big[1]);
assert Genus(X) eq 1;
Xp := ProjectiveClosure(X);
pts := PointSearch(Xp, 20);
E15 := 0;
for pt in pts do
  try
    E15 := MinimalModel(EllipticCurve(Xp, pt));
    break;
  catch e
    continue;
  end try;
end for;
assert Type(E15) eq CrvEll;
printf "X_1(15): %o, Cremona label %o\n", aInvariants(E15), CremonaReference(E15);
T15 := TorsionSubgroup(E15);
r1, r2 := RankBounds(E15);
E15F2 := ChangeRing(E15, GF(2));
printf "rank bounds %o %o, E(Q)_tors = %o, E(F_2) = %o, reduction mod 2 injective: %o\n",
  r1, r2, Invariants(T15), Invariants(AbelianGroup(E15F2)), ReductionMod2IsInjective(E15);
// Claim 6: X_1(15) is 15a8, J_1(15)(Q) = Z/4Z, J_1(15)(F_2) = Z/4Z, and the reduction map is an isomorphism.
assert CremonaReference(E15) eq "15a8";
assert r2 eq 0 and Invariants(T15) eq [4] and Invariants(AbelianGroup(E15F2)) eq [4];
assert ReductionMod2IsInjective(E15);
print "Claim 6 successfully verified";

// --- The other curves in the isogeny class 15a, and X_0(15).
X0 := SmallModularCurve(15);
E0 := MinimalModel(EllipticCurve(X0, PointSearch(X0, 10)[1]));
printf "X_0(15): Cremona label %o\n", CremonaReference(E0);
labels := [ "15a" cat IntegerToString(i) : i in [1..8] ];
for lab in labels do
  EE := EllipticCurve(lab);
  assert IsMinimalModel(EE) and IsIsogenous(EE, E15);
  printf "%o: E(Q)_tors = %o, E(F_2) = %o, reduction mod 2 injective: %o\n", lab,
    Invariants(TorsionSubgroup(EE)), Invariants(AbelianGroup(ChangeRing(EE, GF(2)))), ReductionMod2IsInjective(EE);
end for;
// Claim 7: for the seven curves 15a1, ..., 15a7, the reduction map E(Q)_tors -> E(F_2) is not injective.
// One of them is X_0(15) = 15a1, with E(Q)_tors = Z/2Z x Z/4Z and #E(F_2) = 4.
assert forall{ lab : lab in labels | (lab eq "15a8") eq ReductionMod2IsInjective(EllipticCurve(lab)) };
assert CremonaReference(E0) eq "15a1";
assert Invariants(TorsionSubgroup(E0)) eq [2,4] and #ChangeRing(E0, GF(2)) eq 4;
print "Claim 7 successfully verified";
quit;
