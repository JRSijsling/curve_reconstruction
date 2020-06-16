/***
 *  Testing claimed property of U
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


import "../magma/precomp.m": PrecomputedGamma, GammaFor0;
import "../magma/rosenhain.m": EtaFunction0, EtaFunction, EtaValue, UFromEtaFunction, VectorFromIndex, IndexFromVector, IsEvenVector, LeftActionChar, FindDelta;


prec := 200;
CC := ComplexFieldExtra(prec);
R<x> := PolynomialRing(CC);

//Define curve
f := x^7 + x + 1;
shioda, W := ShiodaInvariants(f);
shioda := WPSNormalize(W, shioda);
X := RiemannSurface(f, 2 : Precision := prec + 20);
P := ChangeRing(X`BigPeriodMatrix, CC);

print "";
print "Finding transformation:";
while true do
    T := RandomSymplecticMatrix(3, 1);
    if &and[ Abs(c) le 2 : c in Eltseq(T) ] then
        break;
    end if;
end while;
print T;
T := Matrix(Integers(), [
[ 0,  0,  1,  0,  0, -1],
[ 2,  1,  1,  0, -1, -1],
[ 0, -2,  0,  1,  0,  0],
[ 2,  0,  0,  0, -1,  1],
[ 0,  1,  1,  0,  0, -1],
[-1, -2,  0,  1,  0,  0]
]);
print "done";
P := P*ChangeRing(T, CC);

P1 := Submatrix(P, 1,1, 3,3); P2 := Submatrix(P, 1,4, 3,3);
tau := P1^(-1)*P2;
tau := ReduceSmallPeriodMatrix(tau);

print "";
print "Check true:";
print IsSmallPeriodMatrix(tau);

print "";
print "In case of no further output everything is A-OK!";

/* Test extended property of U in Theorem 4.2 */
thetas_sq := ThetaSquares(tau);
v0 := FindDelta(thetas_sq)[1];
gamma := PrecomputedGamma(v0);
eta := EtaFunction(gamma);
U := UFromEtaFunction(eta);
Ss := Subsets({ 1..8 });
for S in { S : S in Ss | #S mod 2 eq 0 } do
    v := EtaValue(eta, S);
    theta := thetas_sq[IndexFromVector(3, v)];
    test1 := (Abs(theta) lt 10^-3);
    test2 := #(S sdiff U) ne 4;
    if test1 ne test2 then
        print "PROBLEM";
        print S;
    end if;
end for;
