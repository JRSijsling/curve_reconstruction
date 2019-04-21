/***
 *  Example file
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

prec := 200;
CC := ComplexFieldExtra(prec);
R<x> := PolynomialRing(CC);

//Define curve
f := x^7 + x + 1;
shioda, W := ShiodaInvariants(f);
shioda := WPSNormalizeCC(W, shioda : prec := 30);
X := SE_Curve(f, 2 : Prec := prec + 20);
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
// For characteristic 0:
//T := Matrix(Integers(), [
//[ 0,  0,  1,  0,  0, -1],
//[ 2,  1,  1,  0, -1, -1],
//[ 0, -2,  0,  1,  0,  0],
//[ 2,  0,  0,  0, -1,  1],
//[ 0,  1,  1,  0,  0, -1],
//[-1, -2,  0,  1,  0,  0]
//]);
print "done";
P := P*ChangeRing(T, CC);

P1 := Submatrix(P, 1,1, 3,3); P2 := Submatrix(P, 1,4, 3,3);
tau := P1^(-1)*P2;
tau := ReduceSmallPeriodMatrix(tau);

print "";
print "Check true:";
print IsSmallPeriodMatrix(tau);

print "";
print "Calculating Rosenhain invariants...";
rosens, v0 := RosenhainInvariantsBILV(tau);
print "done";

print "";
print "Vanishing theta characteristic:";
print v0;

print "";
print "Rosenhain invariants:";
print [ ComplexField(5) ! rosen : rosen in rosens ];
