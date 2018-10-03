/***
 *  Example file
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

SetVerbose("EndoFind", 0);

prec := 300;
F := RationalsExtra(prec);
CC := F`CC;

//Define curve
f := x^7 + x + 1;
X := SE_Curve(f, 2 : Prec := prec);
S, W := ShiodaInvariants(f);
P := ChangeRing(X`BigPeriodMatrix, CC);

P1 := Submatrix(P, 1,1, 3,3); P1i := P1^(-1);
P2 := Submatrix(P, 1,4, 3,3);
tau := P1i*P2;

print "Geometric reconstruction over given base:";
Y := ReconstructCurveGeometric(tau, F);
T, W := ShiodaInvariants(Y);

print "Inspect:";
print X; print Y;
print WPSNormalize(W, S) eq WPSNormalize(W, T);

print "";
print "Arithmetic reconstruction over given base:";
Y := ReconstructCurveBase(P, F);
print "Recover correct curve?";
print HyperellipticPolynomials(Y) eq f;

