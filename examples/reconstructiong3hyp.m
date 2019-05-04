/***
 *  Example file
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

SetVerbose("CurveRec", 2);
SetVerbose("EndoFind", 3);

prec := 300;
F := RationalsExtra(prec);
CC := F`CC;

f := x^7 + x + 1;
//D := [-5..5];
//f := x^8 + &+[ Random(D)*x^i : i in [0..7] ];

//Define curve
X := SE_Curve(f, 2 : Prec := prec);
S, W := ShiodaInvariants(f);
P := ChangeRing(X`BigPeriodMatrix, CC);

P1 := Submatrix(P, 1,1, 3,3); P1i := P1^(-1);
P2 := Submatrix(P, 1,4, 3,3);
tau := P1i*P2;

Y := ReconstructCurveGeometric(tau, F : Base := false);
T, W := ShiodaInvariants(Y);
print "";
print "Geometric reconstruction over base:";
print Y;

print "";
print "Inspect invariants:";
print WPSNormalize(W, S) eq WPSNormalize(W, T);
