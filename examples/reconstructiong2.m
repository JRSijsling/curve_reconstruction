/***
 *  Example file
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 1);

prec := 500;
F := RationalsExtra(prec);
CC := F`CC;
R<x> := PolynomialRing(CC);

//Define curve
R<x> := PolynomialRing(F);
f := 2*(x^5 + x + 1);
f := 3*(x^5 + 14*x^2 - x + 1);
D := [-5..5];
repeat
    f := &+[ Random(D)*x^i : i in [0..6] ];
until IsSquarefree(f) and (Degree(f) in [5,6]);
f := -x^6 - 5*x^5 - 3*x^4 + 3*x^2 - 5*x - 2;

print "Input polynomial:";
print f;

X := SE_Curve(f, 2 : Prec := prec);
P := ChangeRing(X`BigPeriodMatrix, CC) / 2;

U := RandomSymplecticMatrix(2, 2);
print "Mess up homology by:";
print U;
P := P*ChangeRing(U, CC);

P1 := Submatrix(P, 1,1, 2,2); P1i := P1^(-1);
P2 := Submatrix(P, 1,3, 2,2);
tau := P1i*P2;

print "Geometric reconstruction:";
Y := ReconstructCurveGeometric(tau, F);
print Y;
print "";

print "Geometric reconstruction over base:";
Y := ReconstructCurveGeometric(tau, F : Base := true);
print Y;
print "";

print "Arithmetic reconstruction over base:";
X := ReconstructCurve(P, F : Base := true);
print X;
print "";

print "Arithmetic reconstruction:";
Y := ReconstructCurve(P, F);
print "Recover correct curve?";
print HyperellipticPolynomials(Y) eq f;
print "";
