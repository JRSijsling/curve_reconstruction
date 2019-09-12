/***
 *  Example file
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

SetVerbose("EndoFind", 1);
SetVerbose("CurveRec", 1);

prec := 300;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);

//Define curve
R<x> := PolynomialRing(F);
f := 2*(x^5 + x + 1);
f := 3*(x^5 + 14*x^2 - x + 1);
D := [-5..5];
repeat
    f := &+[ Random(D)*x^i : i in [0..6] ];
until IsSquarefree(f) and (Degree(f) in [5,6]);
f := -x^6 - 5*x^5 - 3*x^4 + 3*x^2 - 5*x - 2;

X := HyperellipticCurve(f);
P := PeriodMatrix(X);

print "";
print "Curve:";
print X;

print "";
print "Mess up homology by:";
U := RandomSymplecticMatrix(2, 2);
print U;
P := P*ChangeRing(U, BaseRing(P));
tau := SmallPeriodMatrix(P);

print "";
print "Invariant reconstruction:";
I := AlgebraizedInvariants(tau, F);
print I;

print "";
print "Geometric reconstruction:";
X := ReconstructCurveGeometric(tau, F);
print X;

print "";
print "Geometric reconstruction over base:";
X := ReconstructCurveGeometric(tau, F : Base := true);
print X;

print "";
print "Arithmetic reconstruction over base:";
X := ReconstructCurve(P, F : Base := true);
print X;
