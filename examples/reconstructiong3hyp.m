/***
 *  Example file
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

SetVerbose("EndoFind", 1);
SetVerbose("CurveRec", 1);

prec := 500;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);

D := [-5..5];
repeat
  f := x^8 + &+[ Random(D)*x^i : i in [0..7] ];
until Discriminant(f) ne 0;
f := x^7 + x + 1;

X := HyperellipticCurve(f);
P := PeriodMatrix(X);

print "";
print "Curve:";
print X;

print "";
print "Mess up homology by:";
U := RandomSymplecticMatrix(3, 2);
print U;
P := P*ChangeRing(U, BaseRing(P));
tau := SmallPeriodMatrix(P);

print "";
print "Invariant reconstruction:";
I := AlgebraizedInvariants(tau, F : Base := false);
print I;

print "";
print "Geometric reconstruction:";
X := ReconstructCurveGeometric(tau, F);
print X;

print "";
print "Geometric reconstruction over base:";
X := ReconstructCurveGeometric(tau, F : Base := true);
print X;
