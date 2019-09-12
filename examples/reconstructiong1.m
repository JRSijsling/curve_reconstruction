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

f := x^3 + x + 7;
X := HyperellipticCurve(f);
P := PeriodMatrix(X);

print "";
print "Curve:";
print X;

print "";
print "Mess up homology by:";
U := RandomSymplecticMatrix(1, 2);
print U;
P := P*ChangeRing(U, BaseRing(P));
tau := SmallPeriodMatrix(P);

print "";
print "Invariant reconstruction:";
j := AlgebraizedInvariants(tau, F);
print j;

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
