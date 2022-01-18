/***
 *  Example file
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

SetSeed(1);
SetVerbose("EndoFind", 1);
SetVerbose("QuarticRec", 1);
SetVerbose("CurveRec", 2);

prec := 300;
F := RationalsExtra(prec);
R<x,y,z> := PolynomialRing(F, 3);
f := y^3*z - x^4 - z^4;
f := x*y^3 + y^4 + x^3*z - x*y^2*z + y^3*z + x^2*z^2 - 2*x*y*z^2 - x*z^3;
f := x^3*y + x^3*z + x^2*y^2 + 3*x^2*y*z + x^2*z^2 - 4*x*y^3 - 3*x*y^2*z - 3*x*y*z^2 - 4*x*z^3 + 2*y^4 + 3*y^2*z^2 + 2*z^4;

X := PlaneCurve(f);
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

Q := PeriodMatrix(X);
isos := SymplecticIsomorphismsCC(P, Q);

print "";
print "Arithmetic reconstruction over base:";
X := ReconstructCurve(P, F : Base := true);
print X;
