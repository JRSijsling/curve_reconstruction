/***
 *  Example file
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

SetVerbose("EndoFind", 1);
SetVerbose("Reconstruction", 1);
SetVerbose("CurveRec", 1);

// TODO: We should be able to set this lower, the current state of affairs is ridiculous
prec := 200;
F := RationalsExtra(prec);
CC := F`CC;

// Define curve
R<x,y,z> := PolynomialRing(F, 3);

f := x*y^3 + y^4 + x^3*z - x*y^2*z + y^3*z + x^2*z^2 - 2*x*y*z^2 - x*z^3;
f := x^3*y + x^3*z + x^2*y^2 + 3*x^2*y*z + x^2*z^2 - 4*x*y^3 - 3*x*y^2*z - 3*x*y*z^2 - 4*x*z^3 + 2*y^4 + 3*y^2*z^2 + 2*z^4;
//f := y^3*z - x^4 - z^4;

X := PlaneCurve(f);
f := DefiningPolynomial(X);
I, W := DixmierOhnoInvariants(f);
P := PeriodMatrix(X);

/* We see height messing up when not specifying base QQ: */
P1 := Submatrix(P, 1,1, 3,3); P1i := P1^(-1);
P2 := Submatrix(P, 1,4, 3,3);
tau := P1i*P2;

print "";
print "Geometric reconstruction:";
Y := ReconstructCurveGeometric(tau, F : Base := true);
J, W := DixmierOhnoInvariants(Y);
print Y;

print "";
print "Check invariants:";
print WPSNormalize(W, I) eq WPSNormalize(W, J);

Y := ReconstructCurve(P, F : Base := true);
g := R ! DefiningPolynomial(Y);
print "";
print "Arithmetic reconstruction over base:";
print PlaneCurve(g);
assert f eq g;

