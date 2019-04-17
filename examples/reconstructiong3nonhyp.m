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

// Define curve
R<x,y> := PolynomialRing(F, 2); z := 1;
f := x^3*y + x^3*z + x^2*y^2 + 3*x^2*y*z + x^2*z^2 - 4*x*y^3 - 3*x*y^2*z - 3*x*y*z^2 - 4*x*z^3 + 2*y^4 + 3*y^2*z^2 + 2*z^4;
X := PlaneCurve(f);
I, W := DixmierOhnoInvariants(f);
P := PeriodMatrix(X);

P1 := Submatrix(P, 1,1, 3,3); P1i := P1^(-1);
P2 := Submatrix(P, 1,4, 3,3);
tau := P1i*P2;

/* We see height messing up when not specifying base QQ:
print "Geometric reconstruction:";
Y := ReconstructCurveGeometric(tau, F);
J, W := DixmierOhnoInvariants(Y);

print "Inspect:";
print X; print Y;
print WPSNormalize(W, I) eq WPSNormalize(W, J);
*/

print "Geometric reconstruction over base:";
Y := ReconstructCurveGeometric(tau, F : Base := true);
J, W := DixmierOhnoInvariants(Y);

print "Inspect:";
print X; print Y;
print WPSNormalize(W, I) eq WPSNormalize(W, J);

print "Arithmetic reconstruction over base:";
Y := ReconstructCurve(P, F : Base := true);
print Y;
