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
R<x> := PolynomialRing(CC);

// Define curve
f := x^3 + x + 1;
X := SE_Curve(f, 2 : Prec := prec);
P := ChangeRing(X`BigPeriodMatrix, CC);
P1 := Submatrix(P, 1,1, 1,1); P1i := P1^(-1);
P2 := Submatrix(P, 1,2, 1,1);
tau := P1i*P2;

print "Geometric reconstruction over given base:";
X := ReconstructCurveGeometric(tau, F);
print X;
print "";

print "Arithmetic reconstruction over given base:";
X := ReconstructCurveBase(P, F);
print X;
print "";

// BETA
// Define curve
f := x^3 + Sqrt(CC ! 2)*x + Sqrt(CC ! 3) + Sqrt(CC ! 5);
X := SE_Curve(f, 2 : Prec := prec);
P := ChangeRing(X`BigPeriodMatrix, CC);

print "Arithmetic reconstruction allowing for base extension:";
X := ReconstructCurveG1(P, F);
print X;
print "";

