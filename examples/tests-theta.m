/***
 *  Testing speed of theta algorithms
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


import "../precomp.m": PrecomputedGamma, GammaFor0;
import "../rosenhain.m": EtaFunction0, EtaFunction, EtaValue, UFromEtaFunction, VectorFromIndex, IndexFromVector, IsEvenVector, LeftActionChar, ThetaSquares, FindDelta;
import "../fastthetaconstantsgenus3.m": NaiveThetaConstantsGenus3, CalculThetas;


prec := 500;
CC := ComplexField(prec);
R<x> := PolynomialRing(CC);

//Define curve
f := x^7 + x + 1;
shioda, W := ShiodaInvariants(f);
shioda := WPSNormalize(W, shioda);
X := SE_Curve(f, 2 : Prec := prec + 20);
P := ChangeRing(X`BigPeriodMatrix, CC);

print "Finding transformation:";
while true do
    T := RandomSymplecticMatrix(3, 1);
    if &and[ Abs(c) le 2 : c in Eltseq(T) ] then
        break;
    end if;
end while;
print T;
print "done";
//P := P*ChangeRing(T, CC);

P1 := Submatrix(P, 1,1, 3,3); P2 := Submatrix(P, 1,4, 3,3);
tau := P1^(-1)*P2;
tau := ReduceSmallPeriodMatrix(tau);
print "Reduce period matrix:";
print ChangeRing(tau, ComplexField(5));

print "NaiveThetaConstantsGenus3:";
time funds := NaiveThetaConstantsGenus3(tau/2, false);
print "CalculThetas:";
time funds := CalculThetas(tau/2);
