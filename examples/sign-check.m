/***
 *  Example file
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

SetVerbose("CurveRec", 1);

prec := 200;
CC := ComplexFieldExtra(prec);
R<x> := PolynomialRing(CC);

print "";
print "In case of no error of the type \"Wrong sign tuple:\" everything is A-OK!";
print "";
counter := 0;
while true do
    counter +:= 1;
    print "Curve number", counter;
    f := x^7 + &+[ Random([-2..2])*x^i : i in [0..6] ];
    if Abs(Discriminant(f)) gt 10^(-3) then
        shioda, W := ShiodaInvariants(f);
        shioda := WPSNormalizeCC(W, shioda : prec := 30);
        X := SE_Curve(f, 2 : Prec := prec + 20);
        P := ChangeRing(X`BigPeriodMatrix, CC);

        while true do
            T := RandomSymplecticMatrix(3, 1);
            if &and[ Abs(c) le 2 : c in Eltseq(T) ] then
                break;
            end if;
        end while;
        P := P*ChangeRing(T, CC);

        P1 := Submatrix(P, 1,1, 3,3); P2 := Submatrix(P, 1,4, 3,3);
        tau := P1^(-1)*P2;
        tau := ReduceSmallPeriodMatrix(tau);

        rosens, v0 := RosenhainInvariantsBILV(tau);
        CP := CartesianPower([ 1, -1 ], 5);
        tups0 := [ ];
        for tup in CP do
            f_rec := x*(x - 1)*&*[ x - tup[i]*rosens[i] : i in [1..5] ];
            shioda_rec, W := ShiodaInvariants(f_rec);
            shioda_rec := WPSNormalizeCC(W, shioda_rec : prec := 30);
            if Max([ Abs(shioda[i] - shioda_rec[i]) : i in [1..9] ]) lt 10^(-3) then
                Append(~tups0, tup);
            end if;
        end for;

        if #tups0 ne 1 then
            error "Wrong number of sign tuples:", #tups0;
        end if;
        if [ c : c in tups0[1] ] ne [1,1,1,1,1] then
            error "Wrong sign tuple:", tups0[1];
        end if;
    end if;
end while;
