/***
 *  Testing claimed properties of eta
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


import "../magma/precomp.m": PrecomputedGamma, GammaFor0;
import "../magma/rosenhain.m": EtaFunction0, EtaFunction, EtaValue, UFromEtaFunction, VectorFromIndex, IndexFromVector, IsEvenVector, LeftActionChar, ThetaSquares, FindDelta;


function Q1Value(v)
w := Eltseq(v); g := (#w div 2); g := 3;
return &+[ w[i]*w[g + i] : i in [1..g] ];
end function;

function Q2Value(v1, v2)
w1 := Eltseq(v1); w2 := Eltseq(v2); g := (#w1 div 2); g := 3;
return &+[ w1[i]*w2[g + i] + w1[g + i]*w2[i] : i in [1..g] ];
end function;

function EtaValueMod2(eta, S)
return Transpose(Matrix(FiniteField(2), [[ Integers() ! c : c in Eltseq(2*EtaValue(eta, S)) ]]));
end function;


print "In case of no output everything is A-OK!";

/* Zero characteristic */
etas := [ ];
//Append(~etas, EtaFunction(gamma));
v0 := Matrix(Rationals(), [ [ 0 ] : i in [1..6] ]);
gamma := PrecomputedGamma(v0);
eta := EtaFunction(gamma);
//Append(~etas, EtaFunction0());
Append(~etas, eta);

for eta in etas do
    U := UFromEtaFunction(eta);
    T := { 1..8 };

    /* (3), complement */
    Ss := Subsets(T);
    for S in Ss do
        Sc := T diff S;
        if EtaValueMod2(eta, S) ne EtaValueMod2(eta, Sc) then
            print "PROBLEM";
            print S;
        end if;
    end for;

    /* (3), homomorphism */
    for S1 in Ss do
        for S2 in Ss do
            if EtaValueMod2(eta, S1 sdiff S2) ne (EtaValueMod2(eta, S1) + EtaValueMod2(eta, S2))  then
                print "PROBLEM";
                print S1;
                print S2;
            end if;
        end for;
    end for;

    /* (4) */
    for S1 in { S : S in Ss | #S eq 2 } do
        for S2 in { S : S in Ss | #S eq 2 } do
            v1 := EtaValueMod2(eta, S1);
            v2 := EtaValueMod2(eta, S2);
            e2 := Q2Value(v1, v2);
            p := (-1)^(#(S1 meet S2));
            f2 := FiniteField(2) ! ((p - 1) div 2);
            if e2 ne f2 then
                print "PROBLEM";
                print S1;
                print S2;
            end if;
        end for;
    end for;

    /* (5) */
    for S in { S : S in Ss | #S eq 2 } do
        v := EtaValueMod2(eta, S);
        e1 := Q1Value(v);
        p := (-1)^((3 + 1 - #(S sdiff U)) div 2);
        p := Integers() ! p;
        f1 := FiniteField(2) ! ((p - 1) div 2);
        if e1 ne f1 then
            print "PROBLEM";
            print S1;
        end if;
    end for;
end for;
