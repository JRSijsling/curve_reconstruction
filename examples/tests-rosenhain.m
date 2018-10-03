/***
 *  Testing claimed properties of eta
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


import "../precomp.m": PrecomputedGamma, GammaFor0;
import "../rosenhain.m": EtaFunction0, EtaFunction, EtaValue, UFromEtaFunction, VectorFromIndex, IndexFromVector, IsEvenVector, LeftActionChar, ThetaSquares, FindDelta;


function StandardSymplecticMatrix(g)
A := ScalarMatrix(g, 0); B := ScalarMatrix(g, 1); C := -B; D := A;
return VerticalJoin(HorizontalJoin(A, B), HorizontalJoin(C, D));
end function;


function Q1Value(v)
w := Eltseq(v); g := (#w div 2);
return &+[ w[i]*w[g + i] : i in [1..g] ];
end function;


function Q2Value(v1, v2)
w1 := Eltseq(v1); w2 := Eltseq(v2); g := (#w1 div 2);
return &+[ w1[i]*w2[g + i] + w1[g + i]*w2[i] : i in [1..g] ];
end function;


function IsSymplectic1(M)
V := VectorSpace(FiniteField(2), 6);
for v in V do
    q1 := Q1Value(v);
    q2 := Q1Value(v*Transpose(ChangeRing(M, FiniteField(2))));
    if q1 ne q2 then
        return false;
    end if;
end for;
return true;
end function;


function IsSymplectic2(M)
E := StandardSymplecticMatrix(3); F := M*E*Transpose(M);
return &and[ IsEven(Integers() ! c) : c in Eltseq(E - F) ];
end function;


function EtaValueMod2(eta, S)
return Transpose(Matrix(FiniteField(2), [[ Integers() ! c : c in Eltseq(2*EtaValue(eta, S)) ]]));
end function;


print "In case of no output everything is A-OK!";


/* Sanity and compatibility checks follow */
/* Numbering */
n := Random([1..64]);
v := VectorFromIndex(3, n);
m := IndexFromVector(3, v);
if not n eq m then
    print "PROBLEM";
end if;

/* Precomputed gamma should preserve forms */
/* Does not work for 64, but that probably should be so */
vs := [ VectorFromIndex(3, i) : i in [1..63] ];
for v in vs do
    if IsEvenVector(v) then
        gamma := PrecomputedGamma(v);
        test1 := IsSymplectic1(gamma);
        test2 := IsSymplectic2(gamma);
        gamma := Matrix(Integers(), 6,6, [ Integers() ! c : c in Eltseq(gamma) ]);
        test3 := Determinant(ChangeRing(gamma, FiniteField(2))) eq 1;
        if not (test1 and test2 and test3) then
            print "PROBLEM";
            print v;
            print test1;
            print test2;
            print test3;
        end if;
    end if;
end for;

/* Alleged extra property of gamma */
E := StandardSymplecticMatrix(3);
vs := [ VectorFromIndex(3, i) : i in [1..63] ];
v0 := Transpose(Matrix(Rationals(),  [[ 1/2, 1/2, 1/2, 1/2, 0, 1/2 ]]));
for v in vs do
    if IsEvenVector(v) then
        gamma := PrecomputedGamma(v);
        gammav0 := LeftActionChar(gamma, v0);
        if not &and[ IsIntegral(c) : c in Eltseq(v - gammav0) ] then
            print "PROBLEM";
            print v;
        end if;
    end if;
end for;

/* Test if U does the trick */
for n in [1..64] do
    v := VectorFromIndex(3, n);
    if IsEvenVector(v) then
        gamma := PrecomputedGamma(v);
        eta := EtaFunction(gamma);
        U := UFromEtaFunction(eta);
        if not &and[ IsIntegral(c) : c in Eltseq(v - EtaValue(eta, U)) ] then
            print "PROBLEM";
        end if;
    end if;
end for;

/* Zero characteristic */
gamma := GammaFor0();
etas := [ ];
Append(~etas, EtaFunction(gamma));
Append(~etas, EtaFunction0());
n := Random([1..64]); v := VectorFromIndex(3, n); gamma := PrecomputedGamma(v);
Append(~etas, EtaFunction(gamma));
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
