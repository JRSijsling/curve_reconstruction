/***
 *  Main algorithms by Balakrishnan--Ionica--Lauter--Vincent
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


import "fastthetasgenus2.m": ThetaSquaresLabrandeGenus2;
import "fastthetaconstantsgenus3.m": ThetaSquaresLabrandeGenus3;
import "precomp.m": PrecomputedGamma;


intrinsic WPSNormalizeCC(W::SeqEnum, S::SeqEnum : prec := 0) -> SeqEnum
{Better version of WPSNormalize over CC.}
CC := Parent(S[1]);
if prec eq 0 then prec := Precision(CC) / 2; end if;
I := [ i : i in [1..#S] | not Abs(S[i]) lt 10^(-prec) ];
WT := [ W[i] : i in I ]; T := [ S[i] : i in I ];
T0 := WPSNormalize(WT, T); S0 := [ CC ! 0 : s in S ];
for i in [1..#I] do S0[I[i]] := T0[i]; end for;
return S0;
end intrinsic;


/* The first two functions are for compatibility with Labrande's conventions */
/* The vectors returned are column vectors, not because I think that this is so
 * great, but because this respects Magma's conventions */
function VectorFromIndex(g, i)
j := i; v := [ ];
for k:=1 to 2*g do
    Append(~v, (j mod 2)/2);
    j := j div 2;
end for;
return Transpose(Matrix([Reverse(v)]));
end function;


/* Second return value is sign transition from standard value */
function IndexFromVector(g, v);
w := Reverse(Eltseq(v));
i := &+[ ((Integers() ! (2*w[j])) mod 2)*2^(j - 1) : j in [1..(2*g)] ];
if i eq 0 then
    return 2^(2*g);
end if;
v0 := VectorFromIndex(g, i); n := v - v0;
v0 := Eltseq(v0); n := Eltseq(n);
s := &+[ v0[i]*n[g + i] : i in [1..g] ];
if IsIntegral(s) then
    s := 1;
else
    s := -1;
end if;
return i, s;
end function;


/* Checks if a vector is even */
function IsEvenVector(v)
w := Eltseq(v); g := #w div 2;
w1 := w[1..g]; w2 := w[(g + 1)..(2*g)];
return IsZero((Integers() ! (4 * &+[ w1[i]*w2[i] : i in [1..g] ])) mod 2);
end function;


/* Action on vectors */
function LeftActionChar(gamma, v)
g := #Rows(gamma) div 2;
return VectorFromIndex(g, IndexFromVector(g, ChangeRing(gamma, Rationals())*v));
end function;


/* Mumford's eta */
function EtaFunction0();
eta1 := [1/2,0,0,0,0,0]; eta2 := [1/2,0,0,1/2,0,0];
eta3 := [0,1/2,0,1/2,0,0]; eta4 := [0,1/2,0,1/2,1/2,0];
eta5 := [0,0,1/2,1/2,1/2,0]; eta6 := [0,0,1/2,1/2,1/2,1/2];
eta7 := [0,0,0,1/2,1/2,1/2]; eta8 := [0,0,0,0,0,0];
etas := [eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8];
return [ Transpose(Matrix(Rationals(), [eta])) : eta in etas ];
end function;


/* Transform Mumford's eta by gamma */
/* Note that these function are represented by tuples, but we give an evaluation later */
function EtaFunction(gamma)
return [ LeftActionChar(gamma, v) : v in EtaFunction0() ];
end function;


/* Value of eta on a set of branch points */
function EtaValue(eta, S)
g := #Eltseq(eta[1]) div 2;
if #S eq 0 then
    return Transpose(Matrix(Rationals(), [[ 0 : i in [1..2*g] ]]));
end if;
return &+[ eta[i] : i in S ];
end function;


/* Given eta, finds U in Definition 3.3 */
function UFromEtaFunction(eta)
g := #Eltseq(eta[1]) div 2;
return { i : i in [1..#eta] | not IsEvenVector(eta[i]) } join { 2*g + 2 };
end function;


intrinsic ThetaSquares(tau::AlgMatElt : Labrande := true) -> SeqEnum
{Calculate theta null values.}
g := #Rows(tau);
if Labrande and g in [2,3] then
    if g eq 2 then
        FindThetasSq := ThetaSquaresLabrandeGenus2;
    elif g eq 3 then
        FindThetasSq := ThetaSquaresLabrandeGenus3;
    end if;
    CC0 := BaseRing(tau); prec0 := Precision(CC0); RR0 := RealField(CC0);
    eps := RR0 ! (10^(-9.5 * (prec0 div 10)));

    prec := prec0;
    repeat
        CC := ComplexFieldExtra(prec);
        thetas_sq := FindThetasSq(tau);
        prec +:= (prec0 div 5);
        CCnew := ComplexFieldExtra(prec);
        tau := ChangeRing(tau, CCnew);
        thetas_sqnew := FindThetasSq(tau);
        dif := Maximum([ Abs(thetas_sq[i] - thetas_sqnew[i]) : i in [1..#thetas_sq] ]);

        vprint CurveRec, 2 : "";
        vprint CurveRec, 2 : "Precision reached while refining thetas:";
        vprint CurveRec, 2 : RealField(5) ! dif;
    until dif lt eps;
    return ChangeUniverse(thetas_sqnew, CC0);
end if;

g := #Rows(tau);
M0 := ZeroMatrix(Rationals(), g, 1);
return [ Theta(VectorFromIndex(g, i), M0, tau)^2 : i in [1..2^(2*g)] ];
end intrinsic;


intrinsic ThetaValues(tau::AlgMatElt : Labrande := true) -> SeqEnum
{Calculate theta null values.}
CC := BaseRing(tau);
if Labrande then
    thetas_sq := ThetaSquares(tau);
    thetas := [ ];
    g := #Rows(tau);
    M0 := ZeroMatrix(Rationals(), g, 1);
    tausmall := ChangeRing(tau, ComplexField(30));
    for i in [1..2^(2*g)] do
        v := VectorFromIndex(g, i);
        if Abs(thetas_sq[i]) lt CC`epscomp then
            Append(~thetas, 0);
        else
            thetasmall := Theta(v, M0, tausmall);
            thetasqrt := Sqrt(thetas_sq[i]);
            thetasqrts := [ thetasqrt, -thetasqrt ];
            min, ind := Min([ Abs(thetasqrt - thetasmall) : thetasqrt in thetasqrts ]);

            assert min le 10^(-10);
            thetasqrt := thetasqrts[ind];
            Append(~thetas, thetasqrt);
        end if;
    end for;
    return thetas, thetas_sq;
end if;
g := #Rows(tau);
M0 := ZeroMatrix(Rationals(), g, 1);
thetas := [ Theta(VectorFromIndex(g, i), M0, tau) : i in [1..2^(2*g)] ];
thetas_sq := [ theta^2 : theta in thetas ];
return thetas, thetas_sq;
end intrinsic;


/* This finds the unique even zero characteristic */
function FindDelta(thetas : prec := 20)
if prec eq 0 then
    //CC := Parent(thetas[1]); prec := Precision(CC) - 100;
    CC := Parent(thetas[1]); prec := Precision(CC) / 2;
end if;
CP := CartesianPower([0,1/2], 6);
vs := [ Transpose(Matrix([[ c : c in tup ]])) : tup in CP ];
v0s := [ ];
for v in vs do
    theta := thetas[IndexFromVector(3, v)];
    test := (Abs(theta) lt 10^(-prec)) and IsEvenVector(v);
    if test then
        Append(~v0s, v);
    end if;
end for;
return v0s;
end function;


function EtaInnerproduct(eta1, eta2)
eta1 := [ Integers() ! c : c in Eltseq(eta1) ];
eta2 := [ Integers() ! c : c in Eltseq(eta2) ];
g := #eta1 div 2;
return &+[ eta1[i]*eta2[g + i] : i in [1..g] ];
end function;


function LogEpsilon(U, j)
if not j in U then
    return 1;
end if;
return 0;
end function;


/* Calculates sign */
function EpsilonKLM(eta, k, l, m)
U := UFromEtaFunction(eta);
exp := EtaInnerproduct(2*eta[k], 2*(eta[k] + eta[l] + eta[m])) + LogEpsilon(U, k) - 1;
return (-1)^exp;
end function;


/* Theorem 4.5 */
function TakaseQuotient(thetas_sq, eta, k, l, m)
U := UFromEtaFunction(eta);
Bm := { 1..7 }; L := [ bp : bp in (Bm diff { k, l, m }) ];
V := { L[1], L[2] }; W := { L[3], L[4] };
eps := EpsilonKLM(eta, k, l, m);
num1 := thetas_sq[IndexFromVector(3, EtaValue(eta, U sdiff (V join { k, l })))];
num2 := thetas_sq[IndexFromVector(3, EtaValue(eta, U sdiff (W join { k, l })))];
den1 := thetas_sq[IndexFromVector(3, EtaValue(eta, U sdiff (V join { k, m })))];
den2 := thetas_sq[IndexFromVector(3, EtaValue(eta, U sdiff (W join { k, m })))];
return eps*(num1*num2)/(den1*den2);
end function;


intrinsic RosenhainInvariantsBILV(tau::.) -> SeqEnum, .
{Calculates the Rosenhain invariants of the small period matrix tau.}
tau := ReduceSmallPeriodMatrix(tau);
thetas_sq := ThetaSquares(tau); v0s := FindDelta(thetas_sq);
if #v0s ne 1 then
    error "Not right number of even zero characteristics:", v0s;
end if;
v0 := v0s[1];
gamma := PrecomputedGamma(v0); eta := EtaFunction(gamma);
rosens := [ TakaseQuotient(thetas_sq, eta, 1, l, 2) : l in [3..7] ];
return rosens, v0;
end intrinsic;


/* Final function */
function RosenhainInvariantsFromThetaSquares(thetas_sq)
v0s := FindDelta(thetas_sq);
if #v0s ne 1 then
    error "Not right number of even zero characteristics:", v0s;
end if;
v0 := v0s[1];
gamma := PrecomputedGamma(v0); eta := EtaFunction(gamma);
rosens := [ TakaseQuotient(thetas_sq, eta, 1, l, 2) : l in [3..7] ];
return rosens, v0;
end function;


function ShiodaInvariantsFromThetaSquares(thetas_sq)

rosens := RosenhainInvariantsFromThetaSquares(thetas_sq);
CC := Parent(rosens[1]); R<x> := PolynomialRing(CC);
f := x*(x - 1)* &*[ x - rosen : rosen in rosens ];
SCC, W := ShiodaInvariants(f);
return WPSNormalizeCC(W, SCC);

end function;
