/***
 *  Reduction of matrices under the action of the symplectic group
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


function IsSymmetricImproved(M : prec := 0)
/* To avoid stupid behavior of IsSymmetric */
if prec eq 0 then
    CC := BaseRing(M); prec := Precision(CC) - 17;
end if;
CCSmall := ComplexField(prec);
return Max([ Abs(c) : c in Eltseq(M - Transpose(M)) ]) lt 10^(-prec);
return IsSymmetric(ChangeRing(M, CCSmall));
end function;


function IsPositiveDefiniteImproved(M : prec := 0);
/* To avoid stupid behavior of IsPositiveDefinite */
if prec eq 0 then
    RR := BaseRing(M); prec := Precision(RR) - 17;
end if;
RRSmall := RealField(prec);
/* Deal with zero entries that are represented by 10^(-N) */
for i in [1..#Rows(M)] do
    for j in [1..#Rows(M)] do
        if Abs(M[i,j]) lt 10^(-prec) then
            M[i,j] := 0;
        end if;
    end for;
end for;
return IsPositiveDefinite(ChangeRing(M, RRSmall));
end function;


intrinsic IsSmallPeriodMatrix(tau::.) -> BoolElt
{Returns whether tau is (numerically) a small period matrix.}
Imtau := Matrix([ [ Im(c) : c in Eltseq(row) ] : row in Rows(tau) ]);
test := IsSymmetricImproved(tau) and IsPositiveDefiniteImproved(Imtau);
if not test then
    vprint CurveRec : tau;
end if;
return test;
end intrinsic;


intrinsic IsBigPeriodMatrix(P::.) -> BoolElt
{Returns whether P is (numerically) a big period matrix.}
P1 := Submatrix(P, 1,1, 2,2); P1i := P1^(-1);
P2 := Submatrix(P, 1,3, 2,2);
tau := P1i*P2;
return IsSmallPeriodMatrix(tau);
end intrinsic;


intrinsic LeftActionHg(M::AlgMatElt, tau::AlgMatElt) -> AlgMatElt
{Left symplectic action (by some definitions)}

assert IsSmallPeriodMatrix(tau);
g := #Rows(tau);
A := Submatrix(M, 1,  1, g,g); B := Submatrix(M, 1,  g+1, g,g);
C := Submatrix(M, g+1,1, g,g); D := Submatrix(M, g+1,g+1, g,g);
return (A*tau + B)*((C*tau + D)^(-1));

end intrinsic;


function LLLReduceMatrix(tau);
assert IsSmallPeriodMatrix(tau);
tau[2,1] := tau[1,2]; tau[3,1] := tau[1,3]; tau[3,2] := tau[2,3];
Imtau := Matrix([ [ Im(c) : c in Eltseq(row) ] : row in Rows(tau) ]);
_, T := LLLGram(Imtau);
return T*tau*Transpose(T);
end function;


function IntegerReduceMatrix(tau);
tauZZ := Matrix([ [ Round(Re(c)) : c in Eltseq(row) ] : row in Rows(tau) ]);
tauZZ[2,1] := tauZZ[1,2]; tauZZ[3,1] := tauZZ[1,3]; tauZZ[3,2] := tauZZ[2,3];
return tau - ChangeRing(tauZZ, BaseRing(tau));
end function;


function ReduceSmallPeriodMatrixG2Minkowski(tau)
// Dupont; end result is U*tau_orig*Transpose(U)
tau[2,1] := tau[1,2];
assert IsSmallPeriodMatrix(tau);
CC := BaseRing(tau); RR := RealField(CC);

t := true; U := Matrix(CC, [[1,0],[0,1]]);
while t do
    if 2*Abs(Im(tau[1,2])) le Abs(Im(tau[1,1])) then
        if Abs(Im(tau[1,1])) le Abs(Im(tau[2,2])) then
            if Im(tau[1,2]) le 0 then
                T := Matrix(CC, [[1,0],[0,-1]]);
                tau := T*tau*Transpose(T);
                U := T*U;
            end if;
            t := false;
        else
            T := Matrix(CC, [[0,1],[-1,0]]);
            tau := T*tau*Transpose(T);
            U := T*U;
        end if;
    end if;
    q := Round(Im(tau[1,2]) / Im(tau[1,1]));
    T := Matrix(CC, [[1,0],[-q,1]]);
    tau := T*tau*Transpose(T);
    U := T*U;
end while;
tau[2,1] := tau[1,2];
assert IsSmallPeriodMatrix(tau);
a := Im(tau[1,1]); b := Im(tau[1,2]); c := Im(tau[2,2]);
assert c ge a; assert a ge 2*b; assert 2*b ge 0;
return tau, U;

end function;


function ReduceSmallPeriodMatrixG2(tau)
// Dupont; end result is U*tau_orig*Transpose(U)
CC := BaseRing(tau); RR := RealField(CC);

Ms := [
Matrix(CC, [
[1 , 0 , 1 , 0],
[0 , 1 , 0 , 0],
[0 , 0 , 1 , 0],
[0 , 0 , 0 , 1]
]),
Matrix(CC, [
[1 , 0 , 0 , 0],
[0 , 1 , 0 , 1],
[0 , 0 , 1 , 0],
[0 , 0 , 0 , 1]
]),
Matrix(CC, [
[1 , 0 , 0 , 1],
[0 , 1 , 1 , 0],
[0 , 0 , 1 , 0],
[0 , 0 , 0 , 1]
])
];

Ns := [
Matrix(CC, [
[0 , 0 ,-1 , 0],
[0 , 0 , 0 ,-1],
[1 , 0 , 0 , 0],
[0 , 1 , 0 , 0]
]),
Matrix(CC, [
[0 , 0 ,-1 , 0],
[0 , 0 , 0 ,-1],
[1 , 0 , 1 , 0],
[0 , 1 , 0 , 0]
]),
Matrix(CC, [
[0 , 0 ,-1 , 0],
[0 , 0 , 0 ,-1],
[1 , 0 ,-1 , 0],
[0 , 1 , 0 , 0]
]),
Matrix(CC, [
[0 , 0 ,-1 , 0],
[0 , 0 , 0 ,-1],
[1 , 0 , 0 , 0],
[0 , 1 , 0 , 1]
]),
Matrix(CC, [
[0 , 0 ,-1 , 0],
[0 , 0 , 0 ,-1],
[1 , 0 , 0 , 0],
[0 , 1 , 0 ,-1]
]),
Matrix(CC, [
[0 , 0 ,-1 , 0],
[0 , 0 , 0 ,-1],
[1 , 0 , 1 , 0],
[0 , 1 , 0 , 1]
]),
Matrix(CC, [
[0 , 0 ,-1 , 0],
[0 , 0 , 0 ,-1],
[1 , 0 ,-1 , 0],
[0 , 1 , 0 ,-1]
]),
Matrix(CC, [
[0 , 0 ,-1 , 0],
[0 , 0 , 0 ,-1],
[1 , 0 ,-1 , 0],
[0 , 1 , 0 , 1]
]),
Matrix(CC, [
[0 , 0 ,-1 , 0],
[0 , 0 , 0 ,-1],
[1 , 0 , 1 , 0],
[0 , 1 , 0 ,-1]
]),
Matrix(CC, [
[0 , 0 ,-1 , 0],
[0 , 0 , 0 ,-1],
[1 , 0 , 0 , 1],
[0 , 1 , 1 , 0]
]),
Matrix(CC, [
[0 , 0 ,-1 , 0],
[0 , 0 , 0 ,-1],
[1 , 0 , 0 ,-1],
[0 , 1 ,-1 , 0]
]),
Matrix(CC, [
[0 , 0 ,-1 , 0],
[0 , 0 , 0 ,-1],
[1 , 0 , 1 , 1],
[0 , 1 , 1 , 0]
]),
Matrix(CC, [
[0 , 0 ,-1 , 0],
[0 , 0 , 0 ,-1],
[1 , 0 ,-1 ,-1],
[0 , 1 ,-1 , 0]
]),
Matrix(CC, [
[0 , 0 ,-1 , 0],
[0 , 0 , 0 ,-1],
[1 , 0 , 0 , 1],
[0 , 1 , 1 , 1]
]),
Matrix(CC, [
[0 , 0 ,-1 , 0],
[0 , 0 , 0 ,-1],
[1 , 0 , 0 ,-1],
[0 , 1 ,-1 ,-1]
]),
Matrix(CC, [
[1 , 0 ,-1 , 0],
[0 , 1 , 0 ,-1],
[1 , 0 , 0 , 0],
[0 , 0 , 0 , 1]
]),
Matrix(CC, [
[1 , 0 ,-1 , 0],
[0 , 1 , 0 ,-1],
[1 , 0 , 0 , 0],
[0 , 0 , 0 , 1]
]),
Matrix(CC, [
[1 , 0 ,-1 , 0],
[0 , 1 , 0 ,-1],
[0 , 0 , 1 , 0],
[0 , 1 , 0 , 0]
]),
Matrix(CC, [
[1 , 0 , 0 , 0],
[0 , 1 , 0 , 0],
[1 ,-1 , 1 , 0],
[-1, 1 , 0 , 1]
]),
Matrix(CC, [
[-1, 0 , 0 , 0],
[0 ,-1 , 0 , 0],
[1 ,-1 ,-1 , 0],
[-1, 1 , 0 ,-1]
])
];

gamma := IdentityMatrix(CC, 4); taup := tau; t := true;
while t do
    //print taup;
    taup, U := ReduceSmallPeriodMatrixG2Minkowski(taup);
    gamma := BlockMatrix([ [U, 0], [0, Transpose(U^(-1))] ]) * gamma;
    taupjs := [ taup[1,1], taup[2,2], taup[1,2] ];
    for j:=1 to 3 do
        M := Ms[j];
        a := -Round(Re(taupjs[j]));
        taup := LeftActionHg(M^a, taup);
        gamma := M^a * gamma;
    end for;
    t := false;
    for j:=1 to 19 do
        N := Ns[j];
        C := Submatrix(N, 3,1, 2,2);
        D := Submatrix(N, 3,3, 2,2);
        if Abs(Determinant(C*taup + D)) lt 1 then
            t := true;
            taup := LeftActionHg(N, taup);
            gamma := N * gamma;
            //break;
        end if;
    end for;
end while;

taup[2,1] := taup[1,2];
assert IsSmallPeriodMatrix(taup);
return taup, gamma;
end function;


function ReduceSmallPeriodMatrixG3(tau)
N0 := Matrix(Integers(), [
[  0,  0,  0, -1,  0,  0],
[  0,  1,  0,  0,  0,  0],
[  0,  0,  1,  0,  0,  0],
[  1,  0,  0,  0,  0,  0],
[  0,  0,  0,  0,  1,  0],
[  0,  0,  0,  0,  0,  1]
]);
while true do
    tau := LLLReduceMatrix(tau);
    tau := IntegerReduceMatrix(tau);
    if Abs(tau[1,1]) gt 0.99 then
        return tau;
    end if;
    tau := LeftActionHg(N0, tau);
end while;
end function;


intrinsic ReduceSmallPeriodMatrix(tau::.) -> .
{Reduces a small period matrix after Labrande--Thomé.}
assert IsSmallPeriodMatrix(tau);
g := #Rows(tau);
if g eq 2 then
    return ReduceSmallPeriodMatrixG2(tau);
elif g eq 3 then
    return ReduceSmallPeriodMatrixG3(tau);
end if;
return tau;
end intrinsic;
