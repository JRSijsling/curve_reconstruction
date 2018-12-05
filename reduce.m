/***
 *  Reduction of matrices under the action of the symplectic group
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *  using implementations by Marco Streng
 *
 *  See LICENSE.txt for license details.
 */


function IsSymmetricImproved(M : prec := 0)
/* To avoid stupid numerical behavior of IsSymmetric */
if prec eq 0 then
    CC := BaseRing(M); prec := Floor(-Log(CC`epscomp)/Log(10));
end if;
return Max([ Abs(c) : c in Eltseq(M - Transpose(M)) ]) lt 10^(-prec);
end function;


function IsPositiveDefiniteImproved(M : prec := 0);
/* To avoid stupid numerical behavior of IsPositiveDefinite */
if prec eq 0 then
    RR := BaseRing(M); prec := Floor(-Log(RR`epscomp)/Log(10));
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
    vprint CurveRec : "";
    vprint CurveRec : "Not a small period matrix:";
    vprint CurveRec : tau;
end if;
return test;
end intrinsic;


intrinsic SmallPeriodMatrix(P::ModMatFldElt) -> .
{Returns small period matrix associated to P.}
g := #Rows(P);
P1 := Submatrix(P, 1,1,   g,g); P1i := P1^(-1);
P2 := Submatrix(P, 1,g+1, g,g);
tau := P1i*P2;
return tau;
end intrinsic;


intrinsic IsBigPeriodMatrix(P::.) -> BoolElt
{Returns whether P is (numerically) a big period matrix.}
return IsSmallPeriodMatrix(SmallPeriodMatrix(P));
end intrinsic;


intrinsic LeftActionHg(M::AlgMatElt, tau::AlgMatElt) -> AlgMatElt
{Left symplectic action (by some definitions)}
assert IsSmallPeriodMatrix(tau);
g := #Rows(tau);
A := Submatrix(M, 1,  1, g,g); B := Submatrix(M, 1,  g+1, g,g);
C := Submatrix(M, g+1,1, g,g); D := Submatrix(M, g+1,g+1, g,g);
return (A*tau + B)*((C*tau + D)^(-1));
end intrinsic;


function IntegerReduceMatrixG1(tau)
CC := BaseRing(tau);
tauZZ := Matrix(CC, [ [ Ceiling(Re(c) - (1/2)) : c in Eltseq(row) ] : row in Rows(tau) ]);
U := Matrix(CC, [ [ 1, -tauZZ[1,1] ], [ 0, 1 ] ]);
return tau - tauZZ, U;
end function;


function InverseReduceMatrixG1(tau)
CC := BaseRing(tau);
if Abs(tau[1,1]) ge 1 then
    U := IdentityMatrix(CC, 2);
    return tau, U, true;
end if;
U := Matrix(CC, [ [ 0, -1 ], [ 1, 0 ] ]);
return Matrix(CC, [[ -1/tau[1,1] ]]), U, false;
end function;


function ReduceSmallPeriodMatrixG1(tau)
CC := BaseRing(tau); taured := tau; T := IdentityMatrix(CC, 2);
repeat
    taured, U := IntegerReduceMatrixG1(taured);
    T := U*T;

    taured, U, done := InverseReduceMatrixG1(taured);
    T := U*T;
until done;

assert IsSmallPeriodMatrix(taured);
TZZ := Matrix(Integers(), [ [ Round(c) : c in Eltseq(row) ] : row in Rows(T) ]);
assert IsSymplecticMatrix(TZZ);
return taured, T;
end function;


function IntegerReduceMatrixG2(tau)
// Implementation by Marco Streng
CC := BaseRing(tau);
assert Abs(tau[2,1] - tau[1,2]) lt CC`epscomp;
tau[1,2] := tau[2,1];
tauZZ := Matrix(CC, [ [ Ceiling(Re(c) - (1/2)) : c in Eltseq(row) ] : row in Rows(tau) ]);
return tau - tauZZ, -tauZZ;
end function;


function MinkowskiReductionG2RR(M)
// Implementation by Marco Streng
RR := BaseRing(M);
assert Determinant(M) gt 0;
assert Abs(M[1,2] - M[2,1]) le RR`epscomp;
assert M[1,1] gt 0;
Mred := M; T := IdentityMatrix(RR, 2);
done := false;
while not done do
    r := Floor(-Mred[1,2]/Mred[1,1] + (1/2));
    R := Matrix(RR, [[1,0],[r,1]]);
    T := R*T;
    Mred := R*Mred*Transpose(R);
    if Mred[1,1] gt Mred[2,2] then
        U := Matrix(RR, [[0,1],[-1,0]]);
        T := U*T;
        Mred := U*Mred*Transpose(U);
    else
        done := true;
    end if;
end while;
if Mred[1,2] lt 0 then
    U := Matrix(RR, [[1,0],[0,-1]]);
    T := U*T;
    Mred := U*Mred*Transpose(U);
end if;
return Mred, T;
end function;


function MinkowskiReductionG2CC(tau)
// Implementation by Marco Streng
CC := BaseRing(tau); RR := RealField(CC);
assert Abs(tau[2,1] - tau[1,2]) lt CC`epscomp;
tau[1,2] := tau[2,1];
M := Matrix(RR, [ [ Im(c) : c in Eltseq(row) ] : row in Rows(tau) ]);
_, T := MinkowskiReductionG2RR(M);
T := ChangeRing(T, CC);
taured := T*tau*Transpose(T);
return taured, T;
end function;


function GottschlingMatrices()
// Implementation by Marco Streng
Ms := [ ];
Ms cat:= [ Matrix([[0,0,-1,0],[0,1,0,0],[1,0,e,0],[0,0,0,1]]) : e in [-1,0,1] ];
Ms cat:= [ Matrix([[1,0,0,0] ,[0,0,0,-1],[0,0,1,0],[0,1,0,e]]) : e in [-1,0,1] ];
Ms cat:= [ Matrix([[0,0,-1,0],[0,1,0,0],[1,-1,e,0],[0,0,1,1]]) : e in [-2,-1,0,1,2] ];
CP := CartesianPower([-1,0,1], 3);
CP := [ [ c : c in tup ] : tup in CP ];
Ms cat:= [ Matrix([[0,0,-1,0],[0,0,0,-1],[1,0,tup[1],tup[3]],[0,1,tup[3],tup[2]]]) : tup in CP ];
return Ms;
end function;


function GottschlingReduce(tau)
// Implementation by Marco Streng
CC := BaseRing(tau); RR := RealField(CC);
assert Abs(tau[2,1] - tau[1,2]) lt CC`epscomp;
tau[1,2] := tau[2,1];
Ns := [ ChangeRing(N, CC) : N in GottschlingMatrices() ];
impr := RR ! 1; found_impr := false;
for N in Ns do
    C := Submatrix(N, 3,1, 2,2);
    D := Submatrix(N, 3,3, 2,2);
    impr_new := Abs(Determinant(C*tau + D));
    if impr_new lt impr then
        T := N;
        impr := impr_new;
        found_impr := true;
    end if;
end for;
if not found_impr then
    T := IdentityMatrix(CC, 4);
end if;
taured := LeftActionHg(T, tau);
return taured, T, not found_impr;
end function;


function ReduceSmallPeriodMatrixG2(tau)
// Implementation by Marco Streng
CC := BaseRing(tau); taured := tau; T := IdentityMatrix(CC, 4);
repeat
    taured, gamma := MinkowskiReductionG2CC(taured);
    U := BlockMatrix([ [ gamma, 0 ], [ 0, Transpose(gamma^(-1)) ] ]);
    T := U*T;

    taured, gamma := IntegerReduceMatrixG2(taured);
    U := BlockMatrix([ [ 1, gamma ], [ 0, 1 ] ]);
    T := U*T;

    taured, U, done := GottschlingReduce(taured);
    T := U*T;
until done;

assert Abs(taured[2,1] - taured[1,2]) lt CC`epscomp;
taured[2,1] := taured[1,2];
assert IsSmallPeriodMatrix(taured);
TZZ := Matrix(Integers(), [ [ Round(c) : c in Eltseq(row) ] : row in Rows(T) ]);
assert IsSymplecticMatrix(TZZ);
return taured, T;
end function;


function LLLReduceMatrixG3(tau);
assert IsSmallPeriodMatrix(tau);
tau[2,1] := tau[1,2]; tau[3,1] := tau[1,3]; tau[3,2] := tau[2,3];
Imtau := Matrix([ [ Im(c) : c in Eltseq(row) ] : row in Rows(tau) ]);
_, T := LLLGram(Imtau);
return T*tau*Transpose(T);
end function;


function IntegerReduceMatrixG3(tau);
tauZZ := Matrix([ [ Round(Re(c)) : c in Eltseq(row) ] : row in Rows(tau) ]);
tauZZ[2,1] := tauZZ[1,2]; tauZZ[3,1] := tauZZ[1,3]; tauZZ[3,2] := tauZZ[2,3];
return tau - ChangeRing(tauZZ, BaseRing(tau));
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
    tau := LLLReduceMatrixG3(tau);
    tau := IntegerReduceMatrixG3(tau);
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
if g eq 1 then
    return ReduceSmallPeriodMatrixG1(tau);
elif g eq 2 then
    return ReduceSmallPeriodMatrixG2(tau);
elif g eq 3 then
    return ReduceSmallPeriodMatrixG3(tau);
end if;
return tau;
end intrinsic;
