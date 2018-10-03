/***
 *  Reconstruction algorithms
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


import "riemann.m": DixmierOhnoInvariantsFromThetaSquares;
import "rosenhain.m": FindDelta, ShiodaInvariantsFromThetaSquares;


function RationalReconstruction(r);
    r1:=Real(r);
    r2:=Imaginary(r);
    p:=Precision(r2);
    if r2 ne Parent(r2)!0 then
        e:=Log(AbsoluteValue(r2));
    else
        e:=-p;
    end if;
    if -e lt p/2 then
        print "not a real";
        return false;
    end if;
    best:=0;
    i:=p div 10;
    b:=BestApproximation(r1,10^i);
    while b ne best and i le p do
        i:=i+5;
        best:=b;
        b:=BestApproximation(r1,10^i);
        //print b;
    end while;
    if b ne best then
        print "not enough precision";
        return false;
    else
        return b,i;
    end if;
end function;


function ReconstructCurveGeometricG1(tau, K)

assert IsSmallPeriodMatrix(tau);
jCC := jInvariant(tau[1,1]);
j := AlgebraizeElementLLL(jCC, K);
E := EllipticCurveFromjInvariant(j);
if Type(K) eq FldRat then
    E := MinimalModel(E);
end if;
return E;

end function;


function ReconstructCurveGeometricG2(tau, K)
/* Alternative: implement BILV */

assert IsSmallPeriodMatrix(tau);
CC := BaseRing(tau);
P := HorizontalJoin(IdentityMatrix(CC, 2), tau);

/* Reduce small period matrix */
taunew, gamma := ReduceSmallPeriodMatrix(tau);
vprint CurveRec : "Eigenvalues of reduced tau:";
vprint CurveRec : [ ComplexField(5) ! tup[1] : tup in Eigenvalues(taunew) ];

/* Calculate corresponding big period matrix */
A := Transpose(Submatrix(gamma, 1,1, 2,2));
B := Transpose(Submatrix(gamma, 1,3, 2,2));
C := Transpose(Submatrix(gamma, 3,1, 2,2));
D := Transpose(Submatrix(gamma, 3,3, 2,2));
Pnew := P * BlockMatrix([[D, B], [C, A]]);
P1new := Submatrix(Pnew, 1,1, 2,2); P1inew := P1new^(-1);
P2new := Submatrix(Pnew, 1,3, 2,2);

//print ChangeRing(taunew - P1new^(-1)*P2new, ComplexField(5));
//print [ ComplexField(5) ! v[1] : v in Eigenvalues(tau) ];
//print [ ComplexField(5) ! v[1] : v in Eigenvalues(taunew) ];

/* Calculation of theta derivatives at odd two-torsion points */
w1 := (1/2)*taunew*Transpose(Matrix(CC, [[0,1]])) + (1/2)*Transpose(Matrix(CC, [[0,1]]));
w2 := (1/2)*taunew*Transpose(Matrix(CC, [[0,1]])) + (1/2)*Transpose(Matrix(CC, [[1,1]]));
w3 := (1/2)*taunew*Transpose(Matrix(CC, [[1,0]])) + (1/2)*Transpose(Matrix(CC, [[1,0]]));
w4 := (1/2)*taunew*Transpose(Matrix(CC, [[1,0]])) + (1/2)*Transpose(Matrix(CC, [[1,1]]));
w5 := (1/2)*taunew*Transpose(Matrix(CC, [[1,1]])) + (1/2)*Transpose(Matrix(CC, [[0,1]]));
w6 := (1/2)*taunew*Transpose(Matrix(CC, [[1,1]])) + (1/2)*Transpose(Matrix(CC, [[1,0]]));
ws := [ w1, w2, w3, w4, w5, w6 ];
vprint CurveRec : "Calculate theta derivatives...";
theta_derss := [ ThetaDerivatives(taunew, w) : w in ws ];
vprint CurveRec : "done.";
Hs := [ Matrix(CC, [ theta_ders ]) * P1inew : theta_ders in theta_derss ];

/* Determination of ratios = roots */
rats := [ ];
for H in Hs do
    seq := Eltseq(H);
    add := true;
    if Abs(seq[2]) lt Abs(seq[1]) then
        //if Abs(seq[2]/seq[1]) lt 10^(-10) then
        if Abs(seq[2]/seq[1])^2 lt CC`epscomp then
            add := false;
        end if;
    end if;
    if add then
        Append(~rats, -seq[1]/seq[2]);
    end if;
end for;

/* Recover polynomial over CC up to a constant */
RCC := PolynomialRing(CC); R := PolynomialRing(K);
fCC := &*[ RCC.1 - rat : rat in rats ];
vprint CurveRec, 2 : "Reconstructed polynomial up to a constant:";
vprint CurveRec, 2 : fCC;

ICC := IgusaInvariants(fCC); W := [ 2, 4, 6, 8, 10 ];
ICC := WPSNormalizeCC(W, ICC);
Knew := K; Knew`base := K`base; Knew`base_gen := Knew ! K`base_gen; Knew`iota := K`iota;
//Knew`CC := ComplexFieldExtra(Precision(Knew`CC) - 50);
ICC := ChangeUniverse(ICC, Knew`CC);
I := [ AlgebraizeElementLLL(inv, Knew) : inv in ICC ];
g2 := IgusaToG2Invariants(I);
Y := HyperellipticCurveFromG2Invariants(g2);
if Type(Knew) eq FldRat then
    Y := ReducedMinimalWeierstrassModel(Y);
    f, h := HyperellipticPolynomials(Y);
    g := 4*f + h^2;
    coeffs := Coefficients(g);
    d := LCM([ Denominator(coeff) : coeff in coeffs ]);
    coeffs := [ Integers() ! (d*coeff) : coeff in coeffs ];
    e := GCD(coeffs);
    coeffs := [ coeff div e : coeff in coeffs ];
    Y := HyperellipticCurve(Polynomial(coeffs));
    Y := ReducedMinimalWeierstrassModel(Y);
end if;
return Y;

end function;


function ReconstructCurveGeometricG3(tau, K)
/* Alternative: implement Guardia */

/* Calculate thetas and see in which case we are */
assert IsSmallPeriodMatrix(tau);
thetas := ThetaValues(tau);
thetas_sq := [ theta^2 : theta in thetas ];
v0s := FindDelta(thetas_sq);
if #v0s eq 0 then
    error "Uncomment relevant section of reconstruction.m";
    /*
    ICC := DixmierOhnoInvariantsFromThetaSquares(thetas);
    Knew := K; Knew`base := K`base; Knew`base_gen := Knew ! K`base_gen; Knew`iota := K`iota;
    Knew`CC := ComplexFieldExtra(Precision(Knew`CC) - 50);
    ICC := ChangeUniverse(ICC, Knew`CC);
    // TODO: Does not work over extensions right now, the heights involved are too large
    //K, I := NumberFieldExtra(ICC, F);
    I := [ RationalReconstruction(inv) : inv in ICC ];
    return TernaryQuarticFromDixmierOhnoInvariants(I);
    */
elif #v0s eq 1 then
    SCC := ShiodaInvariantsFromThetaSquares(thetas_sq);
    Knew := K; Knew`base := K`base; Knew`base_gen := Knew ! K`base_gen; Knew`iota := K`iota;
    //Knew`CC := ComplexFieldExtra(Precision(Knew`CC) - 50);
    SCC := ChangeUniverse(SCC, Knew`CC);
    S := [ RationalReconstruction(inv) : inv in SCC ];
    return HyperellipticCurveFromShiodaInvariants(S);
else
    error "Too many even characteristics vanish", #v0s;
end if;

end function;


intrinsic ReconstructCurveGeometric(tau::AlgMatElt, K::Fld) -> Crv
{Reconstruct curve from small period matrix tau. The end result will be over the given field K.}

g := #Rows(tau);
if g eq 1 then
    return ReconstructCurveGeometricG1(tau, K);
elif g eq 2 then
    return ReconstructCurveGeometricG2(tau, K);
elif g eq 3 then
    return ReconstructCurveGeometricG3(tau, K);
else
    error "Genus too large!";
end if;

end intrinsic;


function ReconstructCurveBaseG1(P, K)

/* Check small period matrix */
P1 := Submatrix(P, 1,1, 1,1); P1i := P1^(-1);
P2 := Submatrix(P, 1,2, 1,1);
tau := P1i*P2; CC := BaseRing(Parent(tau));
assert IsSmallPeriodMatrix(tau);

/* Classical elliptic functions */
CC := Parent(P[1,1]); RR := RealField(CC);
if Im(P[1,2]/P[1,1]) lt 0 then
    P := [ P[1,2], P[1,1] ];
end if;
g4CC := 120 * (1/P[1,1])^4 * ZetaFunction(RR, 4) * Eisenstein(4, Eltseq(P));
g6CC := 280 * (1/P[1,1])^6 * ZetaFunction(RR, 6) * Eisenstein(6, Eltseq(P));
g4 := AlgebraizeElementLLL(g4CC, K); g6 := AlgebraizeElementLLL(g6CC, K);
R<x> := PolynomialRing(K); f := (4*x^3 - g4*x - g6)/4; h := 0;
X := HyperellipticCurve(f, h);

R<x> := PolynomialRing(CC);
fCC := (4*x^3 - g4CC*x - g6CC)/4; hCC := R ! 0;
Q := ChangeRing(PeriodMatrix([ fCC, hCC ], [ f, h ]), CC);
abs := Max([ Abs(c) : c in Eltseq(P - Q) ]);
assert abs lt CC`epscomp;
return X;

end function;


function ReconstructCurveBaseG2(P, K)

/* Reduce small period matrix */
P1 := Submatrix(P, 1,1, 2,2); P1i := P1^(-1);
P2 := Submatrix(P, 1,3, 2,2);
tau := P1i*P2; CC := BaseRing(Parent(tau));
assert IsSmallPeriodMatrix(tau);

/* Reduce small period matrix */
taunew, gamma := ReduceSmallPeriodMatrix(tau);
vprint CurveRec : "Eigenvalues of reduced tau:";
vprint CurveRec : [ ComplexField(5) ! tup[1] : tup in Eigenvalues(taunew) ];

/* Calculate corresponding big period matrix */
A := Transpose(Submatrix(gamma, 1,1, 2,2));
B := Transpose(Submatrix(gamma, 1,3, 2,2));
C := Transpose(Submatrix(gamma, 3,1, 2,2));
D := Transpose(Submatrix(gamma, 3,3, 2,2));
Pnew := P * BlockMatrix([[D, B], [C, A]]);
P1new := Submatrix(Pnew, 1,1, 2,2); P1inew := P1new^(-1);
P2new := Submatrix(Pnew, 1,3, 2,2);

//print ChangeRing(taunew - P1new^(-1)*P2new, ComplexField(5));
//print [ ComplexField(5) ! v[1] : v in Eigenvalues(tau) ];
//print [ ComplexField(5) ! v[1] : v in Eigenvalues(taunew) ];

/* Calculation of theta derivatives at odd two-torsion points */
w1 := (1/2)*taunew*Transpose(Matrix(CC, [[0,1]])) + (1/2)*Transpose(Matrix(CC, [[0,1]]));
w2 := (1/2)*taunew*Transpose(Matrix(CC, [[0,1]])) + (1/2)*Transpose(Matrix(CC, [[1,1]]));
w3 := (1/2)*taunew*Transpose(Matrix(CC, [[1,0]])) + (1/2)*Transpose(Matrix(CC, [[1,0]]));
w4 := (1/2)*taunew*Transpose(Matrix(CC, [[1,0]])) + (1/2)*Transpose(Matrix(CC, [[1,1]]));
w5 := (1/2)*taunew*Transpose(Matrix(CC, [[1,1]])) + (1/2)*Transpose(Matrix(CC, [[0,1]]));
w6 := (1/2)*taunew*Transpose(Matrix(CC, [[1,1]])) + (1/2)*Transpose(Matrix(CC, [[1,0]]));
ws := [ w1, w2, w3, w4, w5, w6 ];
vprint CurveRec : "Calculate theta derivatives...";
theta_derss := [ ThetaDerivatives(taunew, w) : w in ws ];
vprint CurveRec : "done.";
Hs := [ Matrix(CC, [ theta_ders ]) * P1inew : theta_ders in theta_derss ];

/* Determination of ratios = roots */
rats := [ ];
for H in Hs do
    seq := Eltseq(H);
    add := true;
    if Abs(seq[2]) lt Abs(seq[1]) then
        //if Abs(seq[2]/seq[1]) lt 10^(-10) then
        if Abs(seq[2]/seq[1])^2 lt CC`epscomp then
            add := false;
        end if;
    end if;
    if add then
        Append(~rats, -seq[1]/seq[2]);
    end if;
end for;

/* Recover polynomial over CC up to a constant */
RCC := PolynomialRing(CC); R := PolynomialRing(K);
fCC := &*[ RCC.1 - rat : rat in rats ];
vprint CurveRec, 2 : "Reconstructed polynomial up to a constant:";
vprint CurveRec, 2 : fCC;

/* Identify correct twist */
Y := SE_Curve(fCC, 2 : Prec := Precision(CC));
Q := ChangeRing(Y`BigPeriodMatrix, CC) / 2;
homs := GeometricHomomorphismRepresentationCC(P, Q);
As := [ hom[1] : hom in homs ];
M := Matrix([ [ Re(A[1,1] - A[2,2]), Im(A[1,1] - A[2,2]), Re(A[1,2]), Im(A[1,2]), Re(A[2,1]), Im(A[2,1]) ] : A in As ]);
Ker := IntegralLeftKernel(M);
row := Eltseq(Rows(Ker)[1]);
Lambda := &+[ row[i]*As[i] : i in [1..#As] ];
lam := Lambda[1,1];

/* Recover twisted polynomial over number field */
fCC := lam^2*fCC;
coeffsCC := Coefficients(fCC);
coeffs := [ AlgebraizeElementLLL(coeffCC, K) : coeffCC in coeffsCC ];
f := &+[ coeffs[i]*R.1^(i - 1) : i in [1..#coeffs] ];
Y := HyperellipticCurve(f);

Q := ChangeRing(PeriodMatrix([ fCC ], [ f ]), CC);
/* The next line functions as an assertion */
R := HomologyRepresentation(IdentityMatrix(CC, 2), P, Q);
return Y;

end function;


function ReconstructCurveBaseG3(P, K)

/* Reduce small period matrix */
P1 := Submatrix(P, 1,1, 3,3); P1i := P1^(-1);
P2 := Submatrix(P, 1,4, 3,3);
tau := P1i*P2; CC := BaseRing(Parent(tau));
assert IsSmallPeriodMatrix(tau);
error "Guardia not implemented yet!";
return 0;

end function;


intrinsic ReconstructCurveBase(P::ModMatRngElt, K::Fld) -> Crv
{Reconstruct curve from big period matrix P. The end result will be over the given field K.}

g := #Rows(P);
if g eq 1 then
    return ReconstructCurveBaseG1(P, K);
elif g eq 2 then
    return ReconstructCurveBaseG2(P, K);
elif g eq 3 then
    return ReconstructCurveBaseG3(P, K);
else
    error "Genus too large!";
end if;

end intrinsic;


intrinsic ReconstructCurveG1(P::ModMatRngElt, F::Fld) -> Crv
{Reconstruct curve from period matrix P over field F.}

P := Eltseq(P); CC := Parent(P[1]); RR := RealField(CC);
if Im(P[2]/P[1]) lt 0 then
    P := [ P[2], P[1] ];
end if;
g4CC := 120 * (1/P[1])^4 * ZetaFunction(RR, 4) * Eisenstein(4, P);
g6CC := 280 * (1/P[1])^6 * ZetaFunction(RR, 6) * Eisenstein(6, P);
//F`CC := ComplexFieldExtra(Precision(F`CC) - 50);
K, elts := NumberFieldExtra([ g4CC, g6CC ], F);
g4 := elts[1]; g6 := elts[2];
R<x> := PolynomialRing(K); f := (4*x^3 - g4*x - g6)/4; h := 0;
return HyperellipticCurve(f, h);

end intrinsic;


intrinsic ReconstructCurveG2(P::ModMatRngElt, F::Fld) -> Crv
{Reconstruct curve from period matrix P over field F.}

/* Reduce small period matrix */
P1 := Submatrix(P, 1,1, 2,2); P1i := P1^(-1);
P2 := Submatrix(P, 1,3, 2,2);
tau := P1i*P2; CC := BaseRing(Parent(tau));
assert IsSmallPeriodMatrix(tau);

/* Reduce small period matrix */
taunew, gamma := ReduceSmallPeriodMatrix(tau);
vprint CurveRec : "Eigenvalues of reduced tau:";
vprint CurveRec : [ ComplexField(5) ! tup[1] : tup in Eigenvalues(taunew) ];

/* Calculate corresponding big period matrix */
A := Transpose(Submatrix(gamma, 1,1, 2,2));
B := Transpose(Submatrix(gamma, 1,3, 2,2));
C := Transpose(Submatrix(gamma, 3,1, 2,2));
D := Transpose(Submatrix(gamma, 3,3, 2,2));
Pnew := P * BlockMatrix([[D, B], [C, A]]);
P1new := Submatrix(Pnew, 1,1, 2,2); P1inew := P1new^(-1);
P2new := Submatrix(Pnew, 1,3, 2,2);

//print ChangeRing(taunew - P1new^(-1)*P2new, ComplexField(5));
//print [ ComplexField(5) ! v[1] : v in Eigenvalues(tau) ];
//print [ ComplexField(5) ! v[1] : v in Eigenvalues(taunew) ];

/* Calculation of theta derivatives at odd two-torsion points */
w1 := (1/2)*taunew*Transpose(Matrix(CC, [[0,1]])) + (1/2)*Transpose(Matrix(CC, [[0,1]]));
w2 := (1/2)*taunew*Transpose(Matrix(CC, [[0,1]])) + (1/2)*Transpose(Matrix(CC, [[1,1]]));
w3 := (1/2)*taunew*Transpose(Matrix(CC, [[1,0]])) + (1/2)*Transpose(Matrix(CC, [[1,0]]));
w4 := (1/2)*taunew*Transpose(Matrix(CC, [[1,0]])) + (1/2)*Transpose(Matrix(CC, [[1,1]]));
w5 := (1/2)*taunew*Transpose(Matrix(CC, [[1,1]])) + (1/2)*Transpose(Matrix(CC, [[0,1]]));
w6 := (1/2)*taunew*Transpose(Matrix(CC, [[1,1]])) + (1/2)*Transpose(Matrix(CC, [[1,0]]));
ws := [ w1, w2, w3, w4, w5, w6 ];
vprint CurveRec : "Calculate theta derivatives...";
theta_derss := [ ThetaDerivatives(taunew, w) : w in ws ];
vprint CurveRec : "done.";
Hs := [ Matrix(CC, [ theta_ders ]) * P1inew : theta_ders in theta_derss ];

/* Determination of ratios = roots */
rats := [ ];
for H in Hs do
    seq := Eltseq(H);
    add := true;
    if Abs(seq[2]) lt Abs(seq[1]) then
        //if Abs(seq[2]/seq[1]) lt 10^(-10) then
        if Abs(seq[2]/seq[1])^2 lt CC`epscomp then
            add := false;
        end if;
    end if;
    if add then
        Append(~rats, -seq[1]/seq[2]);
    end if;
end for;

/* Recover polynomial over CC up to a constant */
RCC := PolynomialRing(CC);
fCC := &*[ RCC.1 - rat : rat in rats ];
vprint CurveRec, 2 : "Reconstructed polynomial up to a constant:";
vprint CurveRec, 2 : fCC;

/* Identify correct twist */
Y := SE_Curve(fCC, 2 : Prec := Precision(CC));
Q := ChangeRing(Y`BigPeriodMatrix, CC) / 2;
homs := GeometricHomomorphismRepresentationCC(P, Q);
As := [ hom[1] : hom in homs ];
M := Matrix([ [ Re(A[1,1] - A[2,2]), Im(A[1,1] - A[2,2]), Re(A[1,2]), Im(A[1,2]), Re(A[2,1]), Im(A[2,1]) ] : A in As ]);
Ker := IntegralLeftKernel(M);
row := Eltseq(Rows(Ker)[1]);
Lambda := &+[ row[i]*As[i] : i in [1..#As] ];
lam := Lambda[1,1];

/* Recover twisted polynomial over number field */
fCC := lam^2*fCC; coeffsCC := Coefficients(fCC);
coeffsCC := ChangeUniverse(coeffsCC, F`CC);
K, coeffs := NumberFieldExtra(coeffsCC, F); R := PolynomialRing(K);
f := &+[ coeffs[i]*R.1^(i - 1) : i in [1..#coeffs] ];
Y := HyperellipticCurve(f);

Q := ChangeRing(PeriodMatrix([ fCC ], [ f ]), CC);
/* The next line functions as an assertion */
R := HomologyRepresentation(IdentityMatrix(CC, 2), P, Q);
return Y;

end intrinsic;
