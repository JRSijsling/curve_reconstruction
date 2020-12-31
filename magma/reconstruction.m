/***
 *  Reconstruction algorithms
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

/* The conventions on period matrices are different in these programs, since
 * Guardia uses the current ones. We therefore leave the more elegant
 * conventions of Birkenhake--Lange, but only locally */


import "rosenhain.m": FindDelta;

forward ReconstructCurveGeometricG1;
forward ReconstructCurveGeometricG2;
forward ReconstructCurveGeometricG3;
forward ReconstructCurveG1;
forward ReconstructCurveG2;
forward ReconstructCurveG3;
forward AlgebraizedInvariantsG1;
forward AlgebraizedInvariantsG2;
forward AlgebraizedInvariantsG3;


function TransformForm(f, T : co := true, contra := false)
    R := Parent(f);
    vars := Matrix([ [ mon ] : mon in MonomialsOfDegree(R, 1) ]);
    if (not co) or contra then
        return Evaluate(f, Eltseq(ChangeRing(Transpose(T)^(-1), R) * vars));
    end if;
    return Evaluate(f, Eltseq(ChangeRing(T, R) * vars));
end function;


intrinsic ReconstructCurveGeometric(tau::AlgMatElt, K::Fld : Base := false) -> Crv
{Reconstruct curve from small period matrix tau. The end result will be over an extension of K.}

assert IsSmallPeriodMatrix(tau);
g := #Rows(tau);
if g eq 1 then
    return ReconstructCurveGeometricG1(tau, K : Base := Base);
elif g eq 2 then
    return ReconstructCurveGeometricG2(tau, K : Base := Base);
elif g eq 3 then
    return ReconstructCurveGeometricG3(tau, K : Base := Base);
else
    error "Genus too large!";
end if;

end intrinsic;


intrinsic ReconstructCurve(P::., K::Fld : Base := false) -> Crv
{Reconstruct curve from the big period matrix P. The end result will be over an extension of K.}

assert IsBigPeriodMatrix(P);
g := #Rows(P);
if g eq 1 then
    return ReconstructCurveG1(P, K : Base := Base);
elif g eq 2 then
    return ReconstructCurveG2(P, K : Base := Base);
elif g eq 3 then
    return ReconstructCurveG3(P, K : Base := Base);
else
    error "Genus too large!";
end if;

end intrinsic;


intrinsic AlgebraizedInvariants(tau::., K::Fld : Base := false) -> .
{Reconstruct invariants from the small period matrix tau. The end result will be over an extension of K.}

assert IsSmallPeriodMatrix(tau);
g := #Rows(tau);
if g eq 1 then
    return AlgebraizedInvariantsG1(tau, K : Base := Base);
elif g eq 2 then
    return AlgebraizedInvariantsG2(tau, K : Base := Base);
elif g eq 3 then
    return AlgebraizedInvariantsG3(tau, K : Base := Base);
else
    error "Genus too large!";
end if;

end intrinsic;


function ReconstructCurveGeometricG1(tau, K : Base := false)

jCC := jInvariant(tau[1,1]);
if Base then
    test, j := AlgebraizeElementExtra(jCC, K);
    if not test then
        vprint CurveRec : "";
        vprint CurveRec : "Failed to algebraize j-invariant.";
        return 0, 0, false;
    end if;
    hKL := CanonicalInclusionMap(K, K);
else
    L, Lj, hKL := NumberFieldExtra([ jCC ], K); j := Lj[1];
end if;

E := EllipticCurveFromjInvariant(j);
if Type(BaseRing(E)) eq FldRat then
    E := MinimalModel(E);
end if;
return E, hKL, true;

end function;


function ReconstructCurveG1(P, K : Base := false)
// Reconstruct curve from period matrix P, returned over an extension of the
// base field K.

/* Check small period matrix */
tau := SmallPeriodMatrix(P);
CC := BaseRing(Parent(tau));

/* Reduce small period matrix */
taunew, gamma := ReduceSmallPeriodMatrix(tau);
Imtaunew := Matrix([ [ Im(c) : c in Eltseq(row) ] : row in Rows(taunew) ]);

vprint CurveRec, 2 : "";
vprint CurveRec, 2 : "Eigenvalues of imaginary part of reduced tau:";
vprint CurveRec, 2 : [ ComplexField(5) ! tup[1] : tup in Eigenvalues(Imtaunew) ];

/* Calculate corresponding big period matrix */
A := Submatrix(gamma, 1,1, 1,1);
B := Submatrix(gamma, 1,2, 1,1);
C := Submatrix(gamma, 2,1, 1,1);
D := Submatrix(gamma, 2,2, 1,1);
Pnew := P * Transpose(BlockMatrix([[A, B], [C, D]]));

/* Classical elliptic functions */
CC := Parent(Pnew[1,1]); RR := RealField(CC);
assert Im(Pnew[1,1]/Pnew[1,2]) gt 0;
g4CC := 120 * (1/Pnew[1,2])^4 * ZetaFunction(RR, 4) * Eisenstein(4, Reverse(Eltseq(Pnew)));
g6CC := 280 * (1/Pnew[1,2])^6 * ZetaFunction(RR, 6) * Eisenstein(6, Reverse(Eltseq(Pnew)));

if Base then
    testg4, g4 := AlgebraizeElementExtra(g4CC, K);
    testg6, g6 := AlgebraizeElementExtra(g6CC, K);
    if not (testg4 and testg6) then
        vprint CurveRec : "";
        vprint CurveRec : "Failed to algebraize one of g4 and g6.";
        return 0, 0, false;
    end if;
    L := K; hKL := CanonicalInclusionMap(K, L);
else
    L, elts, hKL := NumberFieldExtra([ g4CC, g6CC ], K);
    g4 := elts[1]; g6 := elts[2];
end if;

R<x> := PolynomialRing(L); f := (4*x^3 - g4*x - g6)/4; h := 0;
X := HyperellipticCurve(f);

R<x> := PolynomialRing(CC);
fCC := (4*x^3 - g4CC*x - g6CC)/4; hCC := R ! 0;
YCC := RiemannSurface(fCC, 2 : Precision := Precision(CC) + 10);
Q := ChangeRing(YCC`BigPeriodMatrix, CC) / 2;

/* The next line functions as an assertion */
vprint CurveRec, 2 : "";
vprint CurveRec, 2 : "Check for isomorphism...";
A := Matrix(CC, [[1]]);
R := HomologyRepresentation(A, P, Q);
vprint CurveRec, 2 : "done.";

return X, hKL, true;

end function;


function AlgebraizedInvariantsG1(tau, K : Base := false)

jCC := jInvariant(tau[1,1]);
if Base then
    test, j := AlgebraizeElementExtra(jCC, K);
    if not test then
        vprint CurveRec : "";
        vprint CurveRec : "Failed to algebraize j-invariant";
        return 0, 0, false;
    end if;
    hKL := CanonicalInclusionMap(K, K);
else
    L, Lj, hKL := NumberFieldExtra([ jCC ], K); j := Lj[1];
end if;
return j, hKL, true;

end function;


function ReconstructCurveGeometricG2(tau, K : Base := false)
/* Alternative: implement variant of BILV */
/* TODO: Add check of not being product of elliptic curves */

CC := BaseRing(tau);
P := HorizontalJoin(tau, IdentityMatrix(CC, 2));

/* Reduce small period matrix */
taunew, gamma := ReduceSmallPeriodMatrix(tau);
Imtaunew := Matrix([ [ Im(c) : c in Eltseq(row) ] : row in Rows(taunew) ]);

vprint CurveRec, 2 : "";
vprint CurveRec, 2 : "Eigenvalues of imaginary part of reduced tau:";
vprint CurveRec, 2 : [ ComplexField(5) ! tup[1] : tup in Eigenvalues(Imtaunew) ];

/* Calculate corresponding big period matrix */
A := Submatrix(gamma, 1,1, 2,2);
B := Submatrix(gamma, 1,3, 2,2);
C := Submatrix(gamma, 3,1, 2,2);
D := Submatrix(gamma, 3,3, 2,2);
Pnew := P * Transpose(BlockMatrix([[A, B], [C, D]]));
P1new := Submatrix(Pnew, 1,1, 2,2);
P2new := Submatrix(Pnew, 1,3, 2,2);
P2inew := P2new^(-1);

/* Calculation of theta derivatives at odd two-torsion points */
w1 := (1/2)*taunew*Transpose(Matrix(CC, [[0,1]])) + (1/2)*Transpose(Matrix(CC, [[0,1]]));
w2 := (1/2)*taunew*Transpose(Matrix(CC, [[0,1]])) + (1/2)*Transpose(Matrix(CC, [[1,1]]));
w3 := (1/2)*taunew*Transpose(Matrix(CC, [[1,0]])) + (1/2)*Transpose(Matrix(CC, [[1,0]]));
w4 := (1/2)*taunew*Transpose(Matrix(CC, [[1,0]])) + (1/2)*Transpose(Matrix(CC, [[1,1]]));
w5 := (1/2)*taunew*Transpose(Matrix(CC, [[1,1]])) + (1/2)*Transpose(Matrix(CC, [[0,1]]));
w6 := (1/2)*taunew*Transpose(Matrix(CC, [[1,1]])) + (1/2)*Transpose(Matrix(CC, [[1,0]]));
ws := [ w1, w2, w3, w4, w5, w6 ];

vprint CurveRec : "";
vprint CurveRec : "Calculating theta derivatives...";
theta_derss := [ ThetaDerivatives(taunew, w) : w in ws ];
vprint CurveRec : "done calculating theta derivatives.";

/* Determination of ratios = roots */
Hs := [ Matrix(CC, [ theta_ders ]) * P2inew : theta_ders in theta_derss ];
rats := [ ];
for H in Hs do
    seq := Eltseq(H);
    add := true;
    if Abs(seq[2]) lt Abs(seq[1]) then
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

ICC := IgusaInvariants(fCC); W := [ 2, 4, 6, 8, 10 ];
ICC := WPSNormalizeCC(W, ICC);
if Base then
    test, I := AlgebraizeElementsExtra(ICC, K);
    if not test then
        vprint CurveRec : "";
        vprint CurveRec : "Failed to algebraize Igusa invariants.";
        return 0, 0, false;
    end if;
    L := K; hKL := CanonicalInclusionMap(K, L);
else
    L, I, hKL := NumberFieldExtra(ICC, K);
end if;

g2 := IgusaToG2Invariants(I);
Y := HyperellipticCurveFromG2Invariants(g2);
if Type(BaseRing(Y)) eq FldRat then
    Y := ReducedMinimalWeierstrassModel(Y);
    f, h := HyperellipticPolynomials(Y);
    g := 4*f + h^2;
    coeffs := Coefficients(g);
    d := LCM([ Denominator(coeff) : coeff in coeffs ]);
    coeffs := [ Integers() ! (d*coeff) : coeff in coeffs ];
    e := GCD(coeffs);
    coeffs := [ coeff div e : coeff in coeffs ];
    Y := HyperellipticCurve(Polynomial(coeffs));
    if Type(K) eq FldRat then
        Y := ReducedMinimalWeierstrassModel(Y);
    end if;
end if;
return Y, hKL, true;

end function;


function ReconstructCurveG2(P, K : Base := false, Dom := [-5..5])
// Reconstruct curve from period matrix P, returned over an extension of the
// base field K.
/* TODO: Add check of not being product of elliptic curves */

/* Reduce small period matrix */
tau := SmallPeriodMatrix(P);
CC := BaseRing(Parent(tau));

/* Reduce small period matrix */
taunew, gamma := ReduceSmallPeriodMatrix(tau);
Imtaunew := Matrix([ [ Im(c) : c in Eltseq(row) ] : row in Rows(taunew) ]);

vprint CurveRec, 2 : "";
vprint CurveRec, 2 : "Eigenvalues of imaginary part of reduced tau:";
vprint CurveRec, 2 : [ ComplexField(5) ! tup[1] : tup in Eigenvalues(Imtaunew) ];

/* Calculate corresponding big period matrix */
A := Submatrix(gamma, 1,1, 2,2);
B := Submatrix(gamma, 1,3, 2,2);
C := Submatrix(gamma, 3,1, 2,2);
D := Submatrix(gamma, 3,3, 2,2);
Pnew := P * Transpose(BlockMatrix([[A, B], [C, D]]));
P1new := Submatrix(Pnew, 1,1, 2,2);
P2new := Submatrix(Pnew, 1,3, 2,2);
P2inew := P2new^(-1);

/* Calculation of theta derivatives at odd two-torsion points */
w1 := (1/2)*taunew*Transpose(Matrix(CC, [[0,1]])) + (1/2)*Transpose(Matrix(CC, [[0,1]]));
w2 := (1/2)*taunew*Transpose(Matrix(CC, [[0,1]])) + (1/2)*Transpose(Matrix(CC, [[1,1]]));
w3 := (1/2)*taunew*Transpose(Matrix(CC, [[1,0]])) + (1/2)*Transpose(Matrix(CC, [[1,0]]));
w4 := (1/2)*taunew*Transpose(Matrix(CC, [[1,0]])) + (1/2)*Transpose(Matrix(CC, [[1,1]]));
w5 := (1/2)*taunew*Transpose(Matrix(CC, [[1,1]])) + (1/2)*Transpose(Matrix(CC, [[0,1]]));
w6 := (1/2)*taunew*Transpose(Matrix(CC, [[1,1]])) + (1/2)*Transpose(Matrix(CC, [[1,0]]));
ws := [ w1, w2, w3, w4, w5, w6 ];

vprint CurveRec : "";
vprint CurveRec : "Calculating theta derivatives...";
theta_derss := [ ThetaDerivatives(taunew, w) : w in ws ];
vprint CurveRec : "done calculating theta derivatives.";

/* Determination of ratios = roots */
Hs := [ Matrix(CC, [ theta_ders ]) * P2inew : theta_ders in theta_derss ];
rats := [ ];
for H in Hs do
    seq := Eltseq(H);
    add := true;
    if Abs(seq[2]) lt Abs(seq[1]) then
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

/* Finding homomorphisms to original matrix */
vprint CurveRec : "";
vprint CurveRec : "Identifying correct twist...";
YCC := RiemannSurface(fCC, 2 : Precision := Precision(CC) + 10);
Q := ChangeRing(YCC`BigPeriodMatrix, CC) / 2;
homs := GeometricHomomorphismRepresentationCC(P, Q);
As := [ hom[1] : hom in homs ]; Rs := [ hom[2] : hom in homs ];
if #As eq 0 then
    error "No geometric homomorphism to original matrix found. Try increasing the precision.";
    //: increase precision or bound in calculating theta derivatives";
end if;

/*  Identify correct twist */
M := Matrix([ [ Re(A[1,1] - A[2,2]), Im(A[1,1] - A[2,2]), Re(A[1,2]), Im(A[1,2]), Re(A[2,1]), Im(A[2,1]) ] : A in As ]);
Ker := IntegralLeftKernel(M); rows := Rows(Ker);
if #rows eq 1 then
    row := Eltseq(rows[1]);
    Lambda := &+[ row[i]*As[i] : i in [1..#As] ];
    lam := Lambda[1,1];

else
    found := false;
    CP := CartesianPower(Dom, #rows);
    for tup in CP do
        row := &+[ tup[i]*rows[i] : i in [1..#rows] ];
        Lambda := &+[ row[i]*As[i] : i in [1..#As] ];
        R := &+[ row[i]*Rs[i] : i in [1..#Rs] ];
        if Abs(Determinant(R)) eq 1 then
            lam := Lambda[1,1];
            found := true;
            break;
        end if;
    end for;

    if not found then
        error "Failed to identify correct twist. Try increasing the precision.";
    end if;
end if;
vprint CurveRec, 2 : "";
vprint CurveRec : "done identifying correct twist.";

/* Recover twisted polynomial over number field */
fCC := lam^2*fCC; coeffsCC := Coefficients(fCC);
coeffsCC := ChangeUniverse(coeffsCC, K`CC);

vprint CurveRec : "";
vprint CurveRec : "Algebraizing...";
if Base then
    test, coeffs := AlgebraizeElementsExtra(coeffsCC, K);
    if not test then
        vprint CurveRec : "Failed to algebraize coefficients.";
        return 0, 0, false;
    end if;
    L := K; hKL := CanonicalInclusionMap(K, L);
else
    L, coeffs, hKL := NumberFieldExtra(coeffsCC, K);
end if;
vprint CurveRec : "done.";

R := PolynomialRing(L);
f := &+[ coeffs[i]*R.1^(i - 1) : i in [1..#coeffs] ];
Y := HyperellipticCurve(f);
YCC := RiemannSurface(fCC, 2 : Precision := Precision(CC) + 10);
Q := ChangeRing(YCC`BigPeriodMatrix, CC) / 2;

/* The next line functions as an assertion */
vprint CurveRec, 2 : "";
vprint CurveRec, 2 : "Check for isomorphism...";
A := IdentityMatrix(CC, 2);
R := HomologyRepresentation(A, P, Q);
vprint CurveRec, 2 : "done.";
return Y, hKL, true;

end function;


function AlgebraizedInvariantsG2(tau, K : Base := false)

CC := BaseRing(tau);
P := HorizontalJoin(tau, IdentityMatrix(CC, 2));

/* Reduce small period matrix */
taunew, gamma := ReduceSmallPeriodMatrix(tau);
Imtaunew := Matrix([ [ Im(c) : c in Eltseq(row) ] : row in Rows(taunew) ]);

vprint CurveRec, 2 : "";
vprint CurveRec, 2 : "Eigenvalues of imaginary part of reduced tau:";
vprint CurveRec, 2 : [ ComplexField(5) ! tup[1] : tup in Eigenvalues(Imtaunew) ];

/* Calculate corresponding big period matrix */
A := Submatrix(gamma, 1,1, 2,2);
B := Submatrix(gamma, 1,3, 2,2);
C := Submatrix(gamma, 3,1, 2,2);
D := Submatrix(gamma, 3,3, 2,2);
Pnew := P * Transpose(BlockMatrix([[A, B], [C, D]]));
P1new := Submatrix(Pnew, 1,1, 2,2);
P2new := Submatrix(Pnew, 1,3, 2,2);
P2inew := P2new^(-1);

/* Calculation of theta derivatives at odd two-torsion points */
w1 := (1/2)*taunew*Transpose(Matrix(CC, [[0,1]])) + (1/2)*Transpose(Matrix(CC, [[0,1]]));
w2 := (1/2)*taunew*Transpose(Matrix(CC, [[0,1]])) + (1/2)*Transpose(Matrix(CC, [[1,1]]));
w3 := (1/2)*taunew*Transpose(Matrix(CC, [[1,0]])) + (1/2)*Transpose(Matrix(CC, [[1,0]]));
w4 := (1/2)*taunew*Transpose(Matrix(CC, [[1,0]])) + (1/2)*Transpose(Matrix(CC, [[1,1]]));
w5 := (1/2)*taunew*Transpose(Matrix(CC, [[1,1]])) + (1/2)*Transpose(Matrix(CC, [[0,1]]));
w6 := (1/2)*taunew*Transpose(Matrix(CC, [[1,1]])) + (1/2)*Transpose(Matrix(CC, [[1,0]]));
ws := [ w1, w2, w3, w4, w5, w6 ];

vprint CurveRec : "";
vprint CurveRec : "Calculating theta derivatives...";
theta_derss := [ ThetaDerivatives(taunew, w) : w in ws ];
vprint CurveRec : "done calculating theta derivatives.";

/* Determination of ratios = roots */
Hs := [ Matrix(CC, [ theta_ders ]) * P2inew : theta_ders in theta_derss ];
rats := [ ];
for H in Hs do
    seq := Eltseq(H);
    add := true;
    if Abs(seq[2]) lt Abs(seq[1]) then
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

ICC := IgusaInvariants(fCC); W := [ 2, 4, 6, 8, 10 ];
ICC := WPSNormalizeCC(W, ICC);
if Base then
    test, I := AlgebraizeElementsExtra(ICC, K);
    if not test then
        vprint CurveRec : "";
        vprint CurveRec : "Failed to algebraize Igusa invariants.";
        return 0, 0, false;
    end if;
    L := K; hKL := CanonicalInclusionMap(K, L);
else
    L, I, hKL := NumberFieldExtra(ICC, K);
end if;
return I, hKL, true;

end function;


function ReconstructCurveGeometricG3(tau, K : Base := Base)

/* Calculate thetas and see in which case we are */
taunew := ReduceSmallPeriodMatrix(tau);
assert IsSmallPeriodMatrix(taunew);
vprint CurveRec, 2 : "";
vprint CurveRec, 2: "Calculating theta values...";
thetas, thetas_sq := ThetaValues(taunew);
vprint CurveRec, 2: "done.";

v0s := FindDelta(thetas_sq);
vprint CurveRec, 2 : "";
vprint CurveRec, 2: "Number of non-zero even theta values:";
vprint CurveRec, 2: #v0s;

if #v0s eq 0 then
    ICC := DixmierOhnoInvariantsFromThetas(thetas);

    if Base then
        test, I := AlgebraizeElementsExtra(ICC, K);
        if not test then
            vprint CurveRec : "";
            vprint CurveRec : "Failed to algebraize Dixmier--Ohno invariants.";
            return 0, 0, false;
        end if;
        L := K; hKL := CanonicalInclusionMap(K, L);

    else
        L, I, hKL := NumberFieldExtra(ICC, K);
    end if;
    f, aut := TernaryQuarticFromDixmierOhnoInvariants(I);
    if #aut eq 0 then
        f := MinimizeC2Quartic(f);
    end if;
    return PlaneCurve(f), hKL, true;

elif #v0s eq 1 then
    ICC := ShiodaInvariantsFromThetaSquares(thetas_sq);

    if Base then
        test, I := AlgebraizeElementsExtra(ICC, K);
        if not test then
            vprint CurveRec : "";
            vprint CurveRec : "Failed to algebraize Shioda invariants.";
            return 0, 0, false;
        end if;
        L := K; hKL := CanonicalInclusionMap(K, L);

    else
        L, I, hKL := NumberFieldExtra(ICC, K);
    end if;
    Y := HyperellipticCurveFromShiodaInvariants(I);
    return Y, hKL, true;

else
    vprint CurveRec : "";
    vprint CurveRec : "Too many even theta characteristics vanish.";
    return 0, 0, false;
end if;

end function;


function ReconstructCurveG3(P, K : Base := Base)
/* Only for plane quartic curves currently, hyperelliptic curves soon to follow */

/* Reduce small period matrix */

tau := SmallPeriodMatrix(P);
CC := BaseRing(Parent(tau));

Y, hKL, test := ReconstructCurveGeometricG3(tau, K : Base := Base);
if not test then
    return 0, 0, false;
end if;
vprint CurveRec, 2 : "";
vprint CurveRec, 2 : "Geometric reconstruction:";
vprint CurveRec, 2 : Y;

g := DefiningPolynomial(Y);
vprint CurveRec, 2 : "";
vprint CurveRec, 2 : "Determining period matrix of geometric reconstruction...";
Q := PeriodMatrix(Y);
if Type(Y) eq CrvHyp then
    vprint CurveRec : "";
    vprint CurveRec : "Arithmetic reconstruction not yet possible for hyperelliptic curves.";
    return 0, 0, false;
end if;
vprint CurveRec, 2 : "done determining period matrix of geometric reconstruction.";

vprint CurveRec, 2 : "";
vprint CurveRec, 2 : "Determining isomorphisms with original period matrix...";
isos := SymplecticIsomorphismsCC(P, Q);
vprint CurveRec, 2 : "done.";

T := isos[1][1];
gCC := EmbedPolynomialExtra(g);
fCC := TransformForm(gCC, T);
CC := BaseRing(Parent(fCC));
coeffs := [ c : c in Coefficients(fCC) | Abs(c) gt CC`epscomp ];
min, ind := Minimum([ Abs(c) : c in coeffs ]);
fCC /:= coeffs[1];

monsCC := Monomials(fCC);
coeffsCC := [ MonomialCoefficient(fCC, monCC) : monCC in monsCC ];
exps := [ Exponents(monCC) : monCC in monsCC ];

if Base then
    test, coeffs := AlgebraizeElementsExtra(coeffsCC, K);
    if not test then
        vprint CurveRec : "";
        vprint CurveRec : "Failed to algebraize coefficients.";
        return 0, 0, false;
    end if;
    L := K; hKL := CanonicalInclusionMap(K, L);
else
    L, coeffs, hKL := NumberFieldExtra(coeffsCC, K);
end if;

R := PolynomialRing(L, 3);
F0 := &+[ coeffs[i]*Monomial(R, exps[i]) : i in [1..#coeffs] ];
X0 := PlaneCurve(F0);
return X0, hKL, true;

end function;


function AlgebraizedInvariantsG3(tau, K : Base := false)

/* Calculate thetas and see in which case we are */
taunew := ReduceSmallPeriodMatrix(tau);
assert IsSmallPeriodMatrix(taunew);
vprint CurveRec, 2 : "";
vprint CurveRec, 2: "Calculating theta values...";
thetas, thetas_sq := ThetaValues(taunew);
vprint CurveRec, 2: "done";

v0s := FindDelta(thetas_sq);
vprint CurveRec, 2 : "";
vprint CurveRec, 2: "Number of non-zero even theta values:";
vprint CurveRec, 2: #v0s;

if #v0s gt 1 then
    vprint CurveRec : "";
    vprint CurveRec : "At least two vanishing even characteristics.";
    return 0, 0, false;
end if;

if #v0s eq 0 then
    ICC := DixmierOhnoInvariantsFromThetas(thetas);
else
    ICC := ShiodaInvariantsFromThetaSquares(thetas_sq);
end if;

if Base then
    test, I := AlgebraizeElementsExtra(ICC, K);
    if not test then
        vprint CurveRec : "";
        vprint CurveRec : "Failed to algebraize invariants.";
        return 0, 0, false;
    end if;
    L := K; hKL := CanonicalInclusionMap(K, L);
else
    L, I, hKL := NumberFieldExtra(ICC, K);
end if;
return I, hKL, true;

end function;
