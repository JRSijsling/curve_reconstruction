/***
 *  Reconstruction algorithms
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic ThetaDerivatives(tau::AlgMatElt, w::ModMatRngElt : B := 30) -> SeqEnum
{Calculates the derivative of the standard theta function at w with respect to all variables z_i.}
/* TODO: Makes this fast after Labrande--Thomé */

g := #Rows(tau); D := [-B..B];
CC<I> := BaseRing(Parent(tau)); pi := Pi(CC);
CP := CartesianPower(D, g);
theta_ders := [ CC ! 0 : i in [1..g] ];
for tup in CP do
    n := Transpose(Matrix(CC, [[ c : c in tup ]]));
    quad := Exp(pi*I*(Transpose(n)*tau*n)[1,1]);
    lin := Exp(2*pi*I*(Transpose(n)*w)[1,1]);
    main := quad*lin;
    for i in [1..g] do
        theta_ders[i] +:= 2*pi*I*tup[i]*main;
    end for;
end for;
return theta_ders;

end intrinsic;
