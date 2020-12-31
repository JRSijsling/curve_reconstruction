/***
 *  Riemann models
 *
 *  Written by Reynald Lercier
 *             and Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *
 *  See LICENSE.txt for license details.
 */


function ModuliFromTheta(thetas);
I:=Parent(thetas[1]).1;
a1:=I*thetas[33]*thetas[5]/(thetas[40]*thetas[12]);
a2:=I*thetas[21]*thetas[49]/(thetas[28]*thetas[56]);
a3:=I*thetas[7]*thetas[35]/(thetas[14]*thetas[42]);
ap1:=I*thetas[5]*thetas[54]/(thetas[27]*thetas[40]);
ap2:=I*thetas[49]*thetas[2]/(thetas[47]*thetas[28]);
ap3:=I*thetas[35]*thetas[16]/(thetas[61]*thetas[14]);
as1:=-thetas[54]*thetas[33]/(thetas[12]*thetas[27]);
as2:=thetas[2]*thetas[21]/(thetas[56]*thetas[47]);
as3:=thetas[16]*thetas[7]/(thetas[42]*thetas[61]);
return [a1,a2,a3,ap1,ap2,ap3,as1,as2,as3];
end function;


function RiemannModelFromModuli(mods);
a1:=mods[1];a2:=mods[2];a3:=mods[3];
ap1:=mods[4];ap2:=mods[5];ap3:=mods[6];
as1:=mods[7];as2:=mods[8];as3:=mods[9];
F:=Parent(a1);
P<x1,x2,x3>:=PolynomialRing(F,3);
k:=1;kp:=1;ks:=1;
M:=Matrix([[1,1,1],[k*a1,k*a2,k*a3],[kp*ap1,kp*ap2,kp*ap3]]);
Mb:=Matrix([[1,1,1],[1/a1,1/a2,1/a3],[1/ap1,1/ap2,1/ap3]]);
U:=-Mb^(-1)*M;
u1:=U[1];
u2:=U[2];
u3:=U[3];
u1:=u1[1]*x1+u1[2]*x2+u1[3]*x3;
u2:=u2[1]*x1+u2[2]*x2+u2[3]*x3;
u3:=u3[1]*x1+u3[2]*x2+u3[3]*x3;
return (x1*u1+x2*u2-x3*u3)^2-4*x1*u1*x2*u2;
end function;


intrinsic DixmierOhnoInvariantsFromThetas(thetas::.) -> .
{Calculate Dixmier--Ohno invariants.}
modsCC := ModuliFromTheta(thetas);
ICC, W := DixmierOhnoInvariants(RiemannModelFromModuli(modsCC));
return WPSNormalizeCC(W, ICC);
end intrinsic;


intrinsic DixmierOhnoInvariants(X::Crv) -> .
{Returns the Dixmier--Ohno invariants of the curve X.}

if not Type(X) eq CrvPln then
    error "Input must be a plane curve";
end if;
if not Degree(DefiningPolynomial(X)) eq 4 then
    error "Input must be a quartic curve";
end if;
return DixmierOhnoInvariants(DefiningPolynomial(X));

end intrinsic;
