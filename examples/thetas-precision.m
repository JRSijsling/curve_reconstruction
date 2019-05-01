SetVerbose("CurveRec", 0);
SetVerbose("EndoFind", 0);

prec := 300;
F := RationalsExtra(prec);
CC := F`CC;

f := x^7 + x + 1;
X := HyperellipticCurve(f);
P := PeriodMatrix(X);
tau := SmallPeriodMatrix(P);
taunew := ReduceSmallPeriodMatrix(tau);
thetas, thetas_sq := ThetaValues(taunew);

taunew0 := taunew;
thetas_sq0 := thetas_sq;

/* Now to higher precision */
prec := 400;
F := RationalsExtra(prec);
CC := F`CC;

f := x^7 + x + 1;
X := HyperellipticCurve(f);
P := PeriodMatrix(X);
tau := SmallPeriodMatrix(P);
taunew := ReduceSmallPeriodMatrix(tau);
thetas, thetas_sq := ThetaValues(taunew);

/* Comparison */
CCSmall := ComplexField(5);
maxtau := CCSmall ! Maximum([ Abs(c) : c in Eltseq(taunew - taunew0) ]);
print "";
print "Difference between period matrices:";
print maxtau;

maxtheta := CCSmall ! Maximum([ Abs(thetas_sq[i] - thetas_sq0[i]) : i in [1..#thetas_sq] ]);
print "";
print "Difference between thetas:";
print maxtheta;
