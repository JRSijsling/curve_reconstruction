Description
--

This repository contains Magma code for reconstructing hyperelliptic curves of genus up to 3 from their period matrices, both geometrically and arithmetically. With some extra work, hyperelliptic curves of arbitrary genus should be feasible. However, the repository is not being actively developed in that direction yet.

Prerequisites
--

An installation of Magma and the dependencies [`edgarcosta/endomorphisms`](https://github.com/edgarcosta/endomorphisms) and [`JRSijsling/quartic`](https://github.com/JRSijsling/quartic).

Installation
--

You can enable the functionality of this package in Magma by attaching the `curve_reconstruction/magma/spec` file with `AttachSpec`. To make this independent of the directory in which you find yourself, and to active this on startup by default, you may want to indicate the relative path in your `~/.magmarc` file, by adding the line
```
AttachSpec("~/Programs/curve_reconstruction/magma/spec");
```

Usage
--

Examples are given in the directory `examples/`.

Verbose comments are enabled by
```
SetVerbose("CurveRec", n);
```
where and `n` is either `1` or `2`. A higher value gives more comments.

Credits
--

This implementation is based on the following works. When using this package, please be aware of the work that you are indirectly applying and please cite it.

For geometric reconstruction in genus 3 (also see their own SageMath implementation at [`christellevincent/genus3`](https://github.com/christellevincent/genus3):

Jennifer Balakrishnan, Sorina Ionica, Kristin Lauter, and Christelle Vincent  
*Constructing genus-3 hyperelliptic Jacobians with CM* (English summary)  
LMS J. Comput. Math. 19 (2016), suppl. A, pp. 283-–300.

For arithmetic reconstruction in genus 2:

Jordi Guàrdia  
*Jacobian Nullwerte and algebraic equations*  
Journal of Algebra 253 (2002) 112–132.

For the fast computation of theta constants needed when geometrically reconstructing in genus 2 or 3:

Hugo Labrande and Emmanuel Thomé  
*Computing theta functions in quasi-linear time in genus 2 and above*  
LMS J. Comput. Math. 19 (2016), suppl. A, pp. 163–-177.

For the calculation of period matrices of plane quartic curves used for arithmetic reconstruction:

Christian Neurohr  
*Efficient integration on Riemann surfaces & applications*  
Ph.D. thesis, Carl-von-Ossietzky-Universität Oldenburg (2018)

In more ad-hoc form, some of the methods in this package were used in

Pinar Kılıçer, Hugo Labrande, Reynald Lercier, Christophe Ritzenthaler, Jeroen Sijsling, Marco Streng  
*Plane quartics over QQ with complex multiplication*  
Acta Arith. 185 (2018), no. 2, 127-156

Finally, this work uses the reduction of genus 2 small period matrices as implemented by Marco Streng and his collaborators in [`mstreng/recip`](https://github.com/mstreng/recip).
