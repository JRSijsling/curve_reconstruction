Description
-----------

This repository contains Magma code for reconstructing hyperelliptic curves of genus up to 3 from their period matrices, both geometrically and arithmetically. With some extra work, plane curves and hyperelliptic curves of arbitrary genus should be feasible. However, the repository is not being actively developed in that direction yet.

Prerequisites
-------------

You need [`https://github.com/edgarcosta/endomorphisms`](https://github.com/edgarcosta/endomorphisms) to install this repository.

Installation
------------

A copy of the computer algebra system Magma is needed to run this code. After being cloned or downloaded, it can be made to run upon startup of Magma by adding the lines

`AttachSpec("[PATH]/spec");  `

to the user's `.magmarc` file, which can typically be found in the home directory). Here `[PATH]` is to be replaced by the directory of the cloned and downloaded repository, so that one could for example have

`AttachSpec("~/Programs/curve_reconstruction/spec");  `

Credits
-------

This implementation is based on the following works. When using this package, please be aware of the work that you are indirectly applying and please cite it.

For geometric reconstruction in genus 3 (also see [their own `SageMath` implementation](`https://github.com/christellevincent/genus3`):

Balakrishnan, Jennifer; Ionica, Sorina; Lauter, Kristin; Vincent, Christelle:  
*Constructing genus-3 hyperelliptic Jacobians with CM.* (English summary)  
LMS J. Comput. Math. 19 (2016), suppl. A, pp. 283-–300.

For arithmetic reconstruction in genus 2:

Guàrdia, Jordi:  
*Jacobian Nullwerte and algebraic equations.*  
Journal of Algebra 253 (2002) 112–132.

For the fast computation of theta constants needed when geometrically reconstructing in genus 2 or 3:

Labrande, Hugo; Thomé, Emmanuel:   
*Computing theta functions in quasi-linear time in genus 2 and above.*  
LMS J. Comput. Math. 19 (2016), suppl. A, pp. 163–-177.

Finally, this work uses the reduction of genus 2 small period matrices as implemented by Marco Streng and his collaborators in [mstreng/recip](https://github.com/mstreng/recip).
