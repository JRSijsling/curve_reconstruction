Description
-----------

This repository contains Magma code for reconstructing hyperelliptic curves of genus up to 3 from their period matrices, both geometrically and arithematically. With some extra work, arbitrary genus should be feasible, though the author is not yet making moves in this direction.

Installation
------------

A copy of the computer algebra system Magma is needed to run this code. After being cloned or downloaded, it can be made to run upon startup of Magma by adding the lines

AttachSpec("[PATH]/spec");  

to the user's .magmarc file, which can typically be found in the home directory). Here [PATH] is to be replaced by the directory of the cloned and downloaded repository, so that one could for example have

AttachSpec("~/Programs/hyperrec/spec");  

Credits
-------

This implementation is based on the following works. Please cite these works, not merely this code, if you have used the code in your work.

Balakrishnan, Jennifer; Ionica, Sorina; Lauter, Kristin; Vincent, Christelle
*Constructing genus-3 hyperelliptic Jacobians with CM.* (English summary)
LMS J. Comput. Math. 19 (2016), suppl. A, pp. 283-–300.

Guàrdia, Jordi
*Jacobian Nullwerte and algebraic equations*
Journal of Algebra 253 (2002) 112–132

Labrande, Hugo; Thomé, Emmanuel
*Computing theta functions in quasi-linear time in genus 2 and above*
LMS J. Comput. Math. 19 (2016), suppl. A, pp. 163–-177.

