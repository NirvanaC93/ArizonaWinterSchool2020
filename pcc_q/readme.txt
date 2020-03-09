////////////////////////////////////////////////// 
// This code is part of the pcc_q MAGMA library //
//                                              //
// copyright (c) 2019 Jan Tuitman               //
//////////////////////////////////////////////////

This is the pcc_q (Point Counting on Curves) MAGMA library v2.18 for computing zeta functions of curves over finite (non-prime) fields.

The code is based on the algorithm presented in my papers:

Counting points on curves using a map to P^1, II, Finite Fields and their Applications 45 (2017), 301--322.
Counting points on curves using a map to P^1, Mathematics of Computation 85 (2016), no. 298, 961--981.

Jan Tuitman
jan.tuitman@wis.kuleuven.be

December 2019.

GETTING STARTED
---------------
To use the code:

1) Put all the files into one and the same directory.

2) Open the file example_q.m and uncomment:

  a) a defining polynomial Q
  b) a line directly above it specifying the cardinality q of the field. 

  or of course add a new line of your own (Q has to be monic in y)!

3) Load "pcc_q.m" and "example_q.m" into MAGMA (in this order).

4) Call chi:=num_zeta(Q,q);

The code will output the numerator chi of the zeta function of the smooth projective model of the curve defined by Q over the finite field with q elements. 

OPTIONAL PARAMETERS
-------------------
N        : p-adic precision for the Frobenius matrix (by default the code will compute a provably correct precision)

verbose  : prints the time and memory spent on the various steps to the screen

exactcoho: carry out the linear algebra for computing a basis for the cohomology to exact precision, which is a lot slower

An example using all of these options: 

chi:=num_zeta(Q,p,n: N:=20,verbose:=true,exactcoho:=true); 

