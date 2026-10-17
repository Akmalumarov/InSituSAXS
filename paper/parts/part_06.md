# Часть 6: с. 12–13

[с. 12]
1358 J. Phys. Chem. B, Vol. 109, No. 4, 2005

Förster et al.
part of the experimentally determined GISAXS pattern.34 Figure
14b shows the calculated diffuse scattering pattern for the same
cylinders with narrow disperse cross-sections and lengths, i.e.,
σR ) 0.05 and σL ) 0.05 to illustrate the periodic oscillations
due to the cylinder form factor.
5. Conclusions

Figure 13. GISAXS out-of-plane scans for a thin (A) and a thick (B)
PS-PI block copolymer film (data taken from ref 34). Solid lines are
fits using eq 30 with the 2D form factor of eq 32.

of the lattice site and particle size distribution of such materials.
Another extension of the model, which is being developed, is
the consideration of orientation order that can be modeled by
appropriate orientational distribution functions.38
4.7. Mesoscopically Structured Thin Films. Figure 13
shows out-of-plane scans of GISAXS patterns measured for a
thin (a) and thick (b) poly(styrene-b-isoprene) (PS-PI) block
copolymer film, which had been spin-coated on top of an oxidecovered Si substrate.34 These scans probe lateral structures in
the film plane. The strong peak in Figure 13a indicates a regular
two-dimensional distribution of domains within the film. The
fit yields an interdomain distance of 168 nm with domain radii
of 65.6. The GISAXS out-of-plane scan of the corresponding
thick PS-PI film (Figure 13b) shows no clear peak, indicating
poor order, with the domains having radii of 48.3 nm. The solid
lines in Figure 13 are fits to eq 30 using the 2D form factor (eq
32), which shows the expected q-3 scaling.
Equations 30 and 32 also allow to calculate the complete twodimensional diffuse scattering pattern of the thin film assuming
cylinders of a length of L ) 100 nm and polydispersity σL )
0.3 (Figure 14a). This pattern compares well with the diffuse

We have outlined a general and versatile method that yields
analytical expressions for the scattering curves of spheres,
cylinders, and lamellae arranged on different lattices. This allows
one to quantitatively fit experimentally measured SAXS and
SANS curves of lipid and block copolymer liquid crystalline
phases, as well as of core/shell nanoparticle superstructures,
ordered nanocomposites, and mesoporous materials. For spheres
arranged in cubic lattices, the description of the scattering curves
compares with the paracrystal model. For lamellae or bilayers,
there is good agreement with the Caillé model. Thus, this
method allows a fast and nearly complete analysis of scattering
curves of mesoscopically ordered materials. The outlined method
can be easily extended to describe the diffusive part of out-ofplane scans of GISAXS and GISANS experiments on mesoscopically ordered thin films.
Acknowledgment. We thank P. Müller-Buschbaum for
providing the original GISAXS data sets and N. Freiberger and
O. Glatter for helpful suggestions. Financial support of the SFB
508 “quantum materials” is gratefully acknowledged.
Appendix
Hypergeometric Functions. The hypergeometric functions
PFQ can be written in terms of series expansions

∞

PFQ[a1,... ,aP,b1,... ,bQ;x] )

∑

(a1)n(a2)n...(aP)n xn

n)0 (b1)n(b2)n...(bQ)n n!

(A1)

where (c)n are the Pochhammer factorials given by

(c)n )

Γ(c + n)
) c(c + 1)(c + 2)...(c + n - 1) (A2)
Γ(c)

In the derivation of the scattering functions, the hypergeometric

Figure 14. (a) Calculated two-dimensional diffuse scattering (GISAXS) for a thin PS-PI block copolymer film (Figure 13a) using the fit parameters
given in Table 5 and assuming a cylinder length of L ) 100 nm with σL ) 0.30. For comparison, the corresponding scattering pattern for narrow
disperse cylinders (σR ) 0.05, σR ) 0.05) is shown in (b) to illustrate the periodic oscillations due to the form factor of the cylinders.

[с. 13]
Scattering Curves of Ordered Mesoscopic Materials

J. Phys. Chem. B, Vol. 109, No. 4, 2005 1359

functions 0F1[b1;x] (eqs 4, 6), 1F2[a1,b1,b2;x] (eqs 5, 7, 20), 2F1[a1,a2,b1;x] (eqs 12, 14), and 3F2[a1,a2,a3,b1,b2;x] (eqs 13, 15,
22) had been considered. In cases where the hypergeometric
functions reduce to simpler functions, these have been compiled
in Table 1. The sums in eq A1 can be easily calculated via
the recurrence relation of the Pochhammer factorial (c)n+1 )
(c + n)(c)n. This enables an efficient implementation of these
functions in fitting algorithms.
In cases where the hypergeometric function reduces to simple
trigonometric functions (see Table 1), the expressions for
2
(q)〉 involve averages for the Schulz-Zimm
〈Fd⊥(q)〉 and 〈F d⊥
distribution which are given by17

〈x2n〉 )

Γ[z + 2n + 1] 2n
z+1 z+2
(2u)2n
u )
2 n 2 n
Γ[z + 1]

(

)(

)

(A3)

e.g., for solid particles as the integrand in eqs 12 and 14 and
the squared integrand in eqs 13 and 15, with R ) 0, as well as
for core/shell particles and in eq 19. For this function, a quickly
converging asymptotic expansion is

[

1F2 a1,b1,b2; -

with

〈sin (ax)x 〉 )
sin [(z + n + 1)arctan(au)]
Γ[z + n + 2]
(A4)
un
Γ[z + 1](z + n + 1)
(1 + a2u2)(z+n+1)/2
Γ[z + n + 1] n cos[(z + n + 1)arctan(au)]
u
〈cos(ax)xn〉 )
Γ[z + 1]
(1 + a2u2)(z+n+1)/2
(A5)

〈cos2(ax)xn〉 )

Γ[z + n + 3]
un
×
Γ[z + 1](z + n + 1)(z + n + 2) 2
cos[(z + n + 1)arctan(2au)]
1(A6)
(1 + 4a2u2)(z+n+1)/2

(

)

Γ[z + n + 1] un
×
Γ[z + 1] 2
cos[(z + n + 1)arctan(2au)]
1+
(A7)
(1 + 4a2u2)(z+n+1)/2

(

)

4

n
Γ[b1]Γ[b2]
x s1+1/2 max
π
1
-1/2
π
cicos x + s1 + + n x-n +
2
2
2
Γ[a1]
i)0
2 -a1 mmax
Γ[b1]Γ[b2]
x
x2 -m
dm (A11)
4
Γ[b1 - a1]Γ[b2 - a1] 4
m)0

n

〈sin 2(ax)xn〉 )

]
( ) ∑ [ ( )]
() ∑ ( )

x2

)

〈sin (ax)cos(ax)x 〉 )
Γ[z + n + 2]
un sin [(z + n + 1)arctan(2au)]
(A8)
Γ[z + 1](z + n + 1) 2
(1 + 4a2u2)(z+n+1)/2

dm )

(a1)m(1 + a1 - b1)m(1 + a1 - b2)m
m!

Here, mmax ) 1 and nmax ) 2 for sufficient accuracy. The first
two terms of dm are d0 ) 1 and d1 ) a1(1 + a1 - b1)(1 + a1
- b2). The other constants are given as

c0 ) 1
c1 )

{

1 1
12
1
+ s1 - s1 +
- 2 - - s 1 + s2
2 3
2
3

(

〈PFQ[a1,... ,aP,b1,... ,bQ; - x2]〉 )
z+1z+2
,
,b1,... ,bQ;(-2u)2 (A9)
P+2FQ a1,... ,aP,
2
2

[

]

which has been used to derive eqs 12, 13, 14, 15, and 22.
Analytic Continuation and Asymptotic Expansion. For
large q, the summation in eq A1 can break down due to rounding
errors. In this case the hypergeometric functions can be
calculated via analytic continuations or asymptotic expansions.
For 2F1[a1,a2,b1;x] in eqs 10 and 14, the following analytic
continuations can be readily applied

x
[
x - 1]
x
) (1 - x) F [a ,b - a ,b ;
(A10)
x - 1]

2F1[a1,a2,b1;x] ) (1 - x)

-a1

2F1 a1,b1 - a2,b2;

-a2
2 1

2

1

1

1

All scattering amplitudes 〈Fd⊥(q)〉 and
can be expressed in terms of the hypergeometric function 1F2[a1,b1,b2;x],
2
〈F d⊥
(q)〉

(

)}

(2c1)2
+
8
s1 3s2
1 3
12
13 1
1
+ s3
s1 +
- s1 +
- s1 + + 4 6 2
2
2
2
2
2
2

c2 )

[( ) ( ) ( ) (

)]

s1 ) a1 - b1 - b2
s2 ) a12 - b12 - b22

n

with u ) jx/(z + 1). For d ) 3, the expression for 〈Fd⊥(q)〉 has
2
previously been derived in ref 23 and for 〈F d⊥
(q)〉 in ref 39.
2n
For the average 〈x 〉, the relation in terms of Pochhammer
factorials as given in eq A3 leads to the general relation

)

s3 ) a13 - b13 - b23
The averages over the 〈xn〉 and 〈cos(ax)xn〉 are given in eqs
A3 and A5. The asymptotic behavior is discussed in ref 17. By
changing from the series expansion (A1) to the asymptotic
expansion (A11) at intermediate q, one calculates the function
over the entire q-range.
Peak Shape Function. The peak shape function L(q) (eq 27)
is related to the Gamma function Γ[z] as

Lhkl (q) )

2

∞

∏

πδ n ) 0

(

1+

γ2νy2

)

(n + ν/2)2

-1

2 || Γ[ν/2 + iγνy] ||
|
|
πδ |
Γ[ν/2]
|
(A12)
2

)

with y ) 4(q - qhkl)2/π2δ2. For ν < 0.001 this function is
practically equal to a Lorentzian, for ν > 1000 practically equal
to a Gaussian. For ν < 1, the product converges quickly and a
number of nmax ≈ 20 products is sufficient for high accuracy.
For ν > 1, convergence becomes slower and values of nmaxf103
are necessary, until for ν > 1000 the Gaussian function becomes
an excellent approximation.
References and Notes
(1) Förster, S.; Plantenberg, T. Angew. Chem., Int. Ed. 2002, 41, 688.
(2) Ying, J. Y.; Mehnert, C. P.; Wong, M. S. Angew. Chem., Int. Ed.
1999, 38, 56.
(3) Seddon, J. M. Biochim. Biophys. Acta 1990, 41, 525.

