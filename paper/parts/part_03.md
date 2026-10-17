# Часть 3: с. 6–7

[с. 6]
1352 J. Phys. Chem. B, Vol. 109, No. 4, 2005

Förster et al.

Figure 5. Peak shapes L(q) (eq 27) for different values of the peak
shape parameter of ν ) 0.001 (gray dotted), 0.1 (black dotted), 0.2
(gray dashed), 0.4 (black dashed), 1(gray solid) and 1000 (black solid)
showing the evolution of the peak shape from Lorentzian to Gaussian.

3.5. Lattice Disorder. For lattices with translational disorder
where the displacements from the ideal lattice points are
statistically independent (static disorder of the first kind), one
can assume a Gaussian lattice point distribution with zero mean
and a relative mean square displacement σ2a. In this case the
lattice factor is given by28

Z(q) ) Z0(q)G(q) + (1 - G(q))

(24)

where the first term, Z0(q)G(q), describes the decay of the
intensity of the Bragg peaks, whereas the second term, 1 G(q), describes the concomitant increase of the diffuse scattering. G(q) is given by

G(q) ) exp[-σ2aja2q2]

(25)

which is known as the Debye-Waller factor for thermal
disorder. Here, ja is the next nearest neighbor distance between
adjacent particles, which for different lattice types is given in
Table 3. In the case of lamellar or bilayer systems, the relative
displacement is related to the Caillé parameter η12,29 as

η ) 8σ2a

(26)

which characterizes thermal displacement fluctuations in lipid
bilayers. If the interlamellar forces are dominated by entropic
undulations, the Caillé parameter has a value of η ) 4/3(1 2R/a)2, where 2R is the bilayer thickness.30 σa in eq 26 is
comparable to the Hosemann g-factor, which describes the
relative standard deviation of the mean distance between
particles in paracrystal theory (static disorder of the second
kind).
3.6. Peak Shape. Lhkl(q) is a normalized peak shape function
for which Lorentzian, Gaussian, pseudo-Voigt, or Pearson
functions are commonly used.31 A general and versatile normalized peak shape function is adopted from a function given by
Burger and Micha32 with

Lhkl(q) )

2

∞

(

∏ 1+

πδn ) 0

γν2

4x2

(n + ν/2)2 π2δ2

)

-1

(27)

Figure 6. Schematic presentation of ideal lattices of monodisperse
homogeneous (A) and core/shell particles (B), as well as distorted
lattices of polydisperse homogeneous (C) and core/shell particles (D).
Cases (C) and (D) are typical for ordered mesoscopic structures and
are described by eq 30.

where x ) (q - qhkl) and

γν ) π1/2

Γ[(ν + 1)/2]
Γ[ν/2]

∞
The peak shape is normalized such that ∫-∞
Lhkl(q)dq ) 1.
Here, qhkl is the peak position and δ the peak width which can
be related to the domain size D of the ordered domain via the
Debye-Scherrer equation D ) 2π/δ. This peak shape function
is convenient since by adjusting the parameter ν the peak shape
can be continuously varied between a Lorentzian (ν f 0) and
a Gaussian (ν f ∞), i.e.,

Lhkl(x) )

{ []
δ/2π
x2 + (δ/2)2

for ν f 0 (Lorentzian)

2
4x2
exp - 2 for ν f ∞ (Gaussian)
πδ
πδ

(28)

where, in case of the Lorentzian, δ is equal to the full width at
half-maximum (fwhm), whereas in case of a Gaussian δ is
related to the standard deviation as σ ) xπ/8δ. Inserting eq 24
into eq 1 yields for the scattered intensity

I(q) ) (b1 - b2)2FNP(q)[1 + β(q)(Z0(q) - 1)G(q)] (29)
which can be rewritten as

I(q) ) (b1 - b2)2FNP(q)S(q)
S(q) ) 1 + β(q)(Z0(q) - 1)G(q)

(30)

where S(q)is the structure factor. This approach considers the
effect of the particles (via their first and second moment of the
particle size distribution) and of the lattice (via the first and
second moment of the distribution of lattice points) on the
scattered intensity I(q). The particle size and lattice point
distribution is schematically depicted in Figure 6. According
to eqs 17 and 25, both effects lead to damping of the higher

[с. 7]
Scattering Curves of Ordered Mesoscopic Materials

J. Phys. Chem. B, Vol. 109, No. 4, 2005 1353

order Bragg reflections by a factor e-(σa2ja2+σR2R2)q2, with additional suppression of Bragg reflections at the zeroes of β(q).
The additivity of both contributions in the argument of the
exponential decay results from the assumption of their statistical
independence.
It is straightforward to consider also other polymer or colloid
structures with their respective form factors P(q) and structure
factors Z0(q) which, e.g., are well compiled in a review by
Pedersen.33
3.7. Grazing Incidence Scattering (GISAXS). Eq 30 also
allows one to describe the 2-dimensional in-plane order of
mesoscopically structured thin films or monolayers. The inplane order can be derived from the diffuse scattering intensity
of grazing incidence small-angle X-ray (GISAXS) or neutron
scattering (GISANS) patterns.34 The total detected intensity is
given by the sum of the specularly reflected intensity, the
intensity scattered due to the surface roughness, and the diffuse
scattering, which is related to the density fluctuations in the
bulk of the thin film. The diffuse scattering is given by35

Id(q) )

k cC

|Ti| |Tf | Idiff (q)
(4πr2)
2

2

∫0 I(q)q dq ) 2π φ(1 - φ)
2

2

(32)

(33)

where q ) xq2xy + q2z , in agreement with ref 35.
Thus, the description of ordered domains using eq 30 is
applicable not only to the scattering of bulk mesoscopic
materials but also to the in-plane scattering of thin films or
monolayers. Such structures are of considerable relevance for
the structuring of surfaces by lithography or self-assembly for
catalysis, sensing, and electrooptical applications. It is straightforward to extend the description to multilayers by considering
the respective reflections of the structure factor S(qxy,qz). By
further consideration of the Fresnel transmission functions T(qz)-

(b1 - b2)2
Vpd2

(34)

which is the total scattering length of the sample. Porod’s
invariant is related to the asymptotic behavior of the scattered
intensity at high q, i.e.,

Q)

π
l (b - b2)2FN lim q4P(q)S(q)
qf∞
4 p 1

(35)

where lp is the Porod length or chord length which is related to
the surface per unit volume as

(31)

Thereby, the in-plane order of thin films as measured by grazing
incidence out-of-plane scans can be described for hexagonal or
square arrangements of cylindrical domains having a homogeneous or core/shell structure. Because of the particular orientation of the cylinders (the long axis being in z-direction, the
cross-section in the xy-plane) the factorization in eq 3 is valid
for any axial ratio of the cylinder, in particular also for disks,
where L , R. The in-plane form factor has an asymptote of
2
lim 〈F 2⊥
(qR)〉 ) 4/π(qR)-3, so that a characteristic q-3 scaling
qf∞
is expected for out-of-plane scans of domains with spherical
cross-sections.
In the case of spherical particles, the form factor contribution
to the diffuse scattering is isotropic and given by
2
(qR)〉
P3(qxy,qz) ) 〈F 3⊥

Q)

∞

lp )

where kc is the vertical component of the momentum vector at
the critical angle of total external reflection, C the illuminated
area, T the Fresnel transmission coefficients of the incident and
reflected waves, and Idiff(q) is the “diffuse structure factor”. With
the film positioned in the xy-plane, the diffuse scattering has
components along qxy (“out of plane of incidence”) and qz
normal to the xy-plane. The scattering amplitude of cylindrical
particles with their long axis normal to the xy-plane is then given
by the cross-sectional contribution F2⊥(qxyR) and the longitudinal
part F2 | (qzL). The diffuse part of the scattering intensity is then
given by Idiff(q) ) IN(q) where the form factor P(q)in eq 30 is
given by
2
(qxy R)〉 〈F 22 | (qzL)〉
P2(qxy,qz) ) 〈F 2⊥

of the thin film, the complete GISAXS or GISANS pattern can
be calculated.
3.8. Integral Scattered Intensity. For two-phase systems,
the factor c in eq 23 ensures that the integral scattered intensity
as characterized by the Porod invariant Q is given by

4φ(1 - φ)
SV

(36)

This relation is rigorous for any microphase-separated system
as long as the interface is sharp. For some materials, e.g.,
polymeric systems, interfaces can be relatively diffuse, so that
this relation may not be applicable. As outlined above, core/
shell particles do not in general have a simple Porod behavior,
so that the following arguments apply only for homogeneous
two-phase systems.
In the case of homogeneous two-phase systems with sharp
interfaces with a given dimensionality d, the surface per unit
volume is given by

SV )

dφ
R

(37)

from which the specific surface area Asp (in m2/g) for a
mesoporous material with wall- or matrix density Fsp can be
calculated as

Asp )

SV
(1 - φ)Fsp

(38)

For q f ∞, the Debye-Waller factor approaches G(q) f 0,
so that S(q) f 1, and the corresponding limit for the scattered
intensity is given by lim I(q) ) Pd∞q-4, with Pd∞ given in
qf∞
Table 1. Combining eqs 34-37 and inserting for the particle
density FN ) n/Vd, i.e., the number of particles per unit cell,
we obtain a closure relation

πRn
(1 - φ)(b1 - b2)2Pd∞ )
dVd

∫0∞ IN(q)q2 dq

(39)

which determines the coefficient c of the lattice factor Z0(q) in
eq 23. In the examples considered below, deviations from c )
1 are only up to a few percent.
3.9. Limiting Behavior for q f 0. In the thermodynamic
limit q f 0, the scattering intensity is related to the fluctuation
in particle number density δFN or the isothermal osmotic
compressibility χT as

lim Ĩ(q) ) 〈δFN2〉 ) kTχT FN
qf0

(40)

where Ĩ(q) ) I(q) (F1 - F2)-2 F-1
N and kT is the thermal energy.
In the limit q f 0 the form factor P(q) (eq 2), the fluctuation

