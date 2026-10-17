# Часть 4: с. 8–9

[с. 8]
1354 J. Phys. Chem. B, Vol. 109, No. 4, 2005

Förster et al.

term β(q) (eq 16), and the Debye-Waller factor G(q) (eq 30)
all approch unity. Due to the q1-d singularity of the lattice factor
Z0(q) (eq 23), the scattering intensity Ĩ(0) will diverge for the
case of cylinders (d ) 2) and spheres (d ) 3), but not for
lamellae (d ) 1). This parallels the behavior of the structure
factor in paracrystal theory, where intensity upturns for qf0
have been reported for spheres (three-dimensional paracrystal)
but not for lamellae (one-dimensional paracrystal).9 In paracrystal theory, the upturn of a scattering curve at qf0 corresponds
to the contribution of the zero-order scattering related to the
finite size of the paracrystal. For the one-dimensional paracrystal, as q approaches zero, the intensity Z(q) continuously
decreases to the limiting value of σ2a and the intensity upturn is
not found. For the three-dimensional paracrystal there is an
2
intensity upturn to a value of σa,eff
, which goes to infinity for a
random orientation of paracrystals with spherically symmetric
distortions of the lattice points.9 In eq 23, the q1-d singularity
similarly occurs as a consequence of assuming a random
orientation of the lattice, so that the angular average of the lattice
factor in reciprocal space yields a factor of (2πq)-1 in the twodimensional, and a factor of (4πq2)-1 in the three-dimensional
case. In eq 23, the zero-order peak is not explicitely considered
because of the q1-d singularity and because its intensity depends
also on the scattering length difference between the lattice inside
and outside, which is a priori not known.
Due to the singularity and since structural features on length
scales d > a, such as packing defects, grain boundaries or
fluctuations decaying on larger length scales are only indirectly
considered via the domain size D, eq 30 is not expected to give
good agreement with experimentally determined scattering
curves in the range of scattering vectors q < 2π/a. However,
for q g 2π/a, this approach describes remarkably well experimentally measured high-resolution scattering curves of mesoscopically ordered materials as will be shown and discussed
below. The use of factorizations (eq 3) and quickly converging
series expansions (eqs 4-7) facilitate computing and allow
efficient implementation into fitting routines. It allows one to
h ,σR,D,ν} from
quickly extract the set of parameters {a,σa,R
measured scattering curves involving spheres, cylinders, and
bilayers or lamellae. In addition, the parameters {Rm,R,F} can
be derived for core/shell systems. Parameters such as the
cylinder length L or the disk radius RD are difficult to determine
because for nearly all mesoscopically ordered materials L,RD
. a. All curves are fitted using a Levenberg-Marquardt
algorithm which quickly converges and allows one to obtain
the relative error of the fit parameters. An executable trial
version of the program can be obtained from the author (http://
www.chemie.uni-hamburg.de/pc/sfoerster/software.html). The
analysis allows also to take into account the contribution of the
background scattering as

Ĩ(q) ) Isample - kIsolvent

(41)

where the constant k is the ratio of the transmission of sample
and solvent, assuming measuring times and primary beam
intensity to be equal. If the constant k is unknown, it can be
treated as an adjustable parameter. If no solvent or background
data are available, a flat q-indpendent background can be
subtracted. Geometrical resolution effects are not taken into
account and are out of the scope of the present paper. They
could be incorporated by convoluting the calculated scattering
curve with the beam profile as outlined ref 33.
We find that the fits to experimentally measured scattering
curves are well comparable to approaches based on paracrystal
theory in the case of cubic micellar structures, or based on the

Figure 7. Calculated scattering curves I(q) (a) and structure factor
S(q) (b) for a hexagonal lattice of cylinders with a ) 50 nm, L ) 500
nm, R ) 15 nm, σR ) 0.1, D ) 500 nm, and mean displacements of
σaja ) 1 (solid), 3 (dashed), 10 (dash-dotted), and 30 (dotted) nm with
Lorentzian line shapes showing the suppression of higher order peaks
with increasing displacement. The inset in (a) shows the characteristic
q-4 envelope in the high-q regime (Porod law).

Caillé approximations in case of lamellar structures, and allow
to nearly quantitatively describe high-resolution synchrotron
SAXS and SANS data of many different mesoscopically ordered
materials.
4. Results and Discussion
4.1. Structure and Form Factor. The influence of both
structure factor S(q) and form factor P(q) on the scattered
intensity I(q) can be observed in Figure 7, where the scattered
intensity I(q) for hexagonally packed cylinders for different
values of the Debye-Waller factor are calculated. The form
factor of the cylinders leads to a q-4 decay of the scattered
intensity (Porod law, see inset in Figure 7a), which considerably
weakens the scattered intensity in the high-q region. In the Porod
regime, the form factor exhibits damped oscillations, whose
positions, periodicity and amplitudes are characteristic for the
particle size, shape, and polydispersity. Generally, the relative
distances between minima increase in the order spheref
cylinderflamellae (see Table 1). Higher order oscillations are
damped due to particle polydispersity.
The structure factor S(q) (Figure 7b) exhibits peaks that
correspond to reflections characteristic for the symmetry and

[с. 9]
Scattering Curves of Ordered Mesoscopic Materials

Figure 8. Synchrotron SAXS curve of Au@SiO2 core/shell particles
that have assembled into core/shell superlayers. There is a set of eight
major oscillations due to scattering of the core layer with smaller
oscillations inbetween at lower q, due to core/shell interferences. The
remarkably small polydispersity of the particles and the respective
superlayer thickness (σR ) 0.03) leads to the large number of form
factor oscillations. The minimum q-value of the instrumental setup is
indicated.

space group of the lattice (see Tables 2, 3, and 4). Different
crystalline lattices have different relative positions of the higher
order reflections. The peak width of the reflections is reciprocally related to the domain size of the ordered domain. The
intensities of the higher-order reflections are damped due to
deviations from the ideal lattice positions (Debye-Waller factor)
and are suppressed at specific q-values due to minima of the
fluctuation term β(q) (see Figure 2). The structure factor
approaches S(q)f1 for q f ∞ and S(q) f ∞ for q f 0.
Most of the structural variables discussed above, e.g., unit
cell dimension and relative displacement, particle radius, and
relative standards deviation, have a specific effect on the
scattering curve, i.e., unit cell dimension a: position of the first
reflection; lattice-type: relative position of higher-order reflections; root-mean-square displacement σa: damping of higherorder reflections; particle radius R: first minimum of oscillations; particle shape (sphere, cylinder, disk): relative position
of higher order oscillations; particle polydispersity σR: damping
of higher order oscillations; ordered domain size D: peak width;
parameter ν: peak shape. This is favorable for the fitting of
scattering curves to the expression of I(q) in eq 30. Since each
fit parameter has a specific effect on the scattering curve, this
leads in most cases to low values of the off-diagonal elements
of the covariance matrix, which indicates reasonable statistical
independence.
Additional features in the scattering curve can arise in the
case of core/shell particles, where additional oscillations from
the particle shell and core/shell interference terms depending
on the total radius Rm, shell/core scattering length ratio F, and
the density profile scaling parameter R, can be observed.
Furthermore, core/shell structures can lead to a suppression of
the scattered intensity at intermediate q, when the contrast
density of core and shell have opposite signs. This will be
discussed below for SAXS curves of hexagonally packed
cylinders in a lyotropic liquid crystalline phase.
4.2. Core/Shell Nanoparticles. Figure 8 shows the high q
part of a synchrotron SAXS curve of Au@SiO2 core/shell
nanoparticles. The nanoparticles have a Au core with diameter

J. Phys. Chem. B, Vol. 109, No. 4, 2005 1355

Figure 9. Comparison of measured and fitted SAXS curves for a PSPI (a, ref 10; b, ref 36) and a lyotropic phase of PI-PEO-9 in
methylcyclohexane (c). The gray lines are best fits to paracrystal theory
(taken from 10 and 36). The black lines are fits to eq 30. The quality
of the description for both models is well comparable. The fit parameters
are summarized in Table 5.

15 nm and a thick SiO2 shell with a total particle diameter of
75 nm. Interestingly, the nanoparticles assemble into layers
giving rise to a scattering curve which can be nearly quantitatively described by a planar core/shell layer with a total diameter
of 2Rm ) 68.2 nm, a core diameter of 2Rc ) 13.7 nm, and an
electron density ratio of shell/core F ) 0.30 with a nearly
homogeneous shell (R ) 0.042). The large number of oscillations indicate that the polydispersity of the layer thickness is
very low, i.e., σR ) 0.030. The eight main oscillations are due
to the high-density core of the layer. The sharp oscillations
between the main oscillations at low q are due to core/shell
interferences. This shows that eqs 18 and 19 give a good
description of the form factor of core/shell structures. The
arrangement of the nanoparticles into layers occurs as a
consequence of the sample preparation (sol/gel reaction under
centrifugation) and is further supported by scanning electron
microscopy.22
4.3. bcc Structures: Comparison to Paracrystal Theory.
Figure 9a shows the SAXS curve of a PS-PI block copolymer
(Mn ) 21.9 × 104 g/mol, 15.4 wt % PI) where spherical PI
domains are arranged on a bcc lattice.10 The first three
reflections, (110), (200), and (211), are clearly observed,
followed by three form factor oscillations at higher q. The
oscillations can well be described yielding the radius and
polydispersity of the PI spheres. The positions, intensities,
widths, and shapes of the reflections are well captured by the
fit to eq 30 (black line) and are comparable or even better than
the description by paracrystal theory (gray line). The results of
paracrystal theory and fits to eq 30 considering particle radius,
standard deviation, mean deviation, and unit cell dimensions
are the same within experimental error and are summarized in
Table 5.
Figure 9b shows a SAXS curve of another PS-PI block
copolymer (Mn ) 17.6 × 104 g/mol, 15.0 wt % PI), where also
in this case spherical PI domains are arranged on a bcc lattice.36
The first three reflections are well observed, followed by three
form factor oscillations at higher q. The gray line is the fit using
paracrystal theory.36 The data can also be well described by
using eq 30 (black line), however, with a slightly poorer fit.
The parameters for the paracrystal fit36 and the parameters

