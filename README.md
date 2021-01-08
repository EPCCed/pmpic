## PMPIC code repository ##

Example usage:

mpiexec -n 2 monc --config=config.mcf

Dependencies:

* MPI

Compilers tested:

* GNU (on laptop and ARCHER)
* Cray (on ARCHER)

Notes on compilation:

* See wiki page on this site

This code simulates the development of a plume using a 
Parcel-in-Cell approach. A detailed technical description is given 
by Dritschel et al 2018 (QJRMS,
https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/qj.3319).

The code is currently constructed so that it uses as simple a set of 
equations as possible. The default test case is a spherical moist 
and warm thermal in a neutrally-stable boundary layer overlaid with 
a stably-stratified atmosphere. The code is currently 
non-dimensional, and uses idealised environmental profiles for 
stratification and environmental humidity. It is also a tool for 
demonstrating the Parcel in Cell method in a meteorological context, 
and shows how intricate representations of cloud properties can be 
efficiently produced using the method.

Mass, specific humidity and liquid-water potential temperature are 
conserved following the motion of each parcel. Parcel vorticity is 
evolved in the dynamical equations. MPIC has the advantages that it 
has an explicit sub-grid representation and numerical dissipation is 
minimal and can be controlled. The number of tuning parameters is 
kept to a minimum.

Effects of latent heating are included by increasing the effective 
buoyancy whenever the parcel specific humudity exceeds the 
height-dependent saturation profile.

Mixing is included in the form of parcel splitting and merging.

The code is written in Fortran, and has been parallelised within the 
eCSE project "A fully Lagrangian dynamical core for the Met Office 
NERC Cloud Model" 

