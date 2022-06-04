# HopTB.jl

HopTB.jl is a tight-binding package written in julia. The package has the ability of dealing with non-orthogonal tight-binding models and aims at both first principle calculations of real materials and model calculations.

For real materials, HopTB.jl currently has interfaces with Wannier90, openMX and FHI-aims. Tight-binding systems are created with these density functional packages and HopTB.jl is a post-processing tool.

For model calculation, HopTB.jl has a similar API as pythtb to construct tight-binding models.

HopTB.jl provides infrastructure for analyzing response function and analyzing band structures. In addition, HopTB.jl contains out-of-box features including
 - Permittivity
 - Drude weight
 - Anomalous Hall effect
 - Spin Hall effect
 - Shift current conductivity
 - Second harmonic generation
 - Symmetrization of tight binding model
 - Intrinsic nonlinear Hall conductivity
 - Berry curvature dipole
 - Second order Drude weight
 - Fermi surface extraction

For more details, see [Documentation](https://hoptb.github.io/HopTB.jl/dev/).
