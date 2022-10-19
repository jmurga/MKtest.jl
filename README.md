# ABC-MK

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://jmurga.github.io/MKtest.jl/dev)  

MKtest.jl is a Julia package including a fast Approximate Bayesian Computation version of the McDonald-Kreitman test (ABC-MK) presented in [Uricchio et al. (2019)](https://doi.org/10.1038/s41559-019-0890-6). The new ABC-MK implementation significantly improves the efficiency of the population genetics inferences. Following [Uricchio et al.(2019)](https://doi.org/10.1038/s41559-019-0890-6), the analytical estimations were used to explore the effect of background selection and selective interference on weakly beneficial alleles. Nonetheless, we developed a more straightforward and computationally efficient ABC-based inference procedure that accounts for the DFE of deleterious and beneficial alleles and partial recombination between selected genomic elements. Our approach estimates $\alpha$, $\alpha_W$, $\alpha_S$, and the Gamma distribution DFE parameters. 

In addition, the package automatizes other MK-like analyses parsing polymorphic and divergence data as well as including several extensions such as [Grapes](https://doi.org/10.1371/journal.pgen.1005774), [aMK](https://doi.org/10.1073/pnas.1220835110), [imputedMK](https://doi.org/10.1093/g3journal/jkac206) or [fwwMK](https://doi.org/10.1038/4151024a).


See the (documentation)[https://jmurga.github.io/MKtest.jl/dev] for details.