# ZipHiC
An R implementation of the ZipHiC as described in "ZipHiC: a novel Bayesian framework to identify enriched interactions and experimental biases in Hi-C data" (Osuntoki et al., 2021). The ZipHiC incorporated the Potts model to account for spatial dependency and to increase the number of components from the standard two components assumed by other methods to three components. The ZipHiC model sources of biases (genomic distances between interacting loci, GC-content, Transposable elements, and Accessibility score of interacting loci) as a regression model in order to understand their effects. Due to the computational intractability of the normalizing constant in the Potts model, ZipHiC used a likelihood-free approach (Approximate Bayesian Computation) (Beaumont et al., 2009) to overcome this problem. Finally, the ZipHiC detects significant interactions using calculated posterior means from the Metropolis-Within-Gibbs sampler.




Mark A. Beaumont, Jean-Marie Cornuet, Jean-Michel Marin, Christian P. Robert, Adaptive approximate Bayesian computation, Biometrika, Volume 96, Issue 4, December 2009, Pages 983–990, https://doi.org/10.1093/biomet/asp052
