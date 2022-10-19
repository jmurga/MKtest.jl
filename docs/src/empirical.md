# Infering the rate and strength of adaptation

We developed an extension of the analytical approximation presented [Uricchio et al. (2019)](https://doi.org/10.1038/s41559-019-0890-6). In such a paper, the authors used the analytical approximation to $\alpha_x$ to explore the effect of BGS on weakly beneficial alleles, but the empirical estimation of $\alpha$ procedure employed was based on computationally intensive forward-in-time simulations. Unfortunately, forward-in-time simulations often require high-performance computing resources and can be prohibitively slow for some models. Here we extend [Uricchio et al. (2019)](https://doi.org/10.1038/s41559-019-0890-6)analytical approximations to develop a simple and computationally efficient ABC-based inference procedure modelling the DFE of deleterious and beneficial alleles and incomplete recombination between selected genomic elements. 

We followed a generic ABC scheme ([Beaumont et al. (2002)](https://doi.org/10.1093/genetics/162.4.2025)) to perform the empirical estimation of $\alpha$. ABC scheme proceeds by first sampling parameter values from prior distributions, simulating the model using these parameter values, calculating informative summary statistics, and comparing the simulated summary statistics to the observed data. Following the standard ABC scheme, we 1) considered $N$ random combinations from prior distributions; 2) solved $N$ independent models to generate informative summary statistics; 3) subset the parameter values producing summary statistics that best match the observed data from an approximate posterior distribution. Additionally, a linear model can be imposed to correct the non-0 distance between the simulated and observed summary statistics.

## Extending analytical estimations
[Uricchio et al. (2019)](https://doi.org/10.1038/s41559-019-0890-6) used the analytical theory to demonstrate the effect of weakly beneficial alleles at $\alpha_{(x)}$. To do that, they solved $\alpha_{(x)}$ through the fixation and polymorphism rates since the locus length ($L$) and the time branch ($T$) estimations at first order $\alpha_{(x)}$ estimation can be canceled

$\mathbb{E}[\alpha_{(x)}] \approx 1 - \frac{LT(\mathbb{E}[d_S])}{LT(\mathbb{E}[d_N])} \frac{LT(\mathbb{E}[p_{S(x)}])}{LT(\mathbb{E}[p_{N(x)}])} \approx 1 - \frac{\mathbb{E}[d_S]}{\mathbb{E}[d_N]} \frac{\mathbb{E}[p_{S(x)}]}{\mathbb{E}[p_{N(x)}]}$

where $\frac{LT}{LT}$ with a constant mutation rate tend to $1$. 

Here, we extended their analytical calculations to generate the summary statistics required at ABC approaches. Thereby we avoided expensive forward simulations. 

Considering fixed values of $T$, $L$, $\mu$, and fixation rates it is possible to solve analytical $\alpha_{(x)}$ performing the multiplication. However, this requires not only branch length estimations but also explicitly locus length selection which highly increases the order of estimations to solve. To avoid performing branch length estimations in our computation, we assumed that the empirically observed number of fixations should be proportional to the length of the evolutionary branch of interest, $T$, the locus length $L$ and mutation ration $\mu$. We take the observed total number of fixations (including both nonsynonymous and synonymous sites) as a proxy for the expected number, and then sample weakly deleterious, neutral, and beneficial substitutions proportional to their relative rates for a fixed set of model parameters. The expected number of substitutions for positively selected substitutions is then


$\mathbb{E}[D] = X \in Poisson\left(\lambda = D_{observed} \times (\mathbb{E}[d_+]+\mathbb{E}[d_-]+\mathbb{E}[d_0])\right)$

To estimate $\alpha_{(x)}$ we draw the expected values considering each fixation category. We sampled the values in two categories according to the relative rates. Following the same procedure:

$\mathbb{E}[D_S] = X \in Poisson\left(\lambda = D_{observed} \times \left[\frac{\mathbb{E}[d_0]}{\mathbb{E}[d_+] + \mathbb{E}[d_-] + \mathbb{E}[d_0]}\right]\right)$

$\mathbb{E}[D_N] = X \in Poisson\left(\lambda = D_{observed} \times \left[\frac{\mathbb{E}[d_+] + \mathbb{E}[d_-]}{\mathbb{E}[d_+] + \mathbb{E}[d_-] + \mathbb{E}[d_0]}\right]\right)$

It should be noted that both sampling variance and process variance affect the number of variable alleles at any particular allele frequency in a sequencing sample. 

To incorporate that variance we made one sample per frequency-bin given the SFS, we sampled a Poisson distributed number of polymorphic alleles at frequency $x$ relative to their rate given the expected frequency spectra. The expected frequency spectra were downsampled using a binomial (with probability of success given by the frequency $\begin{pmatrix} x \\ 2n \end{pmatrix}$ in a sample of $2n$ chromosomes) to account for the sampling variance. Considering the random process at each frequency we made credible intervals for future ABC estimations. Otherwise, we would falsely find higher confidence in our parameter estimates. We draw the expected polymorphic values similarly to fixations, considering the SFS and the expected rates.

$\mathbb{E}[P] = \sum_{x=0}^{x=1} X \in Poisson\left(\lambda = SFS_{(x)\ observed} \times (\mathbb{E}[p_{+(x)}]+\mathbb{E}[p_{-(x)}]+\mathbb{E}[p_{0(x)}])\right)$

It should be noted that we modified the expected SFS and the observed SFS considering each frequency $x$ as the cumulative sum above $x$. This quantities have the same asymptote but are less affected by changing sample size. Therefore, we scaled better the analysis to sample size since most common alleles at frequencies $x$ have very few polymorphic sites in large samples. This way, we finally transformed $\alpha_{(x)}$ to be depending on the previous frequency bin value, which reflects each frequency category over the expected asymptotic shape. Because of that, $\alpha_{(x)}$ analysis is more robust even in low polymorphic populations.

We followed the same procedure described above to subset the expected number of fixation taking into account weakly, strong and deleterious fixation rates categories. 

$\mathbb{E}[D_W] = X \in Poisson\left(\lambda = D_{observed} \times \left[\frac{\mathbb{E}[d_W]}{\mathbb{E}[d_+] + \mathbb{E}[d_-] + \mathbb{E}[d_0]}\right]\right)$

$\mathbb{E}[D_S] = X \in Poisson\left(\lambda = D_{observed} \times \left[\frac{\mathbb{E}[d_S]}{\mathbb{E}[d_+] + \mathbb{E}[d_-] + \mathbb{E}[d_0]}\right]\right)$

In addition, for each model, we retrieve deleterious DFE parameters to perform ABC inference.


# Computational pipeline
The following sections describe a pipeline to estimate the empirical adaptation rate ($\alpha$) using whole-genome protein, Viral Interacting Proteins and Non-Viral Interacting Protein data from the TGP. This analysis is similar to the one performed at [Uricchio et al. (2019)](https://doi.org/10.1038/s41559-019-0890-6). We divided the pipeline into the following step:

 - Expected fixations rates and frequency spectra estimations
 - Summary statistic estimations
 - ABC inference

Please note that we use a relative path called *analysis/* to execute the whole pipeline. This path will be relative to the where Julia or the Command-Line Interface is executed.

The software is prepared to parallelize each pipeline step using Julia using multi-threading computing. You can add the desired number of threads executing `-t` option when running julia (e.g using 8 cores: `julia -t8`). Please make some tests in your machine before executing expensive models. It is almost mandatory to parallelize at least the rate estimations.

The following examples were tested using a laptop with the following hardware:
- Intel i7-7700HQ (8) @ 3.800GHz 
- 16GB RAM DDR4 2400MHz
