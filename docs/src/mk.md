# MK approaches
We included other heuristic MK approaches in our module. All the functions use the formated SFS and divergence data described at [data section](data.md).

## Standard MKT
The standard McDonald & Kreitman test ([McDonald and Kreitman, 1991]) was developed to be applied to protein coding sequences, combining both divergence ($D$) and polymorphism ($P$) sites, and categorizing mutations as synonymous ($P_S$, $D_S$) and non-synonymous ($P_N$, $D_N$). 

If all mutations are either strongly deleterious or neutral, then $Di/D0$ is expected to roughly equal $Pi/P0$. In contrast, if positive selection is operating in the region, adaptive mutations rapidly reach fixation and contribute more to divergence than polymorphism compared to neutral mutations, and then $Di/D0 > Pi/P0$. Assuming that adaptive mutations contribute little to polymorphism but substantially to divergence, the proportion of non-synonymous substitutions is inferred following Smith and Eyre-Walker (2002).

$\alpha = 1 - (\frac{P_N}{P_S}\cdot\frac{D_S}{D_N})$

Please check [`MKtest.standardMK`](@ref) to obtain more info.

```julia
adap = MKtest.parameters(n=661,cutoff=[0.0,1.0])
alpha, sfs, divergence, m= MKtest.parse_sfs(adap, data = "analysis/tgp.txt")

mk = MKtest.standardMK(sfs[1],divergence[1])
```

## Fay, Waycoff and Wu MK extension
[Fay et al. (2002)]([fwwMK](https://doi.org/10.1038/4151024a)) proposed an approach that removes all polymorphisms segregating at a frequency ($j$) below a given threshold (usually $j > 5\%–15\%$). Although there is no consensus about what this threshold should be used, [J. Charlesworth & Eyre-Walker (2008)](https://doi.org/10.1093/molbev/msn005) demonstrated that  estimates are robust using a frequency threshold of 15%, below which most slightly deleterious polymorphisms are found and removed. The estimates are reasonably accurate only when the rate of adaptive evolution is high and the Distribution of Fitness Effects (DFE) of deleterious mutations is leptokurtic.

$\alpha = 1 - (\frac{P_{N (j>15\%)}}{P_{S (j>15\%)}\cdot\frac{D_S}{D_N})$

Please check [`MKtest.fwwMK`](@ref) to obtain more info.

```julia
fww = MKtest.fwwMK(adap,sfs[1],divergence[1],cutoff=0.15)
```

## imputed MKT (in preparation)
The imputed MKT (impMKT) is a modification of the Fay, Waycoff, and Wu MK extension (fwwMK) ([Fay et al. (2002)]([fwwMK](https://doi.org/10.1038/4151024a)) to improve gene-by-gene analyses. The method propose the imputation of slightly deleterious mutations at the SFS rather than removing all variants below a frequency threshold. The imputedMK aims to maximize the information to test the excess of divergence ratio relative to polymorphism at the gene level.

$\alpha$ is estimated as

$\alpha_{imputed} = 1 - \left( \frac{P_{N} - P_{wd}}{P_{S}} \cdot \frac{D_{N}}{D_{S}} \right)$

where $P_{wd}$ is

$P_{wd} \approx P_{wd (j < 15\%)} = P_{N (j<15\%)} - \frac{P_{N (j>15\%) } \cdot P_{S (j<15\%)}}{P_{S (j>15\%)}}$

Please check [`MKtest.imputedMK`](@ref) to obtain more info.

```julia
imputed = MKtest.imputedMK(adap,sfs[1],divergence[1],cutoff=0.15)
```

## Asymptotic MKT
Proposed by Messer and Petrov (2013). This extension is robust to the presence of selective sweeps (genetic hitchhiking) and the segregation of slightly deleterious polymorphisms substitutions (BGS). In this approach, the authors defined $\alpha$ as a function that depends on the SFS of alleles. Therefore, $\alpha$ is estimated in different frequency intervals ($x$), and these values are then adjusted to an exponential function. An exponential fit is suitable as the non-synonymous allele frequency is expected to decay exponentially over the respective levels of synonymous polymorphisms (Messer & Petrov, 2013).

$\alpha$ is estimated as

$\alpha_{fit(x)} = a+b \cdot e^{-cx}$

Please check [`MKtest.aMK`](@ref) to obtain more info.

```julia
amk, ci, model = MKtest.aMK(adap,alpha[1])
```

## Grapes
Grapes is a software ML models of the DFE that assume PRF framework can estimate the expected proportion of adaptive fixations given the inferred DFE from the MKT data. In such approaches, the expected levels of fixations and polymorphism are used to perform likelihood estimates, while considering different evolutionary models. ML estimations of the DFE have varied from the first models inferring constant selection parameters across all loci to including models with continuous distributions of both positive and negative selection coefficients \citep{bierne_genomic_2004, eyre-walker_distribution_2006, boyko_assessing_2008, eyre-walker_estimating_2009, galtier_adaptive_2016, racimo_approximation_2014, galtier_adaptive_2016, tataru_inference_2017, zhen_greater_2021}, correcting the aforementioned assumptions to calculate how many non-adaptive substitutions are expected to become fixed given the empirical DFE. 

Notwithstanding, the newest and more sophisticated ML implementations \citep{galtier_adaptive_2016, tataru_inference_2017} take advantage of that proposed by \cite{eyre-walker_distribution_2006} and \cite{eyre-walker_estimating_2009}, which assumes that non-neutral deleterious mutations arise from a DFE in the form of a Gamma distribution. Nonetheless, unlike \cite{eyre-walker_distribution_2006}, these methods also model the effect of weakly advantageous alleles through an exponentially distributed function, where DFE is a mixture distribution between the gamma and exponential distributions. Therefore, the state-of-the-art methods usually follow a standard population genetic model based on PRF presented in \cite{galtier_adaptive_2016} and \cite{tataru_inference_2017} to later perform ML of the DFE. 

To estimate the expected counts of polymorphism and divergence, the model considers a Wright-Fisher panmictic population of size $N_e$, which diverged in a time $t$ where mutation occurs at a mutation rate $\mu$, per site per generation \citep{galtier_adaptive_2016}. Figure \ref{fig:dfe_prf} illustrate the equations defining the expected polymorphic and divergence counts. Let consider synonymous mutation as neutral ($P_S$) and non-synonymous mutation as selected ($P_N$), from PRF the expected counts given a frequency $i$ is estimated as

\begin{equation}\label{prf_ps}
	\mathbf{E}[P_{S[i]}] = \frac{4N_e \mu L_S}{i}
\end{equation}

\begin{equation}\label{prf_pn}
	\mathbf{E}[P_{N[i]}] = 4N_e \mu L_N \int_{0}^{1} B(i,n,x) H(s,x)dx
\end{equation}

\noindent where $L_N$ and $L_S$ are the total number of synonymous and nonsynonymous sampled alleles, $B(i,n,x)$ is the probability of observing a mutation at frequency $i$ in $n$ sequences when the true allele frequency is $x$

\begin{equation}\label{prf_binom}
	B(i,n,x) = \begin{pmatrix} n \\ i \end{pmatrix} x^{i} (1-x)^{n-i}
\end{equation}

\noindent and $H(s,x)$ is the time that a new semi dominant mutation with selection coefficient $s$ (in the heterozygous) spends between the frequency $x$ and the frequency $x+dx$ from diffusion theory \citep{wright_distribution_1938}

\begin{equation}\label{prf_stationary}
	H(s,x) = \frac{1-e^{-s(1-x)}}{x(1-x)(1-e^{-s})}
\end{equation}

Note that to obtain the expected polymorphic count given the underlying DFE of new mutations equation \eqref{prf_pn} should be integrated over the full DFE. Following \citep{eyre-walker_estimating_2009}, the underlying DFE for new deleterious mutations is defined by expression $\phi$

\begin{equation}\label{prf_gamma}
	\phi(s;a,b) = a^b s^{b-1} \frac{e^{-a s}}{\Gamma(b)}
\end{equation}

\noindent where $a$ and $b$ are scale and shape parameters from the Gamma distribution. Therefore, the expected polymorphic count given a particular DFE is defined as

\begin{equation}\label{prf_total_pn}
	\mathbf{E}[P_{N[i]}] = 2 N_e \mu L_N \int_{-\infty}^{\infty}\int_{0}^{1}B(i,n,x)H(s,x)\phi(s;\alpha,\beta) dx ds
\end{equation}

In addition, the method proposed by \cite{galtier_adaptive_2016} makes the DFE more flexible by considering that the most appropriate way to model the full DFE may not be the classical Gamma distribution over negative alleles. Therefore, Galtier's modeling included two different versions of the Fisher's geometric model, and a model assuming a Beta-shaped distribution of weak effect mutations, instead of a Gamma distribution. However, as explored further in \cite{galtier_adaptive_2016}, \cite{galtier_how_2020} and \cite{rousselle_is_2020}, the shifted negative modeling and the two DFE models based on the Fisher geometric model generally do not perform well in the analysis.

Following \cite{galtier_adaptive_2016} and \cite{tataru_inference_2017} procedure, the expected number of synonymous ($D_S$) and non-synonymous ($D_N$) substitutions can be estimated as

\begin{equation}\label{prf_ds}
	D_S = L_S \mu t
\end{equation}

\begin{equation}\label{prf_dn}
	\mathbf{E}[D_N] = 2 N_e \mu L_N \int_{-\infty}^{\infty} \frac{2s}{1-e^{(-4 N_e s)}} \phi(s) ds
\end{equation}

The proportion of adaptive mutations ($\alpha$) is estimated using the ML inference of DFE parameters and the expected counts of non-adaptive mutations. $\alpha$ therefore is estimated subtracting the non-adaptive substitutions and neutral substitutions from the total observed divergence counts at selected sites.

\begin{equation}\label{prf_adaptive_1}
	\alpha = (d_N - d_{N}^{na}) / d_N)
\end{equation}

\noindent note that $d_{N}^{na}$ is defined following \cite{galtier_adaptive_2016}. The equation decomposition is similar to the equations shown at \cite{tataru_inference_2017} and \cite{eyre-walker_estimating_2009}

\begin{equation}\label{prf_adaptive_2}
	d_{N}^{na} = \frac{2L_N N_e \mu \int_{-\infty}^{s_{adv}} \frac{2s}{1-e^{(-4 Ne s)}} \phi(s) ds}{L_N}
\end{equation}

To estimate the non-adaptive substitutions, the approach proposed by \cite{eyre-walker_estimating_2009} and \cite{tataru_inference_2017} integrate the DFE from $-\infty \to 0$ ($s_{adv} = 0$), taking into account the deleterious fixations and the negative nearly neutral fixations given the DFE. \cite{galtier_adaptive_2016} considers positive $s_{adv}$ values to subtract nearly neutral positive fixations given the $s_{adv}$ threshold too.

\begin{figure}[h]
	\includegraphics[width=1\textwidth]{\dir/introduction/figures/dfe_prf.pdf}
	\caption[Expected polymorphic and divergence sites based on PRF]{Expected polymorphic and divergence sites based on PRF approach. Note that for $s=0$ the equations became similar to equations \eqref{prf_ps} and \eqref{prf_ds}. Adapted from \cite{casillas_molecular_2017}.}
	\label{fig:dfe_prf}
\end{figure}

It is important to emphasize that the shape of the SFS may be biased and differ from the previously expected values. The distortion may mainly be due to the effect of demographic events. To account for such distortions, the different methods usually incorporate nuisance parameters following \citep{eyre-walker_distribution_2006}. The nuisance parameter $r$ modifies each frequency of the SFS individually. Thus, $r$ modifies the effective mutation rate at each frequency $i$, considering the relative mutation rate at $i$ with respect to the mutation rate at singletons \cite{eyre-walker_distribution_2006}. The perturbing parameter $r_i$ considers the same amount of distortion between nonsynonymous and synonymous SFS. Although the assumption is unrealistic, simulations showed that the original correction proposed by \cite{eyre-walker_distribution_2006} is robust enough to take into account demographic effects. However, other models correct for demographic effects by considering explicit changes in the effective population during the modeling process \citep{eyre-walker_estimating_2009,zhen_greater_2021}.

Despite the efforts for modeling the DFE using polymorphic data, \cite{booker_inferring_2020} suggested that inferred parameters from ML approaches must be reviewed. Underlying assumptions regarding the DFE shape, selection coefficient strengths, or population sizes can affect the estimations and capture distinct aspects of the DFE \citep{booker_inferring_2020,zhen_greater_2021}. Such a situation is plausible in several \textit{D. melanogaster} studies. For example, \cite{campos_estimating_2017} and \cite{keightley_inferring_2016} found differences between selection coefficients and probabilities of beneficial alleles. As extensively discussed in \cite{booker_inferring_2020}, it is plausible due to different methodological assumptions or if the DFE for advantageous mutations is bimodal. This latter case should be especially considered for strongly beneficial alleles, which would be undetectable by analyzing the unfolded SFS (uSFS) \citep{booker_inferring_2020}.

\cite{messer_population_2013} compared the performance of their method -the aMKT- and the \texttt{DFE-alpha} method \citep{eyre-walker_estimating_2009} through simulations. They claimed that DFE-alpha correctly estimated  when the model allowed for population size change, but the demography inferred was found to be biased, mainly due to background selection acting at linked sites. Genetic draft leaves signatures in the SFS similar to those observed under a recent population size expansion. The DFE-alpha method inferred systematically a population expansion even though no expansion was set in the simulation . Another limitation of PRF approaches is that it becomes computationally intensive, especially when a change in demographic model is applied. Since DFE-alpha can only consider two population-size changes, it becomes insufficient for capturing the excess of rare variants due to the complex demo