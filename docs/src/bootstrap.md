We included a previously developed bootstrap pipeline to compare adaptation levels between case and non-case (control) genes ([Enard and Petrov 2020](https://doi.org/10.1098/rstb.2019.0575) and [Di et al, 2021](https://doi.org/10.7554/eLife.69026)])

To allow the comparision between case and control genes the bootstrap test uses a straightforward control set-building algorithm that adds control genes to a progressively growing control set, in such a way that the growing control set has the same range of values of confounding factors as the tested set of genes of interest. In addition it is possible to control for a minimum distance between case and control genes. 

For TGP we included the [counfounding dataset](https://raw.githubusercontent.com/jmurga/MKtest.jl/main/data/confounding_factors.txt) used in [Di et al., 2021](https://doi.org/10.7554/eLife.69026)], as well basic [gene coordinates information](https://raw.githubusercontent.com/jmurga/MKtest.jl/main/data/ensembl_gene_coords_v69.bed) need to performt the bootstrap. The counfounding dataset control for the following genomic factors:

 - Average overall expression in 53 GTEx v7 tissues [GTEx Consortium, 2020](https://www.gtexportal.org/home/). We used the log (in base 2) of TPM (Transcripts Per Million).
 - Expression (log base 2 of TPM) in GTEx lymphocytes. Expression in immune tissues may impact the rate of sweeps.
 - Expression (log base 2 of TPM) in GTEx testis. Expression in testis might also impact the rate of sweeps.
 - deCode recombination rates 50 kb and 500 kb: recombination is expected to have a strong impact on iHS and nSL values, with larger, easier to detect sweeps in low recombination regions but also more false positive sweeps signals. The average recombination rates in the gene-centered windows are calculated using the most recent deCode recombination map [Halldorsson et al., 2019](https://doi.org/10.1126/science.aau1043). We use both 50 kb and 500 kb window estimates to account for the effect of varying window sizes on the estimation of this confounding factor (same logic for other factors where we also use both 50 kb and 500 kb windows).
 - GC content is calculated as a percentage per window in 50 kb and 500 kb windows. It is obtained from the USCS Genome Browser.
 - The density of coding sequences in 50 kb and 500 kb windows centered on genes. The density is calculated as the proportion of coding bases respect to the whole length of the window. Coding sequences are Ensembl v99 coding sequences.
 - The density of mammalian phastCons conserved elements [Siepel et al., 2005](https://doi.org/10.1101/gr.3715005) (in 50 kb and 500 k windows), downloaded from the UCSC Genome Browser. We used a threshold considering 10% of the genome as conserved, as it is unlikely that more than 10% of the whole genome is constrained according to previous evidence [Siepel et al., 2005](https://doi.org/10.1101/gr.3715005). Given that each conserved segment had a score, we considered those segments above the 10% threshold as conserved.
 - The density of regulatory elements, as measured by the density of DNASE1 hypersensitive sites (in 50 kb and 500 kb windows) also from the UCSC Genome Browser.
 - The number of protein-protein interactions (PPIs) in the human protein interaction network [Luisi et al., 2015](https://doi.org/10.1093/gbe/evv055). We use the log (base 2) of the number of PPIs.
 - The gene genomic length, that is the distance between the most upstream and the most downstream transcription start sites.
 - The number of gene neighbors in a 50 kb window, and the same number in 500 kb window centered on the focal genes: it is the number of coding genes within 25 kb or within 250 kb.
 - The number of viruses that interact with a specific gene ([Enard and Petrov 2020](https://doi.org/10.1098/rstb.2019.0575)).
 - The proportion of immune genes. The matched control sets have the same proportion of immune genes as disease genes, immune genes being genes annotated with the Gene Ontology terms GO:0002376 (immune system process), GO:0006952 (defense response) and/or GO:0006955 (immune response) as of May 2020 ([Gene Ontology Consortium and Gene Ontology, 2021](https://doi.org/10.1093/nar/gkaa1113)).
 - The average non-synonymous polymorphism rate pN in African populations, and the the synonymous rate pS. We matched pN to build control sets of non-disease genes with the same average amount of strong purifying selection as disease genes. Also, pS can be a proxy for mutation rate and we can build control sets of non-disease genes with similar level of mutation rates.
 - McVicker’s B value which can be used to account for more recent selective constraint ([McVicker et al., 2009](https://doi.org/10.1371/journal.pgen.1000471)).

To run the bootstrap you must to declare a variable containing some basic information about data and bootstrap options using the function ```MKtest.bootstrap_parameters```. To check all the arguments you can access to the function documentation using ```@doc MKtest.bootstrap_parameters```

```julia
b_param = MKtest.bootstrap_parameters(data = "analysis/ensembl_list.txt", annotation = "analysis/ensembl_gene_coords_v69.bed", dist = 1, rep = 3, tol = 0.05, iter = 10, factors = "analysis/confounding_factors.txt", output = "analysis/test")
```
To run the bootstrap just use `MKtest.boostrap` function. It will output one file containing the case genes (`analysis/test_case.txt`) and another containing the control genes (`analysis/test_control.txt`)

```julia
MKtest.bootstrap(b_param)
```