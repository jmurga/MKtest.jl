# Input data

## The gene table

Every method starts from a tab-separated table with one row per gene and no header:

| column | content |
|---|---|
| 1 | gene ID |
| 2 | Pn, the number of nonsynonymous polymorphisms |
| 3 | their derived allele frequencies, comma separated: `,0.0008,0.0015,` |
| 4 | Ps, the number of synonymous polymorphisms |
| 5 | their derived allele frequencies |
| 6 | Dn, nonsynonymous substitutions |
| 7 | Ds, synonymous substitutions |
| 8 | Ln, nonsynonymous sites (optional) |
| 9 | Ls, synonymous sites (optional) |

Frequencies must sit on the grid `i / (sample_size * ploidy)`. Without columns 8 and 9
the site counts are set to 0: alpha is unaffected, but every omega (substitution rate
per site) is NaN and `dfe_ml` cannot be used.

## Bundled datasets (`data/`)

| file | content |
|---|---|
| `MKdata_tgp_may2023.txt` | Human, 1000 Genomes, 661 African individuals; 17,538 genes; 9 columns |
| `tgp.txt` | an older version of the same table, 7 columns (no site counts) |
| `dgn_ral.txt` | *Drosophila melanogaster*, DGN Raleigh: 160 inbred lines (`ploidy = 1`); 12,061 genes; 7 columns |
| `dgn_zi.txt` | *D. melanogaster*, DGN Zambia: 154 inbred lines (`ploidy = 1`); 13,449 genes; 7 columns |
| `ensembl_gene_coords_v109.bed` | Ensembl v109 gene coordinates (`chr start end gene strand`), for control matching |
| `confounding_factors.txt` | per-gene confounding factors, for control matching (see [Case and control sets](case_control.md)) |
| `DNAonly_VIPs.txt.gz` | a case gene set: 1,064 human virus-interacting proteins (VIPs; the DNA-virus set) |
| `DNAonly_nonVIPs.txt.gz` | 1,000 control sets of non-interacting genes matched to it, one set per line |
| `ensembl_list.txt`, `example_bootstrap.txt` | a small case list and 100 control sets, used by the tests |

## Parsing: `parse_sfs`

```julia
sfs, div, alpha_x = MKtest.parse_sfs("data/MKdata_tgp_may2023.txt"; sample_size = 661)
```

returns three values:

- `sfs`: the raw SFS, a `1321 x 3` matrix `[DAC Pn Ps]`
  (`[1.0 36460.0 20983.0; 2.0 20124.0 13555.0; ...]`);
- `div`: `[Dn, Ds, Ln, Ls]` = `[35903.0, 45088.0, 2.19154122e7, 7.8090133e6]`;
- `alpha_x`: `[DAC frequency alpha(x)]`, the MK alpha above each frequency.

Inbred lines are read with `ploidy = 1`:

```julia
zs, zd, _ = MKtest.parse_sfs("data/dgn_zi.txt"; sample_size = 154, ploidy = 1)
# zs is 153 x 3; zd = [297477.0, 249994.0, 0.0, 0.0] (the file has no site counts)
```

### Gene sets

`gene_list` selects genes: a `Vector` of IDs, or a file (plain or `.gz`).

| `gene_list` | `mode` | result |
|---|---|---|
| none | `:case` (default) | one SFS for the whole table |
| none | `:gene` | one SFS per gene, plus the total as the last element |
| one ID per line | `:case` | one SFS for the listed genes |
| one ID per line | `:gene` | one SFS per listed gene, plus their total |
| several IDs per line (`set_1 G1 G2 ...`) | detected automatically | one SFS per line |

```julia
case_sfs, case_div, _ = MKtest.parse_sfs("data/MKdata_tgp_may2023.txt"; sample_size = 661,
                                         gene_list = "data/DNAonly_VIPs.txt.gz")
# case_div = [2176.0, 3112.0, 1.5620831e6, 567018.4]

ctl_sfs, ctl_div, _ = MKtest.parse_sfs("data/MKdata_tgp_may2023.txt"; sample_size = 661,
                                       gene_list = "data/DNAonly_nonVIPs.txt.gz")
# ctl_sfs is a Vector of 1000 SFS, one per control set
```

IDs that are not in the table are skipped silently.

### Frequency cutoff

`cutoff = [lo, hi]` keeps the frequencies in that range (both ends included). The SFS
rows are trimmed, but column 1 keeps the true DAC. A trimmed SFS cannot be projected or
given to `dfe_ml`, which need every bin; MK tests take their own `cutoff` instead.

## Projection and folding

```julia
s20 = MKtest.project(sfs, 20)   # 20 individuals: a 39 x 3 SFS, [1.0 14838.0 12561.0; ...]
MKtest.fold(s20)                # 20 x 3, the minor allele spectrum
```

`project` resamples the SFS to a smaller sample (hypergeometric, as dadi and moments
do); counts are truncated to integers. It needs the full, untrimmed SFS. Both accept a
`Vector` of SFS.

Other helpers: `cumulative_sfs` (tail sums, the input of the alpha(x) curves),
`reduce_sfs(sfs, bins)` (the counts summed within `bins` equal-width frequency bins, each
row at its upper edge `k / bins`) and `bootstrap_data` (Poisson resamples of one SFS).
