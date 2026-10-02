# HCV reference-residue analysis

Data and code for:

> Le DHH, Puenpa J, Nilyanimit P, Honsawek S, Poovorawan Y.
> **Hepatitis C virus resistance a decade after direct-acting antivirals: when the resistance residue is the wild type.**
> *(journal, year, doi — to be added on acceptance)*

Everything needed to reproduce every number and every data panel in that paper.

## What this is

We read the reference amino acid at fourteen resistance-associated positions (ten in NS5A, four in NS5B) in all 93 confirmed HCV genotype and subtype reference genomes, and scored each genotype against its own entry in EASL 2020 Table 7.

Headline results reproduced by this code:

- 112 genotype-by-position cells; 54 are listed by the guideline for that genotype and are therefore evaluable; in 25 of those 54 the residue the guideline lists as resistance is already the reference state.
- Serine 282 is the reference residue in 92 of 93 genomes. The single exception, subtype 4a (Y11604), is a record artefact and is discussed in the paper.
- Genotypes 7 and 8 are listed at no position, so for those two lineages the question cannot be asked.

## Layout

```
data/     input tables, one row per observation, each with its own source tag
R/        analysis and figure scripts
out/      rendered panels (PDF for publication, 600 dpi PNG for checking)
```

| File | What it holds |
|---|---|
| `data/hcv_93_subtype_references.fas` | the 93 reference genomes, unaligned |
| `data/accessions.csv` | the 93 subtype references and their accession numbers |
| `data/fig2_full_residue_matrix.csv` | 93 genomes x 14 positions, the core result |
| `data/parc_ras_easl2020.json` | EASL 2020 Table 7, machine-readable, the scoring standard |
| `data/easl2020_table7_ras_sets.csv` | the same residue sets flattened per panel |
| `data/Table_S2_genotype_position_matrix.csv` | the 112-cell scoring table, Supplementary Table S2 |
| `data/fig2a_residues.csv` | residue counts per position, species level |
| `data/fig2b_genotype_matrix.csv` | the genotype-resolved scoring |
| `data/fig3a_prevalence.csv` | country viraemic prevalence, 128 rows, each tagged with its source |
| `data/fig3a_trial_countries.csv` | the 29 countries that hosted a registration trial site, by trial |
| `data/fig3b_subtypes.csv` | the non-epidemic subtype table |
| `data/fig4c_retreatment.csv` | the retreatment-by-scenario table |

## Reproducing the figures

Requires R 4.3 or later.

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "maps", "ggrepel", "sf", "jsonlite"))
```

From the repository root:

```r
source("R/run_all.R")
```

or one panel at a time:

```r
source("R/fig2a_residue_distribution.R")   # Figure 3A
source("R/fig3a_map.R")                    # Figure 4A
```

## Reproducing the scoring

```r
source("R/score_matrix.R")
```

This reads `data/fig2_full_residue_matrix.csv` and `data/parc_ras_easl2020.json`, applies the scoring rules, rewrites `data/Table_S2_genotype_position_matrix.csv`, and prints the counts reported in the paper:

```
cells                     112
evaluable                 54
listed residue present    25
scored, none listed       29
not scored                58

Figure 3B: cells 40 | evaluable 25 | listed residue present 13
```

Scoring rules, as implemented: fourteen positions, ten in NS5A (24, 28, 30, 31, 32, 38, 58, 62, 92, 93) and four in NS5B (159, 282, 316, 321). Genotype 1 is scored against the union of the 1a and 1b entries; genotypes 2 to 8 against their own entry. A cell is evaluable when the guideline lists at least one substitution at that position for that genotype. NS5B 237 is excluded: it appears in some interpretation rule sets but in neither EASL 2020 Table 7 nor the AASLD-IDSA resistance primer, and no primary phenotypic source for it could be traced.

## Independent check on the extraction

The residues in `fig2_full_residue_matrix.csv` were read in BioEdit by transferring alignment columns from the subtype 6f reference, as described in Supplementary Methods S3. `code/validate_extraction.py` re-derives the same fourteen positions from the raw genomes by a route that shares no step with that analysis: it anchors the polyprotein frame on the NS5B GDD catalytic motif, aligns each polyprotein to the subtype 1a reference with a BLOSUM62 global pairwise alignment, and reads the positions through that alignment. No multiple-sequence alignment is involved.

```
pip install biopython
python code/validate_extraction.py
```

Result:

```
genomes anchored        93 / 93
cells compared          1302
agreement               1298 (99.69%)
disagreement            4
   subtype 4w    FJ025855.1   NS5B-159   primary=X  independent=L
   subtype 6q    EF424625.1   NS5B-159   primary=X  independent=L
   subtype 6u    EU246940.1   NS5A-92    primary=S  independent=A
   subtype 6u    EU246940.1   NS5A-93    primary=T  independent=S
```

Two of the four are ambiguity codes that the primary analysis left unresolved and the independent route called. The remaining two are both in subtype 6u and are a genuine disagreement between the two alignment routes at adjacent positions. Neither affects any count reported in the paper: tyrosine at 93 is unchanged at 50 of 93, and the genotype 6 residue set at position 93 remains I/S/T either way.

`sf` is needed only for the hatching on Figure 4A. Without it that panel still renders; the three seroprevalence-only countries keep their own fill and legend key but are not hatched.

## Method, in short

Coding sequences were aligned as translated nucleotides. Domain coordinates were transferred from the subtype 6f reference DQ835760 column-wise rather than by offset: the alignment column was located once in the anchor, then converted to each genome's own ungapped index before the codon was read. Each codon is therefore read in its own frame and the mapping is unaffected by indels elsewhere.

As a positive control the NS5B GDD catalytic motif was located in every genome. It falls at residue 317 in all 93 without exception.

Full method: Supplementary Methods S3 of the paper.

## Known issues in the source data

**Subtype 4a, Y11604, NS5B 282.** The designated 4a reference encodes threonine where every other 4a sequence examined, including the same isolate re-sequenced as GU814265, carries serine. It is an uncorrected error in a curated record, sitting at exactly the position a resistance analysis reads. We retain it because it is the designated reference and flag it rather than correcting it. Bulk sequence retrieval silently returns serine here; the threonine appears only on a targeted single-record fetch.

**Subtype labels.** An earlier version of `fig2_full_residue_matrix.csv` carried `4w3` and `6xa5`, which are `4w` and `6xa` with an ICTV footnote marker concatenated during extraction from the web table. These are corrected here. The sequences and accessions were always correct.

**ICTV version.** This analysis uses the confirmed list as it stood before the January 2026 revision, which added subtype 5b. Subtype 5b is not included.

## Licence

Data and documentation: CC BY 4.0. Code: MIT. See `LICENSE`.

## Citing

Please cite the paper. To cite this archive specifically, use the Zenodo DOI for the version you used.
