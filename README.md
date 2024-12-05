# Spatially-resolved molecular sex differences at single cell resolution in the adult human hypothalamus

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14285059.svg)](https://doi.org/10.5281/zenodo.14285059)

### Study Design
The hypothalamus is a critical brain area underlying functions with inherent sex differences, such as reproductive physiology, endocrine signaling, and metabolism. In the rodent, these sex-differentiated functions correspond to differences in volume, cell type composition, and gene expression between males and females across individual hypothalamic regions (here, "domains"). The ventromedial hypothalamus (VMH) and arcuate nucleus (ARC) are two hypothalamic regions that influence appetitive/social behaviors and growth/metabolism, respectively. While molecular profiling studies in the rodent hypothalamus have identified specialized cell types with unique transcriptomic signatures, there is a paucity of data describing the molecular architecture of the human HYP, especially in the context of sex-differentiated cell types that drive evolutionarily essential, hypothalamus-mediated behaviors in males and females. 

This study, led by Bernard (Bernie) Mulvey, Kristen Maynard, and Kasper Hansen, profiled the adult human mediobasal hypothalamus (ventromedial nucleus, VMH, and arcuate nucleus, ARC) using  Visium and Xenium spatial transcriptomic platforms (10x Genomics). Atlasing and sex-differential expression efforts using Visium data from 8 donors (4 per sex) were used to guide gene selection for subsequent Xenium assays on adjacent tissue sections from these same 8 donors.

![Experimental Overview](./images/overview.png)

Examining 13 Xenium samples revealed 5 ARC and 4 VMH neuronal populations governing known hypothalamus-specific functions and defined their spatial distributions. Examinining the full dataset of 23 samples, human VMH and ARC demonstrated increases in retinoid pathway gene expression relative to mouse. Sex-DE analyses of VMH and ARC revealed correlated autosomal expression differences, which were localized to *ESR1*- and *TAC3*-expressing neurons in the ARC, and *CRHR2*-expressing neurons in the VMH. VMH- and ARC-residing cell types have a striking number of sex-DE genes linked to sex-biased disorders, including autism, depression, and schizophrenia. By mapping disease associations to hypothalamic regions containing cell types with established roles in mediating sex-divergent physiology and behavior, these data provide insights into mechanistic bases of sex bias in neurodevelopmental and neuropsychiatric disorders.

![Sex DE Analysis with Xenium](./images/Xenium_sex_DE_analysis_schematic.png)
*Sex-DE analysis was performed at the level of Xenium clusters by first subsetting to all cells within the boundaries of the VMH or ARC. Then, cells of each type are tested separately for sex-DE within that domain. Thus, broadly distributed cell clusters/types (e.g., glia) are tested for sex-DE in each domain. Meanwhile, sex-DE testing of domain-specific neuronal clusters (e.g.,* TAC3*-*ESR1 *in ARC) is filtered to the cells given that label and found in their domain (i.e. in the expected anatomic space). This simplified schematic only shows one domain-specific cluster each for ARC and VMH.*

### Data resources

##### Analysis and plotting code
Code for analyses is contained within subdirectories of this repository: 
- Visium: [spatial_HYP](https://github.com/LieberInstitute/spatial_HYP/tree/main/spatial_HYP)
- Xenium: [xenium_HYP](https://github.com/LieberInstitute/spatial_HYP/tree/main/xenium_HYP)
- For manuscript plots: [manuscript_plot_code](https://github.com/LieberInstitute/spatial_HYP/tree/main/manuscript_plot_code)
- For *KISS1*-*ESR1*-*TAC3* smFISH analysis: [here](https://github.com/LieberInstitute/spatial_HYP/tree/main/Br1225_ESR1-TAC3-KISS1_smfish_analysis)

##### Data visualization
Interactive web portals (one for Visium, one for Xenium) have been made available using [Samui](https://samuibrowser.com/) to allow visualization of gene expression and cluster assignments, as well as creation of custom spatial annotations. Documentation is built into these respective visualization tools, which can be found at:

- [Documentation for both browsers](https://github.com/LieberInstitute/spatial_HYP/tree/main/samui_docs)
- [HYP Visium Browser](https://samuibrowser.com/from?url=data.libd.org/samuibrowser/&s=V12D05-348_C1&s=V12D05-348_D1&s=V12D05-350_C1&s=V12D05-350_D1&s=V12D07-075_A1&s=V12D07-075_D1&s=V12Y31-080_A1&s=V13M13-362_A1&s=V13M13-362_D1&s=V13Y24-346_C1)
- [HYP Xenium Browser](https://samuibrowser.com/from?url=data.libd.org/samuibrowser/&s=X36_5459A&s=X36_5459B&s=X36_8667C&s=X86_reg1&s=X86_reg2&s=X86_reg3&s=X97_reg1&s=X97_reg2&s=X97_reg3&s=X99_1225A&s=X99_1225B&s=X99_8741C&s=X99_8741D)

##### Supplemental Data S1-S6 and smFISH microscopy
Supplemental data mentioned in the manuscript, along with smFISH images on Visium/Xenium-adjacent tissue for LAMP5 and KISS1, are available through are available through [http://research.libd.org/globus/](http://research.libd.org/globus/) via the Globus endpoint jhpce#HYP_suppdata. Readme files for these data are also available through the endpoint.

##### Further data availability
FASTQ files from Visium and unfiltered sample-level xeniumranger 1.7 outputs (transcript-level data) compatible with Xenium Explorer software are respectively available through GEO accessions [GSE280316](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280316) and [GSE280460](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280460). Addditional processed data may be added to the Globus endpoint in the process of review and publication. 

### How to Cite
PENDING
