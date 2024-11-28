# Spatially-resolved molecular sex differences at single cell resolution in the adult human hypothalamus

[![DOI](https://zenodo.org)]() PENDING

### Study Design
The hypothalamus is a critical brain area underlying functions with inherent sex differences, such as reproductive physiology, endocrine signaling, and metabolism. In the rodent, these sex-differentiated functions correspond to differences in volume, cell type composition, and gene expression between males and females across individual hypothalamic regions (here, "domains"). The ventromedial hypothalamus (VMH) and arcuate nucleus (ARC) are two hypothalamic regions that influence appetitive/social behaviors and growth/metabolism, respectively. While molecular profiling studies in the rodent hypothalamus have identified specialized cell types with unique transcriptomic signatures, there is a paucity of data describing the molecular architecture of the human HYP, especially in the context of sex-differentiated cell types that drive evolutionarily essential, hypothalamus-mediated behaviors in males and females. 

This study, led by Bernard (Bernie) Mulvey, Kristen Maynard, and Kasper Hansen, profiled the adult human mediobasal hypothalamus (ventromedial nucleus, VMH, and arcuate nucleus, ARC) using  Visium and Xenium spatial transcriptomic platforms (10x Genomics). Atlasing and sex-differential expression efforts using Visium data from 8 donors (4 per sex) were used to guide gene selection for subsequent Xenium assays on adjacent tissue sections from these same 8 donors.

![Experimental Overview](./images/overview.png)

This dataset spanning 23 samples identified 5 ARC and 4 VMH neuronal populations governing known hypothalamus-specific functions and defined their spatial distributions. Compared to rodent VMH and ARC, we found increases in retinoid pathway gene expression in these domains. Sex-DE analysis within VMH and ARC revealed correlated autosomal expression differences, which were localized to *ESR1*- and *TAC3*-expressing neurons in the ARC, and *CRHR2*-expressing neurons in the VMH. VMH- and ARC-residing cell types have a striking number of sex-DE genes linked to sex-biased disorders, including autism, depression, and schizophrenia. By mapping disease associations to hypothalamic regions containing cell types with established roles in mediating sex-divergent physiology and behavior, these data provide insights into mechanistic bases of sex bias in neurodevelopmental and neuropsychiatric disorders.

![Sex DE Analysis with Xenium](./images/Xenium_sex_DE_analysis_schematic.png)
*Sex-DE analysis was performed at the level of Xenium clusters by first subsetting to all cells within the boundaries of the VMH or ARC. Then, cells of each type are tested separately for sex-DE within that domain. Thus, broadly distributed cell clusters/types (e.g., glia) are tested for sex-DE in each domain. Meanwhile, sex-DE testing of domain-specific neuronal clusters (e.g.,* TAC3*-*ESR1 *in ARC) is filtered to the cells given that label and found in their domain (i.e. in the expected anatomic space). This simplified schematic only shows one domain-specific cluster each for ARC and VMH.*

### Data resources

##### Analysis and plotting code
Code for analyses is contained within subdirectories of this repository: for Visium, (spatial_HYP)[https://github.com/LieberInstitute/spatial_HYP/tree/main/spatial_HYP], and for Xenium, (xenium_HYP)[https://github.com/LieberInstitute/spatial_HYP/tree/main/xenium_HYP]. Code for creating the manuscript plots is in [manuscript_plot_code](https://github.com/LieberInstitute/spatial_HYP/tree/main/manuscript_plot_code). Code for analysis of *KISS1*-*ESR1*-*TAC3* smFISH can also be found (here)[https://github.com/LieberInstitute/spatial_HYP/tree/main/Br1225_ESR1-TAC3-KISS1_smfish_analysis].

##### Data visualization
We have created interactive web portals (one for Visium, one for Xenium) using [Samui](https://samuibrowser.com/) that allows for visualization of gene expression and cluster assignments, as well as creation of custom spatial annotations:

- HYP Visium Browser: PENDING
- HYP Xenium Browser: PENDING

##### Supplemental Data S1-S6 and smFISH microscopy
Supplemental data mentioned in the manuscript, along with smFISH images on Visium/Xenium-adjacent tissue for LAMP5 and KISS1, are available through the Globus endpoint PENDING. Readme files for these data are also available through the endpoint.

##### Further data availability
Raw data will be available through GEO, and (additional) processed data through the Globus endpoint above, at the time of publication. 

### How to Cite
PENDING