This is the [Samui](samuibrowser.com) browser for adult human hypothalamus (VMH/ARC) **Xenium** data from Mulvey, et. al. If you're looking for the Visium browser, please go [here](libd.org). 

Information to match sample IDs as listed in browser to those in the manuscript are in a table at the end of this document. 

Using the dropdown menu at the top of this sidebar, you can choose to visualize cell-level log-normalized gene expression, clustering results and domain (VMH/ARC/other) assignments, and QC metrics. You can then use the adjacent text entry field to filter the selected feature (choices will appear if you place your text cursor here). The different choices in the dropdown menu are:
* *Banksy-Cluster*: The 33 cluster labels retained after removing sample/donor-specific clusters.
* *Banksy Cluster Group*: Clusters from above aggregated into broader cell types (e.g., astrocyte clusters all under one "astrocyte" label)
* *Domain of Cell*: The VMH/ARC/other assignment for each cell determined using k nearest-neighbors smoothing as described in the manuscript
* *nCounts*: total number of gene counts for each cell
* *nGenes*: total number of unique genes detected in each cell
* *Cell Area*: area of each xeniumranger 1.7-segmented cell
* *Nucleus Area*: area of each xenium ranger 1.7-segmented nucleus

#### Visualizing two features/annotations simultaneously:
To visualize two separate features at once (one as outlines and the other as fill for each point (cell)), click the "Add Layer" button in the top right of the viewer interface, then choose a second feature. In the list of shown features, the left checkbox toggles whether a feature and filter is displayed as the outline of each point, and the right checkbox will toggle whether that feature and filter is displayed as the fill for each point.

#### Additional resources of note:
* [Analysis and plot code](https://github.com/LieberInstitute/spatial_HYP/)
* Preprint: URL PENDING
* [Globus endpoint](http://research.libd.org/globus/): jhpce#HYP_suppdata

| Samui ID  | Manuscript ID | Sex    | Donor  |
| --------- | ------------- | ------ | ------ |
| X99_1225A | x1225B_M      | Male   | Br1225 |
| X99_1225B | x1225C_M      | Male   | Br1225 |
| X97_reg1  | x1735B_M      | Male   | Br1735 |
| X97_reg2  | x1735C_M      | Male   | Br1735 |
| X36_5459A | x5459B_M      | Male   | Br5459 |
| X36_5459B | x5459C_M      | Male   | Br5459 |
| X86_reg2  | x5993B_F      | Female | Br5993 |
| X86_reg3  | x5993C_F      | Female | Br5993 |
| X86_reg1  | x6197B_M      | Male   | Br6197 |
| X97_reg3  | x6588C_F      | Female | Br6588 |
| X36_8667C | x8667B_F      | Female | Br8667 |
| X99_8741C | x8741C_F      | Female | Br8741 |
| X99_8741D | x8741D_F      | Female | Br8741 |