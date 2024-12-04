This is the [Samui](samuibrowser.com) browser for adult human hypothalamus (VMH/ARC) **Visium** data from Mulvey, et. al. If you're looking for the Xenium browser, please go [here](https://samuibrowser.com/from?url=data.libd.org/samuibrowser/&s=X36_5459A&s=X36_5459B&s=X36_8667C&s=X86_reg1&s=X86_reg2&s=X86_reg3&s=X97_reg1&s=X97_reg2&s=X97_reg3&s=X99_1225A&s=X99_1225B&s=X99_8741C&s=X99_8741D).

Information to match sample IDs as listed in browser are in a table at the end of this document. 

Using the dropdown menu at the top of this sidebar, you can choose to visualize log-normalized gene expression or clustering results (using spatially variable genes as features for dimensionality reduction for clustering with [BayesSpace](https://www.nature.com/articles/s41587-021-00935-2)). You can then use the adjacent text entry field to filter the selected feature (choices will appear if you place your text cursor here). The different choices in the dropdown menu are: 
* *Gene*: normalized log counts in analyzed Visium spots
* *HARMONYlmbna nnsvg10 k15 collapsed*: Clustering at *k*=15, with VMH or ARC clusters collapsed together as "VMH" or "ARC". (This is the clustering used in the main results from the manuscript.)
* *HARMONYlmbna nnsvg10 k15*: *k*=15, but with individual VMH and ARC clusters shown
* *HARMONYlmbna nnsvg10 k20*: *k*=20 with individual VMH and ARC clusters labeled as VMH20.1, ARC20.1..
* *HARMONYlmbna nnsvg10 k31*: *k*=31, with VMH/ARC clusters labeled as above

#### Visualizing two features/annotations simultaneously:
To visualize two separate features at once (one as outlines and the other as fill for each point (cell)), click the "Add Layer" button in the top right of the viewer interface, then choose a second feature. In the list of shown features, the left checkbox toggles whether a feature and filter is displayed as the outline of each point, and the right checkbox will toggle whether that feature and filter is displayed as the fill for each point.

#### Additional resources of note:
* [Analysis and plot code](https://github.com/LieberInstitute/spatial_HYP/)
* Preprint: URL PENDING
* [Globus endpoint](http://research.libd.org/globus/): jhpce#HYP_suppdata


#### Samui/Manuscript Sample ID Conversion

| Samui ID      | Manuscript ID | Sex    | Donor  |
| ------------- | ------------- | ------ | ------ |
| V13Y24-346_C1 | v1225A_M      | Male   | Br1225 |
| V12D07-075_D1 | v1735A_M      | Male   | Br1735 |
| V13M13-362_A1 | v5459A_M      | Male   | Br5459 |
| V12D07-075_A1 | v5993A_F      | Female | Br5993 |
| V12Y31-080_A1 | v6197A_M      | Male   | Br6197 |
| V12D05-348_D1 | v6588A_F      | Female | Br6588 |
| V12D05-350_D1 | v6588B_F      | Female | Br6588 |
| V13M13-362_D1 | v8667A_F      | Female | Br8667 |
| V12D05-348_C1 | v8741A_F      | Female | Br8741 |
| V12D05-350_C1 | v8741B_F      | Female | Br8741 |

