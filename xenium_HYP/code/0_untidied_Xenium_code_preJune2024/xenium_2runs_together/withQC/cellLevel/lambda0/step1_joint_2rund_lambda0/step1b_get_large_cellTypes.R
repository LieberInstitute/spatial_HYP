dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
setwd(dir_out)

cellTypeAnnotations = readRDS("cellTypeAnnotations.rds")
categorized_cell_types <- list(
  GABAergic_neurons = c(
    "cluster 1, Other_GABAergic",
    "cluster 41, GABAergic_unspecified"
  ),
  Microglia = c(
    "cluster 13, Microglia_1",
    "cluster 16, Microglia_(Reactive)",
    "cluster 27, Microglia_2"
  ),
  Astrocytes = c(
    "cluster 17, Mixed_Astrocyte-Arc_POMC",
    "cluster 4, Astrocyte",
    "cluster 7, Astrocyte_2",
    "cluster 30, Astrocyte_3"
  ),
  Oligodendrocytes = c(
    "cluster 2, Oligodendrocyte_Precursor_cells_(OPC)",
    "cluster 11, Oligodendrocyte_3_(OPALIN+)",
    "cluster 6, Oligodendrocyte",
    "cluster 8, Oligodendrocyte",
    "cluster 10, Oligodendrocyte_2_(OPALIN+)",
    "cluster 26, Oligodendrocyte_4_(Optic_Nerve)",
    "cluster 39, Mixed_Oligodendrocytes_and_Supraoptic_nu_(AVP+OXT+)_3"
  ),
  ARC_neurons = c(
    "cluster 12, ARC_1_(POMC)",
    "cluster 33, ARC_3_(TAC3_ESR1)",
    "cluster 25, ARC_2_(GHRH_GAL)"
  ),
  Vascular_cells = c(
    "cluster 9, Vascular_uncertain",
    "cluster 5, Endothelium"
  ),
  Tanycytes = c(
    "cluster 15, Tanycyte_1",
    "cluster 28, Tanycyte_2"
  ),
  VMH_neurons = c(
    "cluster 32, VMH_Lateral_border",
    "cluster 34, VMH_4",
    "cluster 23, VMH_3",
    "cluster 36, VMH_5",
    "cluster 18, VMH_1_(Excitatory)",
    "cluster 19, VMH_2_(Excitatory)",
    "cluster 40, VMH_6_Lateral_bit"
  ),
  Immune_cells = c(
    "cluster 24, Macrophage",
    "cluster 29, Peripheral_immune_NK-or_Tcell"
  ),
  Supraoptic_nucleus_neurons = c(
    "cluster 22, Supraoptic_nu_(AVP+OXT+)_1",
    "cluster 31, Supraoptic_nu_(AVP+OXT+)_2"
  ),
  Other_neurons = "cluster 14, Excitatory_Periventricular_(CRH,_TRH)",
  Discarded_or_unsure_clusters = c(
    "cluster 3, Insufficient_markers",
    "cluster 20, Unsure1",
    "cluster 35, DISCARD",
    "cluster 37, DISCARD",
    "cluster 21, DISCARD",
    "cluster 38, DISCARD"
  )
)

saveRDS(categorized_cell_types, file = "categorized_cell_types.rds")
names_large_clu ="Oligodendrocytes"
vec_large_clus = c()
for(names_large_clu in names(categorized_cell_types)){
  vec_large_clus_tmp = categorized_cell_types[[names_large_clu]]
  vec_large_clus_tmp_ =  rep(names_large_clu, length(vec_large_clus_tmp))
  # names(vec_large_clus_tmp_) = 
  vec_large_clus_tmp = sub(",.*","",vec_large_clus_tmp)
  vec_large_clus_tmp = sub("cluster ","",vec_large_clus_tmp)
  
  names(vec_large_clus_tmp_)  = as.character(vec_large_clus_tmp)
  vec_large_clus = c(vec_large_clus, vec_large_clus_tmp_)
}

saveRDS(vec_large_clus, file = "vec_large_clus.rds")

# > length(categorized_cell_types)
# [1] 12

