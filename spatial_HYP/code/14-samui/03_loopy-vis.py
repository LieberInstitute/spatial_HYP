from pathlib import Path
from pyhere import here
import json
import os
import re
import scanpy as sc

import numpy as np
import pandas as pd
from rasterio import Affine

from loopy.sample import Sample
from loopy.utils.utils import remove_dupes, Url

spot_diameter_m = 55e-6 # 55-micrometer diameter for Visium spot
img_channels = 'rgb'
# default_channels = {'blue': 'DAPI', 'red': 'NeuN'}
default_gene = 'SNAP25'

#   Names of continuous features expected to be columns in the observation 
#   data (colData) of the AnnData
# spe_cont_features = ['PpTau', 'PAbeta']

#   Diagnosis by brain number (not included in the sample_info sheet)
# sample_dx = {
#     'Br3854': 'AD', 'Br3873': 'AD', 'Br3880': 'AD', 'Br3874': 'control'
# }

sample_info_path = here(
    'spatial_HYP','raw-data', 'sample_info', 'demos_by_sampleid.xlsx'
)

spe_path = here("spatial_HYP","processed-data", "14-samui", "spe.h5ad")
# notes_path = str(Path(here('code', '16_samui', 'feature_notes.md')).resolve())
img_path = here('spatial_HYP', 'processed-data', '14-samui', 'rotated_vis_tifs', '{}.tif')
json_path = here(
    'spatial_HYP', 'processed-data', '01_spaceranger', 'all-spaceranger', '{}', 'outs', 'spatial',
    'scalefactors_json.json'
)

out_dir = here('spatial_HYP','processed-data', '14-samui', 'visium', '{}')

################################################################################
#   Read in sample info and clean
################################################################################


sample_info = pd.read_excel(sample_info_path)
sample_info = sample_info.dropna(subset=['brnum', 'Tissue'])

#   Prepend "Br" tp brain number and make a string
sample_info['br_num'] = sample_info['brnum'].astype(str)

#   Add diagnosis using brain number
# sample_info['diagnosis'] = sample_info['br_num'].replace(sample_dx)

#   Fix the experiment number column (use strings of integers), create experiment number for each slide
sample_info['sample_num'] = sample_info.index + 1
sample_info['experiment_num'] = ((sample_info['sample_num'] - 1) // 4  + 1).astype(str)

# create column called sample_id from 'Sample #'
#pattern = '_shk'
sample_info['sample_id'] = sample_info['sample_id']


#sample_info['spaceranger_id'] = sample_info['Slide #'] + '_' + sample_info['Array #']
#sample_info['image_id'] = sample_info['Slide #'] + '_' + sample_info['Array #']
sample_info['spaceranger_id'] = sample_info['sample_id']
sample_info['image_id'] = sample_info['sample_id']

# sample_id_spaceranger = sample_info['Slide #'] + '_' + sample_info['Array #']
# sample_id_image = sample_info['Slide #'] + '_' + sample_info['Array #']

#   Subset all types of IDs to this sample only
sample_id_spaceranger = sample_info['spaceranger_id'].iloc[int(os.environ['SLURM_ARRAY_TASK_ID']) - 1]
sample_id_image = sample_info['image_id'].iloc[int(os.environ['SLURM_ARRAY_TASK_ID']) - 1]
sample_id_samui = sample_info['spaceranger_id'].iloc[int(os.environ['SLURM_ARRAY_TASK_ID']) - 1]


#sample_id_spaceranger = sample_info['spaceranger_id'].iloc[int()]
# #   Update paths for this sample ID
out_dir = Path(str(out_dir).format(sample_id_samui))
json_path = Path(str(json_path).format(sample_id_spaceranger))
img_path = Path(str(img_path).format(sample_id_image))

# #find size of dataframe and assign slide + array info to sample_id_spaceranger
# samplenums = sample_info['sample_num']
# samplerows = list(range(len(sample_info)))
# sample_id_spaceranger = sample_info['Slide #'] + '_' + sample_info['Array #']

# for i in samplerows:
#     sample_id_img = sample_id_spaceranger[i]
    
#     out_dir = Path(str(out_dir).format(sample_id_img))
#     json_path = Path(str(json_path).format(sample_id_img))
#     img_path = Path(str(img_path).format(sample_id_img))

out_dir.mkdir(exist_ok = True)

#   All paths should exist
assert all([x.exists() for x in [out_dir, json_path, img_path]])

################################################################################
#   Read in scale-factors info
################################################################################

#   Read in the spaceranger JSON to calculate meters per pixel for
#   the full-resolution image
with open(json_path, 'r') as f:
    spaceranger_json = json.load(f)

m_per_px = spot_diameter_m / spaceranger_json['spot_diameter_fullres']

################################################################################
#   Gather gene-expression data into a DataFrame to later as a feature
################################################################################

#   Read in AnnData and subset to this sample
spe = sc.read(spe_path)
#path_groups = spe.obs['path_groups'].cat.categories
spe = spe[spe.obs['sample_id'] == sample_id_spaceranger, :]
spe.obs.index.name = "barcode"

#   Convert the sparse gene-expression matrix to pandas DataFrame, with the
#   gene symbols as column names
gene_df = pd.DataFrame(
    spe.X.toarray(),
    index = spe.obs.index,
    columns = spe.var['gene_name']
)

#   Some gene symbols are actually duplicated. Just take the first column in
#   any duplicated cases
gene_df = gene_df.loc[: , ~gene_df.columns.duplicated()].copy()

#   Samui seems to break when using > ~ 5,000 genes. Take just the genes where
#   at least 10% of spots have nonzero counts
#gene_df = gene_df.loc[:, np.sum(gene_df > 0, axis = 0) > (gene_df.shape[0] * 0.1)].copy()

assert default_gene in gene_df.columns, "Default gene not in AnnData"

print('Using {} genes as features.'.format(gene_df.shape[1]))

################################################################################
#   Split 'path_groups' column into binary columns for each of its values
################################################################################

#   Circumvent a Samui bug (https://github.com/chaichontat/samui/issues/84);
#   turn the categorical column 'path_groups' into several numeric columns with
#   just values of 0 and 1
# path_df = pd.DataFrame()
# for path_group in path_groups:
#     path_df[path_group] = (spe.obs['path_groups'] == path_group).astype(int)

################################################################################
#   Use the Samui API to create the importable directory for this sample
################################################################################

this_sample = Sample(name = sample_id_samui, path = out_dir)

this_sample.add_coords(
    spe.obsm['spatial'].rename(
        columns = {'pxl_col_in_fullres': 'x', 'pxl_row_in_fullres': 'y'}
    ),
    name = "coords", mPerPx = m_per_px, size = spot_diameter_m
)

#   Add the IF image for this sample
this_sample.add_image(
    tiff = img_path, channels = img_channels, scale = m_per_px
)

#   Add gene expression results (multiple columns) as a feature
this_sample.add_chunked_feature(
    gene_df, name = "Genes", coordName = "coords", dataType = "quantitative"
)

#   Add additional requested observational columns (colData columns)
this_sample.add_csv_feature(
    spe.obs["HARMONYlmbna_nnsvg10_k15_collapsed"], name = "HARMONYlmbna nnsvg10 k15 collapsed", coordName = "coords",
    dataType = "categorical"
)

#   Add additional requested observational columns (colData columns)
this_sample.add_csv_feature(
    spe.obs["HARMONYlmbna_nnsvg10_k15"], name = "HARMONYlmbna nnsvg10 k15", coordName = "coords",
    dataType = "categorical"
)

#   Add additional requested observational columns (colData columns)
this_sample.add_csv_feature(
    spe.obs["HARMONYlmbna_nnsvg10_k20"], name = "HARMONYlmbna nnsvg10 k20", coordName = "coords",
    dataType = "categorical"
)

#   Add additional requested observational columns (colData columns)
this_sample.add_csv_feature(
    spe.obs["HARMONYlmbna_nnsvg10_k31"], name = "HARMONYlmbna nnsvg10 k31", coordName = "coords",
    dataType = "categorical"
)


#   Add pathology groups
# this_sample.add_csv_feature(
#     path_df, name = "Pathology Group", coordName = "coords",
#     dataType = "quantitative"
# )

this_sample.set_default_feature(group = "Genes", feature = default_gene)

this_sample.write()
