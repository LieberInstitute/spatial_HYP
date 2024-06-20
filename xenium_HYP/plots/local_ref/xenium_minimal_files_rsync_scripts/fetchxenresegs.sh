#!/bin/bash

commonpath=bmulvey@jhpce-transfer01.jhsph.edu:/dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/xenium_HYP/processed-data/01_xeniumranger-resegment1.7_20231019_slides

for outdir in $(cat outdirs.txt)
do
	mkdir nucresegtmp/${outdir}
	for outfile in $(cat outfiles.txt)
	do
		rsync ${commonpath}/${outdir}/outs/${outfile} nucresegtmp/${outdir}
	done
done