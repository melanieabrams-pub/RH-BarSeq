# RH-BarSeq


<b>Pipeline for reciprocal hemizygosity analysis via sequencing, with barcoding (RH-seq with barcoding)</b>


<b>Summary:</b>

The RH-seq with barcoding pipeline found here is described in Abrams and Chuong et al., 2021 (https://doi.org/10.1093/g3journal/jkab412).  Briefly, the analysis uses a modified version of the RBseq pipeline from Coradetti et al. 2018 (v1.1.4) (https://github.com/stcoradetti/RBseq/tree/master/Old_Versions/1.1.4) for the Tn-seq mapping of barcode to transposon insert location, and then to count the barcodes from Barseq competitions; then, the RH-seq analysis is performed with a codebase modified by MAbrams from RH-seq pipeline from Weiss et al., 2018 (https://github.com/weiss19/rh-seq) with modifications from Abrams and Dubin, 2021 (https://github.com/melanieabrams-pub/thermotolerance-loci-across-yeasts) as described in the associated publications.

<b> Steps:</b>

(1) Map the Tn-seq pool with the RBseq_Map_Insertions_v1.1.4_PBa_Jskerker.py and associated metafiles. (See RBseq docs for further details on RBseq metafiles)

(2) Annotate the Tn-seq poolfile with the custom Jupyter notebook annotate_poolfile_for_barseq_v1.2_JSkerker.ipynb and associated metafiles.

(3) Count barcodes with the RBseq pipeline RBseq_Count Barcodes (v 1.1.4, https://github.com/stcoradetti/RBseq/tree/master/Old_Versions/1.1.4) 

(4) Run the entire RH-seq analysis pipeline at once with the optional wrapper script - or run one script at a time.  


<b>Metafiles:</b>

Metafiles and associated files used in the manuscript are provided for convenient replication of results; to do so, adjust file pathways in the metafiles to reflect their new storage directories upon download. 

<b>Environment:</b>

There are several package dependencies for this pipeline.

Prerequisites:
-Python 3.7
-Biopython 1.72 or later
-Anaconda scientific computing environment or the following python libraries: numpy, pandas, matplotlib, json, scipy, and stats models
-NCBI BLAST+ 2.2.30 or later

Most recently tested with the following versions of python packages:
python: 3.7.6
biopython: 1.76
numpy: 1.18.5
pandas: 1.0.5
scipy: 1.5.0
matplotlib: 3.2.2
statsmodels: 0.11.1


<b>Additional credits</b>:

The RBseq codebase from Coradetti et al., 2018, and the RH-seq codebase from Weiss et al., 2018, and the modifications of Abrams and Dubin et al., 2021, form the basis of this pipeline.  This repository includes the jupyter notebook by JSkerker for annotation of the poolFile for RBseq analysis of the yeast Tn-seq in the Sc x Sp hybrid pool used in this study, as well as the GFF3-compliant annotation file used by that jupyter notebook.  Elements of the jupyter notebook are modified from the RBseq annotation script described in Coradetti et al., 2018 (https://github.com/stcoradetti/RBseq)
