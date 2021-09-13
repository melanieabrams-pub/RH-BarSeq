# RH-BarSeq

Coming soon.


Pipeline for RH-seq with barcoding (Reciprocal hemizygosity analysis via sequencing, with barcoding)

The RH-seq v2 pipeline found here is modified by MAbrams from RH-seq pipeline from Weiss et al., 2018 (https://github.com/weiss19/rh-seq) with modifications from Abrams and Dubin, 2020 (https://github.com/melanieabrams-pub/thermotolerance-loci-across-yeasts) as described in the associated publications.

(1) Count barcodes with the RBseq pipeline RBseq_Count Barcodes (v 1.1.4, https://github.com/stcoradetti/RBseq/tree/master/Old_Versions/1.1.4)

(2) Run the entire rh-seq analysis pipeline at once with the optional wrapper script - or run one script at a time.  



This repository also includes the jupyter notebook by JSkerker for annotation of the poolFile for RBseq analysis of the yeast Tn-seq in the Sc x Sp hybrid pool used in this study, as well as the GFF3-compliant annotation file used by that jupyter notebook.  Elements of the jupyter notebook are modified from the RBseq annotation script described in Coradetti et al., 2018 (https://github.com/stcoradetti/RBseq)

Metafiles and associated files used in the manuscript are provided for convenient replication of results; to do so, adjust file pathways in the metafiles to reflect their new storage directories upon download.
