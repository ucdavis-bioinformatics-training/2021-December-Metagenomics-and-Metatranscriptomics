# Metatranscriptomics data analysis
-------------------------------------------------------

This document assumes [preprocessing using HTStream](./01-preproc_htstream_mm.md) has been completed.

The main objectives in metatranscriptomics data analysis is to answer the question:  __what are they doing__. At the same time, it can also answer the question: __who is there__. These two concepts we have seen in the metagenomics data analysis part of this workshop. But we are going to see what is different when we have metatranscriptomics data.

---

<p align = "center">
<img src="overview_figures/P7.png" alt="micribial" width="70%"/>
</p>

---

## Remove/Separate rRNA

Ribosomal RNA is by far the most abundant form of RNA in most cells and can make up to 90% of the sequencing data. They perform essential cellular functions and offer information on a community's structure and have been used in traditional taxonomic profiling of microbial communities. Nonetheless, they do not provide more information on a community. If one is mainly interested in studying the functional units of a microbial community, these rRNA provides no useful information and takes up the majority of computing resources.

<p align = "center">
<img src="metatranscriptome_figures/TvsM.webp" alt="micribial" width="92%"/>
</p>

<p align = "right" style="font-family:Times;font-size:12px;">
Li, etc., Microbiome 7,6 (2019), https://doi.org/10.1186/s40168-019-0618-5
</p>

There are many packages created to perform this task. They all attempt to match the sequencing reads to sequences in rRNA databases. They vary wildly in the requirement in computing resources, but have very similar accuracy.

<p align = "center">
<img src="metatranscriptome_figures/depleterrna.png" alt="micribial" width="92%"/>
</p>

<p align = "right" style="font-family:Times;font-size:12px;">
Li, etc., Microbiome 7,6 (2019), https://doi.org/10.1186/s40168-019-0618-5
</p>


