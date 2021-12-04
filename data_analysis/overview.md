# Overview on metagenomics and metatranscriptomics
-------------------------------------------------------

<p align = "center">
<img src="overview_figures/P1.jpg" alt="micribial" width="60%"/>
</p>

<p align = "right" style="font-family:Times;font-size:12px;">
Vanwonterghem, etc., 2014, Curr. Opin. Biotechnol.
</p>

Microbiomes are ubiquitous in the environment and play essential roles in every aspects of nature and life. Despite its importance, the research on the mcirobial community was severely limited because of the difficulty in cell culturing and separation. In recent years, more effort has been put in studying these microorganisms, facilitated by advances in technology, especially the sequencing technologies, and the significant decrease in the cost.

There are five major areas in microbiom study, metagenomes, metatranscriptomes, metaproteome, microbial metabolomes, and image based cell activity studies. Each of these areas focuses on one aspect of the microbiome community and provides insights into the nature of the microbial communities. In this workshop, we are going to focus on only two areas: metagenomics and metatranscriptomics.


---

## Metagenomics

The field of metagenomics is to study the composition of the microbial communities. There are two approaches available in this field:

* Amplicon sequencing. It only sequences rRNA related genomic materials, such as the most popular 16S target sequencing, and ITS target sequencing for fungal genomes.
  * lower per sample cost
  * using Illumina platform, capable of pooling thousands, or even tens of thousands of barcoded samples/targets per sequening run
  * taxonomy profiling

* __Shotgun-based sequencing. It allows for sequencing of DNA materials from all organisms within a community.__
  * taxonomy profiling
  * metagenomic assembly and binning


Shotgun-base metagenomics approach has the advantage of producing not only the composition of the community, but also provide information on potential functions from all organisms present in a given complex sample. It allows the discovery of previously unknow species even when they could not be classified to known taxa.


<p align = "center">
<img src="overview_figures/P2.webp" alt="micribial" width="60%"/>
</p>

<p align = "right" style="font-family:Times;font-size:12px;">
Quince, C., etc., 2017, Nat. Biotechnol.
</p>

For each of the steps along the metagenomics data analysis using shotgun sequencing, many tools have been developed or adapted from regular genomics data analysis.


<p align = "center">
<img src="overview_figures/P3.png" alt="micribial" width="60%"/>
</p>

---


## Metatranscriptomics

Metatranscriptomics is to use RNASeq technology to profile expressed genes in a microbial community. It allows for the profiling of the composition of the community as well as the activities of the community from gene expression point of view. By mapping the expressed genes to metabolic pathways, it presents a picture of the functional aspect of the microbial community. It distinguishes an active member from a non-active member. It also offers a tool to study the responses that the microbial community has to their changing environmental conditions.


<p align = "center">
<img src="overview_figures/P4.png" alt="micribial" width="60%"/>
</p>

<p align = "right" style="font-family:Times;font-size:12px;">
modified from: Quince, C., etc., 2017, Nat. Biotechnol.
</p>


The analysis of metatranscriptomics borrow some tools from the analysis of metagenomics. For example, for assemblying the metatranscriptome, the assembly tools first designed for metagenomics tools have been used. It also borrows some tools from the regular transcriptomic assembly, such as Trinity, which has been shown to improve the number of annotated genes from the assembled transcripts. There are a few tools designed specifically for metatranscriptomes and take into account the unique features of both transcripts and the complex nature of microbial communities, such as IDBA-MT, TAG.


<p align = "center">
<img src="overview_figures/P5.png" alt="micribial" width="60%"/>
</p>

---

