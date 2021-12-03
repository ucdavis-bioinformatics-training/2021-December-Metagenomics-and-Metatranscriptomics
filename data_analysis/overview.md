# Overview on metagenomics and metatranscriptomics
-------------------------------------------------------

<p align = "center">
<img src="overview_figures/P1.jpg" alt="micribial" width="60%"/>
</p>

<p align = "right" style="font-family:Times;font-size:12px;">
Vanwonterghem, etc., 2014, Curr. Opin. Biotechnol.
</p>

Microbiomes are ubiquitous in the environment and play essential roles in every aspects of nature and life. In recent years, more effort has been put in studying these microorganisms, facilitated by advances in technology, especially the sequencing technologies, and the significant decrease in the cost. There are five major areas in microbiom study, metagenomes, metatranscriptomes, metaproteome, microbial metabolomes, and image based cell activity studies. Each of these areas focuses on one aspect of the microbiome community and provides insights into the nature of the microbial communities. In this workshop, we are going to focus on only two areas: metagenomics and metatranscriptomics.

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

---


----------------------------------------

### Mark duplicates

This step is optional if you have done deduplication in read preprocessing step, as we have done using hts_SuperDeduper. The steps below serves as a helper if you need it in your own data analysis, where you do not do deduplication in the data preprocessing step. 


---

### Base quality score recalibration (BQSR)

<br>

##### <font color='red'> Stop Group Exercise 1: </font>


