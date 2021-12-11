# Metagenomics data analysis
-------------------------------------------------------

This document assumes [preprocessing using HTStream](./01-preproc_htstream_mm.md) has been completed.

The two main objectives in metagenomics data analysis are to answer the questions: __who is there__ and __what can they potentially do__. The first is what we call taxonomy profiling and the second is function profiling.

Taxonomy profiling has traditionally been done using inexpensive 16S rRNA amplicon sequencing. It has limited resolution because of the conservation of the target gene and the length of the amplicon product. It also limits the ability to profile non-bacterial species of a community. In addition, 16S has very limited ability to provide us with the functional capacity of the microbes. Furthermore, [studies](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4523815/) have shown there might be biases introduced during the amplification step to enrich for rRNA gene and in turn introduces biases in quantification fo taxa.

Shotgun metagenomics approach has gained popularity over the recent years in stuying microbial community. It benefits from the constant decrease in sequencing technology and it allows the study of all microorganisms that are difficult to culture. By sequencing the whole DNA materials from a community, it not only permits the discovery of unknown taxa, but also the ability of predicting the functions of microbial members. We will focus only on the shotgun metagenomics data analysis in this workshop.

---

<p align = "center">
<img src="overview_figures/P6.png" alt="micribial" width="70%"/>
</p>

---

## Remove host DNA

Depending on the sample collection source for microbial studies, there are variable amount of host DNA. For example, bovine rumen samples such as the dataset we are using for this workshop, or fecal samples for gut microbiome studies, or tissue biopsies, all have some leve of host DNA remain. Because host genomes are much larger than the microbiome genomes, the sequencing product will have more host reads than microbial reads. For example, if one sequence the same number of human cells as microbial cells, the sequencing would produce over 99% of human reads. Even when microbiome samples are collected from environment/fields, human contamination happens more than we think. Therefore, removing host/contaminant DNA very important. Especially, when the objective is to produce high quality metagenome assembled genomes. It also reduces biases and minimize false positive associations in downstream analysis.

<p align = "center">
<img src="metagenome_figures/hostdna.jpeg" alt="micribial" width="85%"/>
</p>

<p align = "right" style="font-family:Times;font-size:12px;">
Pereira-Marques, etc., 2019, Front. Microbiol. https://doi.org/10.3389/fmicb.2019.01277
</p>

For samples that are suppected to have very high levels of host material, it is worth considering adding a step in the library prep stage to remove host material. I am not going to talk about it because a couple of our lunch sponsors will talk in detail about it. Even with host depletion at the library prep stage, it is important to remove leftover host reads bioinformatically.

[A publication](https://pubmed.ncbi.nlm.nih.gov/32558637/) has done a very thorough study using different software/approaches to detect/remove human reads in microbial datasets. If we are to pick one method, Bowtie2 performs the best overall to remove human reads from the microbial datasets. So, we have chosen to use bowtie2 to remove the host DNA from our dataset.

<p align = "center">
<img src="metagenome_figures/rmhost.gif" alt="micribial" width="85%"/>
</p>

<p align = "center">
<img src="metagenome_figures/rmhost2.gif" alt="micribial" width="85%"/>
</p>

<p align = "center">
<img src="metagenome_figures/rmhost3.gif" alt="micribial" width="85%"/>
</p>

<p align = "right" style="font-family:Times;font-size:12px;">
Bush, etc., 2020, Microb. Genom. https://doi.org/10.1099/mgen.0.000393
</p>

Based on the findings from the above paper, we are going to use Bowtie2 to remove host bovine DNA. The basic idea is to map/align the preprocessed reads to the bovine reference genome. Then we extract the reads that do not align for downstream analysis. Here we are going to sidetrack a little bit to learn two new concepts: __map/align__ vs __assembly__.


### Mapping vs Assembly

Every organism on earth has its own genome sequence, regardless whether we know how it looks like. The main difference between mapping and assembly is whether we know ahead of time the sequence of the genome. When we have the sequence of a genome and we would like to find out from which genomic locations our sequencing data came from, it uses the mapping approach. On the other hand, when we do not have the sequence of a genome and would like to reconstruct it, assembly approach has to be used. Of course, at the very beginning, we have no idea any of the genome sequences look like, so assembly has to be done first, before a regular mapping approach can be used for further research.

**Assembly** seeks to put together the pieces of fragments that have been sequenced to recreate the genome

- The focus is on the reads, how they can be strung together
- The goal is to reconstruct the genome as best as we can
* Large search space to find the overlapping features
* Regions of similarity (aka repeats)
* Sequencing errors
* For metagenomic assemblies, it is furthur complicated by the unequal representation of members of the community, the presence of closely related microorganisms with similar genomes, the presence of several strains of the same microorganism, as well as the insufficient amount of data for minor community members.

<p align = "center">
<img src="metagenome_figures/assembly.png" alt="micribial" width="65%"/>
</p>


**Mapping** (or alignment to a reference) tries to put together the puzzle pieces directly onto an image of the picture._
- The focus is on the puzzle, regions of the puzzle that contain certain characteristics (ex. what feature) that will help you place the piece onto the puzzle.  
- In mapping the question is more, given a small chunk of sequence, where in the genome did this sequence most likely come from.
- The goal then is to find the match(es) with either the “best” edit distance (smallest difference), or all matches with edit distance less than max edit distance. Main issues are:
* Large search space
* Regions of similarity (aka repeats)
* Gaps (INDELS)
* Complexity (RNA, splicing, transcripts)


#### Aligners/Mappers
Many [alignment algorithms](https://en.wikipedia.org/wiki/List_of_sequence_alignment_software
) to choose from. Examples include:
* Aligners that can ’clip’
  * Bowtie2 in local mode
  * bwa-mem
  * minimap2
  * SNAP
* Aligners that does not allow 'clip'
  * Bowtie2 in end-to-end mode
  
#### Reference Genome

Host/contaminant genome sequence fasta files should be identified at the beginning of analysis.

* Host genome fasta files should include all primary chromosomes, unplaced sequences and un-localized sequences, as well as any organelles and alternative haplotypes.
* If you expect contamination, or the presence of additional sequence/genome, add the sequence(s) to the genome fasta file.
---

```bash
cd /share/workshop/meta_workshop/$USER/meta_example
```



