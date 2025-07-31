# Updated Genome Sequencing and Reannotation of the Broad Host Range Rhizobial Symbiont *Sinorhizobium fredii* HH103

**Francisco Fuentes-Romero¹*, Francisco-Javier López-Baena¹, José-María Vinardell¹*, Sebastián Acosta-Jurado¹**

<sup>¹ Department of Microbiology, Faculty of Biology, University of Seville, Spain</sup>  
ffuentesr@us.es, jlopez@us.es, jvinar@us.es, sacosta@us.es

---

### Correspondence

- *Francisco Fuentes-Romero (FF-R)* – ffuentes@us.es | +34 954 55 71 21  
- *José-María Vinardell (J-MV)* – jvinar@us.es | +34 954 55 71 21

---

## Summary

This repository contains the script for the updated sequencing and functional reannotation of *Sinorhizobium fredii* HH103 — a broad host range rhizobial symbiont of agronomic interest.

The updated annotation improves accuracy in coding regions and symbiosis-related genes, enabling further comparative and functional genomics.

This work aims to improve the genome annotation of *S. fredii* HH103 to facilitate the discovery of genes involved in symbiotic nitrogen fixation, host interaction, and plasmid dynamics.


---

## Genome Assembly

The PacBio reads were assembled using **Flye v2.9.3-b1797**. The assembly was performed with the following command:

```bash
flye --pacbio-raw subreads.fastq --genome-size 7m --plasmid -o output_flye -t 24 --scaffold -m 10000
```

---

## Genome Polishing

To improve base-level accuracy of the PacBio assembly, Illumina reads were used for polishing with **Pilon v1.24**.

The Illumina reads were first mapped to the PacBio assembly using **Bowtie2**, generating a sorted and indexed BAM file (`frags.bam`). The polishing step was then performed with the following command:

```bash
pilon --genome assembly_pacbio_HH103.fasta --frags frags.bam
```

---

## Repository Contents

- `scripts/` – Bash/Python scripts used for genome annotation or visualization  
- `annotation_updated/` – Final GFF3 annotation and feature tables  
- `circos_plot/` – Circos configuration and image files  
- `transposases/` – Scripts and output related to transposase gene detection  
- `README.md` – Project overview and citation

---

## Visualizations
<img src="circos_plot/HH103_circos.png" alt="Circos Plot of HH103 Genome" width="400"/>
*Figure 1. Circular representation of the HH103 genome showing key features and plasmids.*

## Accession Numbers

- **BioSample**: [`SAMN47263981`](https://www.ncbi.nlm.nih.gov/biosample/SAMN47263981)  
- **BioProject**: [`PRJNA1233244`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1233244)  
- **SRA Runs**: [`SRR34701253`](https://dataview.ncbi.nlm.nih.gov/object/SRR34701253) ; [`SRR34710307`](https://dataview.ncbi.nlm.nih.gov/object/SRR34710307)

---

## Data Availability

Data submitted to the **NCBI Sequence Read Archive (SRA)**.  
**Release date**: August 31, 2025, or upon publication, whichever is earlier.

---

## Citation

If you use this dataset, please cite as:

> Fuentes-Romero F., López-Baena F.J., Vinardell J.M., Acosta-Jurado S. (2025).  
> *Updated genome sequencing and reannotation of the broad host range rhizobial symbiont Sinorhizobium fredii* HH103.

---

## Contact

For questions or collaboration inquiries, please contact:  
**Francisco Fuentes-Romero** – ffuentes@us.es  
**José-María Vinardell** – jvinar@us.es  
University of Seville
