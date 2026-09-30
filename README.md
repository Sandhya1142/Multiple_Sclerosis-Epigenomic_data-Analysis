**1. OVERVIEW**

An R workflow for analysing DNA methylation data from GEO dataset GSE224455 (Illumina Infinium EPIC array). The pipeline compares brain lesion tissue against normal-appearing white matter (NAWM) and identifies differentially methylated positions (DMPs) using the ChAMP package and limma  

**2. RESEARCH QUESTION ADDRESSED**

Does DNA methylation have a role to play in molecular changes associated with MS Lesions? Are CpGs differentially methylated in Lesion white matter as compared to Normal white matter in MS patients? 

**3. DATASET**

The following dataset has been taken from the GEO dataset (Gene Expression Omnibus) from 9 patients from MS patients consisting of chronically demyelinated MS lesions and matched normal-appearing white matter (NAWM) samples.


**4. ANALYSIS WORKFLOW**
<img width="1040" height="720" alt="Workflow_224455" src="https://github.com/user-attachments/assets/4fcda310-3f18-482a-9f2b-409906d97dae" />

Key Functions used:
1. getGEO(), getGEOSUppFiles():      Download the GEO metadata, raw IDAT files

2. champ.load():                     Loading the IDAT files and detecion p-value filtering

3. champ.QC():                       Quality Control

4. champ.norm()                      Normalization using BMIQ method

5. champ.SVD():                      To inspect the sources of variation

6. champ.runCombat():                To overcome batch effects

7. duplicateCorrelation(), lmFit(), eBayes:  DMP analysis for paired samples

8. probe.features:                   Annotation of DMPs


**5. RESULTS**


1. Total number of probes in the experiment: 8,65,918

2. Number of CpG probes remained after filtering: 699018

3. No. of statistically significant DMP observed (padj < 0.05): 1,32,461

4. Delta_Beta condition: Brain_lesion - Brain_Normal. This implies that a positive DeltaBeta value means higher methylation in Lesion samples and negative value means a higher methylation in Normal white matter samples.

5. No. of probes identified as DMP and with delta_beta > 0 : 62,860 

6. No. of probes identified as DMP and with delta_beta < 0 : 69,601

7. Avergae value of DeltaBeta of DMP values is 0.0040.

8. Median Value of DeltaBeta of DMP values is -0.0152



**Inference**: Among the identified DMPs, 69,601 (52.5%) are hypermethylated in NAWM relative to the lesion, whereas 62,860 (47.5%) are hypermethylated in the lesion relative to NAWM.
The Median value is slightly negative indicating a tendency toward higher methylation in NAWM. However, the average value is close to zero (+0.0040), suggesting that there is no strong overall directional difference in methylation across the DMPs.


9. Annotation results: Number of CpGs beloning to:

- island : 6792
- shore(2Kb) : 14887
- shelf(4Kb) : 6400
- opensea(>4Kb) : 28577
 




**6. TOOLS AND SKILLSET ACQUIRED**
  

DNA methylation analysis 

Epigenomics

R

Bioconductor

ChAMP

Limma

Data Visualization

GitHub
