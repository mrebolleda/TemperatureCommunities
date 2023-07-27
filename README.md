# Data and code for "Predictive microbial community changes across a temperature gradient"

This README file contains the overall structure of all folders and their contents in this repository. However, when more detailed information about data is needed, we made individual README files in the corresponding folders. 

## Information about the paper 

### Abstract: 
A central challenge in community ecology is being able to predict the effects of abiotic factors on community assembly. In particular, microbial communities play a central role in the ecosystem, but we do not understand how changing factors like temperature are going to affect community composition or function. One of the challenges is that we do not understand the mechanistic impacts of temperature on different metabolic strategies, nor how this metabolic plasticity could impact microbial interactions. Dissecting the contribution of environmental factors on microbial interactions in natural ecosystems is hindered by our understanding of microbial physiology and our ability to disentangle interactions from sequencing data. Studying the self-assembly of multiple communities in synthetic environments, here we are able to predict changes in microbial community composition based on metabolic responses of each functional group along a temperature gradient. This research highlights the importance of metabolic plasticity, and metabolic trade-offs in predicting species interactions and community dynamics across abiotic gradients. 

### Authors and contact information: 
We'll be added before publication.

## Overview 
This repository contains data files and R scripts used to analyze the data for the paper. The repository is organized according to different sections in the paper and different kinds of data. The first 3 folders can be run in any order, the last 2 ones require analyses from other folders. 

1. **CommunityAssembly_16s** - Has the data files and code necessary to analyze community composition across communities assembled in 16 different sugars as the only carbon source, and five different temperatures. This corresponds **Figure 1, Figure S7** and most of the first section of results *Community composition changes with temperature in a predictable manner*. 

2. **GrowthAnalysis** - Has the data of optical density over time to measure bacterial growth and determine the thermal performance curves. Growth data is in raw formats from the different plate readers. There are three major steps in the code of this section and the code is numbered accordingly. 
    - The first part is to re-format each plate reader file for the analysis. 
    - The second is fitting growth curves, and analyzing growth data to make **Figure 2** and **Figure S5**. 
    - The third uses the growth parameters obtained in the second section and fits thermal performance curves to fill in the results for *Community assembly at higher temperatures selects isolates with higher Topt* (**Figure 3**,**Figure 4**). 

3. **GR_OA_Estrela2021** - Probably the simplest folder, contains the necessary data and code to make **Figure 5**. Data comes originally from: [Estrela et al. 2021.](https://doi.org/10.1016/j.cels.2021.09.011) and it contains growth rates and organic acid secretions for multiple bacterial strains isolated in communities assembled in glucose at 30C. 

4. **GlucoseEfficiency** - This folder contains data with organic acid secretions and glucose consumption from the fermenter bacterial strains obtained in the communities assembled for this study. Builds from growth data and parameters (GrowthAnalysis folder). It contains the code necessary for **Figure 6** and most of the section *Fermenters face a trade-off between growth and glucose efficiency at higher temperatures*.

5. **PredictionRF** - Builds from data in folders GlucoseEfficiency and CommunityAssembly_16s to calculate predicted respirator to fermenter ratios (R/F) based on our model. It contains the code necessary for section *Changes in the fermenters’ glucose efficiency predict changes in R/F across temperatures*, including **Figure 7**. 

## Detailed repository layout

### CommunityAssembly_16s
Contains five folders: 
- **Code**: Has all the R code necessary for community composition analyses. 
    - **01_AmpliconFullPipeline.Rmd**: Has a detailed pipeline of 16s sequencing analyses, from dereplicating the files to rarefaction and filtering of rarefaction and filtering of the ESV (exact sequence varuiant) table (the matrix of community composition). This code contains the steps to quality filter reads, obtain sequnece variants, remove chimeras, and map reads to taxonomic assignments. Raw sequences are available in the NCBI SRA database (Accession: PRJNA992630). 
    - **02_ESVdiversity.R**: R code with analysis to go from the ESV data obtained in the previous step to barplots of relative abundances across experimental communities, calculation of the R/F ratios (and the linear models to analyse the change in R/F with temperature), as well as calculating dispersion and beta-diversity across community samples. 

- **Data**: Contains relevant data (mostly metadata) for the analyses described above. 
    - **Fermentation_litrev.csv**: File containing the most abundant families, their classification as respirator or fermenter (T is capable of fermentation and F is only respiration), and relevant references. 
    - **TC1_amplicon_metadata.csv**: Contains information about each sequencing sample. Including: SampleID, BarcodeSequence, LinkerPrimerSequence,BarcodePlate,Well (in the original community assembly experiment),Temperature (at which the community was assembled),Carbon (carbon source used as only carbon source during assembly),CarbonID (three letter code for the carbon source),Replicate (number of replicate community from 1-4),Run (samples were sequenced in two runs), Description (short description of sample, same as SampleID). 

- **Output** and **Plots**: The amplicon full pipeline has steps to save intermediate steps (after time-intensive points) as .rds files. Those files get saved in this folder. All of the plots get stored in the Plots directory. 

- **Processed_Data**: This folder is where data processed through the analyses gets saved. These is useful data for new analyses. 
    - **ESV_fulldata_long**: This table has all the ESVs with the sequences, taxonomic assigments, and the metadata for each sample (from TC1_amplicon_metadata.csv). 
    - **ESVtable.csv**: This is the community matrix with samples as rows and ESVs as columns. The numbers in the matrix are the number of reads. 
    - **rarefied_data.csv**: Same as **ESVtable.csv** but rarefied to the minimum number of reads. 
    - **rf_data.csv**: Table with temperature of community assembly (Temp), the carbon source in 3 letter code, replicate, and R/F. 

### GrowthAnalysis
Contains four folders: 
 - **Code**: Folder with R code to process and analyse growth curve data. The analysis takes three steps: 
    1. Open and re-format all growth curves using the **01_ProcessGC_** files. There is one for each temperature (5 files). To format data tables this code calls functions from **GrowthAnalysis/Code/2023-03-06_ProcessGrowthCurves.r**
    2. Merge all processed files and fit growth curves with **02_growth_curves_fit.r**. This code calls 

