# Data and code for "Predictive microbial community changes across a temperature gradient"

This README file contains the overall structure of all folders and their contents in this repository. However, when more detailed information about data is needed, we made individual README files in the corresponding folders. 

## Information about the paper 

### Abstract: 
A central challenge in community ecology is being able to predict the effects of abiotic factors on community assembly. In particular, microbial communities play a central role in the ecosystem, but we do not understand how changing factors like temperature are going to affect community composition or function. One of the challenges is that we do not understand the mechanistic impacts of temperature on different metabolic strategies, nor how this metabolic plasticity could impact microbial interactions. Dissecting the contribution of environmental factors on microbial interactions in natural ecosystems is hindered by our understanding of microbial physiology and our ability to disentangle interactions from sequencing data. Studying the self-assembly of multiple communities in synthetic environments, here we are able to predict changes in microbial community composition based on metabolic responses of each functional group along a temperature gradient. This research highlights the importance of metabolic plasticity, and metabolic trade-offs in predicting species interactions and community dynamics across abiotic gradients. 

### Authors and contact information: 
Corresponding: Maria Rebolleda-Gomez (mreboll1@uci.edu)
Xin Sun (xinsun12@gmail.com)
Jacquelyn Folmar (jackie.folmar@yale.edu)
Ariel Favier (afavier@uci.edu)
Nora Pyenson (nora.pyenson@yale.edu)
Alvaro Sanchez (alvaro.sanchez@cnb.csic.es)


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
    2. Merge all processed files and fit growth curves with **02_growth_curves_fit.r**. This code calls the funtion to fit growth curves using generalized additive models from **growth_curve_fit_function.r**. 
    3. Fit thermal performance curves to the estimated maximum growth rate with **03_thermalperformance.r**. 
In addition, to compare the growth rates across readers we plotted together data from 30C measured with two different plate readers (**Sup_growth30.R**) to make the **Fig S4**. 

- **GrowthCurves**: Contains a folder for each temperature with the raw files from the plate reader for that temperature. 

- **Processed_data**: Folder where data gets saved after processing. First **gc_TCXX.csv** files are the processed growth files for each temperature. They each have a column for the **well** in the plate, the measured OD at 620nm (**abs**), the time (**t**), temperature at which growth was measured(**tmp**), **carbon_source** used for gorwth measurments, and **replicate**. Temperatures 37 and 42 only have one replicate. 

- Aditionally, there is a file on its own: **platemap_gc.csv** is the map of isolates in each plate. It has the well, an internal ID that allows to match each isolate with its 16S Sanger sequence (**SangerID**), the **carbon** source where the community of that isolate assembled, the **temperature** of the source community assembly, the **color** of the strain in Chromgenic media, and finally is the **Fermenter** category of the isolate (their functional group).

### GR_OA_Estrela2021
This directory has three folders:
- One with **Code** to plot organic acid secretions and growth rates from Estrela, et al., 2021 data. 

- The **Data** folder has two files from the paper's [github] (https://github.com/sylestrela/Estrelaetal2021_FunctionalAttractors). The first one (**gr.csv**) contains growth rate data (**gr_max**) for each isolate indicated y its **SangerID**. Data includes growth in many carbon sources (**cs**). We only used the glucose data. The second file (**oa.csv**) contains measurments for different variables (**variable**) including actetate and lactate secreted (the ones relevant for this work), and in **value** is stored the numeric value for this data. Variables were measured at different times, for this work we used 28hrs because it approximates our own time points. Organic acids and in mM. 

- A folder for resulting **Plots**. 

### Glucose efficiency
The code to quantify glucose use efficiency in relation to the consumption of glucose and secretion of organic acids (acetate and lactate) is **Fig6_scrip.R**. 

Relevant **Data** is in the folder as well as other processed data folders (growth rate):
    - **fermenter_flux_data.csv**: 
        - Information about the isolate ("SangerID", "well", "carbon" (from original community), "temperature" (from original community),"color"). 
        - **tmp**: temeprature of measurment
        - **Acetate**: mM concentration of acetate at 24hrs (for all temperatures but 12) or 48 (for 12C). 
        - **Lactate**: mM concentration of lactate at 24hrs (for all temperatures but 12) or 48 (for 12C). 
        - **CorrectedOD**: Measured OD with blanks removed. 
        - **CellNperML**: Cells per mL as estimated from mapping OD and CFUs for some isolates and temperatures and estimating from there the rest. 
        - **Tmatches**: Is tmp==temperature
        - **glucose**: Glucose measured at T2 (ug/uL)
        - **glucose_consumption**: average glucose in control wells - glucose measured at T2 (ug/uL)
        - **Cglucose_consumption**: Total glucose_consumption in mmol-C at T2 standardized by the theoretical mmol-C/L amount of control wells (70 mmol-C/L)
        - **Cglucose_consumption_percell**: Cglucose_consumption / (CellNperML * 1000) glucose consumption per cell at T2, expressed as mmol-C
        - **Cacid_percell**: Cacid/(CellNperML * 1000) organic acid secretion per cell at T2, expressed as mmol-C
        - **acid_glu_ratio**: Cacid/Cglucose_consumption: carbon secreted as organic acids per unit of carbon consumed as glucose
    - **fermenter_tpc_and_secretions.csv**: Is the same file plus growth rate and thermal performance parameters. 

### PredictionRF
The code **predicRF.R** has the calculations to predict R/F from our model and our isolate data. To calculate respirator biomass we used the OD data over time at the community level for the acetate community. The OD data is in the **OD_data** folder. The **carbon_map.csv** has the layout of the plates during the community assembly experiment, providing which wells were acetate communities. Finally, there is a **Plots** folder. 


