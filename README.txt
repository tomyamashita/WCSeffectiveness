This README.txt file was generated on 2024-11-19 by Thomas J. Yamashita

GENERAL INFORMATION

1. Title of Dataset: Predicting species assemblages at wildlife crossing structures using multivariate regression of principal coordinates
	This data is associated with the manuscript, titled Predicting species assemblages at wildlife crossing structures using multivariate regression of principal coordinates, available here: XXXX

2. Author Information
	Thomas J. Yamashita
		Caesar Kleberg Wildlife Research Institute, Texas A&M University - Kingsville
		tjyamashta@gmail.com
		Corresponding Author
	Daniel G. Scognamillo
		Caesar Kleberg Wildlife Research Institute, Texas A&M University - Kingsville
	Kevin W. Ryer
		School of Earth, Environmental, and Marine Sciences, University of Texas Rio Grande Valley
	Richard J. Kline
		School of Earth, Environmental, and Marine Sciences, University of Texas Rio Grande Valley
	Michael E. Tewes
		Caesar Kleberg Wildlife Research Institute, Texas A&M University - Kingsville
	John H. Young Jr. 
		Environmental Affairs Division, Texas Department of Transportation
	Jason V. Lombardi
		Caesar Kleberg Wildlife Research Institute, Texas A&M University - Kingsville

3. Date of Data Collection
	SH 100: September 2016 - May 2019
	FM 106: July 2019 - January 2021
	FM 1847: January 2022 - August 2023

4. Geographic location of data collection: wildlife crossing structures located along SH 100, FM 106, and FM 1847 in Cameron County, Texas, USA

5. Funding Sources: Texas Department of Transportation


DATA & FILE OVERVIEW

1. File List
	CTtable_All_20240301.xlsx: CT Table describing active and total camera trap nights for the duration of the study period.
	EnvData_site_20240419.xlsx: File containing information about experimental units and vegetation data for each site. 
	Species_List_20240108.xlsx: List of all species, including scientific names detected on cameras and those that were included in the analyses. 
	varpart_20240514.xlsx: Primer output of dbRDA analyses used to conduct variation partitioning analysis in R. 
	WCS_effectiveness_Manuscript.R: R script with code used in the analyses

2. Relationship between files: 
	The camera names are the same in the above files and can be used to link the data. 
	The varpart_20240514.xlsx file is separate and provides the Primer output used to conduct variation partitioning analysis.

3. Files not included in this repository. 
	Camera data from this study. This data is available upon request. Please contact Thomas Yamashita


METHODOLOGICAL INFORMATION

1. Description of the methods used for collection/generation and processing of data: 
	Methodology for collection and processing of the data can be found in the manuscript

2. Quality Assurance Procedures: 
	Photographs were processed by trained personnel and went through multiple rounds of checks to ensure accurate detection of all animals 

3. People involved with data collection, processing, and analysis: 
	Thomas J. Yamashita (data collection, data processing, data analysis)
	Kevin W. Ryer (data collection, data processing)
	Daniel G. Scognamillo (data collection)
	Jason V. Lombardi (data collection, data analysis)


DATA SPECIFIC INFORMATION FOR: CTtable_All_20240301.xslx

1. Data Type: Microsoft excel xlsx file

2. Number of Tabs: 3

3. Number of Variables: 28, 16, 34

4. Number of Rows: 101, 80, 63

5. Variable List: 
	ConPeriod: Construction period (late construction or post-construction)
	Camera: Camera Name
	Station: Wildlife crossing side of highway
	Site: Wildlife crossing number
	Session: Numeric construction period
	Highway: Highway
	utm_y: Y coordinate in UTM zone 14N for the camera. This information is redacted and is available upon request.
	utm_x: X coordinate in UTM zone 14N for the camera. This information is redacted and is available upon request.
	Setup_date: Date camera was set up
	Retrieval_date: Date camera was retrieved
	ProblemX_from: (multiple columns): Start date that a camera was inactive. X starts at 1 and increases by 1. Higher numbers indicate that a camera had more problems
	ProblemX_to: (multiple columns): End date that a camera was inactive. Just like above, X starts at 1 and increases by 1 as a camera has more problems


DATA SPECIFIC INFORMATION FOR: EnvData_site_20240419.xlsx

1. Data Type: Microsoft excel xlsx file

2. Number of Tabs: 4

3. Number of Variables: 30, 6, 12, 5

4. Number of Rows: 18, 84, 347, 45

5. Variable List: 
	A full description of the variables in this dataset is provided in the "Metadata" tab in this excel file


DATA SPECIFIC INFORMATION FOR: Species_List_20240108.xlsx

1. Data Type: Microsoft excel xlsx file

2. Number of Tabs: 1

3. Number of Variables: 12

4. Number of Rows: 70

5. Variable List: 
	Species: common name identifier (underscores are used instead of spaces)
	Timelapse: common name used in the Timelapse2 program for classifying species
	Name: full common name of the species
	Class: Taxonomic class
	Order: Taxonomic order
	Family: Taxonomic family
	ScientificName: Full scientific name for the species
	Tag: Identifier for native, exotic, or invasive species. Also included in this are names that represent multiple species (group)
	Individuals: Was the number of individuals calculated for this species
	Interactions: Was interaction type (see manuscript) computed for this species
	Process: Was number of independent events computed for this species
	Analysis: Was this species included in the analysis


DATA SPECIFIC INFORMATION FOR: varpart_20240514.xlsx

1. Data Type: Microsoft excel xlsx file

2. Number of Tabs: 4

3. Number of Variables: 45, 7, 7, 5

4. Number of Rows: 33, 93, 93, 6

5. Variable List: 
	TAB: FiveParts
		This tab shows the workflow for which combination of dbRDAs are necessary for calculating the interactive effects and percent of variation explained by each set of variables
	TAB: Results_cor/Results_nocor:
		Response: Response variable for the dbRDA
		AdjustedR2:  Adjusted R^2 value from the dbRDA for that model
		R2: R^2 value for the dbRDA model
		RSS: Residual sums of squares for the dbRDA model
		No.Groups: Number of sets of predictors used in the dbRDA model
		Selections: Which set of predictors was used in the dbRDA model (see Notes Tab for details)
		RDA: Which model is this. Used to determine the subtractions needed to calculate variation explained
	TAB: Notes
		Labels identifying which set of predictors each number is associated with. 

6. Other Information: 
	The Results_cor nd Results_nocor tabs represent the variation partitioning dbRDAs for accounting for spatial and temporal autocorrelation using dbMEMs (accounting for autocorrelation [Results_cor]) or raw coordinates/months (not accounting for autocorrelation [Results_nocor]).
	Because 5 sets of variables were used in variation partitioning, a custom function needed to be written. See the supplemental information of the manuscript for details and the code. 


DATA SPECIFIC INFORMATION FOR: WCS_effectiveness_Manuscript.R

1. Data Type: R script

5. Other Information: This script provides code for conducting all R-based analyses used in this manuscript. This code will not run properly without the camera data which is not provided in this repository. Camera data is available upon request from Thomas Yamashita.



