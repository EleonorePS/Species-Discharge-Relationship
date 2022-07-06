# Species-Discharge-Relationship

Author: Eleonore Pierrat 

If you use this code, please cite: 
Global modelling of water consumption impacts on freshwater fish biodiversity in LCA, submitted

1. Input.py: set input and output working directories for the geoprocessing files
2. Module.py: utility functions to geoprocess the river basin characteristics
3. SDR_geoprocessing.py: run this code to generate the dataset used to fit the species discharge relationship
    it is necessary to download all the raw data (sources are available in the support information file of the associated publication) and place them in the input  directory.
4. SDR.R: determination of the sdr model fitting the species richness with the characteristics of the river basin (after geoprocessing).
