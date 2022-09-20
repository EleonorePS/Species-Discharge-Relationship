# Species-Discharge-Relationship

Author: Eleonore Pierrat 

If you use this code, please cite: 
Pierrat, E., Barbarossa, V., Núñez, M., Scherer, L., Link, A., Damiani, M., … Dorber, M. (2023). Global water consumption impacts on riverine fish species richness in Life Cycle Assessment. Science of The Total Environment, 854, 158702. https://doi.org/10.1016/j.scitotenv.2022.158702

1. Input.py: set input and output working directories for the geoprocessing files
2. Module.py: utility functions to geoprocess the river basin characteristics
3. SDR_geoprocessing.py: run this code to generate the dataset used to fit the species discharge relationship. It is necessary to download all the raw data manually (sources are available in the support information file of the associated publication) and place them in the input  directory.
4. SDR.R: determination of the sdr model fitting the species richness with the characteristics of the river basin (after geoprocessing).
