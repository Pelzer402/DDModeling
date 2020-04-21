## DDModeling 0.0.1.4
* Several performace improvements (especially to GRID_Import)
* Added ability to tune SIMPLEX parameters in Fit_DDModel (see new argument simplex_coef)
* Added a Compare generic method to DDRep

## DDModeling 0.0.1.3
* Added functionality for fitting with deep learning networks in Fit_DDModel

## DDModeling 0.0.1.2
* Optimization
* Added reference methods to DDFit objects:
  * plot: plots CAF and CDF distributions for a given DDFit object
  * summary: displays eta values (see [here](https://doi.org/10.3758/s13428-020-01366-8)) and booth input and fitted parameters
* Added an import function for GRIDs: Import_GRID
* Added fitting structure customization to Fit_DDModel through a new parameter 'symplex_struc'

## DDModeling 0.0.1.1
* Minor restructuring of some classes in order to increase efficiency
* Added a parameter slot to the DDRep-class
* Added information to several documentations
* Added functionality to the Sim_DDModel function:
  * ability to calculate multiple simulations
  * ability to initialize a simulation with manually choosen parameters

# DDModeling 0.0.1.0  
First package setup, including:
* Full support for three models: [DSTP](https://psycnet.apa.org/buy/2010-14834-002), [DMC](https://www.ncbi.nlm.nih.gov/pubmed/25909766), [SSP](https://psycnet.apa.org/record/2011-23986-003)
* Additional support for customization of each model listed above
* Ability to generate 
  * model predictions
  * model grids
  * CDF/CAF representations from data
* Ability to fit data using a combination of grid-search and downhill simplex (as described [here](https://doi.org/10.3758/s13428-020-01366-8))


