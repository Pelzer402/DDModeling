# DDModeling 0.0.2.0
* Added functionality to perform modelling using the classic drift diffusion model (see [here](https://doi.org/10.1037/0033-295X.85.2.59))
  * set 'model="DDM_classic"' in DDModel
* In addition to the above the lexical decision and repetition memory task were introduced to the package (only for DDM_classic for now)
  * set 'task="RMT_LDT"' in DDModel
* DDRep now has the ability to reshape a given DDrep between different representations
  * Specify 'ddrep' in DDRep (Note: Only reshaping inside an identical DDModel framework is allowed!)
* Added a handy function for extracting scaling using Deep Learning models (see `Scale_DL_Data()`)
* Several performance improvements

## DDModeling 0.0.1.4
* Several performance improvements (especially to GRID_Import)
* Added ability to tune SIMPLEX parameters in Fit_DDModel (see new argument simplex_coef)
* Added a Compare generic method to DDRep

## DDModeling 0.0.1.3
* Added multi thread support to several functions!
  * Defaults to k-1 threads where k is the maximum number of threads available to the system
* Added functionality for fitting with deep learning methods in Fit_DDModel
* Added functionality to easily convert GRIDs into datasets suitable for the training of neural networks using deep learning methods

## DDModeling 0.0.1.2
* Optimizations
* Added reference methods to DDFit objects:
  * plot: plots CAF and CDF distributions for a given DDFit object
  * summary: displays eta values (see [here](https://doi.org/10.3758/s13428-020-01366-8)) and booth input and fitted parameters
* Added an import function for GRIDs: Import_GRID
* Added fitting structure customization to Fit_DDModel through a new parameter 'symplex_struc'

## DDModeling 0.0.1.1
* Minor restructuring of some classes in order to increase efficiency
* Added information to several documentations
* Added functionality to the Sim_DDModel function:
  * ability to calculate multiple simulations
  * ability to initialize a simulation with manually chosen parameters

# DDModeling 0.0.1.0  
First package setup, including:
* Full support for three models: [DSTP](https://psycnet.apa.org/buy/2010-14834-002), [DMC](https://www.ncbi.nlm.nih.gov/pubmed/25909766), [SSP](https://psycnet.apa.org/record/2011-23986-003)
* Additional support for customization of each model listed above
* Ability to generate 
  * model predictions
  * model grids
  * CDF/CAF representations from data
* Ability to fit data using a combination of grid-search and downhill simplex (as described [here](https://doi.org/10.3758/s13428-020-01366-8))


