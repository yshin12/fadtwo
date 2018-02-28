# fadtwo (FActor-Driven TWO-regime regressions)

Last update: 2018-02-28

There are two codes:

  (1) lib_fadtwo.R - the collection of all functions required for estimating the model 
                     model selection is added (2018-02-22)
  
  (2) dat_rz.dat - a data file for Ramey and Zubeira (2017)
  
  (3) rz_fadtwo.R - for estimating the RZ model without model selection
  
  (4) rz_selection_joint.R - for estimating the RZ model with model selection; joint method
  
  (5) rz_selection_iter.R - for estimating the RZ model with model selection; iterative method
  
  (6) rz_dat_original.csv - The original data file from RZ paper; required to generate an updated data for the STATA usage later
  
  (7) Stata/jordagk_llss.do - a STATA file slightly modified from the RZ paper; it uses rz_dat_updated.csv to reflect different economic status from LLSS estimation.
  
  (8) Stata/rz_dat_updated.csv - an updated data file; generated from rz_fadtwo.R file.
  
  
