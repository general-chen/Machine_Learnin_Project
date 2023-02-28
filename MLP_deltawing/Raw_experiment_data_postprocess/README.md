################################################################################################
### This code is correct, but I did not use them for final traning of Machine Learning. This code outputs .h file, while I use pandas Dataframe as final training dataset.
I use matlab to post process the raw data of experiment.
  - individualRuns.m: This script reads the Force and Pressure data collected from Gust Onset experiments. It synchronizes the position, force, and pressure data for each run of the input test case and saves the compiled results, using smoothing.m to smooth data. Phase-averaging is not done in this script.
  - smoothFunc.m: This is a smooth function used in the main code to smooth data.

Output .h files for machine learning (this was not used in the final training of Machine learning process.)
  - shift_average_trim_combine_4test_randoms.m
  - delta_wing_32cases_split_random.txt: it concludes the split k-fold cross-validation information -- the specific cross-validation case name
  
################################################################################################

# I use this code to output .csv file and use pandas_read.csv to transform to Dataframe, and the cross-validation is split in pd.Dataframe, more intuitive, more convenient.
  - output_csv_dataframe_12_28.m: This code outputs the data to .csv, in which, I calculated the 1st and 2nd order derivative of pressure as extra features for Machine Learning, but I find it not help with the improvement of the Machine Learning model, so I did not use the derivatives (although I calculated them).
