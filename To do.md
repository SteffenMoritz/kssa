# kssa 0.0.1

* Added a `NEWS.md` file to track changes to the package.

##############################
# Parameters for kssa function:

#data = target time series including real MD

#start.method = method to perform the first imputation of real MD

#segments = Number of segments of same length in which target time series will be split

#iterations = Number of iterations to simulate new missing data

#metric = metric chosen to evaluate the results (rmse, mase, smape etc)

#maximum Missing Data Siz: User controls the maximum MD size for each segment

#############################
#Homework - search for elegant functions that allow us to split time series and put simulation windows inside the segments

#####################################################################
#### First Fetch upstream in the GH web on my rep, and then pull in R
#####################################################################

#Tasks for netx meeting Feb 15
# 1. Bryan and Felipe works on the functions to select method, impute simulated missing data
  # iterate the process and produce final restults
  
# 2. Steffen works on a list of possible errors or warnings that users should get
when using the package

#Task for next meeting Feb 21
# 1. Brayan solve glitches and deliver a first version of KSS
# 2. During the meeting will test the function for different data and parameters

# Task for next meeting Feb 28
# 1. Brayan do class kssa and apply function to plot on it


# Tasks for Next Meeting Thursday March 3 2022
# 1. Create Documentation files for funcions (Felipe)
# 2. Create shortcut for "all" methods (Bryan)
# 3. Create plot.kssa (comes simple from ggpot) (Bryan)
# 4. Test the function with different time series and find erroes (Steffen, Felipe)



###################################################################
#### New function kssa2 in repositories 
###################################################################

# Commit
# How to connect kssa2 and kssa_plot as a function when we try to build the package
# Actually don't works

