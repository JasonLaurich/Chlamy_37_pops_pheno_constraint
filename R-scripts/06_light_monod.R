# Jason R Laurich
# January 7, 2025

# This script will estimate r for each light condition in the same way that I estimated it for the TPC analysis (04_final_r_estimation)
# I'm going to use a sliding-window approach to identify the exponential phase of the logged linear data, then fit an exponential growth
# curve to un-logged data for that time period for each light level.

# Then I'm going to fit Monod curves to those data (using R2jags) for each population to estimate R* (I*)

############# Packages ########################


############# Upload and examine data #######################


############# Data exploration #####################



############# Loop through all populations ###################


############# Fit Monod curves to data ###################

