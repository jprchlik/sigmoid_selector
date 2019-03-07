Overview
============

In this section I will give a quick overview of how to use the functions in the sigmoid catalog and in what order they work.

The beginning
-------------
First, you need to identify a list of sigmoids you are interested in analyzing and put that information in a csv file.
The csv file should contain at least the following information in the same order.
That is the first parameter listed here should be the first column in the csv file.
Below are a list of columns and a brief explanation of their meaning and how to format them.


    #. ID --> The unique identification (UID) of the sigmoid. The current version of the code suite forces    
    this to be an integer, but you could easily modify the code to suit you.    
    Keep for simplicity keep set this UID early and keep it constant throughout your analysis.    
    #. NOAA --> The NOAA ID for a corresponding active region (AR) if applicable. This value must also be an integer or empty.    
    #. AR   --> The same as NOAA except it is an integer. A value of 0 means there is no corresponding AR with the sigmoid.    
    #. AR_START --> The start time of the magnetic fields in the region, which is not necessarily the AR start time.    
    #. X --> The x-coordinate in arcseconds from disk center at T_best or T_obs. Unfortunately, the coordinates at T_best were    
    not very good in the initial version of the catalog, so for sigmoid IDs greater than 30 the X references T_obs.    
    #. Y --> The y-coordinate in arcseconds from the disk at T_best or T_obs. Check the description of X for more infromation.    
