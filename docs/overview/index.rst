
Overview
============

In this section I will give a quick overview of how to use the functions in the sigmoid catalog and in what order they work.

The csv File 
-------------
First, you need to identify a list of sigmoids you are interested in analyzing and put that information in a csv file.
The csv file should contain at least the following information in the same order.
That is the first parameter listed here should be the first column in the csv file.
Below are a list of columns and a brief explanation of their meaning and how to format them.


    0. UNNAMED --> A dummy index. It is an artifact of exporting the table usings pandas and has no other meaning.

    #. ID --> The unique identification (UID) of the sigmoid. The current version of the code suite forces this to be an integer, but you could easily modify the code to suit you. Keep for simplicity keep set this UID early and keep it constant throughout your analysis.    

    #. NOAA --> The NOAA ID for a corresponding active region (AR) if applicable. This value must also be an integer or empty.    

    #. AR   --> The same as NOAA except it is an integer. A value of 0 means there is no corresponding AR with the sigmoid.    

    #. AR_START --> The start time of the magnetic fields in the region, which is not necessarily the AR start time. The time must be in YYYY-MM-DDTHH:MM:SS.MMM format.

    #. X --> The x-coordinate in arcseconds from disk center at T_best or T_obs. Unfortunately, the coordinates at T_best were not very good in the initial version of the catalog, so for sigmoid IDs greater than 30 the X references T_obs.    

    #. Y --> The y-coordinate in arcseconds from the disk at T_best or T_obs. Check the description of X for more infromation.    

    #. AR_START --> The end time of the magnetic fields in the region, which is not necessarily the AR end time. The time must be in YYYY-MM-DDTHH:MM:SS.MMM format.

    #. SIG_START --> The start time of the sigmoid as seen in Hinode XRT syntopic images. The time must be in YYYY-MM-DDTHH:MM:SS.MMM format.

    #. SIG_END --> The end time of the sigmoid as seen in Hinode XRT syntopic images. The time must be in YYYY-MM-DDTHH:MM:SS.MMM format.

    #. LIFETIME --> The lifetime of the sigmoid in as seen in XRT in days.

    #. TBEST --> The time when the sigmoid looked most sigmoidal in Hinode/XRT observations.
    
    #. TOBS --> The time used to analyze the sigmoid's Full Width at Half Maximum (FWHM).

    #. ORIENTATION --> How the loops look. S = S-Shaped, I = Inverted S-Shaped and H = Horizontal.

    #. HEMISPHERE --> Which hemisphere the sigmoid is located in. N = North and S = South.
  
    #. 171_length --> The length of a 171 Å filament in arcsec during T-Best. If no filament was located the column contains a value of -9999.

    #. 304_length --> The length of a 304 Å filament in arcsec during T-Best. If no filament was located the column contains a value of -9999.

    #. size --> Measurement of the long axis of the sigmoid in arcsec.

    #. tail_length --> The length of the trailing edge of the sigmoid's short axis (arcsec).

    #. lead_length --> The length of the leading edge of the sigmoid's short axis (arcsec).

    #. aspect_ratio -->  Ratio of the long axis to the average short axis.

    #. fwhm --> The full width at half maximum of the sigmoid's core at T-obs (arcsec).
 
    #. height --> The peak of the sigmoid's core at T-obs (DN/s/arcsec^2).

    #. ha_filament --> Is there an Ha filament present during T-Best? Y = Yes, N = No, and ? = Undetermined

    #. ss_at_peak --> The number of sunspots at T-best.

    #. good_mag --> Now meaningless as all magnetic fields are okay.

    #. IAU_ID --> The Internation Astronomical Union (IAU) unique identifier (UID) of the sigmoid at T-best. Not really every used because the X,Y defintion was changed from referencing T-best to T-obs. 

    #. filament_eruption --> Did a filament erupt while the region showed a sigmoidal shape (Y = Yes, N = No, and M = Maybe).

    #. transient_CH --> If there was a Transient Coronal Hole observed during sigmoid lifetime (Y = Yes and N = No, and M = Maybe).

    #. flare_ribbons --> Were there Flare Ribbons observed after a flare? (Y = Yes and N = No).

    #. postflare_loops --> Were there Post-Flare Loops observed after a flare? (Y = Yes and N = No.)

    #. nearby_CH --> If there was a Coronal Hole observed next to the region. (Y = Yes, N = No, and M = Maybe.)



Synoptic Observations
-------------

The first step in this pipeline is to figure out some rough properties of the sigmoid. These rough properties include ID, X, Y, T-best, SIG_START, SIG_END.
From those properties, you need to gather one XRT image of each sigmoid you would like to analyze in some directory. Then you may run sigmoidsize_adv :ref:`sigmoidsize_adv_details` and specify the 
dir keyword to the directory that contains the sigmoids you are interested. sigmoidsize_adv outputs a sav file, which will clobber any save file in the same directory.
As such, you should not put too many sigmoid files in the same directory for sanity's sake. This version broke up the sigmoid's by year.
