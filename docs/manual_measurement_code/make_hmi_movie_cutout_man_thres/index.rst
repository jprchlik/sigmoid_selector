.. _make_hmi_movie_cutout_man_thres:

make_hmi_movie_cutout_man_thres
===============

This code has the same reduction procedure as :ref:`make_hmi_movie_cutout`, 
so I will not redescribed the reduction procedure here.
What this program adds is the ability to dynamically set a threshold detection limit for feature detection.
Generally, a good feature detection threshold is between 0.1 and 0.5, 
but this rule is not hardfast.
Once you confirm a threshold value for a particular sigmoid ID, 
then the program will write the threshold value to manual_edge_dog_threshold.txt.
You should guess a low number in the program first because if the number is too large and cannot detect any features it will error out.

If the program does error out do not worry,
all previous threshold values are already updated in manual_edge_dog_threshold.txt.
As such, the program will restart with the first sigmoid with a threshold value below 0 and skip any subsequent values above 0.
You may exploit this feature if after running the analysis you do not like the threshold value you first selected.
To reset a threshold simply make a threshold value in manual_edge_dog_threshold.txt less than 0 for a sigmoid you want to reanalyze.





NAME:
    make_hmi_movie_cutout_man_thres

PURPOSE
    Creates hmi movie for each sigmoid

CATEGORY:
    Program, sets threshold value

USAGE
    make_hmi_movie_cutout_man_thres,times,hmi_arch='hmi_arch/',rebinv=16

INPUTS
    times      -   A csv file containing times to analyze sigmoid filaments CSV format must be as follows: formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F' readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST, tobs, ORIENTATION, HEMISPHERE, length_171, length_304, length, trail_length, lead_length, aspect_ratio, fwhm, height, format=formats
    hmi_arch   -   Directory to output the HMI files to (Default = 'hmi_arch/')
    cad        -   Cadence to get HMI files in seconds (Default = 30.*60)
    rebinv     -   Rebinning pixel value for median smoothing computationally effciently (Default = 8)

OUTPUTS
    A save file called'manual_edge_dog_threshold.sav' which contains two variables ID and store_thres. ID is the sigmoid ID, while store_thres are the the manually determinated threshold values for each ID. The element number in each array corresponds to the same sigmoid. 
