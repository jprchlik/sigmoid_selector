.. _create_json_for_web:


create_json_for_web
===================

This code creates a JSON file (sigmoid_webpage/json.txt),
which the sigmoid webpage uses to put information in each row,
and is the last code to run in the Sigmoid Catalog suite.
Primarily, the information comes from a .csv file, 
which was created by Patty Bolan in July of 2018 and updated by Jakub Prchlik in December 2018 (SigmiodCatalogAll_filament.csv). 
This .csv file is commonly used as input throughout all programs in this code suite.
If you are running the code with default parameters the call is simply:

 sswidl>create_json_for_web,'SigmoidCatalogAll_filaments.csv'

As stated previously,
running this code creates the JSON file sigmoid_webpage/json.txt.
The heavy lifting the code stems from the need to create links for all movies and plots in the JSON file.
By default it searchs in the "combined_movies/" directory (can be changed by setting the out_arch keyword)
for the directory corresponding to the sigmoid ID for a particular sigmoid.
Once it has the sigmoid ID, 
it links the disk passage HMI movie in the hmi_movie sub-directory,
the sigmoid long AIA/XRT movie in the aia_movie sub-directory,
and all flare movies in the flr_movie sub-directory using relative paths.
From the flare movies, the code counts the number of X, M, C, and B class flares.
Some sigmoids (below a sigmoid ID of 70) are missing flares due to a prolonged data outage at JSOC.
The flare movies were created by a different code, 
so see `make_aia_flare_movies <../make_aia_flare_movies/>`_ for more detailed information on the flare association.
 
All information printed into that file in the code starts with a printf,33 followed by some text.

However, create_json_for_web may take a few additional keyword arguments, which .



Technical Notes
---------------

If for some reason create_json_for_web errors before the free_lun command,
then you will need to run the following command before re-running create_json_for_web:

 sswidl>free_lun,33

This command with free the file unit used in writing json.txt.


NAME:
    create_json_for_web

INPUTS
    times      -   A csv file containing times to analyze sigmoid filaments CSV format must be as follows: formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F' readcol, times, dum, ID, NOAA, AR, AR_START, X, Y, AR_END, SIG_START, SIG_END, lifetime, TBEST, tobs, ORIENTATION, HEMISPHERE,  length_171, length_304, length, trail_length, lead_length, aspect_ratio, fwhm, height, format=formats
    hmi_arch   -  Directory containing HMI/MDI magnetic field movies and observations
    aia_arch   -  Directory containing sigmoid evolution AIA/Hinode-XRT movies
    flr_arch   -  Directory containing AIA flare movies associated with a particular sigmoid
    out_arch   -  Output directory for all movies and plots


OUTPUTS
    A JSON formatted text file for use with the XRT sigmiod webpage at 'sigmoid_webpage/json.txt'.

EXAMPLE
    create_json_for_web,times,hmi_arch='hmi_movie_cutout/',aia_arch='aia_movie/',flr_arch='aia_arch_cutout/',out_arch='combined_movies/'
