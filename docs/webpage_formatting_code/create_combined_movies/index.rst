.. _create_combined_movies:


create_combined_movies
======================

Takes the movies made by make_hmi_movie, aia_movie/make_aia_movie, and make_aia_flare_movies and combines them into one directory based on their IAU indentification number. Currently, creates new directories in combined_movies path then symbolic links to the movies in the other directories.



PROGRAM
    create_combined_movies

USAGE
    create_combined_movies,times,hmi_arch=hmi_arch,aia_arch=aia_arch,flr_arch=flr_arch,out_arch=out_arch


INPUTS
    times      -   A csv file containing times to analyze sigmoid filaments CSV format must be as follows: formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F' readcol, times, dum, ID, NOAA, AR, AR_START, X, Y, AR_END, SIG_START, SIG_END, lifetime, TBEST, tobs, ORIENTATION, HEMISPHERE, length_171, length_304, length, trail_length, lead_length, aspect_ratio, fwhm, height, format=formats
    hmi_arch   -  Directory containing HMI/MDI magnetic field movies and observations
    aia_arch   -  Directory containing sigmoid evolution AIA/Hinode-XRT movies
    flr_arch   -  Directory containing AIA flare movies associated with a particular sigmoid
    out_arch   -  Output directory for all movies and plots

EXAMPLE
    create_combined_movies,times,hmi_arch='hmi_movie_cutout/',aia_arch='aia_movie/',flr_arch='aia_arch_cutout/',out_arch='combined_movies/'

OUTPUTS
    Creates directory structure of symbolic links to make new sigmiod catalog. Combines movies from HMI, sigmoid long movies, and flare movies
