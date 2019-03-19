.. _make_hmi_movie_cutout:

make_hmi_movie_cutout
=====================

Creates a movie of the evolution of a sigmoid in HMI, as well as, measures properties of the sigmoid at a given cadence (36 minutes). First, the code selects files from the local archive created by get_hmi_files_cutout and finds the files that best fit between the sigmoid start and end time with a 36 minute cadence. Then it creates a UID corresponding to the AIU international standard using its position in space and time at Tobs. N.B. this cannot be used as a unique identifier for the sigmoid catalog because the X, Y coordinates at Tbest, which were used early in the catalog are inconsistent with the coordiantes of X, Y at Tobs (Tobs coordinates are better). The Tbest coordinates were human specified, while Tobs comes from the sigmoid select program's automatic region finder. At some point opaque to me the reported X, Y values went from corresponding to Tbest to Tobs without retention of the inferior Tbest X,Y coordinates. The coordinates for the UID are correspond to Tobs.

For the analysis and movie the program will only run with SDO/HMI observations with a quality flag less than 90000. The first step in the analysis is to prep the HMI observation. To do that we rotate the image by 180 degees and select 700x700 pixel region around the position at Tobs. Then the code removes spikes from the image using the HMI prep code from Antonia Savecheva. Next, I create a Gaussian smoothed absolute image from an 8x8 pixel binned image. On that smoothed image I use the IDL procedure edge_dog (Difference of Gaussians) with a radius of 1. and 700./(2*rebinv)-1 , where the 700 comes from the pixel width of the image as long as the sigmoid long axis length is less than 350 pixels. If the long axis lengths is greater than 100 pixels, the widow length is twice the measure sigmoid length plus 250 pixels. The best way I found to set the threshold value was manually unfortunately. For this there is the program make_hmi_movie_cutout_man_thres.pro, which will allow you to set the threshold value at the center of the observations dynamically. After you run that program, come back and run make_hmi_movie_cut with the /man_thres flag set. From the edges of the edge_dog procedure I create an ROI object of whatever region is nearest to the sigmoids position at Tobs accounting for rotation. That ROI objects is then used to measure the unsigned magnetic flux, positive magnetic flux, negative magnetic flux, and magnetic area under the sigmoid in the prepped image.

Once the analysis finishes, the program writes a save file of the form: save,sig_id,out_id,obs_time,obs_qual,tot_ints,pos_ints,neg_ints,pix_area,tot_area,roi_save,phy_save,filename=full_dir+'/'+str_replace(sig_id,':','')+'.sav'
and creates a movie with the ROI overplotted. Both the filenames and the output sigmoid have the unique sigmoid ID numbers.





NAME:
    make_hmi_movie_cutout

PURPOSE
    Creates hmi movie for each sigmoid

CATEGORY:
    Program, movie creation

USAGE
    make_hmi_movie_cutout,times,hmi_arch='hmi_arch/',out_arch='hmi_movie/',rebinv=16,man_thres=man_thres

INPUTS
    times      -   A csv file containing times to analyze sigmoid filaments
CSV format must be as follows:
formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F'
readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
                   length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats
    hmi_arch   -   Directory to output the HMI files to (Default = 'hmi_arch/')
    cad        -   Cadence to get HMI files in seconds (Default = 30.*60)
    out_arch   -   Where to put the output movies and save files (Default = 'hmi_movie_cutout/')
    rebinv     -   Rebinning pixel value for median smoothing computationally effciently (Default = 8)
    man_thres  -   A keyword to use the manually selected threshold value over the one currently selected by fft
    

OUTPUTS
    HMI movie and sav files in hmi_arc
