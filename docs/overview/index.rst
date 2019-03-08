
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
From those properties, you need to gather one XRT image of each sigmoid you would like to analyze in some directory.
The files you put in the directory should be the highest resolution observations of the sigmoid near T-best.
It is possible to use `a python script <https://github.com/jprchlik/cms2_python_helpers>`_ created for CMS2 for this,
but it requires a bit of hacking. Basically, you could comment all other lines in the `download_all <https://github.com/jprchlik/cms2_python_helpers/blob/master/grab_sigmoid_fits_files.py>`_ funciton except for self.get_hinode().

Sigmoid Size Measurements
-------------------------

Next, you may run :ref:`sigmoidsize_adv_details` and specify the 
dir keyword to the directory that contains the sigmoids you are interested. sigmoidsize_adv outputs a sav file, which will clobber any save file in the same directory.
As such, you should not put too many sigmoid files in the same directory for sanity's sake. This version broke up the sigmoid's by year.
From your measurements, you will determine the following parameters in the file csv file: T-obs (and the corresponding X,Y coordinates), size, tail_length, lead_length, aspect_ratio, fwhm, and height.
Those parameters were ingested into the main csv using pandas and readsav in python.

Filament Size Measurements
-------------------------

The next things measured were the filaments observed in 171 Å and 304 Å. The :ref:`filament_selector` reads a csv file with the following IDL command\: readcol,times,ID,RATING,NOAA,AR_START,X,Y,AR_END,SIG_START,SIG_END,TBEST,format='LL,I,A,A,F,F,A,A,A,A',
where all variables were defined above, except rating were defined above. Rating was dropped from this catalog but remained as an artifact in this early program. Unlike the sigmoid measurements there is no need to split the files into year long runs.
Note that this program downloads the required full sun SDO/AIA files as needed.

Flare Association with Sigmoids
-------------------------------

:ref:`flare_cme_sigcat_csv` associates flares with sigmoidal regions based on AR number using the HEK. The program can also associated based on location,
but you should not use location alone for the flare association. I found many flares near the limb, which had poor locations, thus poor associations.
The program takes the final formatted csv file as an argument. You need to run this program before starting to download and create flare movies because
those program require the sav file output by this program. 


Downloading Flare Files
---------------------

After running the flare assocication, you may begin downloading all flares associated with a particular sigmoid. This flare correlation includes anytime the region was on disk, not only when the region is sigmoidal.
Note this flare association will only work for sigmoid's with ARs,
unless you turn on the region association in :ref:`flare_cme_sigcat_csv`. This program uses a hacked version of `rjrlib <http://www.staff.science.uu.nl/~rutte101/rridl/>`_,
which is available at the `github <https://github.com/jprchlik/sigmoid_selector>`_ page for this program suite. Note the downloading program, :ref:`get_aia_files_cutout`, 
will only download files which do not have movies created for them. This step was added because JSOC may be a bit touchy from time to time, and you do not want to redownload all 50 
flares associated with a particular sigmoid. Therefore, I suggest that if downloading fails for any reason to create flare movies, then delete the flare movies that do not 
contain any observations of the flare (discussed more below). Also note, claims to not allow more than 10 connections from a single IP address at a time; however, I can only
download 3 flares at a time before JSOC starts refusing my http requests.
This program will not download flares over the limb because rotating a given coordinate over the limb gives the coordinate -9999,-9999, which is meaningless to a JSOC request.
As such, you may notice the last flare in a few (around 5) sigmoids does not have a movie.
Final hacky note, if you computer does not have wget, then the rjrlib downloading will not work.
You can get around that by creating an alias to wget as "curl -O" in your .*rc file.

Creating Flare movies
--------------------

Use :ref:`make_aia_flare_movies` to create movies for each flare associated with the sigmoid. The program needs a slightly hacked version of aia_mkmovie to work. That hack is needed because the file 
format of the SDO/AIA files from JSOC are not what aia_mkmovie expects. Instead, the hacked aia_mkmovie reads the information from the fits headers. The program tries, and fails to match high cadence
observations of Hinode-XRT observerations. I am unsure why this failure is true because the sigmoid evolution videos do this same task successfully. You should very this program created all the movies
you wanted it to because the movies are how the sigmoid catalog counts the number of flares associated with each region. While this way of counting the flares seems nonintuitive, it is really
useful if you use the region, not AR, flare association in :ref:`flare_cme_sigcat_csv`. As such, you may delete the files that are not associated with a give sigmoid after visiual inspection.
This way of count also means it is better to have a blank movie of a flare than no movie referencing a flare.
