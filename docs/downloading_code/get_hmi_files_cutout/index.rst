.. _get_hmi_files_cutout:

get_hmi_files_cutout
====================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


An IDL code for downloading HMI cutout files from JSOC automatically.
This code needs to run before running make_hmi_movie_cutout.


NAME:
    get_hmi_files_cutout    

PURPOSE
    Download range of hmi files

CATEGORY:
    Program, data gathering

USAGE
    get_hmi_files_cutout,times,hmi_arch='hmi_arch_cutout/'

INPUTS
    times      -   A csv file containing times to analyze sigmoids
    hmi_arch   -   A string containing the location to put the hmi cutout archive (Default = 'hmi_arch_cutout/')
    wave       -   A array of strings containing the wavelengths to download (Default = ['magnetogram'])
    obs        -   IDL anytim value corresponding to when to start looking for HMI instead of MDI observations (Default = anytim('2010-06-27T18:19:00') 
   
 
 

OUTPUTS
    hmi files in hmi_arc
