.. _get_aia_files_cutout:

get_aia_files_cutout
====================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


An IDL code for downloading SDO/AIA cutout files from JSOC automatically. Note that this code will skip downloading files
if flare movies for a particular flare are already created. This features is particularly useful when JSOC disconnects during the middle of a download.



NAME:
    get_aia_files_cutout

PURPOSE
    Download range of aia files

CATEGORY:
    Program, data gathering

USAGE
    get_aia_files_cutout,times,aia_arch='aia_arch_cutout/',wave=['193','304','335'],sel_id=0

INPUTS
    flare_sav  -   A sav file containing flare times and positions
    times      -   A csv file containing times to analyze sigmoids
    aia_arch   -   The directory containing the flare files. Subdirecties exist for each sigmoid's flares by Sigmoid ID
    wave       -   Wavelengths to download for use in the flare movies
    sel_id     -   Only download files associated with a specific sigmoid ID

OUTPUTS
    aia files in aia_arch
