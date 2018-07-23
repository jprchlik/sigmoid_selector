Compliation of sigmoid catalog routines from Sean McKillop and updated by Jakub Prchlik. 
The software suite contains various programs for the XRT sigmoid catalog. 
Unfortunately, the code is a bit less a pipeline than I initially imagined due the impersisence of unique IDs (UIDs). 
I recommend the next time this endeavour is undertaken that UIDs be established from the start.

sigmoidsize_adv,dir=dir
=======================

sigmoidsize_adv performs a few measurements of sigmoids with some user intervention.
The only input to the program is a directory containing any number of sigmoid files. 
First, the program sets the scale for each observation.
Note that XRT images in the catalog are already normalized by exposure time,
but that does not reset the EXPTIME header in the fits file. 
Instead, you must use keyword e_etim, which is in microseconds. 
For the maximum value the program just finds the 99.9 percentile in the image using cgPercentiles.
Then to set the minimum scale value the program first looks for edges in the image using the edge_dog function in idl with
radius1 = 3.6 pixels and radius2 = 15. pixels for 1 threshold value.
Then use the edge_dog solution to create a mask for AR/bright regions in the image. 
If edge_dog solution does not detect any AR/bright regions the program does not apply a mask.
In the masked image we interatively solve for the median background level. 
The program stops searching for the minimum after 95% of the pixels are within 2 sigma of the median or 6 iterations.
Generally, the program finds the minimum level in less than three loops.


Next, the program will ask the user to perform a series of tasks.
First, it asks whether we should analyze the file in the command line. 
Enter 1 if we should skip the file or enter nothing if we should analyze the file.
Only enter 1 if the file does not contain a sigmiod you can analyze.
If you enter 1 the program will give null values for ever parameters and go to the next file.
If you enter nothing the program will ask you to click a number of places on the sigmoid.

Then you will select the long axis of the sigmoid. 
We define the long axis as the longest point on the sigmoid by eye using the following procedure.
First, click the trailing end of the long axis as prompted in the command line. 
Once you click the trailing end, 
you should notice a L appears where you clicked.
The command line prompt change will change to Click leading sigmoid long axis.
Again click the leading edge of the long axis where a L will appear.

After that, we then compute the short axis length of the sigmoid at the trailing and leading ends.
Similar to the long axis, the short axes are defined by inspection.
The program will ask you in order to click the southern trailing short axis, northern trailing
short axis, southern leading short axis, and northern leading short axis.
Everytime you click a new appears.

Again the program asks for more human input. This time it is to put a box around the sigmiod.
The box is important because it helps define the bounding box for the background in the
Full Width Half Maximum (FWHM) calculation, and it helps the program flare_cme_sigcat locate nearby flares.

Finally, the last of the human input. The program asks for the core of the sigmoid to determine the FWHM.
Click first point in the trailing core of the sigmiod, and the second point in the leading core. 
Those two point set the x boundary from which to calculate the FWHM.
Once you select the leading edge of the sigmoid core,
two new windows will pop up in the upper left.
The farthest to the left is an image of the sigmiod rotated by the slope of the line selected in the sigmiod core.
Highlight in the image is the region used to calculate the sigmoids FWHM. 
If you notice any irregularities due to structures not a part of the sigmoid
 reclassify the sigmiod using a new box that exclues the structure.


The calculation of the FWHM is a little tricky,
so I will described it in some detail.
The first problem is how does one compress a 2D sigmoid into a 1D
shape. 
My solution is to rotate the sigmoid about its core and sum pixels along the x-axis for a given rotated y position (you can see the idea in the rotated image).
After you sum the pixels you must now divide by the total pixels used in the sum,
since FWHM is over a small window of the image I sum a rotated mask created from the 
sigmiod selection of the line core and the bounding box. 
I then divide the total counts by the sum of the masked bounding box. 
After summation and normalization, 
I find define the peak value (stored as hght) in the rotated y direction of the sigmiod profile.
Then I break the 1D sigmoid profile into two halves. 
I then find the value closest to the half maximum for each half. 
If a half contains more than 1 value close to the half maximum I use the mean of the values to set the width. 


In order to get a consistent FWHM, 
the program needs to subtract a background flux.
The background finding algorithm uses the edge_dog AR mask of the sigmoid to mask out the bright regions and derive a first median background level.
Then the program begins a loop to further restrict the background flux.
First, it creates a new mask which includes pixels not in the AR mask and pixels less than three sigma (RMS) from the median value.
From the new mask the program computes a new median background level and standard deviation. 
The program exits with 95 percent of the "background" pixel are contained within two sigma or after 6 attempt to find the background value.
The program usually converges in 3 iterations.
Finally, to subtract from the image in physical units and derive the background for the FWHM the pixel size is converted to arcsec.


After all computation are recorded, the program does some classification. First, the program asks for the sigmiod ID number. It will show you
a best guess identification number based on previous identification, but you should verify it is correct before entering. In addtion,
it will return an AR number based on the nearest AR from the HEK +/-7days. Again you should not blindly follow the returned number and 
verify by eye especially if many AR are in the region. 


All the information is saved in an IDL save file in same directory defined in the dir keyword. The stricture of the save file has the following format:


sigdat_mod={sig_id:'',           --> User defined Sigmoid ID   
        NOAA_id:0,               --> NOAA ID of closest AR.   
        filename:'',             --> Fits file name (without directory) of the file used for Sigmoid classification   
        date:'',                 --> Date of observation from the fits header   
        size:0.0,                --> long axis size of the sigmiod in arcseconds   
        sizea:0.0,               --> trailing short axis size of the sigmoid in arcseconds.   
        sizeb:0.0,               --> leading short axis size of the sigmoid in arcseconds 
[early versions duplicated the short axis length](https://github.com/jprchlik/sigmoid_selector/commit/bab8b4ba6db4827b8896cbe12c61558b4bcafa5d#diff-3f0fc4ecafcfd45723c346c14af5de2d).   
        aspect_ratio:0.0,        --> size/((sizea+sizeb)/2.)   
        cx:0.0,                  --> Center of sigmoid in arcseconds from automatic region finding (BETA feature)   
        cy:0.0,                  --> Center of sigmoid in arcseconds from automatic region finding (BETA feature)   
        peri:0.0,                --> Perimeter distance in arcseconds from automatic region finding (BETA feature)   
        area:0.0,                --> Area in square arcseconds from automatic region finding (BETA feature)   
        roi:OBJ_NEW('IDLanROI'), --> Region object from automatic region finding (BETA feature)   
        bkgd:0.0,                --> Background flux in ADU/s/arcsec^2   
        fwhm:0.0,                --> Half width of the sigmiod core in arcseconds   
        hght:0.0,                --> Peak flux of the sigmiod in ADU/s/arcsec^2   
        bboxx:fltarr(5),         --> x-coordinates which make a bounding box for the observation in pixels   
        bboxy:fltarr(5),         --> y-coordinates which make a bounding box for the observation in pixels   
        fwlin1:fltarr(2),        --> x,y coordinates of the trailing end of the sigmiod core in pixels   
        fwlin2:fltarr(2),        --> x,y coordiantes of the leading end of the sigmiod core in pixels   
        longx1:0.0,              --> Trailing Long Axis X coordinate in Pixels       
        longy1:0.0,              --> Leading  Long Axis X coordinate in Pixels   
        longx2:0.0,              --> Trailing Long Axis Y coordinate in Pixels   
        longy2:0.0,              --> Leading  Long Axis Y coordinate in Pixels   
        shrtx1a:0.0,             --> X coordinate of Lower Trailing Sigmoid short axis in pixels       
        shrty1a:0.0,             --> X coordinate of Upper Trailing Sigmoid short axis in pixels   
        shrtx2a:0.0,             --> X coordinate of Lower Leading  Sigmoid short axis in pixels   
        shrty2a:0.0,             --> X coordinate of Upper Leading  Sigmoid short axis in pixels   
        shrtx1b:0.0,             --> Y coordinate of Lower Trailing Sigmoid short axis in pixels   
        shrty1b:0.0,             --> Y coordinate of Upper Trailing Sigmoid short axis in pixels   
        shrtx2b:0.0,             --> Y coordinate of Lower Leading  Sigmoid short axis in pixels   
        shrty2b:0.0}             --> Y coordinate of Upper Leading  Sigmoid short axis in pixels   


flare_cme_sigcat,sigloc,fname=fname,odir=odir
================
flare_cme_sigcat returns a save file with a list of CMEs and flares associated with the sigmiod. In addition,
it computes the time the sigmiod is closest to the central merdian (CM). The only required input is a directory (sigloc) with a 
save file created by sigmoidsize_adv. The output of the program is a save file is in the same directory as
 the sigloc, unless specified. 
The output file outputs 1 row for all sigmiods of a given ID in hte output from sigmoidsize_adv

The program first computes the time sigmoid is nearest to CM using the cx and cy parameters.
Then the program uses the nearest CM point to find the time the sigmoid is at CM.
To solve for the time there is a loop which rotates the cx and cy values.
The program makes a guess of the time offset to CM by using a rough number that the sun
rotates at about 10 arcsecs per hour. 
The using the the guess time it computes the new coordinates. 
If the X coordinate is within 1 arcsec of CM then the program returns that as the time for the 
sigmoid at the CM. If the distance is more than 1 arcsec then the program will try again under
two conditions.
If the new time guess is father than the original then cut down the time guess by a factor of 2
, but if the new guess is closer then the original try another rotation from the update position.
The program usually arrives at a solution after 1 or 2 iterations.   

Next, the program finds all flares associated with a given sigmoid. The program queries the flares
from the HEK for flares identified from the GOES satellites. When the GOES satellites do not return
a position for the flare the program sees if any other coordinates for a flare are available within 
1 minute, so the flares may be correlated with a particular sigmiod. The code correlates flares with
a sigmoid two different ways. First, it uses the AR association of the flare. Then it checks whether
the position of the flare is inside the box defined in sigmoidsize_adv. If either condition is met
the flare is associated with the sigmiod and stored. ~~The output file of the flare correlation has the 
form '("sigmoid_id_",I03,"_",I03,".sav")'~~.
Due to the non unique nature of the IDs the program now outputs in the form '("sigmoid_id_",I04,".sav")'.
This form is base on how the files were analyzed yearly, and not a great unique method for running a pipeline.
With no UIDs the output also now contains sigmoid start and end time to match with the total catalog. 
Currently, there are 4 sigmoids (2 sets of 2) with the same start and end times. 
Therefore, I will need a new unique characteristic of the sigmoids.



test_sig_b = CREATE_STRUCT(
'sigmoid_id',0            , --> User specified sigmiod ID     
'cross_m'   ,'           ', -->  Time flare crossed the meridian    
'sigmd_s','                    ', --> Sigmoid start time from previous sav file
'sigmd_e','                    ', --> Sigmoid end time from previous sav file
'flare_x'   ,fltarr(1000)  , --> FLARE X POSITION arcsec  
'flare_y'   ,fltarr(1000)  , --> FLARE Y POSITION arcsec  
'flare_s'   ,strarr(1000)  , --> FLARE Start time   
'flare_e'   ,strarr(1000)  , --> FLARE End time   
'flare_p'   ,strarr(1000)  , --> FLARE Peak time   
'flare_c'   ,strarr(1000)  , --> FLARE GOES Class   
'cme_x'     ,fltarr(100)  , --> CME X POSITION   
'cme_y'     ,fltarr(100)  , --> CME Y POSITION   
'cme_s'     ,strarr(100)  , --> CME Start time   
'cme_e'     ,strarr(100)  , --> CME End time   
'cme_w'     ,strarr(100)  , --> CME Angular Width   
'cme_v'     ,strarr(100))     --> CME Velocity   




getflux_mod
================



rjrlib
============
Zipped IDL library containing the functions sdo_orderjsoc.pro and sdo_getjsoc.pro, which we exploit to get high cadence cutout images for the catalog. You will need to add these programs to your IDL path.
Full documentation for see http://www.staff.science.uu.nl/~rutte101/rridl/

