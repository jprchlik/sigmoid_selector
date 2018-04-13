Compliation of sigmoid catalog routines from Sean McKillop and updated by Jakub Prchlik

#sigmoidsize_adv,dir=dir
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
two new window will pop up in the upper left.



#flare_cme_sigcat



#getflux_mod
