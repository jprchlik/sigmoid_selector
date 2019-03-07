.. _filament_selector

filament_selector
=================


Program allows a user to manually select the a filament in the core of the sigmoid. To select the filament you need only use the mouse. Points along the filament are selected with a left click. A right click stores that click value and end clicking along the filament. A middle click recenters the observation in case of a bad Tbest X,Y value, which is rather frequent (N.B. Tbest != Tobs, Tobs positions are much better, but this code was developed before Tobs was readily in use). The middle click also reset all clicked points. An immediate right click means there is no filament in the observations, so the code stores a -9999 value.

The output save file has the following format:
fil_d={sig_id:'', --> sigmoid ID 
NOAA_id:0, --> NOAA AR ID
filename:'', --> AIA filename used for the analysis
date:'', --> Date of AIA observation
leng:0.0, --> Length of filament in arcsec
device_arcsecx:0.0, --> Conversion from clicked points to arcsec in X
device_arcsecy:0.0, --> Conversion from clicked points to arcsec in Y
devicex:fltarr(100), --> Clicked points in X
devicey:fltarr(100)} --> Clicked points in Y