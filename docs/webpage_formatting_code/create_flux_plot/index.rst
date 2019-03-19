.. _create_flux_plot:

create_flux_plot
================


Name:
    create_flux_plot

Program
    Creates a flux versuses time plot for the magnetic field under a sigmoid

Usage
    create_flux_plot,filein,outdir,s_stim,e_stim,f_stim=f_stim,f_scls=f_scls,fileout=fileout

Input
    filein -- A save file created by make_hmi_movie.pro
    outdir -- Output directory for the plot
    s_stim -- Start time of observed sigmiod 
    e_stim -- End time of observed sigmiod 
    f_stim -- Time of flare in sigmiod as a string that may pass to anytim
    f_scls -- Class of flare in sigmiod

Output 
    fileout -- Full path to the created file
