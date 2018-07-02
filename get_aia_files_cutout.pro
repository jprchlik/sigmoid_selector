;#############################################################
;
;NAME:
;    get_aia_files    
;
;PURPOSE
;    Download range of aia files
;
;CATEGORY:
;    Program, data gathering
;
;USAGE
;    get_aia_files,times,aia_arch='aia_arch/'
;
;INPUTS
;    times      -   A csv file containing times to analyze sigmoid filaments
;
;OUTPUTS
;    aia files in aia_arc
;
;#############################################################
pro get_aia_files_cutout,flare_sav,aia_arch=aia_arch,wave=wave



;set plot to X Window
set_plot,'X'

;Read in file containing TBEST
;;formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F'
;;readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
;;       length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats
;Set archive directory for download aia files
if keyword_set(aia_arch) then aia_arch = aia_arch else aia_arch = 'aia_arch_cutout/'
aia_arch = aia_arch+'/'
;Set wavelength to download AIA files
if keyword_set(wave) then wave = wave else wave = ['193','304','335']

;restore save files with flare association
restore,flare_sav


;Cadence  
cad = 0.5 ;30 Seconds


;Download aia data for all the best times
for ii=0,n_elements(goodt)-1 do begin

    ;get index for a good time
    i = goodt[ii]
    ;get time range to search over
    ;t1 = tbest[gi]
    ;t2 = anytim(anytim(t1)+24.,/ecs)  
    t1 = sig_start[i]
    t2 = sig_end[i]

    ;Get difference between start and end time in minutes
    diff_t = round((anytim(t2)-anytim(t1))/60.)
    
    ;Get all aia data in date range with 90 minute cadence
    sdo_orderjsoc,t1,diff_t,


    d_cnt = 1
    ;Download files which have the same wavelength as the requested wavelength
    if (d_cnt gt 0)  then s_r = vso_get(s_f,out_dir=aia_arch,/FORCE)


endfor

end
