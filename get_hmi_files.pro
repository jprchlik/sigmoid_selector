;#############################################################
;
;NAME:
;    get_hmi_files    
;
;PURPOSE
;    Download range of hmi files
;
;CATEGORY:
;    Program, data gathering
;
;USAGE
;    get_hmi_files,times,hmi_arch='hmi_arch/'
;
;INPUTS
;    times      -   A csv file containing times to analyze sigmoid filaments.
;                   File will be read with the following call 
;                   readcol,times,ID,RATING,NOAA,AR_START,X,Y,AR_END,SIG_START,SIG_END,TBEST,format='LL,I,A,A,F,F,A,A,A,A'
;    hmi_arch   -   Directory to output the HMI files to (Default = 'hmi_arch/')
;    cad        -   Cadence to get HMI files in seconds (Default = 30.*60)
;
;OUTPUTS
;    HMI files in hmi_arc
;
;#############################################################
pro get_hmi_files,times,hmi_arch=hmi_arch,cad=cad



;set plot to X Window
set_plot,'X'

;Read in file containing TBEST
;readcol,times,dum,ID,RATING,NOAA,AR_START,X,Y,AR_END,SIG_START,SIG_END,TBEST,format='LL,I,A,A,F,F,A,A,A,A',/preserve_null
;Updated with new file format 2018/11/02
;Read in file containing TBEST
formats = 'LL,LL,A,A,A,A,A,A,A,A,F,A,A,A,A,A,F,F,f,F,F'
readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
       length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats,/preserve

;Set archive directory for download aia files
if keyword_set(hmi_arch) then hmi_arch = hmi_arch else hmi_arch = 'hmi_arch/'
hmi_arch = hmi_arch+'/'

;Set cadence for image download
if keyword_set(cad) then cad = cad else cad = 30.*60. ;30 minutes

;Good sigmoid tbest times (i.e. contains time string)
goodt = where(strlen(tbest) eq 23)


;Download HMI data for all the best times
for i=0,n_elements(goodt)-1 do begin

    ;get index for a good time
    gi = goodt[i]
    ;get time range to search over
    ;t1 = tbest[gi]
    ;t2 = anytim(anytim(t1)+24.,/ecs)  
    ;Switched to AR start and AR per Antonia's request 2018/11/02
    t1 = ar_start[gi]
    t2 = ar_end[gi]
    
    ;Get all hmi data in date range with 90 minute cadence
    s_f = vso_search(t1,t2,inst='hmi',provider='jsoc',physobs='los_magnetic_field',sample=cad)


    d_cnt = 1
    ;Download files which have the same wavelength as the requested wavelength
    if (d_cnt gt 0)  then s_r = vso_get(s_f,out_dir=hmi_arch,/FORCE)


endfor

end
