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
;    get_hmi_files,times,hmi_arch='hmi_arch_cutout/'
;
;INPUTS
;    times      -   A csv file containing times to analyze sigmoid filaments
;    hmi_arch   -   A string containing the location to put the 
;                   hmi cutout archive (Default = 'hmi_arch_cutout/')
;    wave       -   A array of strings containing the wavelengths to download
;                   (Default = ['magnetogram'])
;    obs        -   IDL anytim value corresponding to when to start looking for
;                   observations (Default = anytim('2010-06-27T18:19:00') 
;   
; 
; 
;
;OUTPUTS
;    hmi files in hmi_arc
;
;#############################################################
pro get_hmi_files_cutout,times,hmi_arch=hmi_arch,wave=wave,obs=obs



;set plot to X Window
set_plot,'X'

;Read in file containing TBEST
formats = 'LL,LL,A,A,A,A,A,A,A,A,F,A,A,A,A,A,F,F,f,F,F'
readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
       length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats,/preserve
;Set archive directory for download hmi files
if keyword_set(hmi_arch) then hmi_arch = hmi_arch else hmi_arch = 'hmi_arch_cutout/'
hmi_arch = hmi_arch+'/'
;Set wavelength to download hmi files
;if keyword_set(wave) then wave = wave else wave = ['193','304','335']
if keyword_set(wave) then wave = wave else wave = ['magnetogram']
;hmi launch date 
if keyword_set(obs) then obs = obs else obs = anytim('2010-06-27T18:19:00')

;SDO/HMI take over date
sdo_takeover = anytim('2009-04-13T21:48:00')



;Cadence  
cad = '30m' ;30 minutes

email = 'jakub.prchlik@cfa.harvard.edu'
name = ''

;Good sigmoid tbest times (i.e. contains time string)
goodt = where(strlen(tbest) eq 23)



;Download hmi data for all the best times
for ii=0,n_elements(goodt)-1 do begin


    ;get index for a good time
    gi = goodt[ii]
    ;get time range to search over
    ;t1 = tbest[gi]
    ;t2 = anytim(anytim(t1)+24.,/ecs)  
    ;Switched to AR start and AR per Antonia's request 2018/11/02
    t1 = ar_start[gi]
    t2 = ar_end[gi]
    

    ;only get hmi for sigmoids after hmi launch
    if anytim(t1) lt obs then continue

    ;Use MDI if before SDO science data date
    if anytim(t1) lt sdo_takeover then begin
        wave = ['mdi_mag']
    endif else wave = ['magnetogram']
         

    ;Get difference between start and end time in minutes and add 2 hours
    diff_t = round((anytim(t2)-anytim(t1))/60.)

    ;get start time string
    ts = anytim(t1,/hxrbs)
    ts = '20'+ts
    ts = str_replace(ts,', ','T')
    ts = str_replace(ts,'/','-')

    ;rotate sigmoid center to ts
    inp_x = X[gi]
    inp_y = Y[gi]
    inp_t = tobs[gi]
    ;Correct for older sigmoids, which don't have Tobs
    if inp_t eq 0 then inp_t = tbest[gi]
    rot_p = rot_xy(inp_x,inp_y,tstart=inp_t,tend=ts,offlimb=on_limb)

    ;output directory
    full_dir = hmi_arch+strcompress(ID[gi],/remove_all)+'/'
    ;Skip already created directories for now 2018/07/30 J. Prchlik
    if file_test(full_dir) eq 1 then  continue

    ;Add check for limb. If on limb just use y value on limb
    ;Do not do limb flares
    ;if on_limb eq 1 then continue
    
    ;Get all hmi data in date range with 90 minute cadence
    sdo_orderjsoc,strmid(ts,0,16),diff_t,rot_p[0],rot_p[1],email,name,wavs=wave,$
                  xsize=750,ysize=750,cadence=cad,requestidents,requestsizes
    ;Download files
    if requestsizes gt 1 then begin 
        ;make directory if not already created
        if file_test(full_dir) eq 0 then file_mkdir,full_dir 
        ;Wait 3 minute before sending query about files
        wait,3*60
        sdo_getjsoc,requestidents,full_dir
    endif
    ;cd back to base directory because sdo_getjsoc goes down a level
    ;cd,'../'

endfor

end
