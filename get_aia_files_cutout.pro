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
;    flare_sav  -   A sav file containing flare times and positions
;    times      -   A csv file containing times to analyze sigmoid filaments
;
;OUTPUTS
;    aia files in aia_arc
;
;#############################################################
pro get_aia_files_cutout,flare_sav,times,aia_arch=aia_arch,wave=wave



;set plot to X Window
set_plot,'X'

;Read in file containing TBEST
formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F'
readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
       length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats
;Set archive directory for download aia files
if keyword_set(aia_arch) then aia_arch = aia_arch else aia_arch = 'aia_arch_cutout/'
aia_arch = aia_arch+'/'
;Set wavelength to download AIA files
if keyword_set(wave) then wave = wave else wave = ['193','304','335']

;restore save files with flare association (big_str)
restore,flare_sav


;Cadence  
cad = '45s' ;45 Seconds

email = 'jakub.prchlik@cfa.harvard.edu'
name = ''

;Download aia data for all the best times
for ii=0,n_elements(big_str)-1 do begin


    ;get index in sigmoid catalog csv file where tbest are equal
    best_ind = where(ID eq big_str[ii].sigmoid_id,match_count)
    ;get location of nearest sigmoids in catalog csv file
    ;dif_pos = fltarr(n_elements(X))
    ;for k=0,n_elements(X)-1 do begin
    ;    cat_pos = rot_xy(X[k],Y[k],tstart=tbest[k],tend=str_replace(big_str[ii].cross_m,', ','T'))
    ;    dif_pos[k] = sqrt(total((cat_pos-spos)^2))
    ;endfor

    ;;Get index of nearest sigmoid
    ;best_dis = min(dif_pos,best_ind,/Abs)

    if match_count eq 0 then continue ; No sigmoid found

    ;get good flare start times
    poss_t = where((big_str[ii].flare_s ne 'None') and (big_str[ii].flare_s ne ''),count)
 
    ;Exit loop if no flares found
    if count eq 0 then continue

    ;Create directory for output png files
    full_dir = aia_arch+strcompress(big_str[ii].sigmoid_id,/remove_all)+'/'
    if file_test(full_dir) eq 0 then file_mkdir,full_dir

    ;loop over all start times
    for ij = 0,n_elements(poss_t)-1 do begin 
 

        ;get index for a good time
        i = poss_t[ij]

        ;get time range to search over
        t1 = anytim(big_str[ii].flare_s[i])-3600. ; move forward 1 hour
        t2 = anytim(big_str[ii].flare_e[i])

        ;Get difference between start and end time in minutes and add 2 hours
        diff_t = round((t2-t1)/60.)+120

        ;get start time string
        ts = anytim(t1,/hxrbs)
        ts = '20'+ts
        ts = str_replace(ts,', ','T')
        ts = str_replace(ts,'/','-')

        ;rotate sigmoid center to ts
        inp_x = X[best_ind]
        inp_y = Y[best_ind]
        inp_t = tobs[best_ind]
        rot_p = rot_xy(inp_x,inp_y,tstart=inp_t,tend=ts)
        
        stop
        ;Get all aia data in date range with 90 minute cadence
        sdo_orderjsoc,ts,diff_t,rot_p[0],rot_p[1],email,name,wavs=wave,$
                      xsize=750,ysize=750,cadence=cad,requestidents,requestsizes
        ;Download files
        if requestsizes gt 1 then $
            sdo_getjsoc,requestidents,full_dir
        ;cd back to base directory because sdo_getjsoc goes down a level
        ;cd,'../'

    endfor

endfor

end
