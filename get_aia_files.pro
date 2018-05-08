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
pro get_aia_files,times,aia_arch=aia_arch,wave=wave



;set plot to X Window
set_plot,'X'

;Read in file containing TBEST
readcol,times,ID,RATING,NOAA,AR_START,X,Y,AR_END,SIG_START,SIG_END,TBEST,format='LL,I,A,A,F,F,A,A,A,A'
;Set archive directory for download aia files
if keyword_set(aia_arch) then aia_arch = aia_arch else aia_arch = 'aia_arch/'
aia_arch = aia_arch+'/'
;Set wavelength to download AIA files
if keyword_set(wave) then wave = wave else wave = ['193','304','335']


;Good sigmoid tbest times (i.e. contains time string)
goodt = where(strlen(tbest) eq 23)

;Cadence  
cad = 30.*60. ;30 Minutes


;Download aia data for all the best times
for i=0,n_elements(goodt)-1 do begin

    ;get index for a good time
    gi = goodt[i]
    ;get time range to search over
    ;t1 = tbest[gi]
    ;t2 = anytim(anytim(t1)+24.,/ecs)  
    t1 = sig_start[i]
    t2 = sig_end[i]
    
    ;Get all aia data in date range with 90 minute cadence
    s_f = vso_search(t1,t2,inst='aia',provider='jsoc',physobs='intensity',sample=cad,wave=wave)


    d_cnt = 1
    ;Download files which have the same wavelength as the requested wavelength
    if (d_cnt gt 0)  then s_r = vso_get(s_f,out_dir=aia_arch,/FORCE)


endfor

end
