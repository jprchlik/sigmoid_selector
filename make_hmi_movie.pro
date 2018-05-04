;#############################################################
;
;NAME:
;    make_hmi_movie
;
;PURPOSE
;    Creates hmi movie for each sigmoid
;
;CATEGORY:
;    Program, movie creation
;
;USAGE
;    filament_selector,times,hmi_arch='hmi_arch/'
;
;INPUTS
;    times      -   A csv file containing times to analyze sigmoid filaments
;
;OUTPUTS
;    HMI files in hmi_arc
;
;#############################################################

pro make_hmi_movie,times,hmi_arch=hmi_arch
;set plot to X Window
set_plot,'X'

;Read in file containing TBEST
readcol,times,ID,RATING,NOAA,AR_START,X,Y,AR_END,SIG_START,SIG_END,TBEST,format='LL,I,A,A,F,F,A,A,A,A'
;Set archive directory for download aia files
if keyword_set(hmi_arch) then hmi_arch = hmi_arch else hmi_arch = 'hmi_arch/'
hmi_arch = hmi_arch+'/'

;Get list of hmi files
hmi_list = file_search(hmi_arch+"hmi*fits")


;separate file name
fnames = strsplit(hmi_list,'/',/extract) 
fname = strarr(n_elements(fnames))
;Number of columns
splnu = n_elements(fnames[0])

;Get just filenames
for i=0,n_elements(fnames)-1 do fname[i] = fnames[i,splnu-1]

;Remove filename start
times = str_replace(fname,'hmi.m_45s.','')
;Remove fits ending
times = str_replace(times,'_TAI.magnetogram.fits','')
;Switch _ to :
times = str_replace(times,'_',':')

;Replace the first : with a T
pos = strpos(times,':')
strput,times,'T',pos[0]
;Replace the first : with a .
;pos = strpos(times,':',/Reverse_search)
;strput,times,'.',pos[0]


;Get hmi times
hmi_time = double(anytim(times))


;good sigmoid tbest times (i.e. contains time string)
goodt = where(strlen(tbest) eq 23)

;Cadance for image creation in seconds
img_cad = 90.*60.

;Download HMI data for all the best times
for i=0,n_elements(goodt)-1 do begin
    ;get index for a good time
    gi = goodt[i]
    xi = x[i]
    yi = y[i]

    ;get time range to search over
    ;t1 = tbest[gi]
    ;t2 = anytim(anytim(t1)+24.,/ecs)  
    t1 = sig_start[i]
    t2 = sig_end[i]
  
    ;Convert to anytimes
    at1 = anytim(t1)
    at2 = anytim(t2)

    ;Create array of good times at 90 minute cadence
    add_cad  = 1
    counter = 1
    time_arr = [at1]
    
    ;loop until cadence maximum is reached
    while add_cad eq 1 do begin 

        ;Add new cadence to time array 
        new_time = double(counter*img_cad+at1)
        time_arr = [time_arr,new_time]


        ;end if new time greater than endtime
        if new_time gt at2 then add_cad = 0
        counter = counter+1
    endwhile



endfor

end