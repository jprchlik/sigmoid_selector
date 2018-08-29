;###########################################################################
;PROGRAM
;
;CALL
;    sig_id = get_iau_format(tobs,x,y,sig_start,lat=lat,lon=lon)
;
;INPUT
;    tobs - Time of observation for X, Y coordinates
;    X    - coord in HPC arcsec
;    Y    - coord in HPC arcsec
;    sig_start - start time to use a reference for making the IAU coordinate
;    lat  - Carrington latitude (output)
;    lon  - Carrington longitude (output)
;
;
;OUTPUT
;    sig_id - IAU formated Sigmoid ID number  
;
;###########################################################################
function get_iau_format,tobs,x,y,sig_start,lat=lat,lon=lon

;Format input variables
gi = tobs
xi = x
yi = y
t1 = sig_start

;Rotate to start time
rot_p = rot_xy(xi,yi,tstart=gi,tend=t1)

;Create directory for output png files
;Subscribe to IAU standard on output format
iau_time = strsplit(gi,'.',/extract)
iau_time = iau_time[0]

;IAU_format for coordinates
iau_cor = '("L",I03,"C",I03)'

;Get distance to sun
d_sun = get_sun(t1)

WCS_CONV_HPC_HG, rot_p[0], rot_p[1], lon, lat, dsun_obs=d_sun[0], length_units='AU',/carr, /pos_long
iau_pos = string([round(lon),round(lat)],format=iau_cor)

;Create new unique ID with IAU standard
sig_id = 'SOL'+iau_time+iau_pos


return,sig_id
end

;#############################################################################
;
;
;INPUTS
;    times      -   A csv file containing times to analyze sigmoid filaments
;                   CSV format must be as follows:
;                   formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F'
;                   readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
;                   length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats
;
;
;OUTPUTS
;    Creates directory structure of symbolic links to make new sigmiod catalog. 
;    Combines movies from HMI, sigmoid long movies, and flare movies
;
;#############################################################################
pro create_combined_movies,times,hmi_arch=hmi_arch,aia_arch=aia_arch,flr_arch=flr_arch,out_arch=out_arch

;Updates with Patty's new output format 2018/06/13 J. Prchlik
formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F'
readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
       length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats

;Set HMI movie directory 
if keyword_set(hmi_arch) then hmi_arch = hmi_arch else hmi_arch = 'hmi_movie/'
hmi_arch = hmi_arch+'/'

;Set SDO/AIA movie directory 
if keyword_set(aia_arch) then aia_arch = aia_arch else aia_arch = 'aia_movie/'
aia_arch = aia_arch+'/'

;Set flare movie directory 
if keyword_set(flr_arch) then flr_arch = flr_arch else flr_arch = 'aia_arch_cutout/'
flr_arch = flr_arch+'/'

;Set output directory for symbolic links
if not keyword_set(out_arch) then out_arch = 'combined_movies/'
out_arch = out_arch+'/'


;good sigmoid tbest times (i.e. contains time string)
goodt = where(strlen(tobs) eq 23)

;format for output symbolic link directory
out_fmt = '(A30,"/")'

;File to write to
openw,33,'json.txt.test'
printf,'{',33
printf,'}',33

free_lun,33

stop

;Loop over good indices
;for ii=0,n_elements(goodt)-1 do begin
;cut loop short for testing
for ii=0,2 do begin
    ;Set index to value with goodt
    i = goodt[ii]
    ;get index for a good time
    ;Now x, y are measured at tobs instead of tbest 2018/06/14 J. Prchlik
    gi = tobs[i]
    xi = x[i]
    yi = y[i]

    ;if there is no new measurements of the sigmoid just continue
    if gi eq '0' then continue




    ;

    ;Get sigmoid start and end times
    t1 = sig_start[i]
    t2 = sig_end[i]

    sig_id = get_iau_format(gi,xi,yi,t1,lat=lati,lon=loni)

    ;Directory to put the symoblic link
    full_dir = out_arch+string([sig_id],format=out_fmt)
    ;Remove : characters
    full_dir = str_replace(full_dir,':','')
    full_dir = str_replace(full_dir,'-','')

    ;###############################################################################
    ;
    ;HMI Movies
    ;
    ;###############################################################################
    ;get mp4 from hmi directory
    hmi_mov = hmi_dir+sig_id+'_hmi.mp4'

    ;###############################################################################
    ;
    ;SDO/AIA Evolution Movies
    ;
    ;###############################################################################
    ;Create text for link related to sigmoid long mp4 file
    aia_movie = aia_dir+sig_id+'.mp4'

     

    ;###############################################################################
    ;
    ;Flare Movies
    ;
    ;###############################################################################
 

    ;Create directory for output png files
    full_flr = flr_arch+strcompress(ID[i],/remove_all)+'/'

    ;get all flare movies
    flare_files = file_search(full_flr+'*mp4',/full)

    ;flare link
    if n_elements(size(flare_files)) lt 4 then continue

    ;Create text for flare links
    for j=0,n_elements(flare_files)-1 do begin
       ;File name
       file_name = strsplit(flare_files[j],'/',/extract)
       file_name = file_name[n_elements(file_name)-1]


    endfor
    


endfor
end
