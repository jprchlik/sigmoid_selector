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
;PROGRAM
;    print_sig_IAU
;
;USAGE
;    print_sig_IAU,times
;INPUTS
;    times      -   A csv file containing times to analyze sigmoid filaments
;                   CSV format must be as follows:
;                   formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F'
;                   readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
;                   length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats
;
;
;OUTPUTS
;    print the ID side by side with the IAU ID
;#############################################################################
;Updates with Patty's new output format 2018/06/13 J. Prchlik
pro print_sig_IAU,times

formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F'
readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
       length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats,/preserve_null

;Loop over all sigmoids and print the Sig and IAU ID
for i=0,n_elements(ID)-1 do begin

    gi = tobs[i]
    xi = x[i]
    yi = y[i]

    ;if there is no new measurements of the sigmoid just continue
    if gi eq '0' then continue

    ;Get sigmoid start and end times
    t1 = sig_start[i]
    t2 = sig_end[i]

    sig_id = get_iau_format(gi,xi,yi,t1,lat=lati,lon=loni)

    print,ID[i],'   ',sig_id

endfor
end



;#############################################################################
;PROGRAM
;    create_combined_movies

;USAGE
;    create_combined_movies,times,hmi_arch=hmi_arch,aia_arch=aia_arch,flr_arch=flr_arch,out_arch=out_arch
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
       length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats,/preserve_null

;Set HMI movie directory 
if keyword_set(hmi_arch) then hmi_arch = hmi_arch else hmi_arch = 'hmi_movie_cutout/'
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
;Get all sigmoids not just ones with AIA observations 2019/02/22 J. Prchlik
goodt = dindgen(n_elements(tobs))

;format for hmi movie subdirectory
hmi_fmt = '(A30,"/")'

;format for output symbolic link directory
out_fmt = '(A30,"/")'

;Loop over good indices
for ii=0,n_elements(goodt)-1 do begin
    ;Set index to value with goodt
    i = goodt[ii]
    ;get index for a good time
    ;Now x, y are measured at tobs instead of tbest 2018/06/14 J. Prchlik
    gi = tobs[i]
    xi = x[i]
    yi = y[i]

    ;if there is no new measurements of the sigmoid just continue
    ;if gi eq '0' then continue

    ;Get sigmoid start and end times
    t1 = sig_start[i]
    t2 = sig_end[i]

    sig_id = get_iau_format(gi,xi,yi,t1,lat=lati,lon=loni)
    ;Get sigmiod catalog ID  2018/12/03 J. Prchlik
    sig_cid= strcompress(ID[i],/remove_all)

    ;Get sigmoid catalog ID in I03 format 2019/02/07 J. Prchlik
    sig_fid = string(ID[i],format='(I03)')

    ;Directory to put the symoblic link
    ;Now using sigmoid IDs 2018/12/03
    full_dir = out_arch+string([sig_cid]);,format=out_fmt)
    ;Remove : characters
    full_dir = str_replace(full_dir,':','')
    full_dir = str_replace(full_dir,'-','')

    ;make output directory
    hmi_dir = full_dir+'/hmi_movie/'
    if file_test(hmi_dir) eq 0 then file_mkdir,hmi_dir
    aia_dir = full_dir+'/aia_movie/'
    if file_test(aia_dir) eq 0 then file_mkdir,aia_dir
    flr_dir = full_dir+'/flr_movie/'
    if file_test(flr_dir) eq 0 then file_mkdir,flr_dir


    ;###############################################################################
    ;
    ;HMI Movies
    ;
    ;###############################################################################
    ;get files from hmi directory
    ;I used inconsistent start times annoyingly 2018/08/23
    ;Find nearest SOL ID
    ;full_hmi = hmi_arch+string([sig_id],format=out_fmt)
    ;full_hmi = hmi_arch+string([sig_id],format=out_fmt)
    full_hmi = hmi_arch+strcompress(ID[i],/remove_all)+'/'
    ;No longer needed now that I have consistent sigmoid IDs 2018/11/09 J. Prchlik
    ;;;Remove : characters
    ;;full_hmi = str_replace(full_hmi,':','')
    ;;full_hmi = str_replace(full_hmi,'-','')
    ;;;Fix inconsistent ID
    ;;full_hmi = file_search(strmid(full_hmi,0,29)+'*')

    ;;if n_elements(size(full_hmi)) lt 4 then begin
    ;;    print,'Missing HMI '+string([sig_id],format=out_fmt)
    ;;    continue
    ;;endif 

    ;;case 1 of
    ;;    ((n_elements(size(full_hmi)) eq 4) and (n_elements(full_hmi) eq 1)): full_hmi = full_hmi[0]
    ;;    ((n_elements(size(full_hmi)) eq 4) and (n_elements(full_hmi) gt 1)): begin
    ;;        dis = fltarr(n_elements(full_hmi))
    ;;        for j=0,n_elements(full_hmi)-1 do begin
    ;;            cut_one = strsplit(full_hmi[j],'L',/extract) 
    ;;            cut_one = cut_one[n_elements(cut_one)-1]
    ;;            cut_two = fltarr(strsplit(cut_one,'C',/extract))
    ;;            dis[j] = sqrt(total(((cut_two-[loni,lati])*[cos(!dtor*lati),1])^2))
    ;;        endfor
    ;;        best = where(dis eq min(dis))
    ;;        full_hmi = full_hmi[best]
    ;;    end
    ;;endcase

    ;HMI movie
    hmi_mov = file_search(full_hmi+'/*mp4',/FULLY_QUALIFY_PATH ,count=hmi_movie_found)


    ;Only create HMI directory if there are in fact mp4 files
    ;And use unique sigmiod IDs 2018/12/03 J. Prchlik
    if hmi_movie_found ne 0 then begin
        ;Create symbolic link
        ;Removed previously created symbolic link 2018/11/09 J. Prchlik
        ;Switch to file_copy as opposed to symbolic link 2018/12/17 J. Prchlik
        if not file_test(hmi_dir+sig_cid+'_hmi.mp4') then $ 
            ;file_link,hmi_mov[0],hmi_dir+sig_cid+'_hmi.mp4' $
            file_copy,hmi_mov[0],hmi_dir+sig_cid+'_hmi.mp4' $
        else begin
            FILE_DELETE,hmi_dir+sig_cid+'_hmi.mp4'
            ;file_link,hmi_mov[0],hmi_dir+sig_cid+'_hmi.mp4' 
            ;Switch to file_copy as opposed to symbolic link 2018/12/17 J. Prchlik
            file_copy,hmi_mov[0],hmi_dir+sig_cid+'_hmi.mp4' 
        endelse
     endif

    ;###############################################################################
    ;
    ;SDO/AIA Evolution Movies
    ;
    ;###############################################################################
    ;Create link related to sigmoid long mp4 file
    ;Updated with I03 formatted sigmoid ID string 2019/02/07 J. Prchlik
    aia_file = file_search(aia_arch+sig_fid+'.mp4',/fully)
    if ((not file_test(aia_dir+sig_fid+'.mp4')) and (n_elements(size(aia_file)) gt 3)) then begin
        ;Switch to file_copy as opposed to symbolic link 2018/12/17 J. Prchlik
        ;Updated with I03 formatted sigmoid ID string 2019/02/07 J. Prchlik
        ;file_link,aia_file[0],aia_dir+sig_id+'.mp4'
        file_copy,aia_file[0],aia_dir+sig_fid+'.mp4'
    endif
     

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
    if n_elements(size(flare_files)) ge 4 then begin
        ;get array of times and classes for sending to create_flux_plot 2018/09/25 J. Prchlik
        flare_times = strarr(n_elements(flare_files))
        flare_class = strarr(n_elements(flare_files))

        ;Loop and create links
        for j=0,n_elements(flare_files)-1 do begin
           ;File name
           file_name = strsplit(flare_files[j],'/',/extract)
           file_name = file_name[n_elements(file_name)-1]

           ;Added to pass to hmi flux plot
           year  = strmid(file_name,4,4)
           month = strmid(file_name,8,2)
           day   = strmid(file_name,10,2)
           hour  = strmid(file_name,13,2)
           min   = strmid(file_name,15,2)
           class = strmid(file_name,20,2)
           sclass= strmid(file_name,23,1)

           ;Add flare times and class to variables
           flare_times[j] = year+'/'+month+'/'+day+' '+hour+':'+min+':00'
           flare_class[j] = class+'.'+sclass

           ;Create symbolic link
           if not file_test(flr_dir+sig_id+file_name) then $ 
               ;file_link,flare_files[j],flr_dir+sig_id+file_name
                ;Switch to file_copy as opposed to symbolic link 2018/12/17 J. Prchlik
               file_copy,flare_files[j],flr_dir+sig_id+file_name

        endfor
    endif
    
    ;###############################################################################
    ;
    ;HMI Flux evolution plots
    ;
    ;###############################################################################
    ;HMI  save file to use to create a png file
    hmi_png = file_search(full_hmi+'/*sav',/FULLY_QUALIFY_PATH ,count=found_hmi_sav)

    ;If there is no mag. field analysis for the sigmoid yet skip creating a flux evolution plot
    if found_hmi_sav ne 0 then begin
        if n_elements(size(flare_files)) ge 4 then $ 
            create_flux_plot,hmi_png[0],hmi_dir,t1,t2,f_stim=flare_times,f_scls=flare_class,fileout=flux_evo_plot $
        else $
            create_flux_plot,hmi_png[0],hmi_dir,t1,t2,fileout=flux_evo_plot
    endif


endfor
end
