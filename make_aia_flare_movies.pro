;#############################################################
;
;NAME:
;    select_cutout   
;
;PURPOSE
;
;CATEGORY:
;
;USAGE
;    limits = select_cutout(px,py,img_size,img_xmax,img_ymax)   
;
;INPUTS
;
;OUTPUTS
;
;#############################################################
function select_cutout,px,py,img_size,img_xmax,img_ymax

    
    ;compute the min and max x and y pixel ranges
    pxmin = px-(img_size/2.) 
    pxmax = px+(img_size/2.)-1
    pymin = py-(img_size/2.) 
    pymax = py+(img_size/2.)-1 
    
    square = 1
    ;Make sure everything is squared array before returning
    ;Prevents issues at the corners
    while square do begin
        ;make sure pixel values are within image limits
        case 1 of 
            (pxmin lt 0):  begin
                offset = abs(pxmin)
                pxmin = 0
                pxmax = pxmax+offset
            end
            (pymin lt 0):  begin
                offset = abs(pymin)
                pymin = 0
                pymax = pymax+offset
            end
            (pxmax gt img_xmax):  begin
                offset = img_xmax-pxmax
                pxmax = img_xmax-1
                pxmin = pxmin+offset-1
            end
            (pymax gt img_ymax):  begin
                offset = img_ymax-pymax
                pymax = img_ymax-1
                pymin = pymin+offset-1
            end
            else: square= 0
        endcase 
     endwhile
    
    ;force integers
    pxmin = fix(pxmin)
    pxmax = fix(pxmax)
    pymin = fix(pymin)
    pymax = fix(pymax)

    ;return limits for plot
    return,[pxmin,pxmax,pymin,pymax]

end
;#############################################################



;#############################################################
;
;NAME:
;    make_aia_flare_movies   
;
;PURPOSE
;    Make flare movies for sigmoid catalog
;
;CATEGORY:
;    Program, data gathering
;
;USAGE
;    make_aia_flare_movies,times,aia_arch='aia_arch/'
;
;INPUTS
;    flare_sav  -   A sav file containing flare times and positions. This file is created by flare_cme_sigcat.pro.
;                   This program can only run after running get_aia_files_cutout.pro.
;    times      -   A csv file containing sigmoid information
;
;OUTPUTS
;    aia movies
;
;#############################################################
pro make_aia_flare_movies,flare_sav,times,aia_arch=aia_arch,wave=wave,win_w=win_w

;Set window size
if keyword_set(win_w) then win_w = win_w else win_w = 700


;set plot to X Window
set_plot,'X'

;Read in file containing TBEST
formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F'
readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
       length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats,/preserve
;Set archive directory for download aia files
if keyword_set(aia_arch) then aia_arch = aia_arch else aia_arch = 'aia_arch_cutout/'
aia_arch = aia_arch+'/'
;Set wavelength to download AIA files
if keyword_set(wave) then wave = wave else wave = ['193','304','335']

;restore save files with flare association (big_str)
restore,flare_sav


;Cadence  
img_cad = 45. ;45 Seconds


;Download aia data for all the best times
for ii=0,n_elements(big_str)-1 do begin


    ;Switch to sigmoid start and end time because sigmiod ID changed  2018/07/12
    best_ind = where(((SIG_START eq big_str[ii].sigmd_s) and (SIG_END eq big_str[ii].sigmd_e)),match_count)


    ;;Get index of nearest sigmoid
    ;best_dis = min(dif_pos,best_ind,/Abs)

    if match_count eq 0 then continue ; No sigmoid found
    if match_count gt 1 then continue ; Skip dulpicates for now J. Prchlik 2018/07/18

    ;get good flare start times
    poss_t = where((big_str[ii].flare_s ne 'None') and (big_str[ii].flare_s ne ''),count)
 
    ;Exit loop if no flares found
    if count eq 0 then continue

    ;Create directory for output png files
    full_dir = aia_arch+strcompress(ID[best_ind],/remove_all)+'/'
    ;Skip already created directories for now 2018/07/30 J. Prchlik
    
    ;2018/08/06 J. Prchlik gather files for creating flare movie
    ;if file_test(full_dir) eq 0 then file_mkdir,full_dir else continue
    ;Work around for aia_mkmovie 
    ;List all files with a certain wavelength
    ;This is where I have some example files
    l193 = file_search(full_dir+'*193.image.fits')
    l304 = file_search(full_dir+'*304.image.fits')
    l335 = file_search(full_dir+'*335.image.fits')
    lxrt = l335


    ;If no files found just continue
    if ((n_elements(size(l193)) lt 4) and (n_elements(size(l304)) lt 4) and (n_elements(size(l335)) lt 4)) then continue

    ;list of wavelengths for 4 panel movie
    wavs = ['193','304','335']
    ;create pointer array of paths for files
    paths = ptrarr(n_elements(wavs))
    ;create pointer array of for file times
    times = ptrarr(n_elements(wavs))
    ;good sigmoid tbest times (i.e. contains time string)
    goodt = where(strlen(tbest) eq 23)
    
    ;Create new pointers in larger pointer array (Could search and loop)
    paths[0] = ptr_new(l193)
    paths[1] = ptr_new(l304)
    paths[2] = ptr_new(l335)


    ;Convert pointers into times
    for i=0,n_elements(wavs)-1 do begin
        ;extract the string time and date
        temp = *paths[i]
        ;first and last index of string with AIA_YYYYMMDD_HHMMSS_WAVE.fits
        fstr = strsplit(temp[0],'/',length=estr)
        tstr = fstr+estr
        ;number of elements in split 
        nstr = n_elements(tstr)
    
        ;extract string if in sub-directory
        if nstr gt 1 then $
            ;Changed formatting based on different file format from JSOC 2018/08/06
            ;fils = strmid(temp,tstr[nstr-2]+4,15)$
            fils = strmid(temp,tstr[nstr-2]+18,17) $
        ;extract string if in  current directory
        else $
            fils = strmid(temp,4,15)
    
     
        ;Add in time date formatting
        ;Modified with JSOC cutout time format 2018/08/06 J. Prchlik
        any_fmt = strmid(fils,0,4)+'/'+strmid(fils,5,2)+'/'+strmid(fils,8,2)+'T'+strmid(fils,11,2)+':'+strmid(fils,13,2)+":"+strmid(fils,15,2)
    
        ;add times to pointer 
        times[i] = ptr_new(anytim(any_fmt))
    endfor



    ;loop over all start times
    for ij = 0,n_elements(poss_t)-1 do begin 
 
        ;Create sub pointer for each particular time range
        sub_point =  ptrarr(n_elements(wavs))

        ;get index for a good time
        i = poss_t[ij]

    ;get index for a good time
        gi = tbest[ij]
        xi = x[ij]
        yi = y[ij]

        ;get time range to search over
        t1 = anytim(big_str[ii].flare_s[i])-3600. ; move forward 1 hour
        t2 = anytim(big_str[ii].flare_e[i])+3600. ; Add 1 hour after
      
        ;get flare peak time and class
        t_peak = big_str[ii].flare_p[i]
        f_clas = big_str[ii].flare_c[i]

        ;Filename of output movie
        file_out_fmt = 'AIA_'+str_replace(str_replace(t_peak,'-',''),':','')+'_'+str_replace(f_clas,'.','_')+'.mp4'
        ;Do not rerun movie if already created 2018/08/06 J. Prchlik
        if file_test(full_dir+'/'+file_out_fmt) then continue
        

        ;get start time string
        ts = anytim(t1,/hxrbs)
        ts = '20'+ts
        ts = str_replace(ts,', ','T')
        ts = str_replace(ts,'/','-')

        ;get end time string
        te = anytim(t2,/hxrbs)
        te = '20'+te
        te = str_replace(te,', ','T')
        te = str_replace(te,'/','-')

        ;Create array of good times at 90 minute cadence
        add_cad  = 1
        counter = 1
        time_arr = [t1]
        ;loop until cadence maximum is reached
        while add_cad eq 1 do begin 

            ;Add new cadence to time array 
            new_time = double(counter*img_cad+t1)
            time_arr = [time_arr,new_time]


            ;end if new time greater than endtime
            if new_time gt t2 then add_cad = 0
            counter = counter+1
        endwhile

        ;Check to see if the code locates any aia files
        found_fil = 0
        ;Count total number of files
        total_fil = 0
        ;Loop over all pointers and get the best times
        for j=0,n_elements(wavs)-1 do begin

            ;get aia times from pointer
            aia_time = *times[j]
            aia_list = *paths[j]

            ;Create large time array to find the best times for each cadence
            r_ti = double(time_arr ## (-1+dblarr(n_elements(aia_time))))
            p_ti = double(aia_time # (1+dblarr(n_elements(time_arr))))

            ;get time min for matches
            minv = min(r_ti+p_ti,min_loc,/abs,dimension=1)
            ;get min values less than 90 minutes only
            good_min = where(minv lt 50000*img_cad,matches)

            ;leave if there are no good matches
            if matches eq 0 then continue


            ;found aia files to analyze 
            found_fil = 1

            ;clip to only get closest matches
            min_loc = min_loc[good_min]
   
            ;convert 1D indices into 2D indices
            col_i = array_indices(r_ti,min_loc)

            ;get the file indices that match the time
            chk_i = col_i[0,*]


            ;Get matched files
            match_files = aia_list[chk_i]

            ;Get only unique values
            ;Not need J. Prchlik 2018/08/07
            ;match_files=match_files[uniq(match_files)]

            ;Count total number of files
            total_fil = n_elements(match_files)+total_fil

            ;Make sure all wavelengths have the same length in array
            if j eq 0 then begin
                total_l = n_elements(match_files)
            endif else begin
                 if n_elements(match_files) ne total_l then print,'Error Coming in movie creation'
            endelse
        
            ;store file names in pointer 
            sub_point[j] = ptr_new(match_files)
         endfor

        ;leave if no aia files found
         if total_fil lt 4 then continue


        ;prep first aia observation data
        filein = *sub_point[0]
        read_sdo,filein[0],index,data,/uncomp_delete,/noshell

        ;Rotate best point to first time
        ;Plot restricted range
        ;Rotate coordinate to image time
        rot_p = rot_xy(xi,yi,tstart=gi,tend=index(0).date_obs)

        ;Get center pixel value
        pix_x = rot_p[0]/index(0).cdelt1+index(0).crpix1-1
        pix_y = rot_p[1]/index(0).cdelt2+index(0).crpix2-1
 
        ;Get range around pix_x and pix_y values
        cutout = select_cutout(pix_x,pix_y,win_w,index(0).naxis1,index(0).naxis2)
        ;buffer for Y values to include in movie
        y_buff = 284
        
        ;Get XRT data in time range
        ;Query time catalog
        xrt_cat, ts, te, catx, ofiles


        ;If no XRT files found make movie with just AIA and continue
        if n_elements(size(ofiles)) lt 4 then begin
            aia_mkmovie,ts,te,wavs,cadence=1,/multi_panel,path=sub_point,/sequential,ref_times=aia_time,fname=file_out,/delete
            continue
        endif

        ;get X x-ray position at a given time
        x_val = fltarr(n_elements(catx))

        ;Loop and store xvalue for comparison
        for k=0,n_elements(x_val)-1 do begin
           x_pos = rot_xy(xi,yi,tstart=gi,tend=catx(k).date_obs)
           x_val[k] = x_pos[0]

        endfor

        ;Get regions which cover the region
        want = where( $
                     (catx.ec_imty_ eq 'normal') $
                     AND ((catx.ycen+catx.fovy*catx.cdelt1 gt rot_p[1]) AND (catx.ycen-catx.fovy*catx.cdelt1 lt rot_p[1])) $
                     AND (catx.chip_sum le 2) $
                     AND (catx.naxis1 le 513) $
                     AND (catx.naxis2 le 513) $
                     AND ((catx.xcen+catx.fovx*catx.cdelt1 gt x_val) AND (catx.xcen-catx.fovx*catx.cdelt1 lt x_val)))

        ;Find the most common dimensions
        unq_x = catx[UNIQ(catx.naxis1, SORT(catx.naxis1))].naxis1  
        unq_y = catx[UNIQ(catx.naxis2, SORT(catx.naxis2))].naxis2  
        

 
        ;Most common x,y dimension
        x_dim = 0
        y_dim = 0

        ;Loop and find the solution
        check_cnt_x = 0 
        for k=0,n_elements(unq_x)-1 do begin 
            sizer = where(catx[want].naxis1 eq unq_x[k],cnt) 
            if cnt gt check_cnt_x then begin  
                check_cnt_x = cnt 
                x_dim = unq_x[k] 
            endif  
        endfor 

        ;Loop and find the solution
        check_cnt_y = 0  
        for k=0,n_elements(unq_y)-1 do begin  
            sizer = where(catx[want].naxis2 eq unq_y[k],cnt)  
            if cnt gt check_cnt_y then begin   
                check_cnt_y = cnt  
                y_dim = unq_y[k]  
            endif   
        endfor  


        ;Get new solution with the most common values for Naxis1 and 2
        want = where( $
                     (catx.ec_imty_ eq 'normal') $
                     AND ((catx.ycen+catx.fovy*catx.cdelt1 gt rot_p[1]) AND (catx.ycen-catx.fovy*catx.cdelt1 lt rot_p[1])) $
                     AND (catx.chip_sum le 2) $
                     AND (catx.naxis1 le 513) $
                     AND (catx.naxis2 le 513) $
                     AND ((catx.xcen+catx.fovx*catx.cdelt1 gt x_val) AND (catx.xcen-catx.fovx*catx.cdelt1 lt x_val)) $
                     AND ((catx.naxis1 eq x_dim) and (catx.naxis2 eq y_dim)),found_small_xrt)

        if found_small_xrt gt 1 then begin
            ;Force allow for difference size images
            read_xrt, ofiles[want], xrt_index, xrt_data, /quiet
            for i=0, n_elements(xrt_index)-1 do xrt_data[*,*,i] = (xrt_data[*,*,i] - 30) / xrt_index[i].exptime
            xrt_data = alog10(xrt_data >1)
            for i=0, n_elements(xrt_index)-1 do xrt_data[*,*,i] = bytscl(xrt_data[*,*,i], min=median(xrt_data[*,*,i]*0.85))
            xrt_data = byte(xrt_data)

            ;Use aia_mkmovie to make the movie
            ;Need to add ref_time keyword to aia_panel_wrapper and remove aia_prep at line 99 in aia_pane_wrapper 2018/08/06
            ;File with edits located at /home/jprchlik/personaladditions/code/idl/aia_mkmovie_testbed/aia_mkmovie_testbed/aia_mkmovie
            aia_mkmovie,ts,te,wavs,cadence=1,/multi_panel,/sequential,/delete,$
                    path=sub_point,other_index=xrt_index, other_data=xrt_data,ref_times=aia_time,fname=file_out;, $
        endif else aia_mkmovie,ts,te,wavs,cadence=1,/multi_panel,path=sub_point,/sequential,ref_times=aia_time,fname=file_out,/delete

        ;If movie is not made continue because there were no observations in that time period
        if not file_test(file_out) then continue

        ;Move flare movie to new directory
        file_move,file_out,full_dir+'/'+file_out_fmt

    endfor

endfor

end
