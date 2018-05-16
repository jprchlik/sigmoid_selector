
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
;
;NAME:
;    make_aia_movie
;
;PURPOSE
;    Creates aia movie for each sigmoid
;
;CATEGORY:
;    Program, movie creation
;
;USAGE
;    filament_selector,times,aia_arch='aia_arch/'
;
;INPUTS
;    times      -   A csv file containing times to analyze sigmoid filaments
;
;OUTPUTS
;    aia files in aia_arc
;
;#############################################################

pro make_aia_movie,times,aia_arch=aia_arch,out_arch=out_arch
;set plot to Z Window
set_plot,'Z'

;Read in file containing TBEST
readcol,times,ID,RATING,NOAA,AR_START,X,Y,AR_END,SIG_START,SIG_END,TBEST,format='LL,I,A,A,F,F,A,A,A,A'
;Set archive directory for download aia files
if keyword_set(aia_arch) then aia_arch = aia_arch else aia_arch = 'aia_arch/symlinks/'
aia_arch = aia_arch+'/'

;Set archive directory for output png files
if keyword_set(out_arch) then out_arch = out_arch else out_arch = 'aia_movie/'
out_arch = out_arch+'/'

;Work around for aia_mkmovie 
;List all files with a certain wavelength
;This is where I have some example files
l193 = file_search(aia_arch+'*193.fits')
l304 = file_search(aia_arch+'*304.fits')
l335 = file_search(aia_arch+'*335.fits')
lxrt = l335

;list of wavelengths for 4 panel movie
wavs = ['193','304','335']
;create pointer array of paths for files
paths = ptrarr(n_elements(wavs))
;create pointer array of for file times
times = ptrarr(n_elements(wavs))

;Create new pointers in larger pointer array (Could search and loop)
paths[0] = ptr_new(l193)
paths[1] = ptr_new(l304)
paths[2] = ptr_new(l335)


;Convert pointers 
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
        fils = strmid(temp,tstr[nstr-2]+4,15)$
    ;extract string if in  current directory
    else $
        fils = strmid(temp,4,15)

 
    ;Add in time date formatting
    any_fmt = strmid(fils,0,4)+'/'+strmid(fils,4,2)+'/'+strmid(fils,6,2)+'T'+strmid(fils,9,2)+':'+strmid(fils,11,2)+":"+strmid(fils,13,2)

    ;add times to pointer 
    times[i] = ptr_new(anytim(any_fmt))
endfor


;    endif else begin
;        print,"xrt ain't seen nothing yet"
;    endelse
;
;endfor


;;separate file name
;fnames = strsplit(aia_list,'/',/extract) 
;fname = strarr(n_elements(fnames))
;;Number of columns
;splnu = n_elements(fnames[0])


;Cadance for image creation in seconds
img_cad = 30.*60.


;Download aia data for all the best times
for i=0,n_elements(goodt)-1 do begin


   ;Create sub pointer for each particular time range
   sub_point =  ptrarr(n_elements(wavs))


    
    ;get index for a good time
    gi = tbest[i]
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

    ;Create directory for output png files
    full_dir = out_arch+string([id[i]],format=out_fmt)
    if file_test(full_dir) eq 0 then file_mkdir,full_dir

    ;Loop over all pointers and get the best times
    for j=0,n_elements(wavs)-1 do begin

        ;get aia times from pointer
        aia_time = *times[j]
        aia_list = *par[j]

        ;Create large time array to find the best times for each cadence
        r_ti = double(time_arr ## (-1+dblarr(n_elements(aia_time))))
        p_ti = double(aia_time # (1+dblarr(n_elements(time_arr))))

        ;get time min for matches
        minv = min(r_ti+p_ti,min_loc,/abs,dimension=1)
        ;get min values less than 90 minutes only
        good_min = where(minv lt img_cad,matches)

        ;leave if there are no good matches
        if matches eq 0 then continue


        ;clip to only get closest matches
        min_loc = min_loc[good_min]
   
        ;convert 1D indices into 2D indices
        col_i = array_indices(r_ti,min_loc)

        ;get the file indices that match the time
        chk_i = col_i[0,*]


        ;Get matched files
        match_files = aia_list[chk_i]

        ;Get only unique values
        match_files=match_files[uniq(match_files)]
    
        ;store file names in pointer 
        sub_point[j] = ptr_new(match_files)
    endfor

    ;prep first aia observation data
    hmi_prep,*sub_point[0],[0],index,data

    ;Rotate best point to first time
    ;Plot restricted range
    ;Rotate coordinate to image time
    rot_p = rot_xy(xi,yi,tstart=gi,tend=index(0).date_obs)

    ;Get center pixel value
    pix_x = rot_p[0]/index(j).cdelt1+index(j).crpix1-1
    pix_y = rot_p[1]/index(j).cdelt2+index(j).crpix2-1
 
    ;Get range around pix_x and pix_y values
    cutout = select_cutout(pix_x,pix_y,win_w,index(j).naxis1,index(j).naxis2)

    ;Get XRT data in time range
    ;Query time catalog
    xrt_cat, t1, t2, catx, ofiles
    want = where((catx.naxis1 eq 384) AND ((catx.ec_fw2_ eq 'Ti_poly') OR (catx.ec_fw1_ eq 'Al_poly') OR (catx.ec_fw1_ eq 'Be_thin')) AND (catx.ec_imty_ eq 'normal'))
    read_xrt, ofiles[want], xrt_index, xrt_data, /quiet
    for i=0, n_elements(xrt_index)-1 do xrt_data[*,*,i] = (xrt_data[*,*,i] - 30) / xrt_index[i].exptime
    xrt_data = alog10(xrt_data >1)
    for i=0, n_elements(xrt_index)-1 do xrt_data[*,*,i] = bytscl(xrt_data[*,*,i], min=median(xrt_data[*,*,i]*0.85))
    xrt_data = byte(xrt_data)



    ;Use aia_mkmovie to make the movie
    aia_mkmovie,t1,t2,waves(0:2),cadence=1,/no_quality_check,/diff_rot,/multi_panel,cutout=cutout,path=sub_point,other_index=xrt_index, other_data=xrt_data

    stop
  
endfor
end