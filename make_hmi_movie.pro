
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
;Switch . to -
times = str_replace(times,'.','-')

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

;width of image window in pixels
win_w= 500

;title format
title_fmt = '("HMI ID: ",I03," ")'

;Download HMI data for all the best times
for i=0,n_elements(goodt)-1 do begin
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

    ;Create large time array to find the best times for each cadence
    r_ti = double(time_arr ## (-1+dblarr(n_elements(hmi_time))))
    p_ti = double(hmi_time # (1+dblarr(n_elements(time_arr))))

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
    match_files=hmi_list[chk_i]

    ;Get only unique values
    match_files=match_files[uniq(match_files)]
    
  
    ;prep hmi data
    hmi_prep,match_files,findgen(n_elements(match_files)),index,data
    ;hmi_prep,hmi_list[chk_i],[1]findgen(n_elements(chk_i)-1),index,data

    

    ; plot each hmi observation
    for j=0,n_elements(index)-1 do begin

        ;Plot restricted range
        ;Rotate coordinate to image time
        rot_p = rot_xy(xi,yi,tstart=gi,tend=index(j).date_obs)

        ;Get center pixel value
        pix_x = rot_p[0]/index(j).cdelt1+index(j).crpix1-1
        pix_y = rot_p[1]/index(j).cdelt2+index(j).crpix2-1
 
        ;Get range around pix_x and pix_y values
        lims = select_cutout(pix_x,pix_y,win_w,index(j).naxis1,index(j).naxis2)

        ;Store limite seperately
        pxmin = lims[0]
        pxmax = lims[1]
        pymin = lims[2]
        pymax = lims[3]

        print,index(j).date_obs
        ;Get image data to plot
        fimg = data[*,*,j]

        ;Cut image to restricted window
        img = fimg(pxmin:pxmax,pymin:pymax)
        
        plot_image,img,min=-100,max=100,$
            origin=-[(index(j).crpix1-1-pix_x)*index(j).cdelt1, $
            (index(j).crpix2-1-pix_y)*index(j).cdelt2], $
            scale=[index(j).cdelt1,index(j).cdelt2], $
            xtitle='X-postion (arcseconds)', $
            ytitle='Y-position (arcseconds)',$
            title=string([id[i]],format=title_fmt)+index(j).date_obs,$
            xcharsize=1.25, $
            ycharsize=1.25, $
            charsize=2.5
        wait,1
        print,'HERE'

    endfor

endfor

end