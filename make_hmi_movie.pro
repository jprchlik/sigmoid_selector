
;#############################################################
;NAME:
;    gauss_smooth    
;
;USAGE
;    limits = select_cutout(px,py,img_size,img_xmax,img_ymax)   
;
;#############################################################


;#############################################################
;NAME:
;    get_sig    
;
;USAGE
;    sig = get_sig(cube)
;
;DESCRIPTION
;  Compute sigma:
;#############################################################
function get_sig,cube

    
    sig=fltarr(1)
    
    ;Set up dummy variable
    dum=cube
    
    ;Compute histogram
    hist=histogram(dum,binsize=1.,locations=loc)
    w1=where((loc gt -75) and (loc lt 75))
    
    ;Get binned gaussian and wieghts
    hsml=hist(w1)
    vals=loc(w1)
    
    ;Fig gaussian and weights
    gf6=gaussfit(w1,hsml,a6,nterms=6)
    ;Compute FWHM from sigma
    fwhm=2.*sqrt(2.*alog(2))*a6(2)
    
    sig=fwhm/2.
    
    return,sig
end

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

pro make_hmi_movie,times,hmi_arch=hmi_arch,out_arch=out_arch
;set plot to Z Window
set_plot,'Z'

;Read in file containing TBEST
;readcol,times,ID,RATING,NOAA,AR_START,X,Y,AR_END,SIG_START,SIG_END,TBEST,format='LL,I,A,A,F,F,A,A,A,A'
;Updates with Patty's new output format 2018/06/13 J. Prchlik
formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F'
readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
       length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats
;Set archive directory for download aia files
if keyword_set(hmi_arch) then hmi_arch = hmi_arch else hmi_arch = 'hmi_arch/'
hmi_arch = hmi_arch+'/'

;Set archive directory for output png files
if keyword_set(out_arch) then out_arch = out_arch else out_arch = 'hmi_movie/'
out_arch = out_arch+'/'

;Get list of hmi files
hmi_list = file_search(hmi_arch+"hmi*fits")

;Solar radius in cm
phy_rad = 6.957E10 ;cm

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

;These are edgedog radii
 ;radius to scale to select the sigmoid mag field
rad_1 = 3.6
rad_2 = 15.


;good sigmoid tbest times (i.e. contains time string)
goodt = where(strlen(tobs) eq 23)

;Cadance for image creation in seconds
img_cad = 30.*60.

;title format
title_fmt = '("HMI ID: ",A30," ")'

;format for output png file directory
out_fmt = '(A30,"/")'

;IAU_format for coordinates
iau_cor = '("L",I03,"C",I03)'

;Restore sigmiod save file
restore,'Sigmoids2007to2017.sav'


;Rotate images by given angle
rot_mat = [[-1.,0.],[0.,-1.]]

;Download HMI data for all the best times
for ii=0,n_elements(goodt)-1 do begin



    ;Set index to value with goodt
    i = goodt[ii]
    
    ;get index for a good time
    ;Now x, y are measured at tobs instead of tbest 2018/06/14 J. Prchlik
    gi = tobs[i]
    xi = x[i]
    yi = y[i]


    ;if there is no new measurements of the sigmoid just continue
    if gi eq '0' then continue

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
    match_fname = fname[chk_i]

    ;Get only unique values
    match_files=match_files[uniq(match_files)]
    match_fname=match_fname[uniq(match_files)]
    
    ;Check the first file before reading in all files
    mreadfits,match_files[0], index,data

    ;Create directory for output png files
    ;Subscribe to IAU standard on output format
    iau_time = strsplit(gi,'.',/extract)
    iau_time = iau_time[0]

    ;Get Carrington coordinates 
    wcs = fitshead2wcs( index(0) )
    rot_p = rot_xy(xi,yi,tstart=gi,tend=index(0).date_d$obs)
    WCS_CONV_HPC_HG, rot_p[0], rot_p[1], lon, lat, WCS=WCS, /carr, /pos_long
    iau_pos = string([round(lon),round(lat)],format=iau_cor)

    ;Create new unique ID with IAU standard
    sig_id = 'SOL'+iau_time+iau_pos
    full_dir = out_arch+string([sig_id],format=out_fmt)
    ;Remove : characters
    full_dir = str_replace(full_dir,':','')
    full_dir = str_replace(full_dir,'-','')
    if file_test(full_dir) eq 0 then file_mkdir,full_dir
   
    ;If save file exists continue 2018/06/20
    if file_test(full_dir+'/'+str_replace(sig_id,':','')+'.sav') eq 1 then continue
  
    ;prep hmi data
    ;hmi_prep,match_files,findgen(n_elements(match_files)),index,odata
    mreadfits,match_files, index,data
    ;hmi_prep,iindex,data,index,odata
    ;hmi_prep,hmi_list[chk_i],[1]findgen(n_elements(chk_i)-1),index,data


    

    

    ;width of image window in pixels
    win_w= 700
    sc = 3

    ;Get dynamic window to make for larger sigmiods
    ;Get the width of the sigmoid in pixels and add 100
    ;pixels on either side
    sig_p = length[i]/index(0).cdelt1 ; sigmiod length in hmi pixels
    ;Use window width to be at least twice size of the simoid
    tmp_w = 250+2.*sig_p ; Caculate new window width

    ;use the modified window if the sigmoid is larger 
    if tmp_w gt win_w then win_w = round(tmp_w)
    

    ;Set up device
    device,set_resolution=[win_w,win_w],decomposed=0,set_pixel_depth=24
    


    ;Init time, total intensity, neg. intensity, pos. intensity, area
    obs_time = [] ; observation time
    obs_qual = [] ; observation quality
    tot_ints = [] ; Total magnetic field intensity
    pos_ints = [] ; Positive magnetic field intensity
    neg_ints = [] ; Negative magnetic field intensity
    pix_area = [] ; Area of ROI in pixels^2
    tot_area = [] ; Area of ROI in cm^2
    roi_save = [] ; ROI object
    phy_save = [] ; ROI object

    ; plot each hmi observation
    for j=0,n_elements(index)-1 do begin

        ;if image quality greater than 90000 exit
        if index(j).quality gt 90000 then continue
        
        

        ;Plot restricted range
        ;Rotate coordinate to image time
        rot_p = rot_xy(xi,yi,tstart=gi,tend=index(j).date_d$obs)


        ;Get solar radius in arcsec at given time
        sol_rad = index(j).rsun_obs
        ;Get conversion from arcsec to cm
        arc_phy = phy_rad/sol_rad ; cm/arcsec
        ;Area of pixel in cm
        are_pix = arc_phy^2*abs(index(j).cdelt1*index(j).cdelt2)

        ;rotate by 180 degrees
        ;rot_pix = [pix_x,pix_y]#rot_mat
        ;pix_x = rot_pix[0]
        ;pix_y = rot_pix[1]

        ;Center rotated by 180 degrees
        cent_pix = [index(j).crpix1-index(j).naxis1/2.,index(j).crpix2-index(j).naxis2/2.]#rot_mat
        cent_x = cent_pix[0]+index(j).naxis1/2.
        cent_y = cent_pix[1]+index(j).naxis2/2.

        ;rotate delta by 180 degrees
        delt_pix = [index(j).cdelt1,index(j).cdelt2]#rot_mat
        delt_x = -delt_pix[0]
        delt_y = -delt_pix[1]

        ;Get center pixel value
        pix_x = (rot_p[0]/delt_x+cent_x)
        pix_y = (rot_p[1]/delt_y+cent_y) 
        ;Get range around pix_x and pix_y values
        lims = select_cutout(pix_x,pix_y,win_w,index(j).naxis1,index(j).naxis2)

        ;Store limite seperately
        pxmin = lims[0]
        pxmax = lims[1]
        pymin = lims[2]
        pymax = lims[3]

        ;Get image data to plot
        fimg = data[*,*,j]

        ;Rotate the unprepped image
        fimg = rot(fimg,180)

        ;get sigma from image
        img_sig = get_sig(fimg)

        ;set noise thresholds
        zero_lev = 2*img_sig    ; zero level in G
        ;zero_lev2 = 3*sig



        ;Cut image to restricted window
        img = fimg(pxmin:pxmax,pymin:pymax)

        ;Make sure image is the correct size if no loop off left pixels
        img_size = size(img)

        ;Check x
        ;if img_size[1] gt win_w then img = img(0:win_w-1,*)
        ;;Check y
        ;if img_size[2] gt win_w then img = img(*,0:win_w-1)


        ;Create a smoothed version of the image
        simg= img

        

        ;Remove noisy values
        bzzero = where(abs(simg) le  zero_lev)
        simg(bzzero) = 0.0                                ; All pixels below a certain value are set to zero


        ;Despike image
        cormag_p0=fltarr(win_w,win_w)
        cormag_n =fltarr(win_w,win_w)
        
        ;Get positive and negative spikes
        cormag_p0(where(simg ge 0))=simg(where(simg ge 0))
        cormag_n(where(simg lt 0)) =simg(where(simg lt 0))


        ;create spike pixel masks
        ;cormag_dp_p0=nospike(cormag_p0,thre=0.9,bright=0.99,imap=imap_p0)
        cormag_dp_p0=nospike(cormag_p0,thre=0.9,bright=0.85,imap=imap_p0)
        cormag_dp_ni =nospike(cormag_n,thre=0.9,bright=0.85,imap=imap_n)

        ;Remove spikes from image
        ;img(where(imap_p0 eq 1))=0
        ;img(Where(imap_n eq 1))=0
        simg(where(imap_p0 eq 1))=0
        simg(Where(imap_n eq 1))=0

        stop
        ;gaussian smooth image
        gimg = gauss_smooth(abs(simg),30,/edge_truncate)
        ;gimg = abs(simg)



        ;set origin variables correct from center origin to bottom left origin 2018/06/08 J. Prchlik
        org_x = -(cent_x-1-pix_x+win_w/2)*delt_x
        org_y = -(cent_y-1-pix_y+win_w/2)*delt_y

        ;Find the boundaries in the smoothed image
        rad_1 = 1.
        ;rad_2 = 300.
        ;Use the sigmoids measured size +20 pixels to look for features
        ;rad_2 = sig_p+100.
        rad_2 = win_w/2-1
        edge = edge_dog(abs(gimg),radius1=rad_1,radius2=rad_2,threshold=1,zero_crossings=[0,255])
  
        ;Get boundary of created countour
        CONTOUR,edge, LEVEL = 1,  $
               XMARGIN = [0, 0], YMARGIN = [0, 0], $
               /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, $
               XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS;/NODATA

        
        ;Plot image with rotation
        plot_image,img,min=-100,max=100,$
            origin=[org_x, $
            org_y], $
            scale=[delt_x,delt_y], $
            xtitle='X-postion (arcseconds)', $
            ytitle='Y-position (arcseconds)',$
            title=string([sig_id],format=title_fmt)+index(j).date_d$obs,$
            xcharsize=1.50, $
            ycharsize=1.50, $
            xcharthick=1.50, $
            ycharthick=1.50, $
            charsize=2.,$
            charthick=2.;/nosquare


        ;Create new ROI obejct using contour
        loadct,12
        ;search for an object containing the rotated point
        search = 1
        ;index of object
        ind_obj = 0

        ;Use largest area if no good position match
        areas = []
        ;closest area if no good position match
        dists = []

        ;Look through ROIs until you find one with the sigmoid point inside
        while search do begin 
            line = [LINDGEN(PathInfo(ind_obj).N), 0] 
            ;ROI in physical coordinates
            roi_phy = OBJ_NEW('IDLanROI', $
               delt_x*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[0, *]+org_x, $
               delt_y*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[1, *]+org_y) & $
               ;(pathXY(*, pathInfo(0).OFFSET +line ))[0, *], $
               ;(pathXY(*, pathInfo(0).OFFSET +line ))[1, *]) & $

           ;scale roi
           roi_scl = roi_phy
           roi_scl -> Scale,1.2, 1.2
 
           ;check if point is in ROI object
           pnt_chk = roi_phy -> containsPoints(rot_p[0],rot_p[1])

           ;get ROI AREA use largest area if no good point match
           geo_comp = roi_phy.ComputeGeometry(AREA = area, centroid=cent)
           areas = [areas,area]
           ;get ROI closest centriod if no good point match
           distr = sqrt(total((cent[[0,1]]-transpose(rot_p))^2))
           dists = [dists,distr]

           ;if point in object stop while and keep object
           if pnt_chk then search = 0
           ;No region found
           if ind_obj eq n_elements(PathInfo.N)-1 then search = 0 

           ;increment counter
           ind_obj = ind_obj+1
        endwhile

        ;correct for adding 1 to index
        ind_obj = ind_obj-1

        ;If no ROI found save figure and move to next time
        ;Changed to get ROI with largest area if no good point ROI match is found
        if pnt_chk eq 0 then begin
           ;No reason for this with better coordinates 2018/06/14
           ;get ROI with largest area
           ;nearest ROI
           ind_obj = where(dists eq min(dists))  

           ;Do other notmal ROI position stuff
           line = [LINDGEN(PathInfo(ind_obj).N), 0] 
           ;ROI in physical coordinates
           roi_phy = OBJ_NEW('IDLanROI', $
              delt_x*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[0, *]+org_x, $
              delt_y*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[1, *]+org_y) & $
              ;(pathXY(*, pathInfo(0).OFFSET +line ))[0, *], $
          
           ; ;write png file in directory 
           ; TVLCT,r,g,b,/Get
           ; write_png,full_dir+str_replace(match_fname[j],'fits','png'),tvrd(/true),r,g,b
        endif

        


        
        ;ROI in pixel coordinates
        roi_obj = OBJ_NEW('IDLanROI', $
           (pathXY(*, pathInfo(ind_obj).OFFSET + line))[0, *], $
           (pathXY(*, pathInfo(ind_obj).OFFSET + line))[1, *]) & $


        ;Scale up object roi 2018/06/14
        ;Thought about but don't do too much could go wrong
        ;roi_obj -> Scale,1.6, 1.6

        ;DRAW_ROI, roi_obj, COLOR =200,/LINE_FILL
        ;Draw ROI on plot switch to plotting line around ROI instead of ROI in draw_roi
        plots, delt_x*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[0, *]+org_x, delt_y*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[1, *]+org_y,color= 200,thick=3
        loadct,0

        ;write png file in directory 
        TVLCT,r,g,b,/Get
        write_png,full_dir+str_replace(match_fname[j],'fits','png'),tvrd(/true),r,g,b


        ;create mask for image based on edge detection
        ;Switched to img_size because if you can't beat em join em
        maskResult = roi_obj -> ComputeMask(DIMENSIONS = [img_size[1],img_size[2]])            
        IMAGE_STATISTICS, abs(img), MASK = maskResult, $  
                        COUNT = maskArea , data_sum=tot_intensity   
        IMAGE_STATISTICS, img < 0., MASK = maskResult, $  
                        COUNT = neg_maskArea , data_sum=neg_intensity   
        IMAGE_STATISTICS, img > 0., MASK = maskResult, $  
                        COUNT = pos_maskArea , data_sum=pos_intensity   



        ;Save variables
        obs_time = [obs_time,index(j).date_d$obs] ; observation time
        obs_qual = [obs_qual,index(j).quality] ; observation time
        tot_ints = [tot_ints,tot_intensity] ; Total magnetic field intensity
        pos_ints = [pos_ints,pos_intensity] ; Positive magnetic field intensity
        neg_ints = [neg_ints,neg_intensity] ; Negative magnetic field intensity
        pix_area = [pix_area,maskArea] ; Area of ROI in pixels^2
        tot_area = [tot_area,maskArea*are_pix] ; Area of ROI in cm^2
        roi_save = [roi_save,roi_obj] ;ROI object in pixels
        phy_save = [phy_save,roi_phy] ;ROI object in physical coordinates

    endfor


    ;Save output sav file
    out_id = id[i]
    save,sig_id,out_id,obs_time,obs_qual,tot_ints,pos_ints,neg_ints,pix_area,tot_area,roi_save,phy_save,filename=full_dir+'/'+str_replace(sig_id,':','')+'.sav'

    ;create directory for symbolic links
    if file_test(full_dir+'/symlinks/') then file_delete,full_dir+'symlinks/',/recursive
    file_mkdir,full_dir+'/symlinks/'

    ;collect all png files
    png_files = file_search(full_dir+'*png',/FULLY_QUALIFY_PATH,count=png_count)


    ;If no png files found just exit
    if png_count lt 1 then continue 

    ;create symbolic links
    for j=0,n_elements(png_files)-1 do file_link,png_files[j],full_dir+'symlinks/'+string(j,format='(I04,".png")')

    ;Automatically create movie
    ;Where does ffmpeg live?
    if not keyword_set(ffmpeg) then begin 
            spawn, 'which ffmpeg', ffmpeg
            if strpos(ffmpeg, 'Command not found') ne -1 then begin
                    message, 'ERROR: FFMPEG does not appear to be installed.', /informational
                    return
            endif 
            ;ffmpeg = str_replace(ffmpeg, 'opt', 'usr')
    endif

    png_size = string([2*win_w,2*win_w],format='(I04,"x",I04)')
    framerate= strcompress(string(8),/remove_all) ;frames per second
    bitrate = ((24*win_w)^2)*framerate ; number of bits per frame 24 bit colors and win_w^2 image
    bit = trim(bitrate)+'k'
    ;Use aia_ffmpeg for simplicity
    call1 = '-y -f image2 -r '+framerate+' -i ' 
    call2 = '-pix_fmt "yuv420p" -vcodec libx264 -level 41 -crf 18.0 -b '+bit+' -r '+framerate+' '+ $
                    '-bufsize '+bit+' -maxrate '+bit+' -g '+framerate+' -coder 1 -profile main -preset faster ' + $
                    '-qdiff 4 -qcomp 0.7 -flags +loop+mv4 -partitions ' + $
                    '+parti4x4+partp8x8+partb8x8 -subq 7 -me_range 16 -keyint_min 1 -sc_threshold ' + $
                    '40 -i_qfactor 0.71 -rc_eq ''blurCplx^(1-qComp)'' -s '+png_size+' -b_strategy 1 ' + $
                    '-bidir_refine 1 -refs 6 -trellis 1 -x264opts ' + $
                    'keyint='+framerate+':min-keyint=1:bframes=1 -threads 2 '
                
    
    
    ;output file name
    outf = str_replace(sig_id,':','')+"_mag.mp4"
    outf = str_replace(outf,'-','')
    spawn, ffmpeg +' '+ call1 + full_dir+'symlinks/%4d.png'+' ' + call2 + full_dir+outf, result, errResult
endfor

end