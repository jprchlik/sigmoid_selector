
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
    hist=histogram(dum[where(finite(dum))],binsize=1.,locations=loc)
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
;    limits = select_cutout(px,py,img_size,img_xmax,img_ymax,tol=tol)   
;
;INPUTS
;
;OUTPUTS
;
;#############################################################
function select_cutout,px,py,img_size,img_xmax,img_ymax,tol=tol

    ;Number of times it can look to find the solution
    if keyword_set(tol) then tol = tol else tol = 1000

    
    ;compute the min and max x and y pixel ranges
    pxmin = px-(img_size/2.) 
    pxmax = px+(img_size/2.)-1
    pymin = py-(img_size/2.) 
    pymax = py+(img_size/2.)-1 
    
    ;is the image a squared away at the corners
    square = 1
    ;count for the number of loops
    iter = 0
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
        ;exit when looping too long
        ;Return the entire image for cutout
        iter = iter+1
        if iter gt tol then begin
            square=0
            pxmin = -9999.9
            pymin = -9999.9 
            pxmax = -9999.9 
            pymax = -9999.9 
        endif
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
;    make_hmi_movie_man_thres
;
;PURPOSE
;    Creates hmi movie for each sigmoid
;
;CATEGORY:
;    Program, movie creation
;
;USAGE
;    make_hmi_movie_man_thres,times,hmi_arch='hmi_arch/',rebinv=16
;
;INPUTS
;    times      -   A csv file containing times to analyze sigmoid filaments
;                   CSV format must be as follows:
;                   formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F'
;                   readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
;                   length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats
;    hmi_arch   -   Directory to output the HMI files to (Default = 'hmi_arch/')
;    cad        -   Cadence to get HMI files in seconds (Default = 30.*60)
;    rebinv     -   Rebinning pixel value for median smoothing computationally effciently (Default = 8)
;
;OUTPUTS
;    A save file called'manual_edge_dog_threshold.sav' which contains two variables ID and store_thres.
;    ID is the sigmoid ID, while store_thres are the the manually determinated threshold values for each ID.
;    The element number in each array corresponds to the same sigmoid. 
;
;#############################################################

;Running with Rebin = 32 2019/01/08 J. Prchlik
pro make_hmi_movie_cutout_man_thres,times,hmi_arch=hmi_arch,rebinv=rebinv
;set plot to X Window
set_plot,'X'

;Start with BW color table
loadct,0,/silent
;Read in file containing TBEST
;readcol,times,ID,RATING,NOAA,AR_START,X,Y,AR_END,SIG_START,SIG_END,TBEST,format='LL,I,A,A,F,F,A,A,A,A'
;Updates with Patty's new output format 2018/06/13 J. Prchlik
formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F'
readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
       length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats,/preserve_null
;Set archive directory for download aia files
if keyword_set(hmi_arch) then hmi_arch = hmi_arch else hmi_arch = 'hmi_arch_cutout/'
hmi_arch = hmi_arch+'/'

if keyword_set(rebinv) then rebinv = rebinv else rebinv = 16


;Solar radius in cm
phy_rad = 6.957E10 ;cm

;These are edgedog radii
 ;radius to scale to select the sigmoid mag field
rad_1 = 3.6
rad_2 = 15.


;good sigmoid tbest times (i.e. contains time string)
goodt = where(strlen(tobs) eq 23)

;Cadance for image creation in seconds
img_cad = 30.*60.

;title format
title_fmt = '("HMI ID: ",A3," ")'

;format for output png file directory
out_fmt = '(A30,"/")'

;IAU_format for coordinates
iau_cor = '("L",I03,"C",I03)'


;Rotate images by given angle
rot_mat = [[-1.,0.],[0.,-1.]]


;create new threshold varible to store and save for later analysis 2019/02/21 J. Prchlik
store_thres = fltarr(n_elements(tobs))-9999.9

;Name for output save file with the manual threshold values
thres_file = 'manual_edge_dog_threshold.sav'

;restore previous manual thresholds if applicable
if file_test(thres_file) then begin
    ;Manually threshold allow for updating 2019/02/21 J. Prchlik
    restore,thres_file
endif

;Initalize window as an X window
set_plot,'X'

;Download HMI data for all the best times
for ii=0,n_elements(goodt)-1 do begin



    ;Set index to value with goodt
    i = goodt[ii]


    ;skip sigmoid if it is already analyzed with threshold value 2019/02/21 J. Prchlik
    if store_thres[i] gt 0 then continue
    
    ;get index for a good time
    ;Now x, y are measured at tobs instead of tbest 2018/06/14 J. Prchlik
    gi = tobs[i]
    xi = x[i]
    yi = y[i]


    ;if there is no new measurements of the sigmoid just continue
    if gi eq '0' then continue

    ;get time range to search over
    t1 = sig_start[i]
    t2 = sig_end[i]
  
    ;Convert to anytimes
    at1 = anytim(t1)
    at2 = anytim(t2)


    ;Setup output directories
    sig_id = strcompress(ID[i],/remove_all)
   

    ;get list of local files
    match_files = file_search(hmi_arch+'/'+sig_id+"/*hmi*fits",count=file_cnt,/full)

    ;if no files found then continue
    if file_cnt eq 0 then continue
  

    ;Change order of files so first file calculated is in the middle
    ;This way the threshold values are set by an observation "on disk" 2019/01/11 J. Prchlik
    ;Get half way file
    match_moveu = fix(n_elements(match_files)/2)
    ;Get file name of removed file
    first_filen = match_files[match_moveu]
    ;Remove that file from array
    remove,match_moveu,match_files
    ;Add file back in at the beginning
    match_files = [first_filen,match_files]


    ;switch to read_fits because cutout header is missing some information which causes mreadfits to be unhappy 2018/11/05 J. Prchlik
    data = readfits(match_files[0],hdr, exten_no=0, /fpack,/silent) 


    ;width of image window in pixels
    win_w= 700
    sc = 3

    ;Get dynamic window to make for larger sigmiods
    ;Get the width of the sigmoid in pixels and add 100
    ;pixels on either side
    ;switch to fits_read syntax 2018/11/05 J. Prchlik
    sig_p = length[i]/sxpar(hdr,'cdelt1') ; sigmiod length in hmi pixels
    ;Use window width to be at least twice size of the simoid
    tmp_w = 250+2.*sig_p ; Caculate new window width

    ;use the modified window if the sigmoid is larger 
    if tmp_w gt win_w then win_w = round(tmp_w)
    

    
    ;Set up device
    device,decomposed=0;set_pixel_depth=24


    ;Init time, total intensity, neg. intensity, pos. intensity, area
    obs_loca = [] ; location of the observations
    obs_time = [] ; observation time
    obs_qual = [] ; observation quality
    tot_ints = [] ; Total magnetic field intensity
    pos_ints = [] ; Positive magnetic field intensity
    neg_ints = [] ; Negative magnetic field intensity
    pix_area = [] ; Area of ROI in pixels^2
    tot_area = [] ; Area of ROI in cm^2
    roi_save = [] ; ROI object in pixel units
    phy_save = [] ; ROI object in physical units
    elp_save = [] ; ROI object in physical units using a fitting ellipse
    ilp_save = [] ; ROI object in pixel units using a fitting ellipse
    pol_lens = [] ; Stored length of the polarity inversion line 2018/12/07 J. Prchlik


    ;Set Threshold  2018/12/20 Use to set a static threshold for each observations, which keeps ROI at a consisent size
    set_threshold = 1


    ;Get Only plot the middle observation for the threshold value 2019/02/21 J. Prchlik
    j = 0


    ;read in files 1 at a time
    data = readfits(match_files[j],hdr, exten_no=0, /fpack,/silent) 

    ;replace bad pixels
    bad_pix = data le -1.7
    data[bad_pix] = !values.f_nan


    ;if image quality greater than 90000 exit
    ;swtich to fits_read format
    if sxpar(hdr,'quality') gt 90000 then continue


    ;Make sure input fits file is an image
    if n_elements(size(data)) lt 5 then continue
    

    ;Get formatted date
    fmt_dat = str_replace(str_replace(strmid(sxpar(hdr,'T_OBS'),0,19),'.','-'),'_','T')  
    

    ;Plot restricted range
    ;Rotate coordinate to image time
    rot_p = rot_xy(xi,yi,tstart=gi,tend=fmt_dat,offlimb=offlimb)


    ;Get solar radius in arcsec at given time
    sol_rad = sxpar(hdr,'rsun_obs')
    ;Get conversion from arcsec to cm
    arc_phy = phy_rad/sol_rad ; cm/arcsec

    ;Only do the calculations for Magnetic fields less than 50 degrees on the surface
    ;Include if point is off the limb
    ;Switch to 40 degrees per discussion with Antonia 2018/12/06
    if ((sqrt(total(rot_p^2))/sol_rad gt sin(!dtor*40.)) OR (offlimb eq 1)) then continue


    ;calculate the length of sigmoid correcting for projection
    ;Remove projection correction 2018/12/07 because we are going to sit at DC J. Prchlik 
    sig_p = length[i]/sxpar(hdr,'cdelt1');*acos((sqrt(total(rot_p^2))/sol_rad)) ; sigmiod length in hmi pixels


    ;Area of pixel in cm
    ;udated with fits_read syntax
    are_pix = arc_phy^2*abs(sxpar(hdr,'cdelt1')*sxpar(hdr,'cdelt2'))


    ;rotate by 180 degrees
    ;rot_pix = [pix_x,pix_y]#rot_mat
    ;pix_x = rot_pix[0]
    ;pix_y = rot_pix[1]

    ;Center rotated by 180 degrees
    ;Update with fits_ read syntax 2018/11/05 J. Prchlik
    ;cutouts are not rotated by 180
    cent_pix = [sxpar(hdr,'crpix1')-sxpar(hdr,'naxis1')/2.,sxpar(hdr,'crpix2')-sxpar(hdr,'naxis2')/2.];#rot_mat

    cent_x = sxpar(hdr,'crpix1');cent_pix[0]+sxpar(hdr,'naxis1')/2.
    cent_y = sxpar(hdr,'crpix2');cent_pix[1]+sxpar(hdr,'naxis2')/2.

    ;rotate delta by 180 degrees
    ;Update with fits_ read syntax 2018/11/05 J. Prchlik
    ;cutouts are not rotated by 180
    delt_pix = [sxpar(hdr,'cdelt1'),sxpar(hdr,'cdelt2')]#rot_mat
    delt_x =-delt_pix[0];-delt_pix[0] 
    delt_y =-delt_pix[1];-delt_pix[1] 

    ;Get center pixel value
    pix_x = (rot_p[0]/delt_x+cent_x)
    pix_y = (rot_p[1]/delt_y+cent_y) 
    ;Get range around pix_x and pix_y values
    lims = select_cutout(pix_x,pix_y,win_w,sxpar(hdr,'naxis1'),sxpar(hdr,'naxis2'))
 
    ;If the program could not select a cutout coninue
    bad_lims = where(lims lt 0,bad_lim_cnt)
    if bad_lim_cnt gt 0 then continue

    ;Store limite seperately
    pxmin = lims[0]
    pxmax = lims[1]
    pymin = lims[2]
    pymax = lims[3]

    ;Get image data to plot
    fimg = data

    ;Rotate the unprepped image
    ;No longer required for cutouts 2018/11/05 J. Prchlik 
    ;fimg = rot(fimg,180)
    ; Establish error handler. When errors occur, the index of the
    ; error is returned in the variable Error_status:
    ;;;;CATCH, Error_status

    ;;;;;if data is empty continue
;   ;;;;;This statement begins the error handler:
    ;;;;IF Error_status NE 0 THEN BEGIN
    ;;;;   PRINT, 'Error index: ', Error_status
    ;;;;   PRINT, 'Error message: ', !ERROR_STATE.MSG
    ;;;;   ; Handle the error by extending A:
    ;;;;   CATCH, /CANCEL
    ;;;;   ;Make sure to set color table back to default
    ;;;;   loadct,0,/silent
    ;;;;   ;If Error go to the next image
    ;;;;   continue
    ;;;;ENDIF


    ;NaN out off limb pixels because they mess up sigma levels 2019/07/08 J. Prchlik
    ;Confirmed to Fix 2019/01/08 J. Prchlik
    ;This is because there is Extra 0 sigma data from the edges of the HMI image
    ;#######################################################################
    ;NAN out pixels off the limb 2019/01/02 J. Prchlik
    ;This way the distribution of used for spike correction is similar
    ;Get image size
    fimg_size = size(fimg)

    ;fill values outside rsun with median
    fimg_x = dindgen(fimg_size[1]) #  (intarr((fimg_size[1])) +1)-cent_x
    fimg_y = dindgen(fimg_size[2]) ## (intarr((fimg_size[2])) +1)-cent_y
    fimg_x = delt_x*fimg_x
    fimg_y = delt_y*fimg_y
    fimg_r = sqrt(fimg_x^2+fimg_y^2)
    fimg_s =  fimg_r/sol_rad

    ;If you use 0.0 you will skew the sigma levels 2019/01/08
    fimg(where(fimg_s gt 1.)) = !values.F_NAN ;0.0
    ;#######################################################################


    ;get sigma from image
    img_sig = get_sig(fimg)

    ;set noise thresholds
    zero_lev = 2*img_sig    ; zero level in G
    ;zero_lev2 = 3*sig


    ;get size of image
    fimg_size = size(fimg)

    ;Check x
    ;if img_size[1] gt win_w then img = img(0:win_w-1,*)
    ;;Check y
    ;if img_size[2] gt win_w then img = img(*,0:win_w-1)




    ;Remove noisy values
    ;This is level1 data with cutout, so try without rejection of pixels 2018/12/19 J. Prchlik
    bzzero = where(abs(fimg) le  zero_lev)
    fimg(bzzero) = 0.0                                ; All pixels below a certain value are set to zero



    ;Despike image
    fimg_size = size(fimg)
    full_win_x = fimg_size[1]
    full_win_y = fimg_size[2]
    cormag_p0=fltarr(full_win_x,full_win_y)
    cormag_n =fltarr(full_win_x,full_win_y)
    
    ;Get positive and negative spikes
    cormag_p0(where(fimg ge 0))=fimg(where(fimg ge 0))
    cormag_n(where(fimg lt 0)) =fimg(where(fimg lt 0))



    ;create spike pixel masks
    ;Switched back to old way 2018/12/19 J. Prchlik
    ;cormag_dp_p0=nospike(cormag_p0,thre=0.9,bright=0.99,imap=imap_p0)
    ;cormag_dp_p0=nospike(cormag_p0,thre=0.65,bright=0.99,imap=imap_p0,/silent)
    ;Switch to same lower threshold as - ddata 2018/12/21 J. Prchlik
    cormag_dp_p0=nospike(cormag_p0,thre=0.5,bright=0.99,imap=imap_p0,/silent)
    ;cormag_dp_p0=nospike(cormag_p0,thre=0.9,bright=0.85,imap=imap_p0,/silent)
    cormag_dp_ni =nospike(cormag_n,thre=0.9,bright=0.50,imap=imap_n,/silent)
    ;cormag_dp_ni =nospike(cormag_n,thre=0.9,bright=0.85,imap=imap_n,/silent)

    ;Remove spikes from image
    fimg(where(imap_p0 eq 1))=0
    fimg(where(imap_n eq 1))=0
    
    ;Scale down and check the scaling is by an integer
    scale_x = fimg_size[1]/rebinv
    scale_y = fimg_size[2]/rebinv

    ;Accounts for size assuming rebin is not a integer value of the axes size 2018/11/13 J. Prchlik
    ;check if xaxis is equal
    if scale_x*rebinv ne fimg_size[1] then begin
        diff_scale = abs(scale_x*rebinv-fimg_size[1]) 
        fimg= fimg(diff_scale:fimg_size[1]-1,*)
    endif
    ;check if yaxis is equal
    if scale_y*rebinv ne fimg_size[2] then begin
        diff_scale = abs(scale_y*rebinv-fimg_size[2]) 
        fimg= fimg(*,diff_scale:fimg_size[2]-1)
    endif


    ;Get new image size
    fimg_size = size(fimg)

    ;fill values outside rsun with median
    fimg_x = dindgen(fimg_size[1]) #  (intarr((fimg_size[1])) +1)-cent_x
    fimg_y = dindgen(fimg_size[2]) ## (intarr((fimg_size[2])) +1)-cent_y
    fimg_x = delt_x*fimg_x
    fimg_y = delt_y*fimg_y
    fimg_r = sqrt(fimg_x^2+fimg_y^2)
    ;Theta angle for each pixel 2018/11/13 J. Prchlik
    bin_cor_t =  asin(fimg_r/sol_rad)

    ;This must be corrected in Level1 images because you can see a spike in the 
    ;Magnetic Flux near the Edges 2019/01/08 J. Prchlik
    ;Cannot Confirm this, so leaving in
    ;correct pixels by cos(theta)^2
    fimg = fimg/cos(bin_cor_t)^2

    ;zero out all pixels greater than 50 degrees
    ;Only did for exper. Not in final draft 2018/12/03 J. Prchlik
    ;Re implimented for getting sigma from image
    ;Readded limb rejection 2018/12/20 J. Prchlik
    fimg(where(bin_cor_t gt !dtor*65.)) = !values.F_NAN ;0.0

    ;Create a low res smoothed version of the image
    simg = rebin(fimg,fimg_size[1]/rebinv,fimg_size[2]/rebinv)

    ;fill values outside rsun with median
    bin_cor_x = dindgen(fimg_size[1]/rebinv) # (intarr((fimg_size[1]/rebinv)) +1) -cent_x/rebinv
    bin_cor_y = dindgen(fimg_size[2]/rebinv) ## (intarr((fimg_size[1]/rebinv)) +1)-cent_y/rebinv
    bin_cor_x = delt_x*rebinv*bin_cor_x
    bin_cor_y = delt_y*rebinv*bin_cor_y
    bin_cor_r = sqrt(bin_cor_x^2+bin_cor_y^2)


    ;get median value when less than solar radius
    med_sun = median(simg[where(bin_cor_r lt sol_rad)])
    simg[where(bin_cor_r gt sol_rad)] = 0.0
    ;gaussian smooth image
    ;edge zero cuts too far into the image although it is faster (2018/12/18 J. Prchlik)
    ;gimg = gauss_smooth(abs(simg),20*8./rebinv,/edge_trunc,/NAN,kernel=img_kernel)
    ;2D Gaussian function describing the image kernel
    ;Switched to zeroing out edges 2018/12/18 J. Prchlik (it is faster than using gauss_smooth, 2 images per minute compared to 3.5-4)
    ker_val = 20*8./rebinv
    img_kernel = GAUSSIAN_FUNCTION([1,1]*ker_val)
    gimg = CONVOL(abs(simg),img_kernel,total(img_kernel),/edge_trunc,/NAN);,INVALID=255,MISSING=0)
    ;Find the boundaries in the smoothed image
    rad_1 = 1.
    ;rad_2 = 300.
    ;Use the sigmoids measured size +20 pixels to look for features
    ;rad_2 = sig_p+100. 2(half of image width)*8(rebinned pixels)
    rad_2 =  sig_p/(2.*rebinv)-1 ;win_w/(2*rebinv)-1
    ;don't let rad_2 be larger than half of the image
    if rad_2 gt win_w/rebinv/2-1 then rad_2 = win_w/rebinv/2-1



    ;Try manually setting the threshold
    while set_threshold do begin
       

        ;Set the threshold value for edge_dog manually
        READ, thres_val, PROMPT='Threshold Value'
        thres_val = float(thres_val)
     

        ;Don't let threshold value be less than 0
        if thres_val lt 0 then thres_val = 0.01
 

        ;set origin variables correct from center origin to bottom left origin 2018/06/08 J. Prchlik
        org_x = -(cent_x-1-pix_x+win_w/2)*delt_x
        org_y = -(cent_y-1-pix_y+win_w/2)*delt_y

        ;get x,y origin for binned image
        borg_x = -(cent_x/rebinv)*delt_x*rebinv
        borg_y = -(cent_y/rebinv)*delt_y*rebinv
        
        ;Used difference of gaussian to find edges
        ;Use 3 sigma from 0 2018/12/07 J. Prchlik
        edge = edge_dog(abs(gimg),radius1=rad_1,radius2=rad_2,threshold=thres_val,zero_crossings=[0,255])

        ;set threshold value dynamically
        ;Get boundary of created countour
        ;Switched to N levels instead of LEVEL, which just picked the top value, although the image is binary 2018/12/19
        CONTOUR,edge, NLEVEL = 1,  $
               XMARGIN = [0, 0], YMARGIN = [0, 0], $
               /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, $
               XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS;/NODATA

        ;get size of image
        img_size = size(fimg)
        
        loadct,0,/silent
        ;Plot image for testing threshold
        plot_image,fimg,min=-150,max=150,$
            origin=[org_x, $
            org_y], $
            scale=[delt_x,delt_y], $
            xtitle='X-postion (arcseconds)', $
            ytitle='Y-position (arcseconds)',$
            ;Updated title to use fits_read hdr syntax 2018/11/05 J. Prchlik
            title=string([sig_id],format=title_fmt)+fmt_dat,$
            xcharsize =1.20*win_w/700., $
            ycharsize =1.20*win_w/700., $
            xcharthick=1.20*win_w/700., $
            ycharthick=1.20*win_w/700., $
            charsize  =1.00*win_w/700.,$
            charthick =1.00*win_w/700.,Position=[0.2, 0.15, 0.95, 0.90];/nosquare

        ;Get formatted date
        ;fmt_dat = str_replace(str_replace(strmid(sxpar(hdr,'T_OBS'),0,19),'.','-'),'_','T')  
        

        ;Plot restricted range
        ;Rotate coordinate to image time
        ;rot_p = rot_xy(xi,yi,tstart=gi,tend=fmt_dat,offlimb=offlimb)


        ;Get solar radius in arcsec at given time

        print,'##################################################'
        ;Create new ROI obejct using contour
        loadct,12,/silent
        ;Plot rotated sigmoid location
        plots,rot_p[0],rot_p[1],psym=4,symsize=3,thick=3,color=200

        ;Look through ROIs until you find one with the sigmoid point inside
        ind_obj = 0
        search = 1
        while search do begin 
           print,'##################################################'
           print,ind_obj
           line = [LINDGEN(PathInfo(ind_obj).N), 0] 
           ;Draw ROI on plot switch to plotting line around ROI instead of ROI in draw_roi
           plots, rebinv*delt_x*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[0, *]+org_x, rebinv*delt_y*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[1, *]+org_y,color= 200,thick=3

           ;No region found
           if ind_obj eq n_elements(PathInfo.N)-1 then search = 0 

           ;increment counter
           ind_obj = ind_obj+1
        endwhile
        READ, set_threshold, PROMPT='Try a new threshold (1 = Yes, 0 = No): ',format='(I1)'
        print,'##################################################'

    endwhile

    ;Store manual threshold value in an array
    store_thres[i] = thres_val

    ;Manually threshold 2019/02/21 J. Prchlik
    save,ID,store_thres,rebinv,filename=thres_file


endfor

end