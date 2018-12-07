
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
;    make_hmi_movie
;
;PURPOSE
;    Creates hmi movie for each sigmoid
;
;CATEGORY:
;    Program, movie creation
;
;USAGE
;    make_hmi_movie,times,hmi_arch='hmi_arch/',out_arch='hmi_movie/',rebinv=8
;
;INPUTS
;    times      -   A csv file containing times to analyze sigmoid filaments
;                   CSV format must be as follows:
;                   formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F'
;                   readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
;                   length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats
;    hmi_arch   -   Directory to output the HMI files to (Default = 'hmi_arch/')
;    cad        -   Cadence to get HMI files in seconds (Default = 30.*60)
;    out_arch   -   Where to put the output movies and save files (Default = 'hmi_movie_cutout/')
;    rebinv     -   Rebinning pixel value for median smoothing computationally effciently (Default = 8)
;
;OUTPUTS
;    HMI movie and sav files in hmi_arc
;
;#############################################################

pro make_hmi_movie_cutout,times,hmi_arch=hmi_arch,out_arch=out_arch,rebinv=rebinv
;set plot to Z Window
set_plot,'Z'

;Start with BW color table
loadct,0
;Read in file containing TBEST
;readcol,times,ID,RATING,NOAA,AR_START,X,Y,AR_END,SIG_START,SIG_END,TBEST,format='LL,I,A,A,F,F,A,A,A,A'
;Updates with Patty's new output format 2018/06/13 J. Prchlik
formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F'
readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
       length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats,/preserve_null
;Set archive directory for download aia files
if keyword_set(hmi_arch) then hmi_arch = hmi_arch else hmi_arch = 'hmi_arch_cutout/'
hmi_arch = hmi_arch+'/'

;Set archive directory for output png files
if keyword_set(out_arch) then out_arch = out_arch else out_arch = 'hmi_movie_cutout/'
out_arch = out_arch+'/'

if keyword_set(rebinv) then rebinv = rebinv else rebinv = 8


;Solar radius in cm
phy_rad = 6.957E10 ;cm

;No longer need because cutout directory contains the appropriate sigmiod IDs 
;when compared to the final catalog (2018/11/05 J. Prchlik)
;;;;;;Get list of hmi files
;;;;;hmi_list = file_search(hmi_arch+"hmi*fits")
;;;;;;separate file name
;;;;;fnames = strsplit(hmi_list,'/',/extract) 
;;;;;fname = strarr(n_elements(fnames))
;;;;;;Number of columns
;;;;;splnu = n_elements(fnames[0])
;;;;;
;;;;;;Get just filenames
;;;;;for i=0,n_elements(fnames)-1 do fname[i] = fnames[i,splnu-1]
;;;;;
;;;;;;Remove filename start
;;;;;times = str_replace(fname,'hmi.m_45s.','')
;;;;;;Remove fits ending
;;;;;times = str_replace(times,'_TAI.magnetogram.fits','')
;;;;;;Switch _ to :
;;;;;times = str_replace(times,'_',':')
;;;;;;Switch . to -
;;;;;times = str_replace(times,'.','-')
;;;;;
;;;;;;Replace the first : with a T
;;;;;pos = strpos(times,':')
;;;;;strput,times,'T',pos[0]
;;;;;;Replace the first : with a .
;;;;;;pos = strpos(times,':',/Reverse_search)
;;;;;;strput,times,'.',pos[0]
;;;;;
;;;;;
;;;;;;Get hmi times
;;;;;hmi_time = double(anytim(times))

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

;Restore sigmiod save file
;Not needed anymore 2018/07/23 J. Prchlik
;restore,'Sigmoids2007to2017.sav'


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

    ;No longer needed because files are now in directories related to 
    ;there sigmoid IDs not that sigmiod IDs are established
    ;;;;;;Create array of good times at 90 minute cadence
    ;;;;;add_cad  = 1
    ;;;;;counter = 1
    ;;;;;time_arr = [at1]
    ;;;;;
    ;;;;;;loop until cadence maximum is reached
    ;;;;;while add_cad eq 1 do begin 

    ;;;;;    ;Add new cadence to time array 
    ;;;;;    new_time = double(counter*img_cad+at1)
    ;;;;;    time_arr = [time_arr,new_time]


    ;;;;;    ;end if new time greater than endtime
    ;;;;;    if new_time gt at2 then add_cad = 0
    ;;;;;    counter = counter+1
    ;;;;;endwhile

    ;;;;;;Create large time array to find the best times for each cadence
    ;;;;;r_ti = double(time_arr ## (-1+dblarr(n_elements(hmi_time))))
    ;;;;;p_ti = double(hmi_time # (1+dblarr(n_elements(time_arr))))

    ;;;;;;get time min for matches
    ;;;;;minv = min(r_ti+p_ti,min_loc,/abs,dimension=1)
    ;;;;;;get min values less than 90 minutes only
    ;;;;;good_min = where(minv lt img_cad,matches)

    ;;;;;;leave if there are no good matches
    ;;;;;if matches eq 0 then continue



    ;;;;;;clip to only get closest matches
    ;;;;;min_loc = min_loc[good_min]
   
    ;;;;;;convert 1D indices into 2D indices
    ;;;;;col_i = array_indices(r_ti,min_loc)

    ;;;;;;get the file indices that match the time
    ;;;;;chk_i = col_i[0,*]


    ;;;;;;Get matched files
    ;;;;;match_files=hmi_list[chk_i]
    ;;;;;match_fname = fname[chk_i]

    ;;;;;;Get only unique values
    ;;;;;match_files=match_files[uniq(match_files)]
    ;;;;;match_fname=match_fname[uniq(match_files)]
    

    ;Switched to Sigmoid ID directory now that sigmoid IDs
    ;are frozen in (MAKE SURE YOU FORCE EVERYONE TO AGREE TO
    ;IDS first. It will save you a lot of headaches (2018/11/05)
    ;Check the first file before reading in all files
    ;;;;mreadfits,match_files, index,data

    ;;;;;Create directory for output png files
    ;;;;;Subscribe to IAU standard on output format
    ;;;;iau_time = strsplit(gi,'.',/extract)
    ;;;;iau_time = iau_time[0]

    ;;;;;Get Carrington coordinates 
    ;;;;wcs = fitshead2wcs( index(0) )
    ;;;;rot_p = rot_xy(xi,yi,tstart=gi,tend=index(0).date_d$obs)
    ;;;;WCS_CONV_HPC_HG, rot_p[0], rot_p[1], lon, lat, WCS=WCS, /carr, /pos_long
    ;;;;iau_pos = string([round(lon),round(lat)],format=iau_cor)

    ;;;;;Create new unique ID with IAU standard
    ;;;;sig_id = 'SOL'+iau_time+iau_pos
    ;;;;full_dir = out_arch+string([sig_id],format=out_fmt)
    ;;;;;Remove : characters
    ;;;;full_dir = str_replace(full_dir,':','')
    ;;;;full_dir = str_replace(full_dir,'-','')


    ;Setup output directories
    sig_id = strcompress(ID[i],/remove_all)
    full_dir = out_arch+sig_id+'/'
    if file_test(full_dir) eq 0 then file_mkdir,full_dir
   
    ;If save file exists continue 2018/06/20
    if file_test(full_dir+'/'+str_replace(sig_id,':','')+'.sav') eq 1 then continue

    ;get list of local files
    match_files = file_search(hmi_arch+'/'+sig_id+"/*hmi*fits",count=file_cnt,/full)

    ;if no files found then continue
    if file_cnt eq 0 then continue
  
    ;prep hmi data
    ;hmi_prep,match_files,findgen(n_elements(match_files)),index,odata
    ;mreadfits,match_files, index,data
    ;switch to read_fits because cutout header is missing some information which causes mreadfits to be unhappy 2018/11/05 J. Prchlik
    data = readfits(match_files[0],hdr, exten_no=0, /fpack,/silent) 
    ;hmi_prep,iindex,data,index,odata
    ;hmi_prep,hmi_list[chk_i],[1]findgen(n_elements(chk_i)-1),index,data


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
    device,set_resolution=[win_w,win_w],decomposed=0,set_pixel_depth=24
    


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

    ; plot each hmi observation
    for j=0,n_elements(match_files)-1 do begin

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
        CATCH, Error_status

        ;if data is empty continue
;       ;This statement begins the error handler:
        IF Error_status NE 0 THEN BEGIN
           PRINT, 'Error index: ', Error_status
           PRINT, 'Error message: ', !ERROR_STATE.MSG
           ; Handle the error by extending A:
           CATCH, /CANCEL
           ;If Error go to the next image
           continue
        ENDIF

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
        ;cormag_dp_p0=nospike(cormag_p0,thre=0.9,bright=0.99,imap=imap_p0)
        cormag_dp_p0=nospike(cormag_p0,thre=0.65,bright=0.99,imap=imap_p0,/silent)
        ;cormag_dp_p0=nospike(cormag_p0,thre=0.9,bright=0.85,imap=imap_p0)
        cormag_dp_ni =nospike(cormag_n,thre=0.9,bright=0.50,imap=imap_n,/silent)
        ;cormag_dp_ni =nospike(cormag_n,thre=0.9,bright=0.85,imap=imap_n)

        ;Remove spikes from image
        ;img(where(imap_p0 eq 1))=0
        ;img(Where(imap_n eq 1))=0
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

        ;correct pixels by cos(theta)^2
        fimg = fimg/cos(bin_cor_t)^2
        ;zero out all pixels greater than 50 degrees
        ;Only did for exper. Not in final draft 2018/12/03 J. Prchlik
        ;fimg(where(bin_cor_t gt !dtor*50.)) = 0.0

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
        gimg = gauss_smooth(abs(simg),20*8./rebinv,/edge_trunc,/NAN)
        ;Find the boundaries in the smoothed image
        rad_1 = 1.
        ;rad_2 = 300.
        ;Use the sigmoids measured size +20 pixels to look for features
        ;rad_2 = sig_p+100. 2(half of image width)*8(rebinned pixels)
        rad_2 =  sig_p/(2.*rebinv)-1 ;win_w/(2*rebinv)-1
        ;don't let rad_2 be larger than half of the image
        if rad_2 gt win_w/rebinv/2-1 then rad_2 = win_w/rebinv/2-1

        ;2sigma drop threshold assuming rad_2 is an approximation for sigma
        ;Switch to 5 sigma 2018/12/07 J. Prchlik
        sig_cut = 5.0
        thres_val = cgpercentiles(abs(gimg),percentiles=.95)*exp(-(sig_cut)^2/2.)



        ;Get sigma of smoothed image
        ;smt_sig = get_sig(gimg)

        ;Used difference of gaussian to find edges
        ;Use 3 sigma from 0 2018/12/07 J. Prchlik
        edge = edge_dog(abs(gimg),radius1=rad_1,radius2=rad_2,threshold=thres_val,zero_crossings=[0,255])
        ;gimg = abs(simg)
        ;Make gimg back to normal size then cut
        ;nothing is found in image just continue
        if max(edge) lt 1 then continue

        ;Cut image to restricted window
        img = fimg(pxmin:pxmax,pymin:pymax)

        ;Make sure image is the correct size if no loop off left pixels
        img_size = size(img)



        ;set origin variables correct from center origin to bottom left origin 2018/06/08 J. Prchlik
        org_x = -(cent_x-1-pix_x+win_w/2)*delt_x
        org_y = -(cent_y-1-pix_y+win_w/2)*delt_y

        ;get x,y origin for binned image
        borg_x = -(cent_x/rebinv)*delt_x*rebinv
        borg_y = -(cent_y/rebinv)*delt_y*rebinv



        ;Try to calculate polarity inversion line 2018/12/07
        ;gaussian smooth image
        pn_img = gauss_smooth(simg,20*8./rebinv,/edge_trunc,/NAN)
        ;pos_neg = gauss_smooth(simg,20*8./rebinv,/edge_trunc,/NAN)
        edge_pn = edge_dog(pn_img,radius1=rad_1,radius2=rad_2,threshold=0.,zero_crossings=[0,255])


        ;Get +/- Contour close=0 prevents always enclosing a contour
        CONTOUR,edge_pn,LEVEL=1, $
            XMARGIN = [0, 0], YMARGIN = [0, 0], $
            /NOERASE,PATH_XY=cont_pn,PATH_INFO=info_pn, $
            XSTYLE=5,YSTYLE=5,/PATH_DATA_COORDS,closed=0
        
  
        ;Get boundary of created countour
        CONTOUR,edge, LEVEL = 1,  $
               XMARGIN = [0, 0], YMARGIN = [0, 0], $
               /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, $
               XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS;/NODATA

        
        ;Plot image with rotation
        plot_image,img,min=-150,max=150,$
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


        ;Create new ROI obejct using contour
        loadct,12,/silent
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
               rebinv*delt_x*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[0, *]+borg_x, $
               rebinv*delt_y*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[1, *]+borg_y) 
               ;(pathXY(*, pathInfo(0).OFFSET +line ))[0, *], b$
               ;(pathXY(*, pathInfo(0).OFFSET +line ))[1, *]) b& $

           ;scale roi
           ;Not used 2018/06/21 J. Prchlik
           ;roi_scl = roi_phy
           ;roi_scl -> Scale,1.2, 1.2
 
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
              rebinv*delt_x*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[0, *]+borg_x, $
              rebinv*delt_y*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[1, *]+borg_y) & $
              ;(pathXY(*, pathInfo(0).OFFSET +line ))[0, *], $
          
           ; write_png,full_dir+str_replace(match_fname[j],'fits','png'),tvrd(/true),r,g,b
        endif

        

        ;ROI in pixel coordinates
        roi_obj = OBJ_NEW('IDLanROI', $
           rebinv*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[0, *], $
           rebinv*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[1, *]) & $


      

        ;ROI after fittinf for an ellipse
        ;Fit ellipse to roi object
        ;Create ROI in edge coordinate system
        ;; Ellispe is not a good fit 2018/11/07 J. Prchlik
        ;;;roi_img = OBJ_NEW('IDLanROI',(pathXY(*, pathInfo(ind_obj).OFFSET+line))[0,*],(pathXY(*, pathInfo(ind_obj).OFFSET+line))[1,*])
        ;;;;Get where pixels are inside ROIobjec
        ;;;mask = roi_img -> ComputeMask(Dimensions=Size(edge, /Dimensions), Mask_Rule=2)   
        ;;;;Create image where only those pixels have values of 255
        ;;;roiImage  = edge * (mask GT 0)       
        ;;;;Get index location where values are 255 (i.e. the ROIobject)
        ;;;indices = where(roiImage eq 255)
        ;;;;Fit and ellise to the coordinates in the edge image
        ;;;edge_size = size(edge)
        ;;;elp_xy    = fit_ellipse(indices,xsize=edge_size[1],ysize=edge_size[2]) 


        ;;;;create ellipse ROI in physical coordinates
        ;;;roi_elp = OBJ_NEW('IDLanROI', $
        ;;;    rebinv*delt_x*(elp_xy[0, *])+borg_x, $
        ;;;    rebinv*delt_y*(elp_xy[1, *])+borg_y) 
        ;;;;create ellipse ROI in image coordinates
        ;;;roi_ilp = OBJ_NEW('IDLanROI', $
        ;;;    rebinv*(elp_xy[0, *]), $
        ;;;    rebinv*(elp_xy[1, *])) 

        ;Scale up object roi 2018/06/14
        ;Thought about but don't do too much could go wrong
        ;roi_obj -> Scale,1.6, 1.6

        ;DRAW_ROI, roi_obj, COLOR =200,/LINE_FILL
        ;Draw ROI on plot switch to plotting line around ROI instead of ROI in draw_roi
        plots, rebinv*delt_x*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[0, *]+borg_x, rebinv*delt_y*(pathXY(*, pathInfo(ind_obj).OFFSET + line))[1, *]+borg_y,color= 200,thick=3

        ;Do the same thing to find the polarity inversion line
        pn_search = 1
        ind_pn = 0
        pn_dists = []

        ;Set up length of polarity inversion line
        pol_len = 0
 
        ;Look for polarity inversion line contours inside the ROI object
        while pn_search do begin 
            pn_line = [LINDGEN(info_pn(ind_pn ).N), 0] 
            ;Pixel coordinates of the polarity invserion line
            pn_x =  rebinv*(cont_pn(*, info_pn(ind_pn).OFFSET + line))[0, *]
            pn_y =  rebinv*(cont_pn(*, info_pn(ind_pn).OFFSET + line))[1, *]


            ;Find values of the polarity inversion line inside the ROI
            pn_in = roi_obj -> ContainsPoints(pn_x,pn_y)


            ;Get interior points
            int_pnt = where(pn_in eq 1,pn_cnt)

            ;If there are interior points plot them and store the length
            ; Need at least 2 but setting a min of 4
            if pn_cnt gt 3 then begin
                ;Plot positve negative contour 2018/12/07 J. prchlik
                plot_x = delt_x*pn_x[int_pnt]+borg_x
                plot_y = delt_y*pn_y[int_pnt]+borg_y
                ;Get the polarity inversion line length if IDL find more than one add to the previous length
                pol_len = sqrt(ts_diff(plot_x,1)^2+ts_diff(plot_y,1)^2)+pol_len
                ; plot by value with the larger dynamic range
                if max(plot_x)-min(plot_x) gt max(plot_y)-min(plot_y) then $
                    sort_y = sort(plot_x) $
                else  sort_y = sort(plot_y)
                 
                plots, plot_x[sort_y],plot_y[sort_y],color= 100,thick=3
            endif
            ;increment counter
            ind_pn  = ind_pn +1
            ;No region found
            if ind_pn  eq n_elements(info_pn.N)-1 then pn_search = 0 
        endwhile

        ;plots,rebinv*delt_x*(elp_xy[0, *])+borg_x,rebinv*delt_y*(elp_xy[1, *])+borg_y,color=100,thick=3
        plots,rot_p[0],rot_p[1],psym=4,symsize=3,thick=3,color=200
        loadct,0,/silent


        ;Updated format for output file name  
        png_output = strsplit(match_files[j],'/',/extract)
        png_output = png_output[n_elements(png_output)-1]
        png_output = str_replace(png_output,'fits','png')

        ;write png file in directory 
        TVLCT,r,g,b,/Get
        write_png,full_dir+png_output,tvrd(/true),r,g,b,xresolution=600,yresolution=600

        ;create mask for image based on edge detection
        ;Switched to img_size because if you can't beat em join em
        ;Switched to fimg since fimg and simg/gimg have the same FoV
        maskResult = roi_obj -> ComputeMask(DIMENSIONS = [fimg_size[1],fimg_size[2]])            
        IMAGE_STATISTICS, abs(fimg), MASK = maskResult, $  
                        COUNT = maskArea , data_sum=tot_intensity   
        IMAGE_STATISTICS, fimg < 0., MASK = maskResult, $  
                        COUNT = neg_maskArea , data_sum=neg_intensity   
        IMAGE_STATISTICS, fimg > 0., MASK = maskResult, $  
                        COUNT = pos_maskArea , data_sum=pos_intensity   


        ;Save variables
        obs_loca = [obs_loca,[rot_p[0],rot_p[1]]]
        obs_time = [obs_time,fmt_dat] ; observation time
        obs_qual = [obs_qual,sxpar(hdr,'quality')] ; observation time
        tot_ints = [tot_ints,tot_intensity] ; Total magnetic field intensity
        pos_ints = [pos_ints,pos_intensity] ; Positive magnetic field intensity
        neg_ints = [neg_ints,neg_intensity] ; Negative magnetic field intensity
        pix_area = [pix_area,maskArea] ; Area of ROI in pixels^2
        tot_area = [tot_area,maskArea*are_pix] ; Area of ROI in cm^2
        roi_save = [roi_save,roi_obj] ;ROI object in pixels
        phy_save = [phy_save,roi_phy] ;ROI object in physical coordinates
        pol_lens = [pol_lens,pol_len] ; length of the polarity inversion line Add 2018/12/07 J. Prchlik


        ;Cancel out image from memory. Leaking like a sieve
        img = 0
        fimg = 0
        simg = 0 
        bin_cor_x = 0
        bin_cor_y = 0
        bin_cor_r = 0
        roi_obj = 0
        roi_phy = 0
        roi_elp = 0
        roi_ilp = 0
        roi_img = 0


    endfor


    ;Save output sav file
    out_id = id[i]
    save,sig_id,out_id,obs_time,obs_loca,obs_qual,tot_ints,pos_ints,neg_ints,pix_area,tot_area,$
        roi_save,phy_save,pol_lens,filename=full_dir+'/'+str_replace(sig_id,':','')+'.sav'

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

    ;png_size = strcompress(string([win_w,win_w],format='(I6,"x",I6)'),/remove_all)
    ;Use static movie size 2018/12/03 J. Prchlik
    png_size = '700x700'
    framerate= strcompress(string(8),/remove_all) ;frames per second
    bitrate = ((24*win_w)^2)*framerate/1.e6 ; number of bits per frame 24 bit colors and win_w^2 image
    bit = trim(round(bitrate))+'m'
    bit = '8962k'
    ;Use aia_ffmpeg for simplicity
    call1 = '-y -f image2 -r '+framerate+' -i ' 
    call2 = '-pix_fmt "yuv420p" -vcodec libx264 -level 41 -crf 18.0 -b '+bit+' -r '+framerate+' '+ $
                    '-bufsize '+bit+' -maxrate '+bit+' -g '+framerate+' -coder 1 -profile main -preset faster ' + $
                    '-qdiff 4 -qcomp 0.7 -flags +loop+mv4 -partitions ' + $
                    '+parti4x4+partp8x8+partb8x8 -subq 7 -me_range 16 -keyint_min 1 -sc_threshold ' + $
                    '40 -i_qfactor 0.71 -rc_eq ''blurCplx^(1-qComp)'' -s '+png_size+' -b_strategy 1 ' + $
                    '-bidir_refine 1 -refs 6 -trellis 1 -x264opts ' + $
                    'keyint='+framerate+':min-keyint=1:bframes=1 -threads 2 '
                
    
    
    ;Continue in the background after spawning process
    keep_going = ' &'
    ;output file name
    outf = str_replace(sig_id,':','')+"_mag.mp4"
    outf = str_replace(outf,'-','')
    spawn, ffmpeg +' '+ call1 + full_dir+'symlinks/%4d.png'+' ' + call2 + full_dir+outf, result, errResult
    wait,60
   ;Cancel out image from memory. Leaking like a sieve
   index = 0
   data  = 0

endfor

end