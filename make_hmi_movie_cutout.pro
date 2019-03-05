
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
;    make_hmi_movie_cutout,times,hmi_arch='hmi_arch/',out_arch='hmi_movie/',rebinv=16,man_thres=man_thres
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
;    man_thres  -   A keyword to use the manually selected threshold value over the one currently selected by fft
;    
;
;OUTPUTS
;    HMI movie and sav files in hmi_arc
;
;#############################################################

;Running with Rebin = 32 2019/01/08 J. Prchlik
pro make_hmi_movie_cutout,times,hmi_arch=hmi_arch,out_arch=out_arch,rebinv=rebinv,man_thres=man_thres
;set plot to Z Window
set_plot,'Z'

;Start with BW color table
loadct,0
;Read in file containing TBEST
;readcol,times,ID,RATING,NOAA,AR_START,X,Y,AR_END,SIG_START,SIG_END,TBEST,format='LL,I,A,A,F,F,A,A,A,A'
;Updates with Patty's new output format 2018/06/13 J. Prchlik
formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F'
readcol,times,dum,ID_read,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
       length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats,/preserve_null
;Set archive directory for download aia files
if keyword_set(hmi_arch) then hmi_arch = hmi_arch else hmi_arch = 'hmi_arch_cutout/'
hmi_arch = hmi_arch+'/'

;Set archive directory for output png files
if keyword_set(out_arch) then out_arch = out_arch else out_arch = 'hmi_movie_cutout/'
out_arch = out_arch+'/'

if keyword_set(rebinv) then rebinv = rebinv else rebinv = 16


if keyword_set(man_thres) then begin
    ;Name for output save file with the manual threshold values
    ;thres_file = 'manual_edge_dog_threshold.sav'
    ;Switch to text file 2019/02/26 J. Prchlik
    thres_file = 'manual_edge_dog_threshold.txt'
    ;Also updates the rebin value to the value used in make_him_movie_cutout_man_thres
    ;Restores variables ID and store_thres
    ;restore,thres_file
    ;Switch to text file 2019/02/26 J. Prchlik
    readcol,thres_file,thres_id,store_thres,format='A,F'
    ;rename ID to thres_id to avoid any confusion
    ;thres_id = ID
 
;Set man_thres to 0 if keyword is not set
endif else man_thres = 0

;Fix ID variable name to standard ID had to work around in case man_thres is set 2019/02/21 J. Prchlik
ID = ID_read

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
goodt = where(strlen(tbest) eq 23)

;Cadance for image creation in seconds
img_cad = 30.*60.


;format for output png file directory
out_fmt = '(A30,"/")'

;IAU_format for coordinates
iau_cor = '("L",I03,"C",I03)'

;Restore sigmiod save file
;Not needed anymore 2018/07/23 J. Prchlik
;restore,'Sigmoids2007to2017.sav'

;SDO/HMI take over date
sdo_takeover = anytim('2010-06-13T21:48:00')

;Rotate images by given angle
rot_mat = [[-1.,0.],[0.,-1.]]

;Reseting rebinv used for switching from hmi to mdi
rebinv_tmp = rebinv

;Download HMI data for all the best times
for ii=0,n_elements(goodt)-1 do begin


    ;title format
    title_fmt = '("HMI ID: ",A3," ")'

    ;Set index to value with goodt
    i = goodt[ii]
    
    ;Start of observation
    t1 = ar_start[i]
    ;get index for a good time
    ;Now x, y are measured at tobs instead of tbest 2018/06/14 J. Prchlik
    ;Use MDI if before SDO science data date 2019/03/01 J. Prchlik
    if anytim(t1) lt sdo_takeover then begin
        ;tell program it is MDI
        mdi = 1
        gi = tbest[i]
        ;title format
        title_fmt = '("MDI ID: ",A3," ")'

    endif else gi = tobs[i]
    
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
    ;Add mdi to the search 2019/03/01 J. Prchlik
    if mdi eq 1 then match_files = file_search(hmi_arch+'/'+sig_id+"/*mdi*fits",count=file_cnt,/full) $
    else match_files = file_search(hmi_arch+'/'+sig_id+"/*hmi*fits",count=file_cnt,/full)

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


    ;prep hmi data
    ;hmi_prep,match_files,findgen(n_elements(match_files)),index,odata
    ;mreadfits,match_files, index,data
    ;switch to read_fits because cutout header is missing some information which causes mreadfits to be unhappy 2018/11/05 J. Prchlik
    data = readfits(match_files[0],hdr, exten_no=0, /fpack,/silent) 
    ;hmi_prep,iindex,data,index,odata
    ;hmi_prep,hmi_list[chk_i],[1]findgen(n_elements(chk_i)-1),index,data


    ;width of image window in pixels different hmi and mdi
    if mdi then win_w=160 $
    else win_w= 700
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


    ;Do not set threshold automaticall if man_thres is set 2019/02/21 J. Prchlik
    if keyword_set(man_thres) then begin
        ;match index in read sigmoid catalog file to manual thres_hold value
        match_ind = where(thres_id eq ID[i])
        ;Update the threshold value to manually set threshold value
        thres_val = store_thres[match_ind]
        ;If the threshold value is less than 0 (i.e., not set skip this set of sigmoids) 2019/03/01 J. Prchlik
        if thres_val lt 0 then continue

        ;note that the threshold value is already set so no need to compute automatically
        set_threshold = 0

    endif else begin
        ;Set Threshold  2018/12/20 Use to set a static threshold for each observations, which keeps ROI at a consisent size
        set_threshold = 1
    endelse

    ; plot each hmi observation
    for j=0,n_elements(match_files)-1 do begin

        ;Reset rebinv in case scaled by mdi
        rebinv = rebinv_tmp

        ;read in files 1 at a time
        data = readfits(match_files[j],hdr, exten_no=0, /fpack,/silent) 

        ;replace bad pixels
        bad_pix = data le -1.7
        data[bad_pix] = !values.f_nan


        ;if image quality greater than 90000 exit
        ;swtich to fits_read format
        if anytim(t1) gt sdo_takeover then begin
            if sxpar(hdr,'quality') gt 90000 then continue
        endif


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

        ;Get image data to plot
        fimg = data

        ;Get range around pix_x and pix_y values
        ;If larger than size of download just use full window 2019/03/01 J. Prchlik
        size_fimg =size(fimg)
        if ((win_w gt size_fimg[1]) OR (win_w gt size_fimg[2]) OR (mdi eq 1)) then begin
            lims = [0,size_fimg[1]-1,0,size_fimg[2]-1]
            win_w = min(size_fimg(1:2))-1
        endif else lims = select_cutout(pix_x,pix_y,win_w,size_fimg[1],size_fimg[2])
 
        ;If the program could not select a cutout coninue
        bad_lims = where(lims lt 0,bad_lim_cnt)
        if bad_lim_cnt gt 0 then continue

        ;Store limite seperately
        pxmin = lims[0]
        pxmax = lims[1]
        pymin = lims[2]
        pymax = lims[3]


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
           ;Make sure to set color table back to default
           loadct,0
           ;If Error go to the next image
           continue
        ENDIF


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
        ;Try static sigma for MDI 2019/03/04 J. Prchlik
        if mdi eq 1 then $
            img_sig = 25. $
        else $
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
    
        ;If mdi scale down the rebinning 2019/03/01 J. Prchlik
        if mdi then rebinv = rebinv/4    


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


        ;Check the clipped window for rebinning is not smaller than the maximum called index 2019/03/01 J. Prchlik
        if pxmax ge fimg_size[1] then pxmax = fimg_size[1]-1
        if pymax ge fimg_size[2] then pymax = fimg_size[2]-1
 

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
        ;Change scaling for hmi vs. mdi
        if mdi eq 1 then  ker_val = 20*8./rebinv/4 $
        else ker_val = 20*8./rebinv


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



        ;Only do if threshold value for this observation is no set 2018/12/20 J. Prchlik
        if set_threshold then begin
            ;2sigma drop threshold assuming rad_2 is an approximation for sigma
            ;Switch to 5 sigma 2018/12/07 J. Prchlik
            ;sig_cut = 3.0
            ;thres_val = cgpercentiles(abs(gimg),percentiles=.95)*exp(-(sig_cut)^2/2.)
            ;Use a bottom up approach rather than top down 2019/01/08 J. Prchlik
            ;Use computed image sigma multiply that value by the summation of pixels
            ;which gives you the sigma floor value after binnning.
            ;Then multiply the value by the integral of Gaussian Kernel, which we used
            ;to smooth the image 2019/01/09 J. Prchlik (remember img_sig is 2*sigma,
            ; so a 3 sigma detection)
            ;Trying 5 sigma on 2019/01/11
            ;thres_val  = (3./2.)*img_sig*rebinv^2*sqrt(2.*!PI)*ker_val
            ;;gimg_sig = get_sig(gimg)
            ;;;Set threshold value to 1 sigma above the floor
            ;;thres_val = img_sig+min(gimg,/nan)
            ;Tried 1% on 2019/01/08 (also rebin 32). It works okay, but occassionally got too much,
            ;Going to try 5% (roughly 2$\sigma$ if you assume normal dist.) 2019/01/09
            ;thres_val = cgpercentiles(abs(gimg),percentiles=.05);*exp(-(sig_cut)^2/2.)
            ;5% did not work very well 2019/01/09
            ;1% was a bit too inclusive trying 2.5% 2019/01/11 J. Prchlik
            ;thres_val = cgpercentiles(abs(gimg),percentiles=.001);*exp(-(sig_cut)^2/2.)
            ;Use smallest contour level to set min 2019/01/11 J. Prchlik 
            ;DOES NOT WORK and is unpredictable 2019/01/11 J. Prchlik
            ;contour,gimg,path_info=test_info,closed=0,nlevels=100
            ;thres_val = min(test_info.value)
           
            ;Get values of inner pixels to find sigma value for threshold
            size_gimg = size(gimg)
 
            ;Use FFT filter to remove lower level of data 2019/01/11
            fft_gimg = fft(gimg, /center)

            ;Compute the power spectrum of the transform and
            ;apply a log scale.
            pS_gimg = abs(fft_gimg)^2
            lpS_gimg = alog10(pS_gimg)

            ;Set the log power scale maximum to 0
            ;PS0_gimg = lpS_gimg - max(lps_gimg)

            ;Set "noise" threshold
            mask = real_part(lpS_gimg) gt 2.5+min(real_part(lpS_gimg))

            

            ; Apply the mask to the transform to exclude the noise.
            mask_fft_gimg = fft_gimg*mask

            ;invert the transform
            denoised_gimg = real_part(fft(mask_fft_gimg,/inverse,/center))


            ;Set the threshold value to the "noise" floor
            thres_val = min(denoised_gimg,/nan)

            ;Don't let threshold value be less than 0
            if thres_val lt 0 then thres_val = 0.01
 
            ;set min and maximum values of an inner square
            ;Recall the input image will be a square
            ;min_pix = size_gimg[1]/4
            ;max_pix = (size_gimg[1]*3)/4

            ;;Get sigma of inner half of the image and use for the threshold value
            ;thres_val = get_sig(gimg[min_pix:max_pix,min_pix:max_pix])
            

        endif



        ;Get sigma of smoothed image
        ;smt_sig = get_sig(gimg)

        ;Used difference of gaussian to find edges
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
        ;pn_img = gauss_smooth(simg,20*8./rebinv,/edge_zero,/NAN)
        pn_img = convol(simg,img_kernel,total(img_kernel),/edge_trunc,/NAN)
        ;pos_neg = gauss_smooth(simg,20*8./rebinv,/edge_trunc,/NAN)
        ;Trying emboss filter 2018/12/19 J. Prchlik
        ;edge_pn = abs(roberts(pn_img))
        edge_pn = edge_dog(pn_img,threshold=0.,zero_crossings=[0,255])


        ;Get +/- Contour close=0 prevents always enclosing a contour
        CONTOUR,edge_pn,NLEVEL=1, $
            XMARGIN = [0, 0], YMARGIN = [0, 0], $
            /NOERASE,PATH_XY=cont_pn,PATH_INFO=info_pn, $
            XSTYLE=5,YSTYLE=5,/PATH_DATA_COORDS,closed=0
        
  
        ;Get boundary of created countour
        ;Switched to N levels instead of LEVEL, which just picked the top value, although the image is binary 2018/12/19
        CONTOUR,edge, NLEVEL = 1,  $
               XMARGIN = [0, 0], YMARGIN = [0, 0], $
               /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, $
               XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS;/NODATA

        ;Scale up text size for mdi observations
        scale_txt = 1
        if mdi eq 1 then scale_txt = 2
        
        ;Plot image with rotation
        plot_image,img,min=-150,max=150,$
            origin=[org_x, $
            org_y], $
            scale=[delt_x,delt_y], $
            xtitle='X-postion (arcseconds)', $
            ytitle='Y-position (arcseconds)',$
            ;Updated title to use fits_read hdr syntax 2018/11/05 J. Prchlik
            title=string([sig_id],format=title_fmt)+fmt_dat,$
            xcharsize =1.20*scale_txt*win_w/700., $
            ycharsize =1.20*scale_txt*win_w/700., $
            xcharthick=1.20*scale_txt*win_w/700., $
            ycharthick=1.20*scale_txt*win_w/700., $
            charsize  =1.00*scale_txt*win_w/700.,$
            charthick =1.00*scale_txt*win_w/700.,Position=[0.2, 0.15, 0.95, 0.90];/nosquare


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
                 
                
                ;Do not plot PIL it is not a good calculation 2019/02/21 J. Prchlik
                ;plots, plot_x[sort_y],plot_y[sort_y],color= 100,thick=3
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

        ;Change off limb NAN values to 0.0 for summation 2019/01/08
        fimg(where(fimg eq  !values.f_nan)) = 0.0

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
        obs_loca = [[obs_loca],[rot_p[0],rot_p[1]]]
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


       ;If the program get this far then it successfully ran for this observation. Therefore, there is no need to refind the threshold value 
       ;Threfore, keep the threshold static for the remaining observations of this particular sigmoid 2108/12/20 J. Prchlik
       set_threshold = 0


    endfor


    ;Do not save any information if no objects are found 2019/01/11 J. Prchlik
    if obs_time eq !NULL then continue

    ;Resort so middle file is back in sequence
    ; 2019/01/11 J. Prchlik
    real_time = anytim(obs_time)
    sort_time = sort(real_time)

    ;resort all variables as a function of time
    obs_time = obs_time[sort_time]
    obs_loca = obs_loca[*,sort_time]
    obs_qual = obs_qual[sort_time]
    tot_ints = tot_ints[sort_time]
    pos_ints = pos_ints[sort_time]
    neg_ints = neg_ints[sort_time]
    pix_area = pix_area[sort_time]
    tot_area = tot_area[sort_time]
    roi_save = roi_save[sort_time]
    phy_save = phy_save[sort_time]
    pol_lens = pol_lens[sort_time]    


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
   ;Remove Wait 2019/01/08 J. Prchlik
   ; wait,60
    ;Cancel out image from memory. Leaking like a sieve
    index = 0
   data  = 0

   ;stop
endfor

end