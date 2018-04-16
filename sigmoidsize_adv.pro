
;; Use sswidl
;; This version only looks for data in the Scratch disk

;###############################################
;PURPOSE
;    Checks whether input is a string
;
;INPUT
;    input    -   ideally a string
;
;USAGE
;    is_str = checkstring(input)
;
;###############################################
function checkstring,input

    good = 0
    if ((fix(input) eq 0) or (fix(input)) eq 1) then good = 1 

    return,good
end

function set_radii,inp_img,init_rad1=init_rad1,init_rad2=init_rad2
;guesses best input radii and reduce size until area is continous

   print,'Not yet implimented'

end

;
;Compute a series of slopes and intercepts the contain the simgoid
;Usage
;res = comp_limits(inx,iny)
function comp_limits,inx,iny

lx1 = linfit(iny[0:1],inx[0:1])
lx2 = linfit(iny[2:3],inx[2:3])
ly1 = linfit(inx[1:2],iny[1:2])
ly2 = linfit(inx[3:4],iny[3:4])

return,[lx1,lx2,ly1,ly2]
end

;###############################################
;
;PURPOSE
;    Function to calculate the maximum distance and index locations of the maximum difference
;
;INPUT
;    ind     -  a 2D array of indices containing the location of the sigmoid
;    idat    -  a 2D array containing the image for analysis
;    xbox      -       x-coordinates containing the sigmoid (pixels) 
;    ybox      -       y-coordinates containint the sigmoid (pixels)
;
;Usage
;    outp = brute_force_max_dis(inds)
;###############################################
function brute_force_max_dis,inds,idat,xbox,ybox

;break the indices into x and y values
ind_x = inds[0,*]
ind_y = inds[1,*]


;set up variables for maximum distance
max_d  = 0.
max_x1 = 0
max_x2 = 0
max_y1 = 0
max_y2 = 0

;Compute x and y limit functions
lims = comp_limits(xbox,ybox)

;Create array of 1 and 0 for box
xmin = ind_x ge lims[1]*ind_y+lims[0]
xmax = ind_x le lims[3]*ind_y+lims[2]
ymin = ind_y le lims[5]*ind_x+lims[4]
ymax = ind_y ge lims[7]*ind_x+lims[6]
vbox = xmin*xmax*ymin*ymax

;fit a line through all sigmoid points
use_fit = where(vbox eq 1) 
sig_lin = linfit(ind_x[use_fit],ind_y[use_fit],measure_errors=1./idat[use_fit],sigma=sig_unc)

;Loop over all x-values to find the maximum distance for s
for i=0,n_elements(ind_x)-1 do begin

   ;leave right away if given point is outside box
    if vbox[i] eq 0 then continue
   ;temporary diffrance array
   ;Add box cut 2018/02/22 J. Prchlik
   ;difference in x-value
   x_diff = (ind_x[i]-ind_x)
   ;difference in y-value
   y_diff = (ind_y[i]-ind_y)

   ;makesure distance has slope sim. to line regression slope
   s_diff = y_diff/x_diff
   ;only keep distances within 3 sigma of slope
   cm = abs(s_diff-sig_lin[1]) lt 0.75

   print,sig_lin[1],3.*sig_unc[1],mean(s_diff[use_fit]),max(cm),min(cm)
   t_diff = sqrt((x_diff)^2+(y_diff)^2)*vbox*cm

   ;maximum difference for this point and the index location of maximum
   m_diff = max(t_diff)
   i_diff = where(t_diff eq m_diff)

   ;if maximum difference is greater than the current maximum set new points
   if m_diff gt max_d then begin

       ;update maximum differnce
       max_d = m_diff

       ;update initial axis points
       max_x1 = ind_x[i]
       max_y1 = ind_y[i]
       ;update second axis points
       max_x2 = ind_x[i_diff]
       max_y2 = ind_y[i_diff]

   endif

endfor


return,[max_d,max_x1,max_x2,max_y1,max_y2]
end


;###############################################
;
;PURPOSE
;    Function to calculate the maximum distance and index locations of the minimum axis difference
;
;INPUT
;    ind     -  a 2D array of indices containing the location of the sigmoid
;    idat    -  a 2D array containing the image for analysis
;    xbox      -       x-coordinates containing the sigmoid (pixels) 
;    ybox      -       y-coordinates containint the sigmoid (pixels)
;
;Usage
;    outp = brute_force_min_dis(inds,idat,xbox,ybox)
;
;###############################################
function brute_force_min_dis,inds,idat,xbox,ybox

;get x and y indices
ind_x = inds[0,*]
ind_y = inds[1,*]

;set up variables for maximum distance
max_da  = 0.
max_db  = 0.
max_x1a = 0
max_x2a = 0
max_y1a = 0
max_y2a = 0
max_x1b = 0
max_x2b = 0
max_y1b = 0
max_y2b = 0

;Compute x and y limit functions
lims = comp_limits(xbox,ybox)

;Create array of 1 and 0 for box
xmin = ind_x ge lims[1]*ind_y+lims[0]
xmax = ind_x le lims[3]*ind_y+lims[2]
ymin = ind_y le lims[5]*ind_x+lims[4]
ymax = ind_y ge lims[7]*ind_x+lims[6]
vbox = xmin*xmax*ymin*ymax

;fit a line through all sigmoid points
use_fit = where(vbox eq 1) 
sig_lin = linfit(ind_x[use_fit],ind_y[use_fit],measure_errors=1./idat[use_fit])

;get transformed x and y values
trn_y = ind_y-(sig_lin[1]*ind_x+sig_lin[0])
trn_x = ind_x-mean(ind_x[use_fit])

;use inverse of the slope to get tolerance on x pixel deviation allowed
;However set up limits
x_tol = 1./sig_lin[1]
if x_tol lt  1 then x_tol = 1
if x_tol gt 15 then x_tol = 15

;Find the min and maximum values for the sigmoid ends
;Use percentile range to store the results
tmin=cgPercentiles(trn_x[use_fit],percentiles=.10)
tmax=cgPercentiles(trn_x[use_fit],percentiles=.90)


;Loop over all x-values to find the maximum distance for each small axis
for i=0,n_elements(ind_x)-1 do begin

    ;leave right away if given point is outside box
    if vbox[i] eq 0 then continue

    ;find values within x tolerance
    x_pos = abs(ind_x-ind_x[i]) lt x_tol


    ;temporary diffrance array
    t_diff = abs(trn_y-trn_y[i])*vbox*x_pos

    ;maximum difference for this point and the index location of maximum
    m_diff = max(t_diff)
    i_diff = where(t_diff eq m_diff)

    ;if maximum difference is greater than the current maximum set new points
    case 1 of

        (m_diff gt max_da) and (trn_x[i] lt tmin):  begin
            ;update maximum differnce
            max_da = m_diff
            ;update initial axis points
            max_x1a = ind_x[i]
            max_y1a = ind_y[i]
            ;update second axis points
            max_x2a = ind_x[i_diff]
            max_y2a = ind_y[i_diff]
            end
        (m_diff gt max_db) and (trn_x[i] gt tmax):  begin
            ;update maximum differnce
            max_db = m_diff
            ;update initial axis points
            max_x1b = ind_x[i]
            max_y1b = ind_y[i]
            ;update second axis points
            max_x2b = ind_x[i_diff]
            max_y2b = ind_y[i_diff]
            end

         ;do whatever to just continue
         else: v=1
    endcase

endfor


return,[max_da,max_x1a,max_x2a,max_y1a,max_y2a,max_db,max_x1b,max_x2b,max_y1b,max_y2b]
end

;###############################################
;PURPOSE
;    Loop to find the best Edge values
;     Don't use 2018/02/22 J. Prchlik
;
;USAGE
;    result = loop_edge_dog(dat,radius1=radius1,radius2=radius2,threshold=threshold,zero_crossings=zero_crossings)
;
;INPUT 
;;
;###############################################
function loop_edge_dog,dat,radius1=radius1,radius2=radius2,threshold=threshold,zero_crossings=zero_crossings

counter = 0
looper = 1
uponly = 0

while looper gt 0.6 do begin

    ;differen radius tolerence 
    tol = 5

    ;get +/- some resolution in radius for radius 1 make sure radius is always positive
    tradius1 = 0.5*counter*(-1)^counter+radius1
    if ((tradius1 gt 0)  and (uponly eq 0)) then begin 
        radius1 = tradius1 
    endif else begin
        uponly = 1
        radius1 = 0.5+radius1
    endelse

    ;calculate edge dog
    result = edge_dog(dat,radius1=radius1,radius2=radius2,threshold=threshold,zero_crossings=zero_crossings)                 

   
    ;index location of the sigmoid
    sig = where(result gt 254.5,cnts)
    ;create index array for the entire images for the sigmoid
    if cnts gt 0 then ind_loc = array_indices(result,sig) else continue

    xs = transpose(ind_loc[0,*])
    ys = transpose(ind_loc[1,*])

    sortx = sort(xs)

    fxs = xs[sortx]
    fys = ys[sortx]

    ds = sqrt((ts_diff(fxs,1))^2+ts_diff(fys,1)^2)
    if ((mean(ts_diff(ds,1)) lt tol) ) or (counter gt 100) then looper = 0   
    counter = counter+1
endwhile
return,result
end


;###############################################
;PURPOSE
;      Function to get the FWHM of the sigmoid
;
;USAGE
;    fwhm =  real_fwhm(img,xbox,ybox,ax1,ay1,ax2,ay2,sc_x,sc_y,da_x,da_y,bkgd)
;
;INPUT
;     img       -       Input 2d Image with sigmoid
;     xbox      -       x-coordinates containing the sigmoid (pixels) 
;     ybox      -       y-coordinates containint the sigmoid (pixels)
;     ax1       -       trailing fwhm x-point
;     ay1       -       trailing fwhm y-point
;     ax2       -       leading fwhm x-point
;     ay2       -       leading fwhm y-point
;     sc_x      -       device coorinates to pixels x-axis
;     sc_y      -       device coorinates to pixels y-axis
;     da_x      -       arcseconds to pixels x-axis (float)
;     da_y      -       arcseconds to pixels y-axis (float)
;     bkgd      -       background level in #/s/arcsec^2 (float)
;###############################################

function real_fwhm,img,xbox,ybox,ax1,ay1,ax2,ay2,sc_x,sc_y,da_x,da_y,bkgd

sig = where(img gt -9e31)
;create index array for the entire images for the sigmoid
ind_loc = array_indices(img,sig)

;get image shape
shp_img = size(img)

;break into x and y indices
ind_x = ind_loc[0,*]
ind_y = ind_loc[1,*]

;Compute x and y limit functions
lims = comp_limits(xbox,ybox)

;Create array of 1 and 0 for box
xmin = ind_x ge lims[1]*ind_y+lims[0]
xmax = ind_x le lims[3]*ind_y+lims[2]
ymin = ind_y le lims[5]*ind_x+lims[4]
ymax = ind_y ge lims[7]*ind_x+lims[6]

;get slope of FWHM line
fwhm_line = linfit([ax1,ax2]/sc_x,[ay1,ay2]/sc_y)
perp_slop = -1./fwhm_line[1]
perp_int1 = ay1/sc_y-perp_slop*ax1/sc_x
perp_int2 = ay2/sc_y-perp_slop*ax2/sc_x

;get interior limits
imin = ind_y ge perp_slop*ind_x+perp_int1
imax = ind_y le perp_slop*ind_x+perp_int2

;create large array of interior slopes
vbox = xmin*xmax*ymin*ymax*imax*imin

;create new image with stuff outside the box removed
ibox = reform(vbox,shp_img[1],shp_img[2])
;zero out information outside box
m_img = img*vbox

;rotate the angle of the long axis
;first get the rotation angle
rot_rad = atan(ay2-ay1,ax2-ax1)
rot_deg = rot_rad*180./!PI
;rotate image by that angle
rot_img = rot(m_img,rot_deg,1.0,ax1/sc_x,ay1/sc_y)
;rotate full image by that angle
rot_fmg = rot(img,rot_deg,1.0,ax1/sc_x,ay1/sc_y)
;rotate the mask by that angle
rot_ibox = rot(ibox,rot_deg,1.0,ax1/sc_x,ay1/sc_y)



;create rotation maxtrix
rot_mat = [ [ cos(rot_rad), sin(rot_rad)],$
            [-sin(rot_rad), cos(rot_rad)]]
;create new delta keyword parameters from rotation
ccd = [[0.,da_x],$ ; x-coordinates from pixel to arcsec
       [0.,da_y]]  ; y-coordinates from pixel to arcsec

;new coordinates of after the rotation
;Check program test_rotation.pro to understand why this is not need 
; Because IDL does not rotate that way
; standard way only works with square binning
;ccd_new = ccd # rot_mat

;Sum image a long axis
sum_img = total(rot_img,1)
;Sum good pixels along long axis and convert to square arcsec
sum_pix = total(rot_ibox,1)*(da_x*da_y)


;normalize the sum image by pixel area
sum_img = sum_img/sum_pix

;Set sum values to 0 outsize box region
good = finite(sum_img)
lev_img = sum_img-bkgd
max_img = max(lev_img,/nan)
max_arg = where(lev_img eq max_img,cnt_max)
;Indices of the sum
ind_sum = fix(dindgen(n_elements(sum_img)))

;Get the half maximum width value
;first get value of half max
hlf_max = max_img*0.5
;Split into array above and below the maximum
upp_hlf = ind_sum gt mean(max_arg)
low_hlf = ind_sum lt mean(max_arg)

;get the index closest to the half maximum on either side of the max
upp_fun = abs(sum_img*upp_hlf-hlf_max)
upp_ind = where(upp_fun eq min(upp_fun,/nan))
low_fun = abs(sum_img*low_hlf-hlf_max)
low_ind = where(low_fun eq min(low_fun,/nan))

;check sizes of array for more than one index and if so get the average index value
upp_sze = size(upp_ind)
low_sze = size(low_ind)



;Get only 1 center point
if upp_sze[2] gt 1 then upp_ind = mean(upp_ind) else upp_ind = fix(upp_ind[0])
if low_sze[2] gt 1 then low_ind = mean(low_ind) else low_ind = fix(low_ind[0])


;store fwhm value
fwhm = (upp_ind-low_ind)*da_y

;Use the first maximum value
if cnt_max gt 1 then max_arg = fix(max_arg[0]) else max_arg= fix(max_arg[0])

;Get a grid of y-values (i.e. perpendicualr to the Central Sigmoid axis)
;Assume square pixels when binning 2018/02/28 J. Prchlik
xgrid = (findgen(n_elements(sum_img))-max_arg)*da_y ; scale the index by the new delta-y coordinate system

;Plot line core in new window
isize = size(rot_img)
xsize = isize[1]*2
ysize = isize[2]*2

;Value to scale down to if larger than 512x512
bigx = xsize/512.
bigy = ysize/512.

;Scale down tow prevent images bigger than 1024x1024
if xsize gt 512 then xsize = fix(temporary(xsize)/bigx)
if ysize gt 512 then ysize = fix(temporary(ysize)/bigy)

useg = where(finite(alog10(rot_img)))

;Add factor of 4 to really accent selected region
scmin=cgPercentiles(alog10(rot_fmg[useg]+4.*rot_img[useg]),percentiles=.010)
scmax=cgPercentiles(alog10(rot_fmg[useg]+4.*rot_img[useg]),percentiles=.999)
;Plot image
window,2,retain=0,xsize=xsize,ysize=ysize,xpos=0,ypos=1200
;tv,bytscl(rebin(alog10(rot_img),xsize,ysize),min=scmin,max=scmax)
tv,bytscl(rebin(alog10(rot_fmg+4.*rot_img),xsize,ysize),min=scmin,max=scmax)

;Plot rotated image and  Histogram of normalixed counts
window,6,retain=0,xsize=xsize,ysize=ysize,xpos=xsize,ypos=1200
plot,sum_img,xgrid,/device,color=255,ytitle='Distance from Center [``]',xtitle='Counts [#/s/arcsec^2]'
oplot,fltarr(n_elements(xgrid))+0.5*max_img+bkgd,xgrid,color=200
oplot,fltarr(n_elements(xgrid))+bkgd,xgrid,color=200,linestyle=2
oplot,[0,max_img],([upp_ind,upp_ind]-max_arg)*da_y,color=170,linestyle=1
oplot,[0,max_img],([low_ind,low_ind]-max_arg)*da_y,color=170,linestyle=1


return,[fwhm,max_img]
end

;###############################################
;PURPOSE
;    Determine whether an input quantity is valid
;
;INPUT
;    parvl    -  Input quantity (float)
;    error    -  Whether the values is good (1 = good, 0 = bad) 
;    qname    -  The name of the parameter you want to check (string) 
;
;USAGE
;
;###############################################
pro input_quantity,parvl,error=error,qname=qname,qmin=qmin

    error = 0

    print,'Enter paramater for ',qname
    read,parvl

    parvl = float(parvl)

    if parvl lt qmin then error=1

end

;
;NAME:
;    sigmoidsize_adv
;
;PURPOSE
;    Program to compute sigmoid properties of image
;
;CATEGORY:
;    Program, measurements
;
;USAGE
;    sigmoidsize_ave,dir='exammples/'
;
;INPUTS
;    dir        -   Directory containing Sigmoid fits files
;
;OUTPUTS
;    A save file called called sigmoid_sizedata.sav in dir
;    N.B. Currently clobbers existing save file
pro sigmoidsize_adv,dir=dir
;set up plot as X plot
set_plot,'X'
;scratch_path='/Volumes/Scratch/Users/ehanson/XRT_fits/'
if keyword_set(dir) then scratch_path = dir else  scratch_path='examples/'
fits_files=find_file(scratch_path+'*.fits')
nfiles=n_elements(fits_files)
if (nfiles eq 1) then begin
  if (fits_files eq '') then begin
     message,'No .fits files in the path: '+scratch_path,/cont
     message,'Please verify the path and try again.'
  endif
endif
mreadfits,fits_files,index0
;help,index0,/st
;stop

sigdat_mod={sig_id:'',           $
        NOAA_id:0,               $
        filename:'',             $
        date:'',                 $
        size:0.0,                $
        sizea:0.0,               $
        sizeb:0.0,               $
        aspect_ratio:0.0,        $
        cx:0.0,                  $ 
        cy:0.0,                  $
        peri:0.0,                $
        area:0.0,                $
        roi:OBJ_NEW('IDLanROI'), $
        bkgd:0.0,                $
        fwhm:0.0,                $
        hght:0.0,                $
        bboxx:fltarr(5),         $
        bboxy:fltarr(5),         $
        fwlin1:fltarr(2),         $
        fwlin2:fltarr(2),         $
        longx1:0.0,              $
        longy1:0.0,              $
        longx2:0.0,              $
        longy2:0.0,              $
        shrtx1a:0.0,             $
        shrty1a:0.0,             $
        shrtx2a:0.0,             $
        shrty2a:0.0,             $
        shrtx1b:0.0,             $
        shrty1b:0.0,             $
        shrtx2b:0.0,             $
        shrty2b:0.0}
;Create new sigmoids array
sigmoids=replicate(sigdat_mod,nfiles)
;Init assume start counting from 0
start = 0

;Check if save file already exists
f_chck = file_test(scratch_path+'sigmoid_sizedata.sav')

;restart from last if file already exists
if f_chck eq 1 then begin
    restore,scratch_path+'sigmoid_sizedata.sav'
    ;check to see filled sigmoid values
    comp = where(sigmoids.filename ne '',cnt_chck)
    start = cnt_chck
endif

print,'##########################'
start_fmt = '("Starting at index ",I6)'
print,string([start],format=start_fmt)
print,'##########################'
for xx=start,nfiles-1 do begin
  mreadfits,fits_files[xx],index1,data1,/verbose

  ;Pick default window Size
  xwdw_size=1024
  ywdw_size=1024

  ;Get axis size for the image
  img_xsize=index1[0].naxis1
  img_ysize=index1[0].naxis2

  ;create new window size based on the image size if not square
  if img_xsize ne img_ysize then begin
      x_rat = xwdw_size/img_xsize
      y_rat = ywdw_size/img_ysize

      ;if X image size is smaller make x window size smaller
      if img_xsize lt img_ysize then xwdw_size = y_rat*img_xsize $
      ;If Y image sie is smaller make y window size smaller
      else ywdw_size = x_rat*img_ysize 
  endif

  ;Make sure scaling is an integer number 
  ; 2018/03/01 J. Prchlik
  xwdw_size = round(xwdw_size/float(img_xsize))*img_xsize
  ywdw_size = round(ywdw_size/float(img_ysize))*img_ysize


  ;scale down images larger than 1024 to 1024
  bigx = img_xsize/1024.
  bigy = img_ysize/1024.

  ;Scale down tow prevent images bigger than 1024x1024
  if xwdw_size gt 1024 then xwdw_size = fix(temporary(img_xsize)/bigx)
  if ywdw_size gt 1024 then ywdw_size = fix(temporary(img_ysize)/bigy)

  ;normalize data by exposure time
  ;use actually normalization exp. time and convert from microseconds
  data1 = temporary(float(data1))/float(index1[0].e_etim)*1.e6

  filesplit=strsplit(fits_files[xx],'/',/extract)
  flnm=filesplit[n_elements(filesplit)-1]
  sigmoids[xx].filename=flnm
  sigmoids[xx].date=index1[0].DATE_OBS ;strmid(flnm,3,15)
  print,'Current file: '+sigmoids[xx].filename
  print,'Current date: '+sigmoids[xx].date
  loadct,3
  scmin=0.1
  useg = where(finite(alog10(data1)))
  scmin=cgPercentiles(alog10(data1[useg]),percentiles=.010)
  scmax=cgPercentiles(alog10(data1[useg]),percentiles=.999)
  imgok=0
  expire=0
  rscl=0
  rescale_image=''
  while not(imgok) do begin
     window,5,xs=xwdw_size,ys=ywdw_size



     ;Calculate edgedog right away and use to set the minimum for plotting and image background subtraction
     ;radius to scale to select the sigmoid
     ;rad_2 = abs(float(px3-px1))/abs(float(py3-py1))
     rad_1 = 3.6
     ;if rad_2 lt 1 then rad_2 = 1./rad_2
     ;rad_2 = rad_1*rad_2
     rad_2 = 15.

     ;Scale to find image edges 2018/02/22 J. Prchlik
     ;result = edge_dog(data1,radius1=rad_1,radius2=rad_2,threshold=1,zero_crossings=[0,255])                 
     ;Changed to 0 1 for mask 2018/02/23 J. Prchlik
     ;Multiplied exptime time back in because edge_dog is abs number dependent for some reason
     ; Removed exposure time because XRT L1 is already normalized
     result = edge_dog(data1,radius1=rad_1,radius2=rad_2,threshold=1,zero_crossings=[0,255])                 

     ;Create mask to remove ARs
     armask = abs(result-255)/255

     ;if AR mask is all 0 then set the entire mask array to 1
     if total(armask) eq 0 then armask = armask+1

     ;Get the median where ar does exists
     med_bkgd = median(data1[where(armask)])
     med_stdd = stddev(data1[where(armask)]-med_bkgd) 

     ;loop a couple times to find the best background
     looper = 1
     counter= 0

     while looper eq 1 do begin
         ;cut out points more than 3 sigma from the median and are in ARs
         good = where(armask and abs(data1-med_bkgd) lt 3.*med_stdd,cnt_good)

         ; Calculate the new backgroudn and standard deviation
         datat = data1[good]
         med_bkgd = median(datat)
         med_stdd = stddev(datat-med_bkgd) 

         ;if 95% percent of points are within 2 standard deviations compute the background level
         sigm = abs(datat-med_bkgd) lt 2.*med_stdd
         pert = float(total(sigm))/float(n_elements(datat)) 
         if pert gt 0.94 then looper = 0

         ;Quit after 6 tries
         if counter gt 5 then looper = 0 
         counter = counter+1
     endwhile

 
     ; guess scmin is the calculated background level
     if med_bkgd gt 0 then scmin = alog10(med_bkgd)

      

     image = bytscl(rebin(alog10(data1),xwdw_size,ywdw_size),min=scmin,max=scmax)
     tv,image
     ;show filter and time 
     estr = strcompress(index1[0].exptime)+'s'
     ;filter wheels
     fw1 = index1[0].EC_FW1_
     fw2 = index1[0].EC_FW2_
     plt_lab = fw1+'/'+fw2+'/'+estr
     xyouts,10,10,plt_lab,/device,charsize=2.0,charthick=2.0



     ;overlay image filter
     ;result = edge_dog(data1,radius1=6.0,radius2=20,threshold=15,zero_crossings=[0,255])
     ;radius help isolate the sigmoid
     ;Moved to after polygon 2018/02/28 J. Prchlik
     ;result = loop_edge_dog(data1,radius1=3.0,radius2=15.0,threshold=1,zero_crossings=[0,255])                 
     ;tv,bytscl(rebin(result,xwdw_size,ywdw_size))

     ;get the location of the maximum axis
     ;max_axis = brute_force_max_dis(ind_loc)

     ;oplot,[max_axis[1],max_axis[2]],[max_axis[3],max_axis[4]],color=200,thick=3
     ;min_axis = brute_force_min_dis(ind_loc,max_axis[1],max_axis[3],max_axis[2],max_axis[4])
     ;oplot,[min_axis[1],min_axis[2]],[min_axis[3],min_axis[4]],color=255,thick=3
     ;oplot,[min_axis[6],min_axis[7]],[min_axis[8],min_axis[9]],color=255,thick=3


     print,'Current image settings:'
     print,'Lower bound: ',scmin
     print,'Upper bound: ',scmax
     print,'Do you want to rescale the image? Enter 1 for yes, 0 for no.'
     read,rescale_image
     if not(checkstring(rescale_image)) then begin
       message,'WARNING: You must enter 0 or 1.',/cont
       expire=expire+1
     endif else begin
       rscl=fix(rescale_image)
       if ((rscl lt 0) OR (rscl gt 1)) then begin
         message,'WARNING: You must enter 0 or 1',/cont
         expire=expire+1
       endif else begin
         if (rscl) then begin
           oldmin=scmin
           oldmax=scmax
           input_quantity,scmin,error=err,qname='lower bound',qmin=-3.0
           if (err) then begin
             message,'Invalid inputs. The lower bound will not be changed.',/cont
             scmin=oldmin
           endif
           input_quantity,scmax,error=err,qname='upper bound',qmin=scmin
           if (err) then begin
             message,'Invalid inputs. The upper bound will not be changed.',/cont
             scmax=oldmax
           endif
           expire=0
         endif else begin
           expire=0
           imgok=1
         endelse
        endelse
      endelse
      if (expire gt 3) then begin
        message,'Too many consecutive invalid inputs.',/cont
        message,'The program will now quit automatically.'
      endif
  endwhile
;;  arcsec_per_pixel=0.6  ;; Conversion for AIA data
  arcsec_per_pixel=1.0286  ;; Conversion for XRT data
  arcsec_per_pixel_x=index1[0].cdelt1  ;; any data with reasonable headers
  arcsec_per_pixel_y=index1[0].cdelt2  ;; any data with reasonable headers
  ;Remove binscale because it is not need with proper cdelt keywords
  ;binscale=index1[0].chip_sum
  img_xsize=index1[0].naxis1
  img_ysize=index1[0].naxis2
  devicex_to_imgx=float(xwdw_size)/(img_xsize);*binscale)
  devicey_to_imgy=float(ywdw_size)/(img_ysize);*binscale)
  arcsec_per_devicex=arcsec_per_pixel_x/devicex_to_imgx
  arcsec_per_devicey=arcsec_per_pixel_y/devicey_to_imgy

  ;convert pixel to absolute x,y values
  xr0 = index1[0].crval1
  yr0 = index1[0].crval2
  xp0 = index1[0].crpix1
  yp0 = index1[0].crpix2
  ;offset if the lower left cornern so 
  xoff = xr0-arcsec_per_pixel_x*xp0
  yoff = yr0-arcsec_per_pixel_y*yp0

  skp=''
  skipthis=0
  print,''
  print,'Now that you have seen it...'
  print,'Do you want to use this file or skip it?'
  print,'Enter 1 to SKIP THIS FILE, or any other key to continue measuring the sigmoid.'
  read,skp
  if (checkstring(skp)) then begin
     skipthis=fix(skp)
     if (skipthis) then done=1 else done=0
  endif else done=0
  while not(done) do begin
     ;Changed to create polygon around Sigmoid
     ;Change back to manual can work with auto later
     ; 2018/02/23 J. Prchlik
     print,''
     print,"Click Trailing Sigmoid Long Axis (Point 1 of 2)"
     cursor,px1,py1,/down,/device
     xyouts,px1,py1,'L',/device,alignment=0.5
     print,"Click Leading Sigmoid Long Axis (Point 2 of 2)"
     cursor,px2,py2,/down,/device
     xyouts,px2,py2,'L',/device,alignment=0.5
     print,"Click Sigmoid Southern Trailing Short Axis (Point 1 of 2)"
     cursor,px3,py3,/down,/device
     xyouts,px3,py3,'X',/device,alignment=0.5
     print,"Click Sigmoid Northern Trailing Short Axis (Point 2 of 2)"
     cursor,px4,py4,/down,/device
     xyouts,px4,py4,'X',/device,alignment=0.5
     print,"Click Sigmoid Southern Leading Short Axis (Point 1 of 2)"
     cursor,px5,py5,/down,/device
     xyouts,px5,py5,'X',/device,alignment=0.5
     print,"Click Sigmoid Northern Leading Short Axis (Point 1 of 2)"
     cursor,px6,py6,/down,/device
     xyouts,px6,py6,'X',/device,alignment=0.5

     ;Changed to polygon around sigmoid
     print,"Click Lower Left Bounding Box (Point 1 of 4)"
     cursor,bx1,by1,/down,/device
     xyouts,bx1,by1,'B',/device,alignment=0.5
     print,"Click Upper Left Bounding Box (Point 2 of 4)"
     cursor,bx2,by2,/down,/device
     xyouts,bx2,by2,'B',/device,alignment=0.5
     print,"Click Upper Right Bounding Box (Point 3 of 4)"
     cursor,bx3,by3,/down,/device
     xyouts,bx3,by3,'B',/device,alignment=0.5
     print,"Click Lower Right Bounding Box (Point 4 of 4)"
     cursor,bx4,by4,/down,/device
     xyouts,bx4,by4,'B',/device,alignment=0.5

     ;Create sigmoid region where to calculate the FWHM
     print,''
     print,"Click Bottom FWHM Point (Point 1 of 2)"
     cursor,sx1,sy1,/down,/device
     xyouts,sx1,sy1,'F',/device,alignment=0.5
     print,"Click Top FWHM Point (Point 2 of 2)"
     cursor,sx2,sy2,/down,/device
     xyouts,sx2,sy2,'F',/device,alignment=0.5

     ;Create line tracing the sigmoid
     xvals = [bx1,bx2,bx3,bx4,bx1];,bx5,bx6,bx7] Changed to bolygon J. brchlik 2018/02/22
     yvals = [by1,by2,by3,by4,by1];,by5,by6,by7] Changed to bolygon J. brchlik 2018/02/22
     ;Changed to polygon 2018/02/22 J. Prchlik
     ;samp = 100
     ;xgrid = findgen(samp)*(max(xvals)-min(xvals))/float(samp)+min(xvals)

     ;result = roberts(data1)                 
     ;result = sobel(data1)                 
     ;result = prewitt(data1)                 
     ;result = shift_diff(data1)                 
     ;result = edge_dog(data1)                 
     ;result = laplacian(data1)                 
     ;result = emboss(data1)

     ;adjust the input by the markup scale
     scale_x = float(xwdw_size)/float(img_xsize)
     scale_y = float(ywdw_size)/float(img_ysize)
     adjxv = xvals/scale_x
     adjyv = yvals/scale_y

     ;store the shape of the array so you can reform the indices after creating mask
     ar_sp = size(result)

     ;index location of the bright regions
     ; Get indices for full array 2018/02/23 J. Prchlik
     sig = where(result gt -.50)
     ;create index array for the entire images for the sigmoid
     ind_loc = array_indices(result,sig)

     ;break into x and y indices
     ind_x = ind_loc[0,*]
     ind_y = ind_loc[1,*]
     ;Compute x and y limit functions
     lims = comp_limits(adjxv,adjyv)
     ;Create array of 1 and 0 for box
     xmin = ind_x ge lims[1]*ind_y+lims[0]
     xmax = ind_x le lims[3]*ind_y+lims[2]
     ymin = ind_y le lims[5]*ind_x+lims[4]
     ymax = ind_y ge lims[7]*ind_x+lims[6]
     vbox = xmin*xmax*ymin*ymax
     ibox = reform(vbox,ar_sp[1],ar_sp[2])

     ;Idea stems from /Applications/exelis/idl81/help/pdf/image.pdf Chapter 6 working with ROIs
     ;Programmatically Defining ROIs 2018/02/23 J. Prchlik
     ;Create a mask which highlights interior sigmoid region and is inside bounding box
     mask = ibox*result


     ;if mask is empty do not try to autocompute AR properties
     ; Added 2018/03/05 J. Prchlik
     if total(mask) gt 0 then begin

         ;Get boundary of created countour
         CONTOUR, mask, LEVEL = 1,  $
                XMARGIN = [0, 0], YMARGIN = [0, 0], $
                /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, $
                XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS



         ;Create new ROI obejct using contour
          line = [LINDGEN(PathInfo(0).N), 0] & $
          roi_obj = OBJ_NEW('IDLanROI', $
             (scale_x*pathXY(*, pathInfo(0).OFFSET + line))[0, *], $
             (scale_y*pathXY(*, pathInfo(0).OFFSET + line))[1, *]) & $
          ;Draw ROI on plot
          DRAW_ROI, roi_obj, COLOR = 0,/device,/LINE_FILL


         ;Create image object
         ;img_obj = OBJ_NEW( 'IDLanROI',ind_x*vbox,ind_y*vbox)  

         ;Compute image stats from region
         ;area = pixel area
         ;cent = centriod
         ;peri = perimeter distance in pixels
         ;SPATIAL OFFSET is the zero point scale is the input pixel to arcsecond conversion
         ;SPATIAL scale is the input pixel to arcsecond conversion
         geo_comp = roi_obj.ComputeGeometry(AREA = area,$
                          CENTROID=cent, PERIMETER=peri,$ 
                          SPATIAL_SCALE=[arcsec_per_pixel_x,arcsec_per_pixel_y])


         ;calibrated centriod
         ;Use FWHM for centriod for now
         ;cal_cent = (cent/[scale_x,scale_y,1]-[xp0,yp0,0])*[arcsec_per_pixel_x,arcsec_per_pixel_y,1]+[xr0,yr0,0]
         cal_cent = (([sx1+sx2,sy1+sy2]/2./[scale_x,scale_y])-[xp0,yp0])*[arcsec_per_pixel_x,arcsec_per_pixel_y]+[xr0,yr0]
     endif else begin

         area = -9999.9
         peri = -9999.9
         ;Use FWHM points to compute center points
         cal_cent = (([sx1+sx2,sy1+sy2]/2./[scale_x,scale_y])-[xp0,yp0])*[arcsec_per_pixel_x,arcsec_per_pixel_y]+[xr0,yr0]

     endelse


     ;Compute the background into per arcsec^2
     med_bkgd = med_bkgd/(arcsec_per_pixel_x*arcsec_per_pixel_y)
  
     ;Plot the edge image
     ;Put in and remove for testing purposes 2018/02/28 J. Prchlik
     ;test_plot = [[[image]],[[image]],[[bytscl(rebin(result,xwdw_size,ywdw_size))]]]
     ;Not using Test plot for now
     ;tv,test_plot,true=3

     ;get the location of the maximum axis
     ;Do after defining polygon 2018/02/22 J. Prchlik 
     ;Remove auto max length 2018/02/23 J. Prchlik
     ;max_axis = brute_force_max_dis(ind_loc,data1,adjxv,adjyv)

     ;plots,[max_axis[1],max_axis[2]]*scale_x,[max_axis[3],max_axis[4]]*scale_y,color=200,thick=3,/device

     ;Get values for the minimum axis
     ;Do after defining polygon 2018/02/22 J. Prchlik 
     ;min_axis = brute_force_min_dis(ind_loc,data1,adjxv,adjyv)
     ;Remove automatic min values 2018/02/23 J. Prchlik
     ;plots,[min_axis[1],min_axis[2]]*scale_x,[min_axis[3],min_axis[4]]*scale_y,color=240,thick=3,/device
     ;plots,[min_axis[6],min_axis[7]]*scale_x,[min_axis[8],min_axis[9]]*scale_y,color=240,thick=3,/device

     ;interpolate quadratic
     ;Changed to polygon 2018/02/22 J. Prchlik
     ;py_inp = interpol(yvals,xvals,xgrid,/Quadratic)

     ;cgplot,xvals,yvals,color=255,thick=3,linestyle=1,/device,/noerase,xrange=[0,xwdw_size],yrange=[0,ywdw_size],xstyle=1,ystyle=1,/overplot
     ;plot quadratic interpolation
     ;Changed to polygon 2018/02/22 J. Prchlik
     ;plots,xgrid,py_inp,color=255,thick=5,linestyle=0,/device
     ;Keep Bounding box plot and add axis plots
     plots,xvals,yvals,color=255,thick=5,linestyle=0,/device
     plots,[px1,px2],[py1,py2],color=255,thick=2,linestyle=2,/device
     plots,[px3,px4],[py3,py4],color=255,thick=2,linestyle=2,/device
     plots,[px5,px6],[py5,py6],color=255,thick=2,linestyle=2,/device
 
 
     ;Get the FWHM Value J. Prchlik 2018/02/23
     outp = real_fwhm(data1,adjxv,adjyv,sx1,sy1,sx2,sy2,scale_x,scale_y,arcsec_per_pixel_x,arcsec_per_pixel_y,med_bkgd)
     fwhm = outp[0]
     hght = outp[1]

     ;temp quick answers
     ; Updated back with manual parameters 2018/02/23 J. Prchlik
     lx1 = px1
     lx2 = px2
     ly1 = py1
     ly2 = py2
     ;Trailing is a leading is b
     sx1a = px3
     sx2a = px4
     sy1a = py3
     sy2a = py4
     sx1b = px5
     sx2b = px6
     sy1b = py5
     sy2b = py6

     ;Get long Axis information
     lx_dev=abs(lx1-lx2)
     ly_dev=abs(ly1-ly2)
     lx_arc=arcsec_per_devicex*lx_dev
     ly_arc=arcsec_per_devicey*ly_dev
     long_axis_xy=sqrt((lx_dev^2)+(ly_dev^2))
     long_axis_arc=sqrt((lx_arc^2)+(ly_arc^2))

     ;Get short Axis information Trailing
     sxa_dev=abs(sx1a-sx2a)
     sya_dev=abs(sy1a-sy2a)
     sxa_arc=arcsec_per_devicex*sxa_dev
     sya_arc=arcsec_per_devicey*sya_dev
     short_axisa_xy =sqrt((sxa_dev^2)+(sya_dev^2))
     short_axisa_arc=sqrt((sxa_arc^2)+(sya_arc^2))

     ;Get short Axis information Trailing
     sxb_dev=abs(sx1b-sx2b)
     syb_dev=abs(sy1b-sy2b)
     sxb_arc=arcsec_per_devicex*sxb_dev
     syb_arc=arcsec_per_devicey*syb_dev
     short_axisb_xy =sqrt((sxb_dev^2)+(syb_dev^2))
     short_axisb_arc=sqrt((sxb_arc^2)+(syb_arc^2))

     ;Get center pixel informaiton
     ;Get center pixel informaiton
     print,''
     print,'The long axis size in DEVICE-units is:'
     print,strcompress(string(long_axis_xy),/remove_all)
     print,'The short axis size in DEVICE-units is:'
     print,strcompress(string(short_axisa_xy),/remove_all),' ',strcompress(string(short_axisb_xy),/remove_all)
     print,'The aspect ratio is:'
     print,strcompress(string(long_axis_xy/(short_axisa_xy+short_axisb_xy)*2.),/remove_all)

     print,''
     print,'The long axis size in ARCSEC is:'
     print,strcompress(string(long_axis_arc),/remove_all)
     print,'The short axis size in ARCSEC is:'
     print,strcompress(string(short_axisa_arc),/remove_all),' ',strcompress(string(short_axisb_arc),/remove_all)
     print,'The aspect ratio should be the same as above:'
     print,strcompress(string(long_axis_arc/(short_axisa_arc+short_axisb_arc)*2.),/remove_all)
     print,string([fwhm,hght,area],format='("The FWHM =",F6.2,"arcsec, Height = ",F11.1,"#/s/arcsec^2, Area = ",F15.1,"arcsec^2")')
     print,strcompress(string(long_axis_arc/(short_axisa_arc+short_axisb_arc)*2.),/remove_all)

     print,string(cal_cent[0:1],format='("x = ",F6.1," y = ",F6.1)')

     continue=''
     print,''
     print,'Shall we try this same file again?'
     print,'Enter 1 for yes, any other key for no.'
     read,continue
     if (continue ne '1') then done=1 else begin
       window,5,xs=xwdw_size,ys=ywdw_size
       tv,image
       xyouts,10,10,plt_lab,/device,charsize=2.0,charthick=2.0
     endelse
  endwhile

  print,''
  print,'You have just completed this file: '+sigmoids[xx].filename
  print,'for this date and time: '+sigmoids[xx].date

  if not(skipthis) then begin
     sigmoids[xx].longx1=lx1/scale_x  ;Trailing Long Axis X coordinate in Pixels
     sigmoids[xx].longx2=lx2/scale_x  ;Leading  Long Axis X coordinate in Pixels
     sigmoids[xx].longy1=ly1/scale_y  ;Trailing Long Axis Y coordinate in Pixels
     sigmoids[xx].longy2=ly2/scale_y  ;Leading  Long Axis Y coordinate in Pixels
     sigmoids[xx].size=long_axis_arc  ;Length of long axis in Arcsec
     sigmoids[xx].sizea=short_axisa_arc ;Length of trailing short axis in Arcsec
     sigmoids[xx].sizeb=short_axisb_arc ;Length of leading short axis in Arcsec
     sigmoids[xx].aspect_ratio=long_axis_arc/(short_axisa_arc+short_axisb_arc)*2. ;Aspect ratio
     sigmoids[xx].bkgd=med_bkgd ;Median background in #/s/arcsec^2        
     sigmoids[xx].fwhm=fwhm     ; Sigmoid FWHM in Arcsec    
     sigmoids[xx].hght=hght     ; Sigmoid Height in #/s/arcsec^2    
     sigmoids[xx].area=area     ; Sigmoid Area in arcsec^   
     sigmoids[xx].peri=peri     ; Sigmoid Perimeter in arcsec 
     sigmoids[xx].roi=roi_obj   ; ROI object of sigmoid 
     sigmoids[xx].cx=float(cal_cent[0]) ; Centriod point of sigmoid in X in arcsec
     sigmoids[xx].cy=float(cal_cent[1]) ; Centriod point of sigmoid in Y in arcsec
     sigmoids[xx].fwlin1=[sx1,sy1]/[scale_x,scale_y] ;FWHM line trailing points in pixels        
     sigmoids[xx].fwlin2=[sx2,sy2]/[scale_x,scale_y] ;FWHM line leading  points in pixels        
     sigmoids[xx].bboxx=xvals/scale_x  ;bounding box x-values of sigmoid in pixels     
     sigmoids[xx].bboxy=yvals/scale_y  ;bounding box y-values of sigmoid in pixels     
     sigmoids[xx].shrtx1a=sx1a/scale_x  ;X coordinate of Lower Trailing Sigmoid short axis in pixels    
     sigmoids[xx].shrtx2a=sx2a/scale_x  ;X coordinate of Upper Trailing Sigmoid short axis in pixels    
     sigmoids[xx].shrtx1b=sx1b/scale_x  ;X coordinate of Lower Leading  Sigmoid short axis in pixels    
     sigmoids[xx].shrtx2b=sx2b/scale_x  ;X coordinate of Upper Leading  Sigmoid short axis in pixels
     sigmoids[xx].shrty1a=sy1a/scale_y  ;Y coordinate of Lower Trailing Sigmoid short axis in pixels    
     sigmoids[xx].shrty2a=sy2a/scale_y  ;Y coordinate of Upper Trailing Sigmoid short axis in pixels    
     sigmoids[xx].shrty1b=sy1b/scale_y  ;Y coordinate of Lower Leading  Sigmoid short axis in pixels    
     sigmoids[xx].shrty2b=sy2b/scale_y  ;Y coordinate of Upper Leading  Sigmoid short axis in pixels
     sig_id=''
     noaa_id=''


      
      ;Guess the sigmoid ID based on rotation
      ;2018/03/05 J. Prchlik
      ;Get previous sigmoids
      comp = where(sigmoids.sig_id ne '',cnt_chck)

     print,'What is the ID of the sigmoid you just identified?'
     if cnt_chck gt 0 then begin 

         ;get important values from sigmoid structure
         pos_x = sigmoids[comp].cx
         pos_y = sigmoids[comp].cy
         pos_t = sigmoids[comp].date
         pos_id= sigmoids[comp].sig_id
         sig_guess = -100
         
         ;Create large array for storing rotation
         rot_s = fltarr([2,cnt_chck])
         ;Rotate position to obs time
         for j=0, cnt_chck-1 do $
             rot_s[*,j] = rot_xy(pos_x[j], pos_y[j], tstart=pos_t[j], tend=index1[0].DATE_OBS)

         ;Store x,y in separate array
         rot_x = rot_s[0,*]
         rot_y = rot_s[1,*]
  
         ;get the closest ar after rotation
         dis_m = sqrt((rot_x-cal_cent[0])^2+(rot_y-cal_cent[1])^2)
         ;give 50 arcsec buffer and needs to be min
         min_i = where((dis_m lt 100) and (dis_m eq min(dis_m,/nan)),min_cnt)
  
         ;Print rotated answer if found
         if min_cnt gt 0 then begin 
             sig_guess = fix(pos_id[min_i])
             print,string(sig_guess,format='("Best Guess = ",I5)')
         endif
     endif

     read,sig_id
     sigmoids[xx].sig_id=sig_id

     ;Turn header obs time into anytime seconds
     obs_tim = anytim(index1[0].DATE_OBS)
     ;get +/- 7 Days
     offset = 7.*24.*3600. ; days to seconds
     obs_tim_s = anytim(obs_tim-offset,/ex)
     obs_tim_e = anytim(obs_tim+offset,/ex)

     ;covert to strings to pass to query
     dat_fmt = '(I04,"-",I02,"-",I02)'
     obs_str_s = string([obs_tim_s[6],obs_tim_s[5],obs_tim_s[4]],format=dat_fmt)
     obs_str_e = string([obs_tim_e[6],obs_tim_e[5],obs_tim_e[4]],format=dat_fmt)

     ;Added best guess of NOAA number
     query=ssw_her_make_query(obs_str_s,obs_str_e,/ar,x1=cal_cent[0],x2=cal_cent[1])
     her=ssw_her_query(query,/str) 
     if n_elements(size(her)) gt 3 then begin 
         ; Get time, postion and name
         pos_x = her.ar.required.event_coord1
         pos_y = her.ar.required.event_coord2
         pos_t = her.ar.required.EVENT_STARTTIME
         ar_guess = her.ar.optional.AR_NOAANUM


         rot_p = fltarr([2,n_elements(ar_guess)])

         ;Rotate position to obs time
         for j=0, n_elements(ar_guess)-1 do $
             rot_p[*,j] = rot_xy(pos_x[j], pos_y[j], tstart=pos_t[j], tend=index1[0].DATE_OBS)

        
         ;Store x,y in separate array
         rot_x = rot_p[0,*]
         rot_y = rot_p[1,*]
  
         ;get the closest ar after rotation
         dis_m = sqrt((rot_x-cal_cent[0])^2+(rot_y-cal_cent[1])^2)
         min_i = where(dis_m eq min(dis_m[where(ar_guess gt 0)],/nan))
         ar_guess = fix(ar_guess[min_i])
     endif

     print,'What is the NOAA active region number of the sigmoid?'
     print,'(If there is no NOAA #, just hit enter.)'
     if n_elements(size(her)) gt 3 then print,string(ar_guess,format='("Best Guess = ",I5)')
     read,noaa_id
     sigmoids[xx].NOAA_id=noaa_id
     save,filename=scratch_path+'sigmoid_sizedata.sav',sigmoids
  endif else begin
     sigmoids[xx].sig_id='file skipped'
  endelse
endfor

;stop
;
;return
end
