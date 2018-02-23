
;; Use sswidl
;; This version only looks for data in the Scratch disk

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
;Function to calculate the maximum distance and index locations of the maximum difference
;
;ind is a 2D array of indices containing the location of the sigmoid
;
;Usage
;outp = brute_force_max_dis(inds)
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
;Function to calculate the maximum distance and index locations of the minimum axis difference
;
;ind is a 2D array of indices containing the location of the sigmoid
;nx1,ny1 are the left x,y coordinates of the maximum distance in the sigmoid
;nx1,ny1 are the right x,y coordinates of the maximum distance in the sigmoid
;
;Usage
;outp = brute_force_min_dis(inds,idat,xbox,ybox)
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

;Loop to find the best Edge values
; Don't use 2018/02/22 J. Prchlik
;
;
;;
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

pro sigmoidsize_adv
;scratch_path='/Volumes/Scratch/Users/ehanson/XRT_fits/'
scratch_path='examples/'
xwdw_size=1024
ywdw_size=1024
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

sigdat_mod={sig_id:'',        $
        NOAA_id:0,        $
        filename:'',      $
        date:'',          $
        size:0.0,         $
        sizea:0.0,        $
        sizeb:0.0,        $
        aspect_ratio:0.0, $
        fwhm:0.0,         $
        bboxx:fltarr(5),  $
        bboxy:fltarr(5),  $
        longx1:0.0,       $
        longy1:0.0,       $
        longx2:0.0,       $
        longy2:0.0,       $
        shrtx1a:0.0,      $
        shrty1a:0.0,      $
        shrtx2a:0.0,      $
        shrty2a:0.0,      $
        shrtx1b:0.0,      $
        shrty1b:0.0,      $
        shrtx2b:0.0,      $
        shrty2b:0.0}
sigmoids=replicate(sigdat_mod,nfiles)


for xx=0,nfiles-1 do begin
  mreadfits,fits_files[xx],index1,data1,/verbose
  filesplit=strsplit(fits_files[xx],'/',/extract)
  flnm=filesplit[n_elements(filesplit)-1]
  sigmoids[xx].filename=flnm
  sigmoids[xx].date=strmid(flnm,3,15)
  print,'Current file: '+sigmoids[xx].filename
  print,'Current date: '+sigmoids[xx].date
  loadct,3
  scmin=0.1
  scmin=cgPercentiles(data1,percentiles=.005)
  scmax=cgPercentiles(data1,percentiles=.999)
  imgok=0
  expire=0
  rscl=0
  rescale_image=''
  while not(imgok) do begin
     window,5,xs=xwdw_size,ys=ywdw_size

     image = bytscl(rebin(data1,xwdw_size,ywdw_size),min=scmin,max=scmax)
     tv,image


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
           input_quantity,scmin,error=err,qname='lower bound',qmin=0
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
  binscale=index1[0].chip_sum
  img_xsize=index1[0].naxis1
  img_ysize=index1[0].naxis2
  devicex_to_imgx=float(xwdw_size)/(img_xsize*binscale)
  devicey_to_imgy=float(ywdw_size)/(img_ysize*binscale)
  arcsec_per_devicex=arcsec_per_pixel/devicex_to_imgx
  arcsec_per_devicey=arcsec_per_pixel/devicey_to_imgy

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
     xyouts,bx1,by1,'X',/device,alignment=0.5
     print,"Click Upper Left Bounding Box (Point 2 of 4)"
     cursor,bx2,by2,/down,/device
     xyouts,bx2,by2,'X',/device,alignment=0.5
     print,"Click Upper Right Bounding Box (Point 3 of 4)"
     cursor,bx3,by3,/down,/device
     xyouts,bx3,by3,'X',/device,alignment=0.5
     print,"Click Lower Right Bounding Box (Point 4 of 4)"
     cursor,bx4,by4,/down,/device
     xyouts,bx4,by4,'X',/device,alignment=0.5


     ;Create line tracing the sigmoid
     xvals = [bx1,bx2,bx3,bx4,bx1];,bx5,bx6,bx7] Changed to bolygon J. brchlik 2018/02/22
     yvals = [by1,by2,by3,by4,by1];,by5,by6,by7] Changed to bolygon J. brchlik 2018/02/22
     ;Changed to polygon 2018/02/22 J. Prchlik
     ;samp = 100
     ;xgrid = findgen(samp)*(max(xvals)-min(xvals))/float(samp)+min(xvals)

     ;radius to scale to select the sigmoid
     ;rad_2 = abs(float(px3-px1))/abs(float(py3-py1))
     rad_1 = 3.6
     ;if rad_2 lt 1 then rad_2 = 1./rad_2
     ;rad_2 = rad_1*rad_2
     rad_2 = 15.

     ;Scale to find image edges 2018/02/22 J. Prchlik
     ;result = edge_dog(data1,radius1=rad_1,radius2=rad_2,threshold=1,zero_crossings=[0,255])                 
     ;Changed to 0 1 for mask 2018/02/23 J. Prchlik
     result = edge_dog(data1,radius1=rad_1,radius2=rad_2,threshold=1,zero_crossings=[0,255])                 

      
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
      DRAW_ROI, roi_obj, COLOR = 80,/device


     ;Create image object
     ;img_obj = OBJ_NEW( 'IDLanROI',ind_x*vbox,ind_y*vbox)  

     ;Compute image stats from region
     ;area = pixel area
     ;cent = centriod
     ;peri = perimeter distance in pixels
     geo_comp = roi_obj.ComputeGeometry(AREA = area,$
                      CENTROID=cent, PERIMETER=peri,$ 
                      SPATIAL_OFFSET=soff, SPATIAL_SCALE=sscl)




  
     ;Plot the edge image
     ;Put in and remove for testing purposes 2018/02/28 J. Prchlik
     test_plot = [[[image]],[[image]],[[bytscl(rebin(result,xwdw_size,ywdw_size))]]]
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
     sxb_dev=abs(sx1a-sx2a)
     syb_dev=abs(sy1a-sy2a)
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

     continue=''
     print,''
     print,'Shall we try this same file again?'
     print,'Enter 1 for yes, any other key for no.'
     read,continue
     if (continue ne '1') then done=1 else begin
       tv,bytscl(rebin(data1,xwdw_size,ywdw_size),min=scmin,max=scmax)
     endelse
  endwhile

  print,''
  print,'You have just completed this file: '+sigmoids[xx].filename
  print,'for this date and time: '+sigmoids[xx].date

  if not(skipthis) then begin
     sigmoids[xx].longx1=lx1
     sigmoids[xx].longy1=ly1
     sigmoids[xx].longx2=lx2
     sigmoids[xx].longy2=ly2
     sigmoids[xx].size=long_axis_arc
     sigmoids[xx].sizea=short_axisa_arc
     sigmoids[xx].sizeb=short_axisb_arc
     sigmoids[xx].aspect_ratio=long_axis_arc/(short_axisa_arc+short_axisb_arc)*2.
     sigmoids[xx].fwhm=0.0         
     sigmoids[xx].bboxx=xvals       
     sigmoids[xx].bboxy=yvals       
     sigmoids[xx].shrtx1a=sx1a      
     sigmoids[xx].shrty1a=sy1a      
     sigmoids[xx].shrtx2a=sx2a      
     sigmoids[xx].shrty2a=sy2a      
     sigmoids[xx].shrtx1b=sx1b      
     sigmoids[xx].shrty1b=sy1b      
     sigmoids[xx].shrtx2b=sx2b
     sigmoids[xx].shrty2b=sy2b
     sig_id=''
     noaa_id=''
     print,'What is the ID of the sigmoid you just identified?'
     read,sig_id
     sigmoids[xx].sig_id=sig_id

     ;Added best guess of NOAA number

     print,'What is the NOAA active region number of the sigmoid?'
     print,'(If there is no NOAA #, just hit enter.)'
     read,noaa_id
     sigmoids[xx].NOAA_id=noaa_id
     save,filename=scratch_path+'sigmoid_sizedata.sav',sigmoids
  endif else begin
     sigmoids[xx].sig_id='file skipped'
  endelse
endfor

stop

return
end
