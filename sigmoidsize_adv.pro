
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

;###############################################
;
;Function to calculate the maximum distance and index locations of the maximum difference
;
;ind is a 2D array of indices containing the location of the sigmoid
;
;Usage
;outp = brute_force_max_dis(inds)
;###############################################
function brute_force_max_dis,inds

;break the indices into x and y values
ind_x = inds[0,*]
ind_y = inds[1,*]


;set up variables for maximum distance
max_d  = 0.
max_x1 = 0
max_x2 = 0
max_y1 = 0
max_y2 = 0


;Loop over all x-values to find the maximum distance for s
for i=0,n_elements(ind_x)-1 do begin
   ;temporary diffrance array
   t_diff = sqrt((ind_x[i]-ind_x)^2+(ind_y[i]-ind_y)^2)

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
;outp = brute_force_min_dis(inds,nx1,ny1,nx2,ny2)
;###############################################
function brute_force_min_dis,inds,nx1,ny1,nx2,ny2

;get unique x values
ind_x = inds[0,*]
unq_x = sort(ind_x)
unq_x = ind_x[unq_x]
unq_x = unq_x[uniq(unq_x)]

;break the indices into x and y values and correct for maximum distance slope
ind_y = inds[1,*]
m1 = (ny2-ny1)/(nx2-nx1)
b1 = ny1-nx1*m1

;rotate cooridinates to be in x direction with the x center at 0
ind_cx = ind_x-(nx2+nx1)/2.
ind_cy = ind_x*m1+b1

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


;Loop over all x-values to find the maximum distance for each small axis
for i=0,n_elements(unq_x)-1 do begin

   
   ;get yvalues were x values equal to some value
   x_use = where(ind_x eq unq_x[i])

   ;get transformed xvalue
   x_trn = unq_x[i]-(nx2+nx1)/2.

   ;interested yvalues=
   ind_iy = ind_cy[x_use]

   ;temporary diffrance array
   t_diff = abs(max(ind_iy)-min(ind_iy))

   ;maximum difference for this point and the index location of maximum
   m_diff = max(t_diff)
   i_diff = where(t_diff eq m_diff)

   ;if maximum difference is greater than the current maximum set new points
   case 1 of

       (m_diff gt max_da) and (x_trn lt 0):  begin
           ;update maximum differnce
           max_da = m_diff
           ;update initial axis points
           max_x1a = unq_x[i]
           max_y1a = min(ind_iy)
           ;update second axis points
           max_x2a = unq_x[i]
           max_y2a = max(ind_iy)
           end
       (m_diff gt max_db) and (x_trn gt 0):  begin
           ;update maximum differnce
           max_db = m_diff
           ;update initial axis points
           max_x1b = unq_x[i]
           max_y1b = min(ind_iy)
           ;update second axis points
           max_x2b = unq_x[i]
           max_y2b = max(ind_iy)
           end

        ;do whatever to just continue
        else: v=1
   endcase

endfor


return,[max_da,max_x1a,max_x2a,max_y1a,max_y2a,max_db,max_x1b,max_x2b,max_y1b,max_y2b]
end

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
   print,counter,radius1,mean(ts_diff(ds,1))
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

sigdat={sig_id:'',        $
        NOAA_id:0,        $
        filename:'',      $
        date:'',          $
        size:0.0,         $
        aspect_ratio:0.0, $
        longx1:0.0,       $
        longy1:0.0,       $
        longx2:0.0,       $
        longy2:0.0,       $
        shrtx1:0.0,       $
        shrty1:0.0,       $
        shrtx2:0.0,       $
        shrty2:0.0}
sigmoids=replicate(sigdat,nfiles)


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
     tv,bytscl(rebin(data1,xwdw_size,ywdw_size),min=scmin,max=scmax)


     ;overlay image filter
     ;result = edge_dog(data1,radius1=6.0,radius2=20,threshold=15,zero_crossings=[0,255])
     ;radius help isolate the sigmoid
     result = loop_edge_dog(data1,radius1=3.0,radius2=15.0,threshold=1,zero_crossings=[0,255])                 
     ;tv,bytscl(rebin(result,xwdw_size,ywdw_size))

     ;index location of the sigmoid
     sig = where(result gt 254.5)
     ;create index array for the entire images for the sigmoid
     ind_loc = array_indices(result,sig)

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
     print,''
     print,"Click Along Sigmoid Axis (Point 1 of 7)"
     cursor,px1,py1,/down,/device
     xyouts,px1,py1,'X',/device,alignment=0.5
     print,"Click Along Sigmoid Axis (Point 2 of 7)"
     cursor,px2,py2,/down,/device
     xyouts,px2,py2,'X',/device,alignment=0.5
     print,"Click Along Sigmoid Axis (Point 3 of 7)"
     cursor,px3,py3,/down,/device
     xyouts,px3,py3,'X',/device,alignment=0.5
     print,"Click Along Sigmoid Axis (Point 4 of 7)"
     cursor,px4,py4,/down,/device
     xyouts,px4,py4,'X',/device,alignment=0.5
     print,"Click Along Sigmoid Axis (Point 5 of 7)"
     cursor,px5,py5,/down,/device
     xyouts,px5,py5,'X',/device,alignment=0.5
     print,"Click Along Sigmoid Axis (Point 6 of 7)"
     cursor,px6,py6,/down,/device
     xyouts,px6,py6,'X',/device,alignment=0.5
     print,"Click Along Sigmoid Axis (Point 7 of 7)"
     cursor,px7,py7,/down,/device
     xyouts,px7,py7,'X',/device,alignment=0.5


     ;Create line tracing the sigmoid
     xvals = [px1,px2,px3,px4,px5,px6,px7]
     yvals = [py1,py2,py3,py4,py5,py6,py7]
     samp = 100
     xgrid = findgen(samp)*(max(xvals)-min(xvals))/float(samp)+min(xvals)
  
     ;interpolate quadratic
     py_inp = interpol(yvals,xvals,xgrid,/Quadratic)

     print,px1,px6,min(xgrid),max(xgrid)
     print,py1,py6,min(py_inp),max(py_inp)

     ;cgplot,xvals,yvals,color=255,thick=3,linestyle=1,/device,/noerase,xrange=[0,xwdw_size],yrange=[0,ywdw_size],xstyle=1,ystyle=1,/overplot
     ;plot quadratic interpolation
     plots,xgrid,py_inp,color=255,thick=5,linestyle=0,/device
     
     ;temp quick answers
     lx1 = px1
     lx2 = px6
     ly1 = py1
     ly2 = py6
     sx1 = px1
     sx2 = px2
     sy1 = py1
     sy2 = py2

     ;Get long Axis information
     lx_dev=abs(lx1-lx2)
     ly_dev=abs(ly1-ly2)
     lx_arc=arcsec_per_devicex*lx_dev
     ly_arc=arcsec_per_devicey*ly_dev
     long_axis_xy=sqrt((lx_dev^2)+(ly_dev^2))
     long_axis_arc=sqrt((lx_arc^2)+(ly_arc^2))

     ;Get short Axis information
     sx_dev=abs(sx1-sx2)
     sy_dev=abs(sy1-sy2)
     sx_arc=arcsec_per_devicex*sx_dev
     sy_arc=arcsec_per_devicey*sy_dev
     short_axis_xy=sqrt((sx_dev^2)+(sy_dev^2))
     short_axis_arc=sqrt((sx_arc^2)+(sy_arc^2))

     ;Get center pixel informaiton

     print,''
     print,'The long axis size in DEVICE-units is:'
     print,strcompress(string(long_axis_xy),/remove_all)
     print,'The short axis size in DEVICE-units is:'
     print,strcompress(string(short_axis_xy),/remove_all)
     print,'The aspect ratio is:'
     print,strcompress(string(long_axis_xy/short_axis_xy),/remove_all)

     print,''
     print,'The long axis size in ARCSEC is:'
     print,strcompress(string(long_axis_arc),/remove_all)
     print,'The short axis size in ARCSEC is:'
     print,strcompress(string(short_axis_arc),/remove_all)
     print,'The aspect ratio should be the same as above:'
     print,strcompress(string(long_axis_arc/short_axis_arc),/remove_all)

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
     sigmoids[xx].shrtx1=sx1
     sigmoids[xx].shrty1=sy1
     sigmoids[xx].shrtx2=sx2
     sigmoids[xx].shrty2=sy2
     sigmoids[xx].size=long_axis_arc
     sigmoids[xx].aspect_ratio=long_axis_arc/short_axis_arc
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
