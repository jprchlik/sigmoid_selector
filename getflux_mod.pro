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


pro get_hmi_file,date


end


pro getflux_mod,sigloc,fname=fname,odir=odir,line_file=line_file

;; Use sswidl
;; The line_file keyword should contain the filename of the IDL .sav file
;;   where you have stored the times at which vertical lines should be
;;   drawn on the plot. This file should be stored in the same directory
;;   where you intend to save the plot.


;; Directory containing HMI images
dir='/Volumes/Scratch/Users/ehanson/hmi_files/'
files=find_file(dir+'*.fits')
nfiles=n_elements(files)

filedates=lonarr(nfiles)
fdate_str=strarr(nfiles)
filetimes=lonarr(nfiles)
ftime_str=strarr(nfiles)
nchardir=strlen(dir)

for xx=0,nfiles-1 do begin
  fdate_str[xx]=strmid(files[xx],nchardir+3,8)
  filedates[xx]=long(fdate_str[xx])
  ftime_str[xx]=strmid(files[xx],nchardir+12,6)
  filetimes[xx]=long(ftime_str[xx])
endfor


;; UI to get sigmoid data

sig_id=''
noaa_id=''
print,''
print,'Please enter the ID of the sigmoid you want to work with.'
read,sig_id
print,'Please enter the NOAA ID# of the active region.'
print,'If there is no NOAA number associated with the AR, just hit <enter>.'
read,noaa_id
if (noaa_id ne '') then noaa_id='_AR'+noaa_id

print,''
print,'Now you must enter the date and time where you want to START measuring'
print,'the flux of the sigmoid.'
input_date,sigsd,error=err,minyr=2010
if (err) then message,'The program cannot proceed without a start date. Quitting...'
input_timeofday,sigst,error=err
if (err) then message,'The program cannot proceed without a start time. Quitting...'
sdat_int=long(sigsd.yearstr+sigsd.mnthstr+sigsd.daystr)
stim_int=long(sigst.hrstr+sigst.minstr+sigst.secstr)

print,''
print,'Now you must enter the date and time where you want to STOP measuring'
print,'the flux of the sigmoid.'
input_date,siged,error=err,minyr=sigsd.year
if (err) then message,'The program cannot proceed without an end date. Quitting...'
input_timeofday,siget,error=err
if (err) then message,'The program cannot proceed without an end time. Quitting...'
edat_int=long(siged.yearstr+siged.mnthstr+siged.daystr)
etim_int=long(siget.hrstr+siget.minstr+siget.secstr)

input_savepath,savepath,error=err, $
  defaultpath='/Volumes/Scratch/Users/ehanson/Sigmoid_data_by_region/'
if (err) then message,'The program needs a valid path to save plots to. Quitting...'


;; Identify files within desired time range

good=intarr(nfiles)
good[*]=0
for xx=0,nfiles-1 do begin
  startok=0
  endok=0
  if (filedates[xx] gt sdat_int) then begin
     startok=1
  endif else if (filedates[xx] eq sdat_int) then begin
     if (filetimes[xx] ge stim_int) then startok=1
  endif
  if (filedates[xx] lt edat_int) then begin
     endok=1
  endif else if (filedates[xx] eq edat_int) then begin
     if (filetimes[xx] le etim_int) then endok=1
  endif
  if ((startok) AND (endok)) then good[xx]=1
endfor
if (total(good) lt 1) then message,'No HMI files fall within your desired range. Quitting...'
keep=where(good eq 1)
files=files[keep]
nfiles=n_elements(files)
filedates=filedates[keep]
filetimes=filetimes[keep]
fdate_str=fdate_str[keep]
ftime_str=ftime_str[keep]


;; Continue execution of flux measurement code

read_sdo,files,index,data
;help,index,/st


sz=size(data)
n=sz(3)

posflux=dblarr(n)
negflux=dblarr(n)
r=dblarr(n)
crossed=intarr(nfiles)

for kk=0,nfiles-1 do begin
  sun=get_sun(index[kk].date_obs)
  r[kk]=sun[1]/index[0].cdelt1
  area=(index.cdelt1/206265.*sun[0]*149598000.*10^5.)^2.
  window,0,xs=1024,ys=1024
  tv,bytscl(congrid(data(*,*,kk),1024,1024),min=-200,max=200)
  print,''
  print,'MAGNETOGRAM  # '+strcompress(string(kk+1),/remove_all)+ $
     '   out of total   '+strcompress(string(nfiles),/remove_all)
  print,fdate_str[kk]+'     '+ftime_str[kk]

  print,'CLICK LOWER LEFT'
  cursor,a1,b1,/down,/device
  print,'CLICK UPPER RIGHT'
  cursor,a2,b2,/down,/device

  window,1,xs=(a2-a1)*4.+1.,ys=(b2-b1)*4.+1.
  jump1:
  tv,bytscl(data(4*a1:4*a2,4*b1:4*b2,kk),min=-200,max=200)

  ;; you can change the number of points where you click here
  xx=fltarr(10)
  yy=fltarr(10)
  print,'CLICK ON 10 POINTS AROUND THE REGION'
  for ii=0,9 do begin
     cursor,x,y,/down,/device
     xx[ii]=x
     yy[ii]=y
     print,'YOU HAVE CLICKED AT POINT #',ii+1
     if (ii gt 0) then begin 
       plots,xx[ii],yy[ii],/device,thick=5,color=256,/continue
     endif else begin  
       plots,xx[ii],yy[ii],/device,thick=5,color=256
     endelse
  endfor
  ;; and here +1

  x1=fltarr(11)
  y1=fltarr(11)
  x1[0:9]=xx
  x1[10]=xx[0]
  y1[0:9]=yy
  y1[10]=yy[0]

  centerx=0.5*(max(x1)+min(x1))
  if (centerx le index[kk].xcen) then begin
     crossed[kk]=0
  endif else begin
     crossed[kk]=1
  endelse

  plots,x1,y1,/device,thick=5,color=256

  p=''
  print,'Are you satisfied with the contour? Y/N'
  print,'(Enter B if this data file is bad.)'
  read,p
  if (strlowcase(p) eq 'n') then goto,jump1


  in=intarr((a2-a1)*4.+1.,(b2-b1)*4.+1.)
  for ii=0,4*(a2-a1) do begin
     for jj=0,4*(b2-b1) do begin
       in[ii,jj]=inside(ii,jj,x1,y1)
     endfor
  endfor

  posfluxm=dblarr((a2-a1)*4.+1.,(b2-b1)*4.+1.)
  negfluxm=dblarr((a2-a1)*4.+1.,(b2-b1)*4.+1.)

  if (strlowcase(p) ne 'b') then begin

     for ii=(4*a1),(4*a2)-1 do begin
       for jj=(4*b1),(4*b2)-1 do begin
         if (in[ii-(4*a1),jj-(4*b1)] eq 1) then begin
           if (data[ii,jj,kk] gt 0.) then begin
             l=(ii-(index[kk].crpix1))^2.+(jj-(index[kk].crpix2))^2.
             posfluxm[ii-(4*a1),jj-(4*b1)]=data[ii,jj,kk]/(1.-l/r[kk]^2.)*area
           endif else begin
             l=(ii-(index[kk].crpix1))^2.+(jj-(index[kk].crpix2))^2.
             negfluxm[ii-(4*a1),jj-(4*b1)]=data[ii,jj,kk]/(1.-l/r[kk]^2.)*area
           endelse
         endif
       endfor
     endfor
  endif else begin
     posfluxm[*]=0.0
     negfluxm[*]=0.0
  endelse

  posflux[kk]=total(posfluxm)
  negflux[kk]=total(negfluxm)

  tempfile=sig_id+'_'+fdate_str[0]+'_'+ftime_str[0]+'_peakflux_data_temp.sav'
  save,filename=savepath+tempfile, $
     posflux,negflux,files,sigsd,sigst,siged,siget,filedates,filetimes, $
     fdate_str,ftime_str,sig_id,noaa_id,kk
  message,tempfile+'   has been saved to   '+savepath,/cont
endfor

tempfile=sig_id+'_'+fdate_str[0]+'_'+ftime_str[0]+'_'+fdate_str[nfiles-1]+ $
  '_'+ftime_str[nfiles-1]+'_peakflux_data_basics.sav'
save,filename=savepath+tempfile, $
  posflux,negflux,files,sigsd,sigst,siged,siget,filedates,filetimes, $
  fdate_str,ftime_str,sig_id,noaa_id,kk
message,tempfile+'   has been saved to   '+savepath,/cont

compute_hours,fdate_str,ftime_str,day_of_year,hrs_past_mdnite,hrs_from_first

fluxfile=savepath+sig_id+'_'+fdate_str[0]+'_'+ftime_str[0]+'_'+fdate_str[nfiles-1]+ $
  '_'+ftime_str[nfiles-1]+'_peakflux_data_for_plotting.sav'
save,filename=fluxfile, $
  posflux,negflux,files,sigsd,sigst,siged,siget,filedates,filetimes, $
  sig_id,noaa_id,day_of_year,hrs_past_mdnite,hrs_from_first,crossed
message,fluxfile+'   has been saved!',/cont

;fluxplot,fluxfile,line_file=line_file,savepath=savepath


stop
end
