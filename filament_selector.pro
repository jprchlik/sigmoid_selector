

;#############################################################
;
;NAME:
;    filament_selector
;
;PURPOSE
;    Program to calculate length of sigmiod filament in AIA images
;
;CATEGORY:
;    Program, measurements
;
;USAGE
;    filament_selector,times,wavelnth=wavelnth,aia_arch='aia_arch/'
;
;INPUTS
;    times      -   A csv file containing times to analyze sigmoid filaments
;
;OUTPUTS
;    A save file called called filament_sizedata.sav
;
;#############################################################
pro filament_selector,times,wavelnth=wavelnth,aia_arch=aia_arch,outf=outf


;Read in file containing TBEST
readcol,times,ID,RATING,NOAA,AR_START,X,Y,AR_END,SIG_START,SIG_END,TBEST,format='LL,I,A,A,F,F,A,A,A,A'
;Set archive directory for download aia files
if keyword_set(aia_arch) then aia_arch = aia_arch else aia_arch = 'aia_arch/'
aia_arch = aia_arch+'/'

;if output file not set just use the same as the input file replaceing csv with sav
if not(keyword_set(outf)) then outf = str_replace(times,'.csv','.sav')

;wavelengths to retrieve for analysis
if not(keyword_set(wavelnth)) then wavelnth = [171,304]

;Good sigmoid tbest times
goodt = where(strlen(tbest) eq 23)

;Download AIA data for all the best times
for i=0,n_elements(goodt)-1 do begin

    ;get index for a good time
    gi = goodt[i]
    ;get time range to search over
    t1 = tbest[gi]
    t2 = anytim(anytim(t1)+20.,/ecs)  
    
    ;Get all aia data in wavelength range
    s_f = vso_search(t1,t2,min_wav=strcompress(min(wavelnth),/remove_all),max_wav=strcompress(max(wavelnth),/remove_all),unit_wav='Angstrom',provider='jsoc')


    ; Leave in there are no matches
    if n_elements(size(s_f)) le 3 then continue

    ;Create large arrays to match indices with the same wavelenth as given in wavelnth array
    s_fi = fix(s_f.wave.min ## (-1+fltarr(n_elements(wavelnth))))
    a_fi = fix(wavelnth # (1+fltarr(n_elements(s_f.wave.min))))
    ;Will produce zeros where wavelengths are the same as the requested ones
    t_fi = a_fi+s_fi

    ;Get where values are 0
    good_wav = where(t_fi eq 0)

    ;get two columns of indices where values are 0
    col_i = array_indices(t_fi,good_wav)

    ;get the indices that match wavelenth
    mat_i = col_i[1,*]

    ;Check if files already exist locally
    ;wavelength string
    w_str = strcompress(fix(s_f[mat_i].wave.min),/remove_all)
    ;time string
    t_str = str_replace(s_f[mat_i].time.start,':','_')
    ;searching string to see if files exist
    s_str = file_test(aia_arch+'aia.lev1.'+w_str+'A_'+t_str+'*.image_lev1.fits')

    ;files to get for download 
    d_ind = where(s_str eq 0,d_cnt)
    ;Download files which have the same wavelength as the requested wavelength
    if (d_cnt gt 0)  then s_r = vso_get(s_f[mat_i[d_ind]],out_dir=aia_arch)


endfor

;Loop over all good tbest sigmoid times
for i=0,n_elements(goodt)-1 do begin
    initialized=0
    done=0
    while not(done) do begin
    
      no_aia=0
      fits_path=aia_path
      fits_files=find_file(fits_path+'*_'+wavelnth+'.fits')
      nfiles=n_elements(fits_files)
      if ((nfiles eq 1) AND (fits_files[0] eq '')) then begin
         no_aia=1
      endif else begin
         xwdw_size=1024
         ywdw_size=1024
      endelse
      if (no_aia) then begin
         message,'There are no AIA files to match your date and hour.'
      endif else begin
         shortfiles=strarr(nfiles)
         filetimes=fltarr(nfiles)
         diff=fltarr(nfiles)
         for xx=0,nfiles-1 do begin
           shortfiles[xx]=strmid(fits_files[xx],strlen(fits_path)+12,6)
           filetimes[xx]=fix(shortfiles[xx])
           diff[xx]=abs(filetimes[xx]-fix(sigtime.hrstr+sigtime.minstr+sigtime.secstr))
         endfor
         least=min(diff,bestidx)
         best_file=fits_files[bestidx]
    ;     aia_prep,best_file,[0],index1,data1,/verbose
         aia_prep,fits_files,[bestidx],index1,data1,/verbose
         aia_lct,r,g,b,wavelnth=index1[0].wavelnth,/load
         window,5,xs=xwdw_size,ys=ywdw_size
         tv,bytscl(rebin(data1,xwdw_size,ywdw_size),min=5,max=700)
         print,''
         print,'Hopefully this file produces a good image.'
         print,"If not, don't worry, because you'll have an opportunity to"
         print,'do this again in a minute or so.'
         arcsec_per_pixel=index1[0].cdelt1 ;; Conversion for AIA data assuming square
         binscale=index1[0].chipsum
         img_xsize=index1[0].naxis1
         img_ysize=index1[0].naxis2
         devicex_to_imgx=float(xwdw_size)/(img_xsize*binscale)
         devicey_to_imgy=float(ywdw_size)/(img_ysize*binscale)
         arcsec_per_devicex=arcsec_per_pixel/devicex_to_imgx
         arcsec_per_devicey=arcsec_per_pixel/devicey_to_imgy
    
         print,''
         print,"Click on the filament"
         print,'Click now to designate Endpoint 1.'
         cursor,lx1,ly1,/down,/device
         xyouts,lx1,ly1,'.',/device,alignment=0.5
         print,'Click now to designate Endpoint 2.'
         cursor,lx2,ly2,/down,/device
         xyouts,lx2,ly2,'.',/device,alignment=0.5
    
         print,''
         print,"Now click on the endpoints of the sigmoid's shortest axis."
         print,'Click now to designate Endpoint 1.'
         cursor,sx1,sy1,/down,/device
         xyouts,sx1,sy1,'.',/device,alignment=0.5
         print,'Click now to designate Endpoint 2.'
         cursor,sx2,sy2,/down,/device
         xyouts,sx2,sy2,'.',/device,alignment=0.5
    
         lx_dev=abs(lx1-lx2)
         ly_dev=abs(ly1-ly2)
         lx_arc=arcsec_per_devicex*lx_dev
         ly_arc=arcsec_per_devicey*ly_dev
         long_axis_xy=sqrt((lx_dev^2)+(ly_dev^2))
         long_axis_arc=sqrt((lx_arc^2)+(ly_arc^2))
    
         sx_dev=abs(sx1-sx2)
         sy_dev=abs(sy1-sy2)
         sx_arc=arcsec_per_devicex*sx_dev
         sy_arc=arcsec_per_devicey*sy_dev
         short_axis_xy=sqrt((sx_dev^2)+(sy_dev^2))
         short_axis_arc=sqrt((sx_arc^2)+(sy_arc^2))
    
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
      endelse
    
    
      cont=''
      print,''
      print,'Shall we try again?'
      print,'Enter 1 for yes, any other key for no.'
      read,continue
      if (cont ne '1') then done=1
    endwhile
    
endfor
    
return
end
