





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




fil_d={sig_id:'',                $
        NOAA_id:0,               $
        filename:'',             $
        date:'',                 $
        leng:0.0,                $
        device_arcsecx:0.0,      $
        device_arcsecy:0.0,      $
        devicex:fltarr(10),      $
        devicey:fltarr(10)}


;set plot to X Window
set_plot,'X'

;Read in file containing TBEST
readcol,times,ID,RATING,NOAA,AR_START,X,Y,AR_END,SIG_START,SIG_END,TBEST,format='LL,I,A,A,F,F,A,A,A,A'
;Set archive directory for download aia files
if keyword_set(aia_arch) then aia_arch = aia_arch else aia_arch = 'aia_arch/'
aia_arch = aia_arch+'/'

;if output file not set just use the same as the input file replaceing csv with sav
if not(keyword_set(outf)) then outf = str_replace(times,'.csv','.sav')

;wavelengths to retrieve for analysis
if not(keyword_set(wavelnth)) then wavelnth = [171,304]

;Good sigmoid tbest times (i.e. contains time string)
goodt = where(strlen(tbest) eq 23)
;for testing purposes
goodt = [goodt[100]]

;matched files
mat_f = strarr(n_elements(wavelnth))

;Download AIA data for all the best times
for i=0,n_elements(goodt)-1 do begin

    ;get index for a good time
    gi = goodt[i]
    ;get time range to search over
    t1 = tbest[gi]
    t2 = anytim(anytim(t1)+12.,/ecs)  
    
    ;Get all aia data in wavelength range
    s_f = vso_search(t1,t2,min_wav=strcompress(min(wavelnth),/remove_all),max_wav=strcompress(max(wavelnth),/remove_all),unit_wav='Angstrom',provider='jsoc')


    ; Leave in there are no matches
    if n_elements(size(s_f)) le 3 then mat_f = [[mat_f],['None'+strarr(n_elements(wavelnth))]]

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

    ;output files
    o_str = file_search(aia_arch+'aia.lev1.'+w_str+'A_'+t_str+'*.image_lev1.fits')
    mat_f = [[mat_f],[o_str]] 

endfor

;set window size
wind_size = 1024


;Create new filaments array
filament=replicate(fil_d,n_elements(goodt))

;Init assume start counting from 0
start = 0

;Check if save file already exists
f_chck = file_test(outf)

;restart from last if file already exists
if f_chck eq 1 then begin
    restore,outf
    ;check to see filled sigmoid values
    comp = where(filament.sig_id ne '',cnt_chck)
    start = cnt_chck
endif


;Loop over all good tbest sigmoid times
for i=start,n_elements(goodt)-1 do begin
    
    ;get index for a good time
    gi = goodt[i]

    no_aia=0
    ;Get fits files from matched array
    fits_files=mat_f[*,i+1]
    ;Count number of files
    nfiles=n_elements(fits_files)

    ;We have AIA Data
    if ((nfiles eq 1) AND (fits_files[0] ne 'None')) then begin
       no_aia=1
    endif else begin
       xwdw_size=wind_size
       ywdw_size=wind_size
    endelse
    if (no_aia) then begin
       print,'There are no matching AIA files'
       continue
    endif


    ;Prep all wavelength images
    aia_prep,fits_files,findgen(nfiles),index1,data1,/verbose

    ;loop over all wavelengths for each time
    for j=0,nfiles-1 do begin
        initialized=0
        done=0
        ;Try to select filament
        while not(done) do begin
             ;Load color table for wavelength
             aia_lct,r,g,b,wavelnth=index1[j].wavelnth,/load

             ;Get AIA data header infromation
             arcsec_per_pixel=index1[j].cdelt1 ;; Conversion for AIA data assuming square
             ;binscale  = index1[j].chipsum
             binscale  = 1 ; Assume all full res for now
             img_xsize = index1[j].naxis1
             img_ysize = index1[j].naxis2
             img_xphys = index1[j].crval1
             img_yphys = index1[j].crval2
             img_xcntr = index1[j].crpix1
             img_ycntr = index1[j].crpix2


             ;convert x,y best coordinates into pixel values
             px = (x[gi]+img_xphys)/arcsec_per_pixel+img_xcntr
             py = (y[gi]+img_yphys)/arcsec_per_pixel+img_ycntr

             ;compute the min and max x and y pixel ranges
             pxmin = px-(wind_size/2.) 
             pxmax = px+(wind_size/2.)-1
             pymin = py-(wind_size/2.) 
             pymax = py+(wind_size/2.)-1 
          
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
                 (pxmax gt img_xsize):  begin
                     offset = img_xsize-pxmax
                     pxmax = img_xsize-1
                     pxmin = pxmin+offset
                 end
                 (pymax gt img_ysize):  begin
                     offset = img_ysize-pymax
                     pymax = img_ysize-1
                     pymin = pymin+offset
                 end
                 else: square= 0
             endcase 

             ;force integers
             pxmin = fix(pxmin)
             pxmax = fix(pxmax)
             pymin = fix(pymin)
             pymax = fix(pymax)

             ;device to physical cooridinates
             devicex_to_imgx=float(xwdw_size)/(img_xsize*binscale)
             devicey_to_imgy=float(ywdw_size)/(img_ysize*binscale)
             arcsec_per_devicex=arcsec_per_pixel/devicex_to_imgx
             arcsec_per_devicey=arcsec_per_pixel/devicey_to_imgy
        
             ;Plot image
             window,5,xs=xwdw_size,ys=ywdw_size
             tv,bytscl(rebin(data1(pxmin:pxmax,pymin:pymax),xwdw_size,ywdw_size),min=5,max=700)
             ;print,''
             ;print,'Hopefully this file produces a good image.'
             ;print,"If not, don't worry, because you'll have an opportunity to"
             ;print,'do this again in a minute or so.'

             lx = []
             ly = []

             ;allow clicks along filament
             click = 1
             while click do begin

                 print,''
                 print,"Click continuously along the filament (right click to end)"
                 cursor,lx1,ly1,/down,/device
                 xyouts,lx1,ly1,'x',/device,alignment=0.5
                
                 if !MOUSE.button eq 4  then click = 0
                 
                 ;Add lx and ly values to array 
                 lx = [lx,lx1]
                 ly = [ly,ly1]
             endwhile

             plots,lx,ly,color=255,/device,linestyle=1,thick=2

             ;get running difference of x values
             dx = ts_diff(lx,1)
             dy = ts_diff(ly,1)
        
             ;total distance in device coordinates
             lr_dev=sqrt(total(dx^2+dy^2))
             ;total distance in physical cooridinates       
             lr_phy=sqrt(total((arcsec_per_devicex*dx)^2+(arcsec_per_devicey*dy)^2))

             print,''
             print,'The axis size in DEVICE-units is:'
             print,strcompress(string(lr_dev),/remove_all)
        
             print,''
             print,'The axis size in ARCSEC is:'
             print,strcompress(string(lr_phy),/remove_all)
        
        
          cont=''
          print,''
          print,'Shall we try again?'
          print,'Enter 1 for yes, any other key for no.'
          read,cont
          if (cont ne '1') then done=1
        endwhile
readcol,times,ID,RATING,NOAA,AR_START,X,Y,AR_END,SIG_START,SIG_END,TBEST,format='LL,I,A,A,F,F,A,A,A,A'

        ;Update information in save file
        filament[i].sig_id=ID[gi] 
        filament[i].NOAA_id=NOAA[gi]        
        filament[i].filename=fits_files[j]      
        filament[i].date=tbest[gi]   
        filament[i].leng=lr_phy         
        filament[i].device_arcsecx=arcsec_per_devicex
        filament[i].device_arcsecy=arcsec_per_devicey
        filament[i].devicex=lx
        filament[i].devicey=ly
 
        ;save sav file
        save,filename=outf,filament

    endfor
    
endfor
    
end
