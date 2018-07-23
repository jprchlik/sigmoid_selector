
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
;                   The format for the file is as follows:
;                   readcol,times,ID,RATING,NOAA,AR_START,X,Y,AR_END,SIG_START,SIG_END,TBEST,format='LL,I,A,A,F,F,A,A,A,A'
;    wavelnth   -   An array of wavelengths to look at filaments (Default = [171,304])
;    aia_arch   -   Archive to put and read in AIA files (Default = aia_arch)    
;    outf       -   Name of output file containing save information (Default = times file name with the csv replaced by a sav extension)
;                   fil_d={sig_id:'',            -->  Sigmiod ID (not positive that Patty kept this constant as I expected)     
;                           NOAA_id:0,           -->  NOAA AR ID 
;                           filename:'',         -->  AIA filename used for the analysis        
;                           date:'',             -->  Date of AIA observation
;                           leng:0.0,            -->  Length of filament in arcsec     
;                           device_arcsecx:0.0,  -->  Conversion from clicked points to arcsec in X    
;                           device_arcsecy:0.0,  -->  Conversion from clicked points to arcsec in Y    
;                           devicex:fltarr(100), -->  Clicked points in X   
;                           devicey:fltarr(100)} -->  Clicked points in Y
;                   
;
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
        devicex:fltarr(100),      $
        devicey:fltarr(100)}


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

;matched files
mat_f = strarr(n_elements(wavelnth))

;Download AIA data for all the best times
for i=0,n_elements(goodt)-1 do begin

    ;get index for a good time
    gi = goodt[i]
    ;get time range to search over
    t1 = tbest[gi]
    t2 = anytim(anytim(t1)+24.,/ecs)  
    
    ;Get all aia data in wavelength range
    s_f = vso_search(t1,t2,min_wav=strcompress(min(wavelnth),/remove_all),max_wav=strcompress(max(wavelnth),/remove_all),unit_wav='Angstrom',provider='jsoc')


    ; Leave in there are no matches
    ; But fill with none values
    if n_elements(size(s_f)) le 3 then begin
         mat_f = [[mat_f],['None'+strarr(n_elements(wavelnth))]]
         continue
    endif

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
    chk_i = col_i[1,*]

    ;Only get unique values for each wavelength
    uni_i = UNIQ(s_f[chk_i].wave.min, SORT(s_f[chk_i].wave.min))
    mat_i = chk_i[uni_i]

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

    ;Make sure all wavelengths  files are found
    if n_elements(o_str) lt n_elements(wavelnth) then mat_f = [[mat_f],$
        ['None'+strarr(n_elements(wavelnth)-n_elements(o_str)),o_str]] $
    else mat_f = [[mat_f],[o_str]] 

endfor

;set window size
wind_size = 512 ;1024
scale = 2


;Create new filaments array
sfilament=replicate(fil_d,2*n_elements(goodt))

;Init assume start counting from 0
start = 0

;Check if save file already exists
f_chck = file_test(outf)

;restart from last if file already exists
if f_chck eq 1 then begin
    restore,outf
    ;check to see filled sigmoid values
    if isa(nfilament) then  begin
        comp = where(nfilament.sig_id ne '',cnt_chck)
        start = cnt_chck/2
    endif
    ;check for old formatted sigmiod filaments
    if isa(filament) then  begin
        comp = where(filament.sig_id ne '',cnt_chck)
        save,filament,filename=outf+'.backup'
        start = 0
    endif
endif

;cludge to make this earlier mistake work
nfilament = sfilament

;Loop over all good tbest sigmoid times
for i=start,n_elements(goodt)-1 do begin
    
    ;get index for a good time
    gi = goodt[i]

    no_aia=0
    ;Get fits files from matched array
    fits_files=mat_f[*,i+1]
    ;Count number of files
    gfiles=where(fits_files ne 'None',nfiles)

    ;We have AIA Data
    if (nfiles eq 0) then begin
       no_aia=1
    endif else begin
       xwdw_size=wind_size*scale
       ywdw_size=wind_size*scale
    endelse
    if (no_aia) then begin
       print,'There are no matching AIA files'
       continue
    endif


    ;Only return good fits files
    fits_files = fits_files[gfiles]

    ;Prep all wavelength images
    aia_prep,fits_files,findgen(n_elements(fits_files)),index1,data1,/verbose

    ;loop over all wavelengths for each time
    for j=0,nfiles-1 do begin
        initialized=0
        done=0
        ;check to see if file is already analyzed 
        ;if so update parameters in new save file and return to start of loop
        if isa(filament) then begin
            print,fits_files[j]
            analyzed =  where((fits_files[j] eq filament.filename) and (ID[gi] eq filament.sig_id),cnt)
            if cnt eq 1 then begin
                nfilament[2*i+j].sig_id          = filament[analyzed].sig_id 
                nfilament[2*i+j].NOAA_id         = filament[analyzed].NOAA_id        
                nfilament[2*i+j].filename        = filament[analyzed].filename
                nfilament[2*i+j].date            = filament[analyzed].date
                nfilament[2*i+j].leng            = filament[analyzed].leng
                nfilament[2*i+j].device_arcsecx  = filament[analyzed].device_arcsecx
                nfilament[2*i+j].device_arcsecy  = filament[analyzed].device_arcsecy
                nfilament[2*i+j].devicex         = filament[analyzed].devicex
                nfilament[2*i+j].devicey         = filament[analyzed].devicey
                save,filename=outf,nfilament
                print,'HERE'
                continue
            endif
        endif
      
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

            ; Get coorinate limits for image plot
            limits = select_cutout(px,py,wind_size,img_xsize,img_ysize)   
            pxmin = limits[0]
            pxmax = limits[1]
            pymin = limits[2]
            pymax = limits[3]


            ;device to physical cooridinates
            devicex_to_imgx=float(xwdw_size)/(wind_size)
            devicey_to_imgy=float(ywdw_size)/(wind_size)
            arcsec_per_devicex=arcsec_per_pixel/devicex_to_imgx
            arcsec_per_devicey=arcsec_per_pixel/devicey_to_imgy
   
            ;Get image limits
            lim = cgPercentiles(data1(pxmin:pxmax,pymin:pymax),percentiles=[0.01,0.99])
        
            ;Plot image
            window,5,xs=xwdw_size,ys=ywdw_size
            img = bytscl(rebin(asinh(double(data1(pxmin:pxmax,pymin:pymax))),xwdw_size,ywdw_size),min=asinh(double(lim[0])),max=asinh(double(lim[1])))
            tv,img
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
                print,"Right click once if no filament"
                print,"Center click once to recenter window"
                cursor,lx1,ly1,/down,/device
                xyouts,lx1,ly1,'x',/device,alignment=0.5
               
                ;Right click to end
                if !MOUSE.button eq 4 then click = 0

                ;Middle click to recenter
                if !MOUSE.button eq 2 then begin
                   ; Get coorinate limits for image plot from center click
                   ;correct for the previous center location 
                   px = lx1/devicex_to_imgx+temporary(pxmin)
                   py = ly1/devicey_to_imgy+temporary(pymin)

                   ;Set limits for cutout
                   limits = select_cutout(px,py,wind_size,img_xsize,img_ysize)   
                   pxmin = limits[0]
                   pxmax = limits[1]
                   pymin = limits[2]
                   pymax = limits[3]

                   ;Get image limits
                   lim = cgPercentiles(data1(pxmin:pxmax,pymin:pymax),percentiles=[0.01,0.99])
        
                   ;Plot image
                   window,5,xs=xwdw_size,ys=ywdw_size
                   img = bytscl(rebin(asinh(double(data1(pxmin:pxmax,pymin:pymax))),xwdw_size,ywdw_size),min=asinh(double(lim[0])),max=asinh(double(lim[1])))
                   tv,img
 
                   ;reset filament length arrays
                   lx = []
                   ly = []

                   ;Exit current iteration without saving clicks
                   continue

                endif
                
                ;Add lx and ly values to array 
                lx = [lx,lx1]
                ly = [ly,ly1]
            endwhile

            ;If no filament create fillvalues
            if n_elements(ly) eq 1 then begin
                lr_phy = -9999.0 


            ;If a filament is present solve length
            endif else begin

                ;plot path along route
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
            endelse
        
            cont=''
            print,''
            print,'Shall we try again?'
            print,'Enter 1 for yes, any other key for no.'
            read,cont
            if (cont ne '1') then done=1
        endwhile

        ;Update information in save file
        nfilament[2*i+j].sig_id=ID[gi] 
        nfilament[2*i+j].NOAA_id=NOAA[gi]        
        nfilament[2*i+j].filename=fits_files[j]      
        nfilament[2*i+j].date=tbest[gi]   
        nfilament[2*i+j].leng=lr_phy         
        nfilament[2*i+j].device_arcsecx=arcsec_per_devicex
        nfilament[2*i+j].device_arcsecy=arcsec_per_devicey
        nfilament[2*i+j].devicex=lx
        nfilament[2*i+j].devicey=ly
 
        ;save sav file
        save,filename=outf,nfilament

    endfor
    
endfor
    
end
