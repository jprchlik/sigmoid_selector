;=====================================
;Function
;    Creates a flux versuses time plot for the magnetic field under a sigmoid
;
;Usage
;    create_flux_plot,filein
;
;Input
;    filein -- A save file created by make_hmi_movie.pro
;    outdir -- Output directory for the plot
;    s_stim -- Start time of observed sigmiod 
;    e_stim -- End time of observed sigmiod 
;    f_stim -- Time of flare in sigmiod as a string that may pass to anytim
;    f_scls -- Class of flare in sigmiod
;
;Ouput 
;    fileout -- Full path to the created file
;
;=====================================

pro create_flux_plot,filein,outdir,s_stim,e_stim,f_stim=f_stim,f_scls=f_scls,fileout=fileout

;Setup plotting parameters
set_plot,'Z'
device,decomposed=0,set_pixel_depth=24,set_resolution=[1200,750]
loadct,12

;Add / to outdirectory
outdir = outdir+'/'

;Restore save file created by make_hmi_movie
restore,filein

;get anytime format of dates
time = anytim(obs_time)


;Get the filename part of the plot only
fileonly = strsplit(filein,'/',/extract)
fileonly = fileonly[n_elements(fileonly)-1]

;output png file
fileout = outdir+str_replace(fileonly,'sav','png')

;set up date format
dummy = LABEL_DATE(DATE_FORMAT=["%M-%D"])


;Store X,Y positions of sigmoid center
cx = 4096./2.
cy = 4096./2.
delt0 = 0.504297 
x = fltarr(n_elements(roi_save))
y = fltarr(n_elements(roi_save))

;get coordinates from ROI objects
for i=0,n_elements(roi_save)-1 do begin
    test = roi_save[i]->ComputeGeometry(centroid = cent) 
    x[i] = (cent[0]-cx)*delt0
    y[i] = (cent[1]-cy)*delt0
endfor

;compute distance from D.C
r = sqrt(x^2+y^2)

;get maximum allows radius using 60 deg cut off angle
sun_par = get_sun(time[0])
r_max = sun_par[1]*sin(50.*!dtor)

;Remove magnetic field measurements outsided alloted radius
good_par = where(r le r_max,good_cnt)


;create plot if there are good measurements
if good_cnt gt 0 then begin

    ;Get constance to convert from Gauss to Maxwell
    cont = (delt0/sun_par[1]*6.955E10)^2

    ;remove near limb observing times
    pos_ints = pos_ints[good_par] * cont * 1E-21
    neg_ints = neg_ints[good_par] * cont * 1E-21
    time     = time[good_par]

    ;combined pos and negative arrays
    com_ints = [pos_ints,-neg_ints]
     
    ;set limits for plot
    y_min = min(com_ints)
    y_max = max(com_ints)
    y_rng = y_max-y_min
    x_min = min(time)
    x_max = max(time)
    x_rng = x_max-x_min
    
    ;Set up plot
    utplot,[0,0],[0,0],'1-jan-79',ytitle="Flux [x10!U21!N Mx]",$
                XSTYLE=1,$;timerange=['24-aug-16,05:59:00','24-aug-16,8:00:00'],$
                xrange=[min(time)-3*3600.,max(time)+3*3600.],YSTYLE=1,$
                /nodata,yrange=[y_min-.1*y_rng,y_max+.1*y_rng],background=cgColor('white'),color=0,$
                charthick=3,charsize=2.5,xminor=12,xtitle='Time [UTC]' ;yrange=[80,120]
    
    ;Add time evolution of magnetic field strength
    oplot,time,pos_ints,linestyle=3,color=200,thick=7
    oplot,time,-neg_ints,linestyle=0,color=0,thick=7
    

    ;Add flares onto plot 
    if (keyword_set(f_stim)) and (keyword_set(f_scls)) then begin


        ;counter to vary the label y-coordinate
        y_off = 0
       
        ;Use line ID plot to plot the time
        for k=0,n_elements(f_stim)-1 do begin
            ; skip out of range flares
            if ((anytim(f_stim[k]) lt x_min) or (anytim(f_stim[k]) gt x_max)) then continue

            ;Plot flares
            oplot,[0.,0.]+anytim(f_stim[k]),[y_min-.1*y_rng,y_max+.1*y_rng],linestyle=2,thick=4,color=100
            xyouts,anytim(f_stim[k])-0.005*x_rng,y_min+y_rng*y_off/10.,f_scls[k],ORIENTATION=90,color=100, CHARSIZE=2, CHARTHICK=4
           
            ;incremenet y-counter by 1
            y_off +=1
            ;do not vary the offset by more than a factor of 10
            if y_off gt 10 then y_off = 0
        endfor
    endif


    ;Add polygon shading
    ;POLYFILL, X [, Y [, Z]] [, IMAGE_COORD=array] [, /IMAGE_INTERP] [, /LINE_FILL] [, PATTERN=array] [, SPACING=centimeters] [, TRANSPARENT=value]
    POLYFILL,anytim([s_stim,s_stim,e_stim,e_stim]),[y_min-.1*y_rng,y_max+.1*y_rng,y_max+.1*y_rng,y_min-.1*y_rng], $
          color=65, TRANSPARENT=0.5,/LINE_FILL,ORIENTATION=45

    write_png,fileout,tvrd(/true)


  

endif

    
end