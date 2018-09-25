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
;
;Ouput 
;    fileout -- Full path to the created file
;
;=====================================

pro create_flux_plot,filein,outdir,fileout=fileout

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

;output png file
fileout = outdir+str_replace(filein,'sav','png')

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
r_max = sun_par[1]*sin(60.*!dtor)

;Remove magnetic field measurements outsided alloted radius
good_par = where(r le r_max,good_cnt)


;create plot if there are good measurements
if good_cnt gt 0 then begin


    ;remove near limb observing times
    pos_ints = pos_ints[good_par]
    neg_ints = neg_ints[good_par]
    time     = time[good_par]

    ;combined pos and negative arrays
    com_ints = [pos_ints,-neg_ints]
     
    ;set limits for plot
    y_min = min(com_ints)
    y_max = max(com_ints)
    y_rng = y_max-y_min
    
    utplot,[0,0],[0,0],'1-jan-79',ytitle="Mag. Field Strength [Gauss]",$
                XSTYLE=1,$;timerange=['24-aug-16,05:59:00','24-aug-16,8:00:00'],$
                xrange=[min(time)-3*3600.,max(time)+3*3600.],YSTYLE=1,$
                /nodata,yrange=[y_min-.1*y_rng,y_max+.1*y_rng],background=cgColor('white'),color=0,$
                charthick=3,charsize=2.5,xminor=12,xtitle='Time [UTC]' ;yrange=[80,120]
    
    ;Add time evolution of magnetic field strength
    oplot,time,pos_ints,linestyle=3,color=200,thick=3
    oplot,time,-neg_ints,linestyle=0,color=0,thick=3
    
    write_png,fileout,tvrd(/true)
endif
    
end