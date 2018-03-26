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

function get_sigmoid_flares,obs_tim_s,obs_tim_e,obs_time_c,xbox,ybox,cx,cy,arnum,$
                            ffl_x,ffl_y,ffl_ts,ffl_te,ffl_tp,ffl_mx,ffl_cl,cme=cme

    ;Added best guess of NOAA number
    if keyword_set(cme) then begin 
        query=ssw_her_make_query(obs_tim_s,obs_tim_e,/ce)
    endif else begin
        ;search_goes = '&sparam0=OBS_Instrument&op0==&value0=GOES'
        query=ssw_her_make_query(obs_tim_s,obs_tim_e,/fl,search_array=['obs_observatory=GOES'],result_limit=5000)
        
        ;Add GOES class Instrument rescriction to query
        ;query = query;search_goes
        ;Get positions from AIA data 2018/03/26 J. Prchlik
        q_aia=ssw_her_make_query(obs_tim_s,obs_tim_e,/fl,search_array=['obs_observatory=SDO','obs_channelid=171'],result_limit=5000)
        aia=ssw_her_query(q_aia)
    endelse

    her=ssw_her_query(query) 
    if n_elements(size(her)) gt 3 then begin 
        ; Get time, postion and name
        if keyword_set(cme) then begin
            fl_x  = her.ce.event_coord1
            fl_y  = her.ce.event_coord2
            fl_u  = her.fl.EVENT_COORDUNIT
            fl_ts = her.ce.EVENT_STARTTIhME
            fl_te = her.ce.EVENT_ENDTIME
            fl_tp = her.ce.EVENT_PEAKTIME
            fl_mx = her.ce.FL_PEAKFLUX
            fl_cl = her.ce.FL_GOESCLS

        endif else begin
            fl_x  = her.fl.event_coord1
            fl_y  = her.fl.event_coord2
            fl_u  = her.fl.EVENT_COORDUNIT 
            fl_ts = her.fl.EVENT_STARTTIME
            fl_te = her.fl.EVENT_ENDTIME
            fl_tp = her.fl.EVENT_PEAKTIME
            fl_mx = her.fl.FL_PEAKFLUX
            fl_cl = her.fl.FL_GOESCLS

        endelse

        ;Rotation array
        rot_p = fltarr([2,n_elements(fl_x)])


        ;Rotate position to obs time
        for j=0, n_elements(fl_x)-1 do begin


            ;Update GOES flux with AIA coordinates if CME not set
            if keyword_set(cme) eq 0 then begin
                asdo_t=double(anytim(aia.fl.event_peaktime))
                goes_t=double(anytim(her.fl[j].event_peaktime))
                ;Get minimum time index
                getmin=min(abs(goes_t-asdo_t),imin)
  
                ;Only update pointing GOES pointing with AIA pointing if within a minute
                if min(abs(goes_t-asdo_t)) lt 1.*60. then begin
                    fl_x[j] = round(float(aia.fl[imin].event_coord1))
                    fl_y[j] = round(float(aia.fl[imin].event_coord2))
                    fl_u[j] = aia.fl[imin].event_coordunit
                endif

            endif

            ;Set up coordinate information to pass to rotation
            hpc_x = fl_x[j]
            hpc_y = fl_y[j]

            ;Get coordinates WCS value to send to coordinate coversion 
            WCS_CONV_FIND_DSUN, DSUN, RSUN, WCS=WCS,DATE_OBS=fl_tp[j] 
            ;If a heliographic projection convert to arcsec for consistency
            if fl_u[j] eq 'degrees' then WCS_CONV_HG_HPC,fl_x[j],fl_y[j],hpc_x,hpc_y,wcs=WCS ,/arcseconds
 
            ;reset fl coordiantes to be in arcsec
            fl_x[j] = round(float(hpc_x)) 
            fl_y[j] = round(float(hpc_y)) 

            ;rotate to observed central sigmoid time
            rot_p[*,j] = rot_xy(hpc_x, hpc_y, tstart=fl_tp[j], tend=obs_time_c)
        endfor

       
        ;Store x,y in separate array
        rot_x = rot_p[0,*]
        rot_y = rot_p[1,*]
 
        ;get flares inside box after rotation
        ;Compute x and y limit functions
        lims = comp_limits(xbox,ybox)
        
        ;Create array of 1 and 0 for box
        xmin = rot_x ge lims[1]*rot_y+lims[0]
        xmax = rot_x le lims[3]*rot_y+lims[2]
        ymin = rot_y le lims[5]*rot_x+lims[4]
        ymax = rot_y ge lims[7]*rot_x+lims[6]
        ;Add NOAA number deg along with inside box
        vbox = where(((xmin*xmax*ymin*ymax) or (fix(her.fl.ar_noaanum) eq fix(arnum))),cnt)

        
        ;Get the flares inside the sigmoid box
        if cnt gt 0 then begin
            ffl_x  = fl_x [vbox]
            ffl_y  = fl_y [vbox]
            ffl_ts = fl_ts[vbox]
            ffl_te = fl_te[vbox]
            ffl_tp = fl_tp[vbox]
            ffl_mx = fl_mx[vbox]
            ffl_cl = fl_cl[vbox]
        endif else begin
            ffl_x  = -1.e30
            ffl_y  = -1.e30
            ffl_ts = 'None'
            ffl_te = 'None'
            ffl_tp = 'None'
            ffl_mx = -1.e30
            ffl_cl = 'None'
 

        endelse

    endif else begin
       ffl_x  = -1.e30
       ffl_y  = -1.e30
       ffl_ts = 'None'
       ffl_te = 'None'
       ffl_tp = 'None'
       ffl_mx = -1.e30
       ffl_cl = 'None'


    endelse

    return,0 ;[[ffl_x],[ffl_y],[ffl_ts],[ffl_te],[ffl_tp],[ffl_mx],[ffl_cl]]
end

;sigloc is the directory location of save files and sigmoid analysis files
pro flare_cme_sigcat, sigloc,fname=fname
; Give program active region lifetime start and end times to get out flares and associated cmes

sigloc = sigloc+'/'

if keyword_set(fname) then fname = fname else fname = 'sigmoid_sizedata.sav'
;Get save file in directory
restore,sigloc+fname ;structure name is sigmoids

;Sort sigmioid IDs
sort_id = sort(sigmoids.sig_id)
ssig_id = sigmoids[sort_id].sig_id

;get unqiue sigmiod ids
uniq_id = uniq(ssig_id)

;get the sigmoid number for unique values
usig_id = ssig_id[uniq_id]

;loop over all sigmoid ids
for i=0,n_elements(usig_id)-1 do begin

    ;create variable for sigmoid id 
    sig_id = usig_id[i]
   
    ;find where sigmoid ids match the input value
    this_sig = sigmoids.sig_id eq sig_id


    ;Get maximum and minimum analysis data
    min_date = min(sigmoids[where(this_sig)].DATE)
    max_date = max(sigmoids[where(this_sig)].DATE)

    ;find where sigmoid nearest to center
    cntr_sig = abs(this_sig-1)*1e31+sigmoids.cx
    cntr_idx = where(cntr_sig eq min(cntr_sig),count)
    
    ;Get the time at central meridian
    cnt_date = sigmoids[cntr_idx].DATE
    cnt_x = sigmoids[cntr_idx].cx
    cnt_y = sigmoids[cntr_idx].cy
    ;Get the fits_header file information for the nearest to center sigmoid
    if count eq 1 then hdr = headfits(sigloc+sigmoids[cntr_idx].filename) $
    else continue

    ;get crval, crpix, cddelta keywords
    cdelt1 = sxpar(hdr,'cdelt1')
    cdelt2 = sxpar(hdr,'cdelt2')
    crval1 = sxpar(hdr,'crval1')
    crval2 = sxpar(hdr,'crval2')
    crpix1 = sxpar(hdr,'crpix1')
    crpix2 = sxpar(hdr,'crpix2')


    ;convert box values into physical coordinates
    xvals = cdelt1*(sigmoids[cntr_idx].bboxx-crpix1)+crval1
    yvals = cdelt2*(sigmoids[cntr_idx].bboxy-crpix2)+crval2


    ;first guess of rotation time to center
    fg = (0-cnt_x)/(10.)*3600. ;distance from center in arcsec and guess 10arcsec per hour from center
    ;find the closest value to when sigmoid is at the central meridian
    cpos = [cnt_x,cnt_y]
    spos = rot_xy(cnt_x,cnt_y,fg ,date=cnt_date)
    new_dat = anytim(cnt_date)+fg


    loop = 1
    counter = 0
    while loop eq 1 do begin
        case 1 of  
           ( abs(cpos[0]) gt abs(spos[0])): begin
               fg = (0-spos[0])/(10.)*3600.
               new_dat = new_dat+fg
               spos = rot_xy(spos[0],spos[1],fg ,date=new_dat)
           end
           ( abs(cpos[0]) lt abs(spos[0])): begin
               fg = (0-cpos[0])/(10.)*3600./2.
               spos = rot_xy(cpos[0],cpos[1],fg ,date=new_dat)
           end
        endcase
        
        if ((abs(spos[0]) lt 1.0) or (counter gt 100)) then loop = 0
        counter = counter+1
    endwhile
    
    ;Get +/- 7 days from disk center for sigmoid
    rng = 7.*3600.*24.
    min_date = anytim(new_dat-rng,/yymmdd)
    max_date = anytim(new_dat+rng,/yymmdd)
    ;Get flares around sigmoid Use +/- 7 days around disk center
    ;Comment out since LOCKHEED MARTIN is down for the day
    outvals = get_sigmoid_flares(min_date,max_date,cnt_date,xvals,yvals,cnt_x,cnt_y,igmoids[cntr_idx].NOAA_ID,
                                 ffl_x,ffl_y,ffl_ts,ffl_te,ffl_tp,ffl_mx,ffl_cl)

    ;Print test information
    print,'#############################################################'
    print,'#############################################################'
    print,sig_id,',',min_date,',',cnt_date,',',max_date
    for m=0,n_elements(ffl_x)-1 do print,ffl_x[m],' ,',ffl_y[m],' ,',ffl_ts[m],' ,',ffl_te[m],' ,',ffl_tp[m],' ,',ffl_mx[m],' ,',ffl_cl[m]

    print,counter,' ,Date = 20',anytim(new_dat,/yymmdd), '(x,y) = ',string(spos,format='(F6.1)'),' AR = ',string(sigmoids[cntr_idx].NOAA_ID,format='(I7)')
    print,'#############################################################'
    print,'#############################################################'

endfor


;if count eq 0 then print, 'No Flare Events Found For This Region Within the Time Specified'
;
;if count ne 0 then begin
;   selsize = size(sel, /n_elements)
;   fl_st = strarr(selsize)
;   fl_en = strarr(selsize)
;   fl_pk = strarr(selsize)
;   fl_cl = strarr(selsize)
;   for i=0, selsize-1 do begin
;      fl_st(i) = events(sel(i))._begin
;      fl_en(i) = events(sel(i))._end
;      fl_pk(i) = events(sel(i)).max
;      fl_cl(i) = events(sel(i)).class
;
;   endfor
;
;
;endif

; Now need to take that information and look for the CME links




;stop
end
