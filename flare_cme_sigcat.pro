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

function get_sigmoid_flares,obs_tim_s,obs_tim_e,obs_time_c,xbox,ybox,cx,cy,cme=cme

    ;Added best guess of NOAA number
    if keyword_set(cme) then begin 
        query=ssw_her_make_query(obs_tim_s,obs_tim_e,/ce,x1=cy-100.,x2=cy+100.)
    endif else begin
        search_goes = '&sparam0=OBS_Instrument&op0==&value0=GOES'
        query=ssw_her_make_query(obs_tim_s,obs_tim_e,/fl,y1=cy-100.,y2=cy+100.)
        
        ;Add GOES class Instrument rescriction to query
         query = query+search_goes
    endelse

    her=ssw_her_query(query,/str) 
    if n_elements(size(her)) gt 3 then begin 
        ; Get time, postion and name
        if keyword_set(cme) then begin
            fl_x  = her.ce.required.event_coord1
            fl_y  = her.ce.required.event_coord2
            fl_ts = her.ce.required.EVENT_STARTTIME
            fl_te = her.ce.required.EVENT_ENDTIME
            fl_tp = her.ce.required.EVENT_PEAKTIME
            fl_mx = her.ce.optional.FL_PEAKFLUX
            fl_cl = her.ce.optional.FL_GOESCLS

        endif else begin
            fl_x  = her.fl.required.event_coord1
            fl_y  = her.fl.required.event_coord2
            fl_ts = her.fl.required.EVENT_STARTTIME
            fl_te = her.fl.required.EVENT_ENDTIME
            fl_tp = her.fl.required.EVENT_PEAKTIME
            fl_mx = her.fl.optional.FL_PEAKFLUX
            fl_cl = her.fl.optional.FL_GOESCLS

        endelse

        rot_p = fltarr([2,n_elements(fl_x)])

        ;Rotate position to obs time
        for j=0, n_elements(ar_guess)-1 do $
            rot_p[*,j] = rot_xy(fl_x[j], fl_y[j], tstart=fl_tp[j], tend=obs_time_c)

       
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
        vbox = xmin*xmax*ymin*ymax
        
       ;Get the flares inside the sigmoid box
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
       ffl_ts = -1.e30
       ffl_te = -1.e30
       ffl_tp = -1.e30
       ffl_mx = -1.e30
       ffl_cl = -1.e30


    endelse

    jlkjlk = jjlkjkl
    return,[[ffl_x],[ffl_y],[ffl_ts],[ffl_te],[ffl_tp],[ffl_mx],[ffl_cl]]
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

    ;Get flares around sigmoid
    ;Comment out since LOCKHEED MARTIN is down for the day
    outvals = get_sigmoid_flares(min_date,max_date,cnt_date,xvals,yvals,cnt_x,cnt_y)

    print,'#############################################################'
    print,sig_id,',',min_date,',',cnt_date,',',max_date
    print,sigmoids[where(this_sig)].sig_id
    for m=0,n_elements(outvals[*,0])-1 do print,outvals[m,*]

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
    print,counter,' ,Date = 20',anytim(new_dat,/yymmdd), '(x,y) = ',string(spos,format='(F6.1)'),' AR = ',string(sigmoids[cntr_idx].NOAA_ID,format='(I7)')
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
