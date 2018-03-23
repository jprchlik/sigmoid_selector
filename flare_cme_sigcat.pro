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
    if keyword_set(cme) then $
        query=ssw_her_make_query(obs_tim_s,obs_tim_e,/ce,x1=cx,x2=cy)
    else $
        query=ssw_her_make_query(obs_tim_s,obs_tim_e,/fl,x1=cx,x2=cy)

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

        rot_p = fltarr([2,n_elements(ar_guess)])

        ;Rotate position to obs time
        for j=0, n_elements(ar_guess)-1 do $
            rot_p[*,j] = rot_xy(fl_x[j], fl_y[j], tstart=fl_t[j], tend=obs_time_c)

       
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
    min_date = min(sigmoids[this_sig].DATE)
    max_date = max(sigmoids[this_sig].DATE)

    ;find where sigmoid nearest to center
    cntr_sig = abs(this_sig-1)*1e31+sigmoids.cx
    cntr_idx = where(cntr_sig eq min(cntr_sig))
    
    ;Get the time at central meridian
    cnt_date = sigmoids[cntr_idx].DATE
    cnt_x = sigmoids[cntr_idx].cx
    cnt_y = sigmoids[cntr_idx].cy
    ;Get the fits_header file information for the nearest to center sigmoid
    hdr = headfits(sigloc+sigmoids[cntr_idx].filename)

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
    outvals = get_sigmoid_flares(min_date,max_date,cnt_date,xvals,yvals,cnt_x,cnt_y)

    print,sig_id,min_date,cnt_date,max_date
    for m=0,n_elements(outvals[*,0])-1 do print,outvals[m,*]


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




stop
end
