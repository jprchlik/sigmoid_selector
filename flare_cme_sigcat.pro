

pro compile_structure

    test_sig_b = CREATE_STRUCT('sigmoid_id',0,$
                             'cross_m','                    ',$; Time flare crossed the meridian 
                             'sigmd_s','                    ',$; Sigmoid start time 
                             'sigmd_e','                    ',$; Sigmoid end time 
                             'flare_x',fltarr(1000),$;FLARE X POSITION
                             'flare_y',fltarr(1000),$;FLARE Y POSITION
                             'flare_s',strarr(1000),$;FLARE Start time
                             'flare_e',strarr(1000),$;FLARE End time
                             'flare_p',strarr(1000),$;FLARE Peak time
                             'flare_c',strarr(1000),$;FLARE GOES Class
                             'cme_x',fltarr(100),$;CME X POSITION
                             'cme_y',fltarr(100),$;CME Y POSITION
                             'cme_s',strarr(100),$;CME Start time
                             'cme_e',strarr(100),$;CME End time
                             'cme_w',strarr(100),$;CME Angular Width
                             'cme_v',strarr(100)) ;CME Velocity



end

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
        ;Get LASCO CMES
        query=ssw_her_make_query(obs_tim_s,obs_tim_e,/ce,search_array=['obs_observatory=LASCO','FRM_NAME=CACTus+(Computer+Aided+CME+Tracking)'])
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

            ;Only get CME size greater than 45 Deg
            cme_clip = where(her.ce.CME_ANGULARWIDTH gt 45)
            ffl_x  = her.ce[cme_clip].HPC_X
            ffl_y  = her.ce[cme_clip].HPC_Y
            ffl_u  = her.ce[cme_clip].EVENT_COORDUNIT
            ffl_ts = her.ce[cme_clip].EVENT_STARTTIME
            ffl_te = her.ce[cme_clip].EVENT_ENDTIME
            ffl_tp = her.ce[cme_clip].EVENT_STARTTIME
            ffl_mx = her.ce[cme_clip].CME_ANGULARWIDTH
            ffl_cl = her.ce[cme_clip].CME_RADIALLINVEL

        endif else begin
            fl_x  = her.fl.hpc_x
            fl_y  = her.fl.hpc_y
            fl_u  = her.fl.EVENT_COORDUNIT 
            fl_ts = her.fl.EVENT_STARTTIME
            fl_te = her.fl.EVENT_ENDTIME
            fl_tp = her.fl.EVENT_PEAKTIME
            fl_mx = her.fl.FL_PEAKFLUX
            fl_cl = her.fl.FL_GOESCLS


            ;Moved inside flare loop only 2018/03/29 J. Prchlik
            ;Rotation array
            rot_p = fltarr([2,n_elements(fl_x)])


            ;Rotate position to obs time
            for j=0, n_elements(fl_x)-1 do begin

                ;Update GOES flux with AIA coordinates if AIA flare found
                if n_elements(aia) gt 3 then begin
                    asdo_t=double(anytim(aia.fl.event_peaktime))
                    goes_t=double(anytim(her.fl[j].event_peaktime))
                    ;Get minimum time index
                    getmin=min(abs(goes_t-asdo_t),imin)

                    ;Use GOES information as degfault
                    fl_u[j] = 'arcseconds'
  
                    ;Only update pointing GOES pointing with AIA pointing if within a minute 
                    ; and GOES position is 0
                    if ((min(abs(goes_t-asdo_t)) lt 1.*60.) and (fl_x[j] eq 0))  then begin
                        fl_x[j] = round(float(aia.fl[imin].event_coord1))
                        fl_y[j] = round(float(aia.fl[imin].event_coord2))
                        fl_u[j] = 'degrees'
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
            ;change to using ROI J. Prchlik 2018/07/11
            ;xmin = rot_x ge lims[1]*rot_y+lims[0]
            ;xmax = rot_x le lims[3]*rot_y+lims[2]
            ;ymin = rot_y le lims[5]*rot_x+lims[4]
            ;ymax = rot_y ge lims[7]*rot_x+lims[6]

            ;compute sigmoid roi box in physical coordinates
            roi_phy = OBJ_NEW('IDLanROI',xbox,ybox)

            ;check if flare point is in ROI object
            pnt_chk = roi_phy -> containsPoints(rot_x,rot_y)

            ;Add NOAA number deg along with inside ROI box
            vbox = where((pnt_chk and (fix(her.fl.ar_noaanum) eq 0)) or (fix(her.fl.ar_noaanum) eq fix(arnum)),cnt)

            
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


;####################################################
;FUNCTION
;Returns full path for file in local XRT archive given XRT filename
;
;USAGE
;  res = return_full_path(xrtf)
;####################################################
function return_full_path,xrtf

    ;parse fname to get the full path
    fyear = strmid(xrtf,6,4)
    fmont = strmid(xrtf,10,2)
    fdate = strmid(xrtf,12,2)
    fhour = 'H'+strmid(xrtf,15,2)+'00'

    res = '/'+fyear+'/'+fmont+'/'+fdate+'/'+fhour+'/'
    return,res
end

;####################################################################
;
;NAME:
;    flare_cme_sigcat
;
;PURPOSE:
;    Program to determine CME and flare in the same time and space range 
;    as the sigmoid. It uses the sigmoid box and AR number given in 
;    sigmoidsize_adv.
;
;USAGE:
;    flare_cme_sigcat, sigloc,fname=fname,odir=odir
;
;INPUTS:
;    sigloc     - Directory containting a sav file exported from sigmoidsize_adv
;    year       - String of year of observation considering cosistent format (e.g. '2017')
;    odir       - Output directory of save file created by flare_cme_sigcat (Default = sigloc)
;
;OUTPUTS:
;    An IDL save file with the format '("sigmoid_id_",I03,"_",I03,".sav")', 
;    where I03 are the first and last sigmoid ids in the list 
;    (e.g. sigmoid_id_001_010.sav)
;
;
;####################################################################
pro flare_cme_sigcat, sigloc,year,odir=odir,sigcat=sigcat
; Give program active region lifetime start and end times to get out flares and associated cmes

sigloc = sigloc+'/'

;Output directory
if keyword_set(odir) then odir = odir else odir = sigloc
;Text file with sigmoid catalog
if keyword_set(sigcat) then sigcat = sigcat else sigcat = 'SigmoidCatalogAll_filament.csv' 

;directory location in the hinode archive
xrt_arch = '/archive/hinode/xrt/level1/'


;Filename for sigmoid position and size
fname = "sigmoid_sizedata"+year+".sav"
;Get save file in directory
restore,sigloc+fname ;structure name is sigmoids

;Read in csv file to match with save file
cname = "Sigmoids"+year+".csv"
readcol,sigloc+cname,format='I,I,A,A,F,F,A,A,A,A',real_sig_id,rating,noaa,ar_start,b_x,b_y,AR_END,SIG_START,SIG_END,TBEST
;Switch to final sigmoid catalog
;Does not work because X,Y coordinates are no good with Tbest or Tobs in catalog
;readcol,sigcat,format='I,I,A,A,F,F,A,A,A,A,A',dum,real_sig_id,NOAA,AR,AR_START,b_X,b_Y,AR_END,SIG_START,SIG_END,lifetime,TBEST_old,TBEST,/preserve

;convert NA to zeors in noaa
noaa[where(noaa eq 'NA')] = '0'
noaa = fix(noaa)

;Sort sigmioid IDs
sort_id = sort(float(sigmoids.sig_id))
ssig_id = sigmoids[sort_id].sig_id

;get unqiue sigmiod ids
uniq_id = uniq(ssig_id)

;get the sigmoid number for unique values
usig_id = ssig_id[uniq_id]



;loop over all sigmoid ids
for i=0,n_elements(usig_id)-1 do begin

    ;create variable for sigmoid id 
    sig_id = usig_id[i]

    ;leave if sig id is a long string (i.e. file skipped)
    type_sig = strlen(sig_id)
    if type_sig ge 4 then continue
   
    ;find where sigmoid ids match the input value
    this_sig = sigmoids.sig_id eq sig_id


    ;Get maximum and minimum analysis data
    min_date = min(sigmoids[where(this_sig)].DATE)
    max_date = max(sigmoids[where(this_sig)].DATE)

    ;find where sigmoid nearest to center
    ;Fixed using the absolure value of x insteand of just X 
    cntr_sig = abs(this_sig-1)*1e7+abs(sigmoids.cx)
    cntr_idx = where(cntr_sig eq min(cntr_sig),count)



    
    
    ;Get the fits_header file information for the nearest to center sigmoid
    if count eq 1 then fname = sigmoids[cntr_idx].filename
    if count gt 1 then begin 
        cntr_idx = cntr_idx[0]
        fname = sigmoids[cntr_idx].filename
    endif
 
    ;Leave if no files found
    if count eq 0 then continue


    ;Make sure filename is in correct format
    if strlen(fname) ne 28 then begin
        ;split string into comps
        comps = strsplit(fname,'.',/extract)

        ;fix bad ms in filename
        if strlen(comps[1]) gt 1 then begin
            comps[1] = strmid(comps[1],0,1)
            fname = strjoin(comps,'.')
        endif

    endif

    ;Get full path to file
    full_path  = return_full_path(fname)

    ;get header information
    hdr = headfits(xrt_arch+full_path+fname,ERRMSG=ERRMSG)

    ;fix for a bad filename 
    if strlen(ERRMSG) gt 1 then begin
        stop
    end 

    ;Get the time at central meridian
    cnt_date = sigmoids[cntr_idx].DATE
    cnt_x = sigmoids[cntr_idx].cx
    cnt_y = sigmoids[cntr_idx].cy


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

    ;Less than 2 arsecond from center just call it the center 2018/07/03 (prevents MORINT error)
    if abs(cnt_x) lt 2 then begin
        loop = 0 
        new_dat = anytim(cnt_date)
        spos = cpos
    endif else begin
        spos = rot_xy(cnt_x,cnt_y,fg ,date=cnt_date)
         new_dat = anytim(cnt_date)+fg
        ;Loop to find merdian crossing
         loop = 1
         counter = 0
    endelse


    ;Loop to find where sigmoid crosses central merdian
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
    outvals = get_sigmoid_flares(min_date,max_date,cnt_date,xvals,yvals,cnt_x,cnt_y,sigmoids[cntr_idx].NOAA_ID,$
                                 ffl_x,ffl_y,ffl_ts,ffl_te,ffl_tp,ffl_mx,ffl_cl)

    ;Get CME values
    outvals = get_sigmoid_flares(min_date,max_date,cnt_date,xvals,yvals,cnt_x,cnt_y,sigmoids[cntr_idx].NOAA_ID,$
                                 cme_x,cme_y,cme_ts,cme_te,cme_tp,cme_da,cme_vl,/cme)


    ;Create dummy arrays to fill in structure
    cross_m = '20'+str_replace(anytim(new_dat,/yymmdd),', ','T');merdian crossing time
    flare_x = fltarr(1000)-1.e31;FLARE X POSITION
    flare_y = fltarr(1000)-1.e31;FLARE Y POSITION
    flare_u = strarr(1000);FLARE POSITION Units
    flare_s = strarr(1000);FLARE Start time
    flare_e = strarr(1000);FLARE End time
    flare_p = strarr(1000);FLARE Peak time
    flare_c = strarr(1000);FLARE GOES Class
    cmevl_x = fltarr(100)-1.e31;CME X POSITION
    cmevl_y = fltarr(100)-1.e31;CME Y POSITION
    cmevl_u = strarr(100);CME POSITION Units
    cmevl_s = strarr(100);CME Start time
    cmevl_e = strarr(100);CME End time
    cmevl_w = strarr(100);CME Angular Width
    cmevl_v = strarr(100) ;CME Velocity


    ;Match Sigmiod in save file with sigmiod in text file
    test_ar = where(this_sig)
    best_ind = where(fix(noaa) eq fix(sigmoids[test_ar[0]].NOAA_ID), count_match)

    ;If the Sigmiod has no AR number try matching based on position
    if ((fix(sigmoids[test_ar[0]].NOAA_ID) eq 0) or (count_match ne 1)) then begin

    ;get location of nearest sigmoids in catalog csv file if no AR number exists
        ;compute sigmoid roi box in physical coordinates to time nearest D.C.
        roi_phy = OBJ_NEW('IDLanROI',xvals,yvals)


        ;Add NOAA number deg along with inside ROI box
        dif_pos_x = fltarr(n_elements(b_x))-9999.9
        dif_pos_y = fltarr(n_elements(b_x))-9999.9
        dif_pos_r = fltarr(n_elements(b_x))-9999.9

        ;Rotate point to best observed time with bounding box for each sigmiod in catalog
        for k=0,n_elements(b_x)-1 do begin 
           ;Do not allow sigmoid to go over the limb for comp (Give 1.5 week time pad
           dif_pos_t = abs((anytim(tbest[k])-anytim(cnt_date))/(3600.*24.*7.))

           ;remove values when over 1.5 weeks away
           if dif_pos_t gt 1.5 then continue

           ;store rotated values
           cat_pos = rot_xy(b_x[k],b_y[k],tstart=tbest[k],tend=cnt_date,offlimb=testlimb,error=testerror) 
           dif_pos_x[k] = cat_pos[0]
           dif_pos_y[k] = cat_pos[1] 
           dif_pos_r[k] = sqrt((cat_pos[0]-cnt_x)^2+(cat_pos[1]-cnt_y)^2)
        endfor

        ;check if flare point is in ROI object
        pnt_chk = roi_phy -> containsPoints(dif_pos_x,dif_pos_y)

        ;Get index of nearest sigmoid
        best_ind = where(pnt_chk,pnt_cnt)

        print,'NO NOAA NUM'
        print,best_ind
        print,fix(real_sig_id[best_ind[0]])
        ;If no matching sigmoid found find closest value in catalog
        if pnt_cnt ne 1 then begin
           best_ind = where(dif_pos_r eq min(dif_pos_r),cnt_min)
           ;Give up if no unique solution found
           if cnt_min ne 1 then continue
        endif
    endif
 

    ;Print test information
    ;print,'#############################################################'
    ;print,'#############################################################'
    ;print,sig_id,',',min_date,',',cnt_date,',',max_date
    ;print,'Flares'
    ;Create fixed width arrays containing flare info
    for m=0,n_elements(ffl_x)-1 do begin
        flare_x[m] = ffl_x[m]
        flare_y[m] = ffl_y[m]
        flare_s[m] = ffl_ts[m]
        flare_e[m] = ffl_te[m]
        flare_p[m] = ffl_tp[m]
        flare_c[m] = ffl_cl[m]
    endfor

    ;print,counter,' ,Date = 20',anytim(new_dat,/yymmdd), '(x,y) = ',string(spos,format='(F6.1)'),' AR = ',string(sigmoids[cntr_idx].NOAA_ID,format='(I7)')

    ;Create fixed width arrays containing CME info
    for m=0,n_elements(cme_x)-1 do begin
        cmevl_x[m] =cme_x[m]
        cmevl_y[m] =cme_y[m]
        cmevl_s[m] =cme_ts[m]
        cmevl_e[m] =cme_te[m]
        cmevl_w[m] =cme_da[m]
        cmevl_v[m] =cme_vl[m]
     endfor
     ;stop
     ;print,'#############################################################'
     ;print,'#############################################################'
     ;Create single row in structure
     tmp = {test_sig_b,$
         sigmoid_id:fix(real_sig_id[best_ind[0]]),$ ;Use index in sav file to call real ID in csv file input to make save file
         sigmd_s:SIG_START[best_ind[0]],$ ;Use index in sav file to call sigmiod start time in csv file input to make save file
         sigmd_e:SIG_END[best_ind[0]],$ ;Use index in sav file to call sigmoid end time in csv file input to make save file
         cross_m:cross_m,$; Time flare crossed the meridian 
         flare_x:flare_x,$;FLARE X POSITION
         flare_y:flare_y,$;FLARE Y POSITION
         flare_s:flare_s,$;FLARE Start time
         flare_e:flare_e,$;FLARE End time
         flare_p:flare_p,$;FLARE Peak time
         flare_c:flare_c,$;FLARE GOES Class
         cme_x:cmevl_x,$;CME X POSITION
         cme_y:cmevl_y,$;CME Y POSITION
         cme_s:cmevl_s,$;CME Start time
         cme_e:cmevl_e,$;CME End time
         cme_w:cmevl_w,$;CME Angular Width
         cme_v:cmevl_v} ;CME Velocity


     ;Create large structure
     ;Added check if variable exists instead of just the first index for more versitility
     if isa(big_str) eq 0 then big_str = tmp $
     else big_str = [big_str,tmp] 
endfor

;Save large structure to file and encopass the range of simoid IDs in save file
outf = '("sigmoid_id_",I04,".sav")'
save,big_str,filename=string([fix(year)],format=outf)


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




end
