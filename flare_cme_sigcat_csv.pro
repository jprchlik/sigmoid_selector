

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

            endfor

       
            
            ;Add NOAA number deg along with inside ROI box
            vbox = where((fix(her.fl.ar_noaanum) eq fix(arnum)),cnt)

            
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
;    flare_cme_sigcat_csv
;
;PURPOSE:
;    Program to determine CME and flare in the same time and space range 
;    as the sigmoid. It uses the sigmoid box and AR number given in 
;    sigmoidsize_adv.
;
;USAGE:
;    flare_cme_sigcat_csv,sigcat=sigcat
;
;INPUTS:
;    sigcat     - Directory containting a sav file exported from sigmoidsize_adv
;
;OUTPUTS:
;    An IDL save file with the format '("sigmoid_id_",I03,"_",I03,".sav")', 
;    where I03 are the first and last sigmoid ids in the list 
;    (e.g. sigmoid_id_001_010.sav)
;
;
;####################################################################
pro flare_cme_sigcat_csv,sigcat=sigcat


;Text file with sigmoid catalog
if keyword_set(sigcat) then sigcat = sigcat else sigcat = 'SigmoidCatalogAll_filament.csv' 


;Switch to final sigmoid catalog
;Does not work because X,Y coordinates are no good with Tbest or Tobs in catalog
;They are now as fixed by J. Prchlik 2019/02/20
readcol,sigcat,format='I,I,A,A,A,A,A,A,A,A,A',dum,real_sig_id,NOAA,AR,AR_START,b_X,b_Y,AR_END,SIG_START,SIG_END,lifetime,TBEST_old,TBEST,/preserve

;convert NA to zeors in noaa
noaa[where(noaa eq 'NA')] = '0'
noaa = fix(noaa)


;loop over all sigmoid ids
for i=0,n_elements(noaa)-1 do begin

    ;create variable for sigmoid id 
    sig_id = real_sig_id[i]

    ;leave if sig id is a long string (i.e. file skipped)
    ;type_sig = strlen(sig_id)
    ;if type_sig ge 4 then continue
   
    ;find where sigmoid ids match the input value
    ;this_sig = sigmoids.sig_id eq sig_id


    ;Get maximum and minimum analysis data
    ;The AR start and end dates are not ideal for searching so expanding search range 2019/03/04 J. Prchlik
    min_date = anytim(AR_START[i])-3600.*24.*7
    max_date = anytim(AR_END[i])+3600.*24.*7


    ;Comment out if LOCKHEED MARTIN is down for the day
    outvals = get_sigmoid_flares(min_date,max_date,cnt_date,xvals,yvals,cnt_x,cnt_y,noaa[i],$
                                 ffl_x,ffl_y,ffl_ts,ffl_te,ffl_tp,ffl_mx,ffl_cl)


    ;Create dummy arrays to fill in structure
    cross_m = '20'+str_replace(anytim(TBEST[i],/yymmdd),', ','T');merdian crossing time
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


    ;Create fixed width arrays containing flare info
    for m=0,n_elements(ffl_x)-1 do begin
        flare_x[m] = ffl_x[m]
        flare_y[m] = ffl_y[m]
        flare_s[m] = ffl_ts[m]
        flare_e[m] = ffl_te[m]
        flare_p[m] = ffl_tp[m]
        flare_c[m] = ffl_cl[m]
    endfor


     ;Create single row in structure
     tmp = {test_sig_b,$
         sigmoid_id:fix(real_sig_id[i]),$ ;Use index in sav file to call real ID in csv file input to make save file
         sigmd_s:SIG_START[i],$ ;Use index in sav file to call sigmiod start time in csv file input to make save file
         sigmd_e:SIG_END[i],$ ;Use index in sav file to call sigmoid end time in csv file input to make save file
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
outf = '("sigmoid_id_",I03,"_",I03,".sav")'
save,big_str,filename=string([fix(min(real_sig_id)),fix(max(real_sig_id))],format=outf)


end
