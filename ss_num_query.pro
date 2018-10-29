pro ss_num_query,times



    ;Updates with Patty's new output format 2018/06/13 J. Prchlik
    formats = 'LL,LL,A,A,A,F,F,A,A,A,F,A,A,A,A,A,F,F,f,F,F'
    ;formats = 'A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A'
    readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
           length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats,/preserve_null


    ;string print format
    str_fmt = '(I4,3x,I6,3x,I4)'

    ;Loop over all sigmiods
    for i=0,n_elements(dum)-1 do begin

        ;if ID[i] lt 30 then continue

        ;get date of best sigmiod 
        obs_tim_s = strmid(TBEST[i],0,11)+'00:00:00'
        obs_tim_e = strmid(TBEST[i],0,11)+'23:59:59'
        ;Get AR number as a string for sigmiod
        if ar[i] lt 1 then continue ;skip sigmoids without an AR
        ar_val = ar[i]
        ar_string = strcompress(ar_val,/remove_all)

        query=ssw_her_make_query(obs_tim_s,obs_tim_e,/AR,search_array=['FRM_NAME=NOAA SWPC Observer','AR_NOAANUM='+ar_string])
        her=ssw_her_query(query)

        ;if nothing found just continue
        if typename(her) eq 'INT' then continue

        print,string([ID[i],AR[i],her.ar.AR_NUMSPOTS],format=str_fmt)
    endfor

end