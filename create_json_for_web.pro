;###########################################################################
;PROGRAM
;
;CALL
;    sig_id = get_iau_format(tobs,x,y,sig_start,lat=lat,lon=lon)
;
;INPUT
;    tobs - Time of observation for X, Y coordinates
;    X    - coord in HPC arcsec
;    Y    - coord in HPC arcsec
;    sig_start - start time to use a reference for making the IAU coordinate
;    lat  - Carrington latitude (output)
;    lon  - Carrington longitude (output)
;
;
;OUTPUT
;    sig_id - IAU formated Sigmoid ID number  
;
;###########################################################################
function get_iau_format,tobs,x,y,sig_start,lat=lat,lon=lon

;Format input variables
gi = tobs
xi = x
yi = y
t1 = sig_start

;Rotate to start time
rot_p = rot_xy(xi,yi,tstart=gi,tend=t1)

;Create directory for output png files
;Subscribe to IAU standard on output format
iau_time = strsplit(gi,'.',/extract)
iau_time = iau_time[0]

;IAU_format for coordinates
iau_cor = '("L",I03,"C",I03)'

;Get distance to sun
d_sun = get_sun(t1)

WCS_CONV_HPC_HG, rot_p[0], rot_p[1], lon, lat, dsun_obs=d_sun[0], length_units='AU',/carr, /pos_long
iau_pos = string([round(lon),round(lat)],format=iau_cor)

;Create new unique ID with IAU standard
sig_id = 'SOL'+iau_time+iau_pos


return,sig_id
end

;#############################################################################
;
;
;INPUTS
;    times      -   A csv file containing times to analyze sigmoid filaments
;                   CSV format must be as follows:
;                   formats = 'LL,LL,A,A,A,A,F,F,A,A,F,A,A,A,A,A,F,F,f,F,F'
;                   readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
;                   length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats
;
;
;OUTPUTS
;    A JSON formatted text file for use with the XRT sigmiod webpage
;
;#############################################################################
pro create_json_for_web,times,hmi_arch=hmi_arch,aia_arch=aia_arch,flr_arch=flr_arch,out_arch=out_arch

;Updates with Patty's new output format 2018/06/13 J. Prchlik
formats = 'LL,LL,A,A,A,F,F,A,A,A,F,A,A,A,A,A,F,F,f,F,F'
;formats = 'A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A'
readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
       length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats,/preserve_null

;Set HMI movie directory 
if keyword_set(hmi_arch) then hmi_arch = hmi_arch else hmi_arch = 'hmi_movie/'
hmi_arch = hmi_arch+'/'

;Set SDO/AIA movie directory 
if keyword_set(aia_arch) then aia_arch = aia_arch else aia_arch = 'aia_movie/'
aia_arch = aia_arch+'/'

;Set flare movie directory 
if keyword_set(flr_arch) then flr_arch = flr_arch else flr_arch = 'aia_arch_cutout/'
flr_arch = flr_arch+'/'

;Set output directory for symbolic links
if not keyword_set(out_arch) then out_arch = 'combined_movies/'
out_arch = out_arch+'/'


;good sigmoid tbest times (i.e. contains time string)
goodt = where(strlen(tobs) eq 23)

;format for output symbolic link directory
out_fmt = '(A30,"/")'


;Select button at the start of each row
select_button = '<img class=\"details_icon\" style=\"cursor:pointer;\" title=\"Click for details\" alt=\"Open Details\" src=\"./scripts/DataTables/examples/examples_support/details_open.png\"></img>",'

;3 indents for the start of each column in a raw
ind_3 = '            "'
select_button = ind_3+select_button



;Table containing flare string
flare_str_fmt = '("B: ",I4," C: ",I4," M: ",I4," X: ",I4)'

;File to write to
openw,33,'sigmoid_webpage/json.txt'
printf,33,'{'
printf,33,'    "aaData": ['

max_it = n_elements(goodt)-1 ;maximum iterator
;max_it = 2 ;maximum iterator
;Loop over good indices
;for ii=0,n_elements(goodt)-1 do begin
;cut loop short for testing
for ii=0,max_it do begin

     ;start bracket for storing data in a "row"
     printf,33,'        ['

    ;Set index to value with goodt
    i = goodt[ii]
    ;get index for a good time
    ;Now x, y are measured at tobs instead of tbest 2018/06/14 J. Prchlik
    gi = tobs[i]
    xi = x[i]
    yi = y[i]

    ;if there is no new measurements of the sigmoid just continue
    if gi eq '0' then continue




   
   
 

    ;Get sigmoid start and end times
    t1 = sig_start[i]
    t2 = sig_end[i]

    sig_id = get_iau_format(gi,xi,yi,t1,lat=lati,lon=loni)

    ;Directory to put the symoblic link
    full_dir = out_arch+string([sig_id],format=out_fmt)
    ;Remove : characters
    full_dir = str_replace(full_dir,':','')
    full_dir = str_replace(full_dir,'-','')
    ;make output directory
    hmi_dir = full_dir+'/hmi_movie/'
    aia_dir = full_dir+'/aia_movie/'
    flr_dir = full_dir+'/flr_movie/'

    ;###############################################################################
    ;
    ;HMI Movies
    ;
    ;###############################################################################
    ;get mp4 from hmi directory
    hmi_mov = hmi_dir+sig_id+'_hmi.mp4'

    ;hmi flux plot
    ;###############################################################################
    ;
    ;HMI Flux
    ;
    ;###############################################################################
    hmi_png = file_search(hmi_dir+'*.png',/full)

   




    ;###############################################################################
    ;
    ;SDO/AIA Evolution Movies
    ;
    ;###############################################################################
    ;Create text for link related to sigmoid long mp4 file
    aia_mov = aia_dir+sig_id+'.mp4'

     

    ;###############################################################################
    ;
    ;Flare Movies
    ;
    ;###############################################################################
 

    ;Create directory for output png files
    ;full_flr = flr_arch+strcompress(ID[i],/remove_all)+'/'

    ;get all flare movies
    flare_files = file_search(flr_dir+'*mp4',/full,count=flare_cnt)

    ;flare link
    ;if n_elements(size(flare_files)) lt 4 then continue

    ;Get counters for flares
    b_cnt = 0
    c_cnt = 0
    m_cnt = 0
    x_cnt = 0
    
    
    

    ;initial flare text
    mat_flr   = '<b>Flares from this region: </b></A>'
    ;Create text for flare links
    if flare_cnt gt 0 then begin
        ;get array of times and classes for sending to create_flux_plot 2018/09/25 J. Prchlik
        flare_times = strarr(n_elements(flare_files))
        flare_class = strarr(n_elements(flare_files))
  
        for j=0,n_elements(flare_files)-1 do begin
           ;File name
           file_name = strsplit(flare_files[j],'/',/extract)
           file_name = file_name[n_elements(file_name)-1]
           year  = strmid(file_name,34,4)
           month = strmid(file_name,38,2)
           day   = strmid(file_name,40,2)
           hour  = strmid(file_name,43,2)
           min   = strmid(file_name,45,2)
           class = strmid(file_name,50,2)
           sclass= strmid(file_name,53,1)

           ;Create line with flare class and time
           add_flare = '<br><A HREF=\"'+flare_files[j]+'\">'+class+'.'+sclass+' '+month+'/'+day+' '+hour+':'+min+'</A>'

           ;add to matched flare text
           mat_flr = mat_flr+add_flare

           ;Add flare times and class to variables
           flare_times[j] = year+'/'+month+'/'+day+' '+hour+':'+min+':00'
           flare_class[j] = class+'.'+sclass
  
         ;Add one to counter if a flare class is found
         case strmid(class,0,1) of 
             'B': b_cnt += 1
             'C': c_cnt += 1
             'M': m_cnt += 1
             'X': x_cnt += 1
         endcase

        endfor
     endif
    ;add table break
    mat_flr = mat_flr+'</div>'


    ;Flare string for table
    flare_str = string([b_cnt,c_cnt,m_cnt,x_cnt],format=flare_str_fmt)





    ;readcol,times,dum,ID,NOAA,AR,AR_START,X,Y,AR_END,SIG_START,SIG_END,lifetime,TBEST,tobs,ORIENTATION,HEMISPHERE, $
    ;       length_171,length_304,length,trail_length,lead_length,aspect_ratio,fwhm,height,format=formats
    ;Text to include in image information
    ar_sstr  = '<div align=\"left\" class=\"block\"><b>AR Start: </b>'+str_replace(strmid(ar_start[i],0,16),'T',' ')
    ar_estr   = '<br><b>AR End: </b>'+str_replace(strmid(ar_end[i],0,16),'T',' ')
    sig_sstr = '<br><b>Sigmoid Start: </b>'+str_replace(strmid(sig_start[i],0,16),'T',' ')
    sig_estr   = '<br><b>Sigmoid End: </b>'+str_replace(strmid(sig_end[i],0,16),'T',' ')
    mag_mov   = '<br><b>Movies: </b><A HREF=\"../'+hmi_mov+'\">Magnetogram</A>'
    sdo_mov   = '    <A HREF=\"../'+aia_mov+'\">SDO/AIA and Hinode/XRT</A>'
    flr_cat   = '<br><A HREF=\"http://xrt.cfa.harvard.edu/flare_catalog/all_full.html?search='+string(NOAA[i],format='(I5)')+'\" target=\"_blank\">'

    ;Series of still images to add to the sigmoid catalog linking to solar monitor
    
    ;base link to webpage
    base_solar_m = '\"https://solarmonitor.org/data/'


    ;Commented out to test other things 2018/09/25 J. Prchlik
    ;get the solarmonitory links 
    ;get_solarmonitor_links_dirty,str_replace(strmid(TBEST[i],0,10),'-',''),solar_links

    ;Add solar monitor links in text format
    solmon_img = ''
    for k=0,n_elements(solar_links)-1 do begin
        add_img = '<div align=\"middle\" class=\"block\"><A HREF=\"'+solar_links[k]+'\" rel=\"lightbox\">SRC=\"'+solar_links[k]+'\" height=\"200\" width=\"200\"></A></div>'
        solmon_img = solmon_img+add_img
    endfor 

    ;<A HREF=\"./images/001/001_euv.png\" rel=\"lightbox\">SRC=\"./images/001/001_euv.png\" height=\"200\" width=\"200\"></A></div><div align=\"middle\" class=\"block\"><A HREF=\"./images/001/001_mag.png\" rel=\"lightbox\">SRC=\"./images/001/001_mag.png\" height=\"200\" width=\"200\"></A></div><div align=\"middle\" class=\"block\"><A HREF=\"./images/001/001_halpha.png\" rel=\"lightbox\">SRC=\"./images/001/001_halpha.png\" height=\"200\" width=\"200\"></A></div>"'


    ;Combine all infromation into one table
    image_info = ind_3+ar_sstr+ar_estr+sig_sstr+sig_estr+mag_mov+sdo_mov+flr_cat+mat_flr+solmon_img+'"'
;    image_info = ind_3+'<div align=\"left\" class=\"block\"><b>AR Start: </b>02/26 18:03<br><b>AR End: </b>03/07 06:14<br><b>Sigmoid Start: </b>02/26 11:17<br><b>Sigmoid End: </b>02/27 11:40<br><b>Movies: </b><A HREF=\"./images/001/001_mag.mov\">Magnetogram</A><br><A HREF=\"http://xrt.cfa.harvard.edu/flare_catalog/full.html?search=2007+10944\" target=\"_blank\"><b>Flares from this region: </b></A><br>B2.5 03/02 05:29</div><div align=\"middle\" class=\"block\"><A HREF=\"./images/001/001_sigmoid.png\" rel=\"lightbox\">SRC=\"./images/001/001_sigmoid.png\" height=\"200\" width=\"200\"></A></div><div align=\"middle\" class=\"block\"><A HREF=\"./images/001/001_euv.png\" rel=\"lightbox\">SRC=\"./images/001/001_euv.png\" height=\"200\" width=\"200\"></A></div><div align=\"middle\" class=\"block\"><A HREF=\"./images/001/001_mag.png\" rel=\"lightbox\">SRC=\"./images/001/001_mag.png\" height=\"200\" width=\"200\"></A></div><div align=\"middle\" class=\"block\"><A HREF=\"./images/001/001_halpha.png\" rel=\"lightbox\">SRC=\"./images/001/001_halpha.png\" height=\"200\" width=\"200\"></A></div>"'
    
    ;print sigmoid information 
    printf,33,select_button
    printf,33,ind_3+string(ID[i],format='(I03)')          +'",';Sigmoid ID
    printf,33,ind_3+'0'            +'",'         ;Rating
    printf,33,ind_3+string(NOAA[i],format='(I5)')        +'",' ;NOAA AR number
    printf,33,ind_3+strmid(TBEST[i],0,19)       +'",' ;Best observing time in XRT images
    printf,33,ind_3+string(X[i],format='(F6.1)')           +'",' ;X location at time T_obs
    printf,33,ind_3+string(Y[i],format='(F6.1)')           +'",' ;Y location at time T_obs
    printf,33,ind_3+string(length[i]      ,format='(F6.1)')+'",' ;length of the sigmoid in arcsec
    printf,33,ind_3+string(aspect_ratio[i],format='(F6.1)')+'",' ;aspect ratio of sigmoid
    printf,33,ind_3+orientation[i] +'",' ;Orientation of sigmoid
    printf,33,ind_3+hemisphere[i]  +'",' ;Hemisphere of the Sigmoid
    printf,33,ind_3+string(length_171[i],format='(F9.1)')  +'",' ;does the sigmoid have a EUV filament
    printf,33,ind_3+string(length_304[i],format='(F9.1)')  +'",' ;Does the sigmoid have an H alpha filament
    printf,33,ind_3+'5'            +'",' ;Number of sunspot in AR (get from HEK)
    printf,33,ind_3+'N'            +'",' ;Magnetic Flux emergence 
    printf,33,ind_3+'N'            +'",' ;Magnetic Flux cancellation 
    printf,33,ind_3+'<A HREF=\"'+hmi_png+'\" >Flux Plot</A>'+'",' ;Get flux values for sigmoid
    printf,33,ind_3+flare_str      +'",' ;flare table string
    printf,33,ind_3+'0'            +'",' ;Number of CMEs
    printf,33,ind_3+'-'            +'",' ;Filament eruption
    printf,33,ind_3+'-'            +'",' ;Transient CH
    printf,33,ind_3+'-'            +'",' ;Flare Ribbons
    printf,33,ind_3+'-'            +'",' ;Post-Flare Loops
    printf,33,ind_3+'-'            +'",'  ;Nearby CH
    printf,33,ind_3+'-'            +'",'  ;Nearby AR
    printf,33,image_info
    ;end bracket for storing data in a "row"
    if ii eq  max_it then printf,33,'        ]' $
    else printf,33,'        ],'


endfor


;close the file
printf,33,'    ]'
printf,33,'}'
free_lun,33
end
