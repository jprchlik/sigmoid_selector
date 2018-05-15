;####################################################
;
;Description
;    Create formatted file formatted links required by aia_mkmovie
;####################################################
pro create_formatted_links,aia_arch=aia_arch,date_arr=date_arr,time_arr=time_arr,wave_arr=wave_arr,out_arch=out_arch,ext=ext,outfmt=outfmt

;set default extension
if keyword_set(ext) then ext = ext else ext = '.fits'

;set default file output format
if keyword_set(outfmt) then outfmt = outfmt else outfmt = '("AIA",I8,"_",I06,"_",I04,".fits")'

;Default to aia_arch if aia_arch not set
if keyword_set(aia_arch) then aia_arch = aia_arch else aia_arch = 'aia_arch/'
aia_arch = aia_arch+'/'

;Default to out_arch if aia_arch/symlinks/ not set
if keyword_set(out_arch) then out_arch = out_arch else out_arch = aia_arch+'symlinks/'
out_arch = out_arch+'/'

;Default 18-25 for date
if keyword_set(date_arr) then date_arr = date_arr else date_arr = [17,25]

;Default 27-32 for time
if keyword_set(time_arr) then time_arr = time_arr else time_arr = [26,32]

;Default 35-37 for wavelength
if keyword_set(wave_arr) then wave_arr = wave_arr else wave_arr = [34,37]


;get fits files in aia_arch
fils = file_search(aia_arch+'*'+ext,/full)

;create new symbolic links formatted for aia_mkmovie
for i=0L,n_elements(fils)-1 do begin

    ;get the root file name
    fname = strsplit(fils[i],'/',/extract)
    fname = fname[n_elements(fname)-1]

    ;format file for output
    date_out = strmid(fname,date_arr[0],date_arr[1]-date_arr[0])
    time_out = strmid(fname,time_arr[0],time_arr[1]-time_arr[0])
    wave_out = strmid(fname,wave_arr[0],wave_arr[1]-wave_arr[0])

    ;output file name
    oname = out_arch+string(long([date_out,time_out,wave_out]),format=outfmt)
   
    ;only continue if filename is in the proper format
    if strlen(strcompress(oname,/remove_all)) ne strlen(oname) then continue

    ;if symbolic link already exits continue
    if file_test(oname,/symlink) then continue


    ;remove empty links
    if file_test(oname,/dangling_symlink) then file_delete,oname
    
    ;Create symbolic link
    file_link,fils[i],oname

end


end