;############################################################
;DESCRIPTION
;    Using the reference to an image on the home page of 
;    solarmonitor.org find the full disk png file path
;    stored one solarmonitor.org
;
;
;USAGE
;   link_path = get_png_file(ref_text,base)
;
;INPUT
;    ref_text -- Line in base solarmonitor page referencing 
;    the wavelength/observation you want as php query
;    
;    base -- base file directory (https:/solarmonitor.org)
;
;OUTPUT
;    link_path -- The link for a specific full sun image on 
;    solarmonitor
;
;############################################################
function get_png_file_link,ref_text,base

    print,ref_text
    link_arry = strsplit(ref_text,'"',/extract) 
    link_arry = link_arry[1] ;Get the element containing the php query 
    ;Query the webpage to get Full disk png file
    spawn,"wget -S -O - '"+base+'/'+link_arry+"'",full_page
    ;Get line with png file name
    full_file = full_page[where(STRMATCH(full_page,'*png*fulldisk*') EQ 1,full_cnt)]
    ;split array of full file
    full_splt = strsplit(full_file)
    ;get the location with the png file
    full_pngf = str_replace(full_splt[1],'src=','')

    link_path = base+'/'+full_pngf

    return,link_path
end

;############################################################
;DESCRIPTION
;    Returns links to full Sun images on Solarmonitor.org for
;    the following observations: LOS Magnetic Fields, X-ray,
;    H alpha, and SDO/AIA 171 or SWAP 174
;
;USAGE
;    get_solarmonitor_links_dirty,date,links
;
;INPUT
;    date -- Date formated in YYYYMMDD to get solarmonitor
;    full sun images
;
;OUTPUT
;    links -- An array of links to the solar monitor website
;
;
;
;############################################################
pro get_solarmonitor_links_dirty,date,links



    base = 'https://solarmonitor.org'
    spawn,"wget -S -O - '"+base+"?date="+date+"'",base_page

    ;Get magnetogram page
    mag_cut = base_page[where(STRMATCH(base_page,'*full_disk*MAG*indexnum=1*',/FOLD_CASE) EQ 1,mag_cnt)]
    ;Get XRT page
    xrt_cut = base_page[where(STRMATCH(base_page,'*full_disk*hxrt*indexnum=1*',/FOLD_CASE) EQ 1,xrt_cnt)]
    ;Get Halpha page
    hal_cut = base_page[where(STRMATCH(base_page,'*full_disk*halph*indexnum=1*',/FOLD_CASE) EQ 1,hal_cnt)]
    ;Get 171 Fe line page
    fel_cut = base_page[where(STRMATCH(base_page,'*full_disk*00171*indexnum=1*',/FOLD_CASE) EQ 1,fel_cnt)]
    ;if no 171 use SWAP 174
    if fel_cnt eq 0 then fel_cut = base_page[where(STRMATCH(base_page,'*full_disk*00174*indexnum=1*',/FOLD_CASE) EQ 1,fel_cnt)]


    links = []

    if mag_cnt gt 0 then begin $
       mag_link = get_png_file_link(mag_cut,base)
       links = [mag_link]
    endif
    if xrt_cnt gt 0 then begin $
       xrt_link = get_png_file_link(xrt_cut,base)
       links = [links,xrt_link]
    endif
    if hal_cnt gt 0 then begin $
       hal_link = get_png_file_link(hal_cut,base)
       links = [links,hal_link]
    endif
    if fel_cnt gt 0 then begin $
       fel_link = get_png_file_link(fel_cut,base)
       links = [links,fel_link]
    endif


end