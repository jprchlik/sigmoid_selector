;====================================================
;USAGE
;    combine_sigmoid_structures,arch=arch,subdir=subdir, $
;                              ,syear=syearc,eyear=eyear
;
;DESCRIPTION 
;    Combines seperate sigmoid property save files with
;    keyword sigmoid into a larger structure containing
;    all sigmiod properties
;
;====================================================
pro combine_sigmoid_structures,arch=arch,subdir=subdir,syear=syear,eyear=eyear

;Set up variables if specified, but none needed
if keyword_set(arch) then arch=arch else arch = '/data/tardigrade/pbolan/CMS2/sigmoid_selector'
arch = arch+'/'
if keyword_set(subdir) then subdir = subdir else subdir = '("SigmoidCatalog",I4)'
if keyword_set(syear) then syear = syear else syear = 2007
if keyword_set(eyear) then eyear = eyear else eyear = 2017


;loop over all years
for i=syear,eyear do begin

    ;Restore sigmoid save file
    restore,arch+string([i],format=subdir)+'/sigmoid_sizedata.sav'

    ;if it is the first year init the large object
    if i eq syear then all_sig = sigmoids $ 
    else all_sig = [all_sig,sigmoids]

endfor

save,all_sig,filename=string([syear,eyear],format='("sigmoid_sizedata_",I4,"_to_",I4,".sav")')

end