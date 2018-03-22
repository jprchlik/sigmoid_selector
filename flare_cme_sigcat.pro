pro flare_cme_sigcat, date1, date2, arnum

; Give program active region lifetime start and end times to get out flares and associated cmes

; Call the flare program 
date1 = date1
date2 = date2
; example dates:
; '2012/05/03 00:30', '2012/05/14 16:30'

; database line fore when running program outside of desktop

database = '/Users/smckillo/Desktop/Sig_Cat/xra_database.sav'

events = noaa_xra_query(date1, date2, database = database, /quiet) 

arnum = strmid(arnum, 1, 4)

sel = where(events.reg eq arnum, count)

if count eq 0 then print, 'No Flare Events Found For This Region Within the Time Specified'

if count ne 0 then begin
   selsize = size(sel, /n_elements)
   fl_st = strarr(selsize)
   fl_en = strarr(selsize)
   fl_pk = strarr(selsize)
   fl_cl = strarr(selsize)
   for i=0, selsize-1 do begin
      fl_st(i) = events(sel(i))._begin
      fl_en(i) = events(sel(i))._end
      fl_pk(i) = events(sel(i)).max
      fl_cl(i) = events(sel(i)).class

   endfor


endif

; Now need to take that information and look for the CME links




stop
end
