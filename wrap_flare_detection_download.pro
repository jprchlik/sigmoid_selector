pro wrap_flare_detection_download,years=years


   ; years = ['2012','2013','2014','2015','2016','2017','2018']
    if not keyword_set(years) then years = ['2012','2013','2014','2015','2016','2017']

    for i=0,n_elements(years)-1 do begin
        year = years[i]
        ;Associate flares with sigmoids
;        flare_cme_sigcat,'/data/tardigrade/pbolan/2018/',year,odir='./'
        ;Download cutouts around flares
        get_aia_files_cutout,'sigmoid_id_'+year+'.sav','SigmoidCatalogAll.csv'
    endfor
end