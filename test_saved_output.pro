pro test_saved_output
set_plot,'X'

restore,'examples2/sigmoid_sizedata.sav'
loadct,3

nfiles=n_elements(sigmoids)

for xx=0,nfiles-1 do begin
    mreadfits,'examples2/'+sigmoids[xx].filename,index1,data1,/verbose

    ;Get axis size for the image
    xwdw_size=index1[0].naxis1
    ywdw_size=index1[0].naxis2

    ;Scale image max and min cheaply
    useg = where(finite(alog10(data1)))
    scmin=cgPercentiles(alog10(data1[useg]),percentiles=.010)
    scmax=cgPercentiles(alog10(data1[useg]),percentiles=.999)

    ;Set up window
    window,5,xs=xwdw_size,ys=ywdw_size


    image = bytscl(alog10(data1),min=scmin,max=scmax)
    tv,image

     xyouts,sigmoids[xx].longx1,sigmoids[xx].longy1,'L',/device,alignment=0.5
     xyouts,sigmoids[xx].longx2,sigmoids[xx].longy2,'L',/device,alignment=0.5

     xyouts,sigmoids[xx].shrtx1a,sigmoids[xx].shrty1a,'X',/device,alignment=0.5
     xyouts,sigmoids[xx].shrtx2a,sigmoids[xx].shrty2a,'X',/device,alignment=0.5
     xyouts,sigmoids[xx].shrtx1b,sigmoids[xx].shrty1b,'X',/device,alignment=0.5
     xyouts,sigmoids[xx].shrtx2b,sigmoids[xx].shrty2b,'X',/device,alignment=0.5

     ;Changed to polygon around sigmoid
     plots,sigmoids[xx].bboxx,sigmoids[xx].bboxy,color=255,/device

     ;Create sigmoid region where to calculate the FWHM
     xyouts,sigmoids[xx].fwlin1[0],sigmoids[xx].fwlin1[1],'F',/device,alignment=0.5,color=0
     xyouts,sigmoids[xx].fwlin2[0],sigmoids[xx].fwlin2[1],'F',/device,alignment=0.5,color=0
     wait,1.0



endfor

end
