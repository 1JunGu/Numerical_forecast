begin
    a = addfile("interpolated_hgt.nc","r")
    hgt = a->hgt
    lon = a->lon
    lat = a->lat
    hgt@_FillValue = -999

    wks = gsn_open_wks("png", "interpolation_hgt using Fortran")
    gsn_define_colormap(wks, "WhiteBlueGreenYellowRed")

    res = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@cnFillOn = True
    res@cnLinesOn = False
    res@gsnAddCyclic = False

    res@gsnSpreadColors = True

    res@mpLimitMode                 = "LatLon"
    res@mpMinLatF                 = 15.
    res@mpMinLonF                 = 75.
    res@mpMaxLatF                 = 55.
    res@mpMaxLonF                 = 135.
    res@mpOutlineBoundarySets = "National"
    res@mpDataBaseVersion     = "MediumRes"
    res@mpDataSetName         = "Earth..4"
    res@mpOutlineSpecifiers   = (/"China:States","Taiwan"/)

    res@lbBoxEndCapStyle      = "TriangleBothEnds"
    res@pmLabelBarHeightF     = 0.06
    res@pmLabelBarOrthogonalPosF = 0.11
    res@tiMainString = "interpolation use Fortran"
    plot = gsn_csm_contour_map_ce(wks,hgt,res)

    draw(plot)
    frame(wks)
end


    

