begin
    data = asciiread("../input/ex4.txt", (/177,5/), "float")
    data@_FillValue = -999

    lat = fspan(15,55,41)
    lat@units ="degrees_north"
    
    lon = fspan(75,135,121)
    lon@units = "degrees_east"

    index_nomiss = ind(.not.ismissing(data(:,2)))
    ind_num=dimsizes(index_nomiss)
    ;print(ind_num)
    data_nomiss = new((/ind_num,3/), "float")
    data_nomiss(:,0)=data(index_nomiss,0)
    data_nomiss(:,1)=data(index_nomiss,1)
    data_nomiss(:,2)=data(index_nomiss,2)
;   Exponent quals 0.5
    dssetp("exp", 0.5)
    grid = dsgrid2(data_nomiss(:,0), data_nomiss(:,1), data_nomiss(:,2), lon, lat)

    grid!0 = "lon"  ;I think it is very important. To name and identify dimension and assign value to dimension.
    grid!1 = "lat"

    grid&lon = lon
    grid&lat = lat
    ;printVarSummary(grid) 
    ;print(grid)

    wks = gsn_open_wks("png", "dsgrid2 using NCL")
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

    res@gsnLeftString = "Geopotential Height"
    res@gsnRightString = "m"
    res@tiMainString   = "dsgrid2 use NCL"

    res@lbBoxEndCapStyle      = "TriangleBothEnds"
    res@pmLabelBarHeightF     = 0.06
    res@pmLabelBarOrthogonalPosF = 0.11
    if(all(ismissing(grid))) then
   	    print("Your data_nomiss is all missing. Cannot create plot.")
    else
    	map = gsn_csm_contour_map_ce(wks, grid(lat|:,lon|:), res)
    end if

    resMaker = True
    resMaker@gsMarkerColor = "red"
    resMaker@gsMarkerIndex = 16
    resMaker@gsMarkerSizeF = .004

    markers = gsn_add_polymarker(wks, map, data_nomiss(:,0), data_nomiss(:,1), resMaker)

    draw(map)
    frame(wks)
end


    

