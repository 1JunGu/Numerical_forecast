begin
	f =addfile("../output/hgt.850hPa.199807.nc","r")
	U = f->Ug
	V = f->Vg
	hgt = f->hgt

	wks =gsn_open_wks("png","../pics/geostrophic_winds")
;;creat plots
	res		      		  = True
	res@gsnFrame          = False
	res@gsnDraw	          = False
    res@gsnMaximize       = True
	res@tmXTOn            = False
    res@tmYROn            = False
	res@tmXBLabelFontHeightF = 0.008
	res@tmYLLabelFontHeightF = 0.008
	res@vpHeightF		  = 0.5
	res@vpWidthF          = 0.4
    res@gsnLeftString     = ""
    res@gsnRightString    = ""
	;res@tiMainString      = "geostrophic winds and mean geopotential height"
	res@tiMainString      = "ex2"
	res@tiMainFontHeightF = 0.015
	;lres@gsnSpreadColors   = True
	
;;set map
	mpres		      		= res
 ;       mpres@mpDataSetName         	= "Earth..4"
 ;       mpres@mpDataBaseVersion     	= "MediumRes"
 ;       mpres@mpOutlineOn       	= True          
 ;       mpres@mpOutlineSpecifiers   	= (/"China:States","Taiwan"/)
 ;       mpres@mpOutlineBoundarySets 	= "National"
	;	mpres@mpFillOn				= True
	;	;mpres@mpFillDrawOrder	   	= "PostDraw"
	;set area
        mpres@mpLimitMode       	= "LatLon"
        mpres@mpMinLatF         	= 5.
        mpres@mpMinLonF         	= 60.
        mpres@mpMaxLatF         	= 60.
        mpres@mpMaxLonF         	= 140.

;;set contour
	cnres				= res
	;cnres@cnFillDrawOrder		= "PreDraw"
	;cnres@cnFillOn 			= True
	;cnres@cnLinesOn 		= False
	;;res@gsnScalarContour     = True 
	;cnres@cnLevelSelectionMode 	= "ManualLevels"
	;cnres@cnLevelSpacingF	  	= 5.
	;cnres@cnMinLevelValF	  	= 1440.
	;cnres@cnMaxLevelValF	  	= 1540.
	;cnres@pmLabelBarDisplayMode		= "Always"
	;cnres@lbOrientation         	= "Vertical"
	;cnres@lbBoxEndCapStyle      	= "TriangleBothEnds"
	;cnres@pmLabelBarWidthF      	= 0.03
	;cnres@pmLabelBarHeightF			= 0.4
	;cnres@lbLabelFontHeightF          = 0.005
	;cnres@pmLabelBarOrthogonalPosF 	= -0.04
	;cnres@pmLabelBarParallelPosF	= 0.2
	;;subset a color map
	;cmap				= read_colormap_file("BlRe")
	;cnres@cnFillPalette		= cmap
	cnres@gsnLeftString		= "monyhly mean geopotential height /m"
	cnres@gsnLeftStringFontHeightF  = 0.015 
        cnres@gsnAddCyclic      	= False
	
;;set vector
	vcres				= res
	vcres@vcGlyphStyle         	= "CurlyVector"    
	vcres@vcRefMagnitudeF      	= 10.0             
	vcres@vcRefLengthF         	= 0.04             
	vcres@vcLineArrowColor     	= "black"          
	vcres@vcMinDistanceF       	= 0.02              
	vcres@vcVectorDrawOrder    	= "PostDraw" 
	vcres@gsnRightString         	= "geostrophic wind"	
	vcres@gsnRightStringFontHeightF = 0.015 
    vcres@gsnAddCyclic      	= False

	vcres@vcRefAnnoOn               = True
    vcres@vcRefAnnoString1          = "10"
    vcres@vcRefAnnoSide             = "Top"
    vcres@vcRefAnnoString2On        = False
    vcres@vcRefAnnoPerimOn          = False
    vcres@vcRefAnnoOrthogonalPosF   = -0.12
    vcres@vcRefAnnoParallelPosF     = 1
    vcres@vcRefAnnoBackgroundColor  = "Red"
	
;;plot
	map				= gsn_csm_map(wks,mpres)
	contour				= gsn_csm_contour(wks,hgt,cnres)
	vector				= gsn_csm_vector(wks,U,V,vcres)

;;overlay
	overlay(map,contour)
	overlay(map,vector)
	;print(U)

	draw(map)
	frame(wks)
end
