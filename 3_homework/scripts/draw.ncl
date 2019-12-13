begin
	;namelists = systemfunc("ls -tr /lustre/data/jungu/Senior_autumn/Numerical_forecast/3_homework/output/time_step=180/*.nc")
	;namelists = systemfunc("ls -tr /lustre/data/jungu/Senior_autumn/Numerical_forecast/3_homework/output/time_step=180/*.nc")
	namelists = systemfunc("ls -tr /lustre/data/jungu/Senior_autumn/Numerica_forecast/3_homework/output/time_step=120/*.nc")
	number_nc = dimsizes(namelists)
	print(namelists)
	system("mkdir ../pics/time_step=120")
	do i=0,number_nc
		;name="/lustre/data/jungu/Senior_autumn/homework_zrj/3_homework/output/time_step=200/zuv_0.nc"
		name = namelists(i)
		;print(name)
		f = addfile(name,"r");;;
		U = f->u 
		V = f->v 
		hgt = f->z 
		str_array=str_split(name,"/")
		;print(str_array)
		;system("mkdir ../pics/"+str_array(7))
		str_array2=str_split(str_array(8),".")
		nc_name = str_array2(0)
		;i=1	
		;name_str ="zuv_" + sprinti("%0.4i",i)
		;print(name_str)
		str_array3 = str_split(nc_name,"_")
		sn_str  = str_array3(1)
		sn_int  = stringtoint(sn_str)
		;print(sn_int)
		name_str="zuv_" + sprinti("%0.3i",sn_int)
		print(name_str)
		;system("rm -r "+name_str+".png")
		wks =gsn_open_wks("png","../pics/"+str_array(7)+"/"+name_str);;;
;;creat plots
	res		      = True
	res@gsnFrame          = False
	res@gsnDraw	      = False
   	res@gsnMaximize       = True
	res@tmXTOn            = False
    	res@tmYROn            = False
	res@tmXBLabelFontHeightF = 0.008
	res@tmYLLabelFontHeightF = 0.008
	res@vpHeightF		  = 0.5
	res@vpWidthF          = 0.4
    	res@gsnLeftString     = ""
    	res@gsnRightString    = ""
	;res@tiMainString      = "geostrophic winds and mean geopotential height sn="+sn_str
	res@tiMainString      = sn_str
	res@tiMainFontHeightF = 0.015
	;lres@gsnSpreadColors   = True
	
;;set map
	mpres		      		= res
    	;mpres@mpDataSetName         	= "Earth..4"
    	;mpres@mpDataBaseVersion     	= "MediumRes"
    	;mpres@mpOutlineOn       	= True          
    	;mpres@mpOutlineSpecifiers   	= (/"China:States","Taiwan"/)
    	;mpres@mpOutlineBoundarySets 	= "National"
    	;mpres@mpFillOn				= True
	;mpres@mpFillDrawOrder	   	= "PostDraw"
	;set area
        mpres@mpLimitMode       	= "LatLon"
        mpres@mpMinLatF         	= 10.
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
	;cnres@gsnLeftString		= "monyhly mean geopotential height /m"
	cnres@gsnLeftStringFontHeightF  = 0.015 
    	cnres@gsnAddCyclic      	= False
	
;;set vector
	vcres				= res
	vcres@vcGlyphStyle         	= "CurlyVector"    
	vcres@vcRefMagnitudeF      	= 10.0             
	vcres@vcRefLengthF         	= 0.03             
	vcres@vcLineArrowColor     	= "Blue"          
	vcres@vcMinDistanceF       	= 0.03              
	vcres@vcVectorDrawOrder    	= "PostDraw" 
	;vcres@gsnRightString         	= "geostrophic wind"	
	vcres@gsnRightStringFontHeightF = 0.015 
    	vcres@gsnAddCyclic      	= False

	vcres@vcRefAnnoOn               = True
   	vcres@vcRefAnnoString1          = "10"
   	vcres@vcRefAnnoSide             = "Top"
   	vcres@vcRefAnnoString2On        = False
   	vcres@vcRefAnnoPerimOn          = False
   	vcres@vcRefAnnoOrthogonalPosF   = -0.14
   	vcres@vcRefAnnoParallelPosF     = 1
   	vcres@vcRefAnnoBackgroundColor  = "white"
	
;;plot
	map				= gsn_csm_map(wks,mpres)
	contour				= gsn_csm_contour(wks,hgt,cnres)
	vector				= gsn_csm_vector(wks,U,V,vcres)

;;overlay
	overlay(map,contour)
	overlay(map,vector)

	draw(map)
	frame(wks)
	end do
end
