;---Path to file served by OPeNDAP server at NOAA
begin
  	url      = "http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/"
  	filename = "hgt.mon.mean.nc"

  	exists = isfilepresent(url+filename)
  	if(.not.exists) then
    		print("OPeNDAP test unsuccessful.")
    		print("Either the file doesn't exist, or NCL doesnot have OPeNDAP cabilities on this system.")
  	else
    		f = addfile(url + filename,"r")
    		;vnames = getfilevarnames(f)
    		;print(vnames)
  	end if
	time = cd_calendar(f->time,-1);convers a mixed Julian/Gregorian date to a UT-referenced date
	time_ind = ind(time .eq. 199807);seek the time index of 199807

	hgt = f->hgt(time_ind, {850}, {5:60}, {60:140});use normal subscripting and coordinate subscripting to find
						       ;level=850hPa,lat from 5N to 60N,lon from 60E to 140E
	hgt = hgt * hgt@scale_factor + hgt@add_offset;calculate true value of hgt

	wks = gsn_open_wks("png","hgt_199807_850hPa")
	
	cnres			= True
	cnres@gsnMaximize 	= True
	cnres@gsnAddCyclic	= False

    	cnres@mpOutlineOn       = True        	
	cnres@mpOutlineBoundarySets = "National"
    	cnres@mpDataBaseVersion     = "MediumRes"
    	cnres@mpDataSetName         = "Earth..4"
    	cnres@mpOutlineSpecifiers   = "China:States"
    	
        ;cnres@mpGridAndLimbOn       = True
    	;cnres@mpGridLineDashPattern = 2
	
	cnres@mpLimitMode	= "LatLon"
	cnres@mpMinLatF		= 5.
	cnres@mpMinLonF		= 60.
	cnres@mpMaxLatF		= 60.
	cnres@mpMaxLonF		= 140.
	
	plot = gsn_csm_contour_map_ce(wks,hgt,cnres)
	
	system("rm -rf hgt.850hPa.199807.nc")
	fc = addfile("hgt.850hPa.199807.nc","c")
	fc->hgt = hgt 
	delete(fc)
	
	txt = True
	txt@fout = "hgt.850hPa.199807.txt"
	write_matrix(hgt,"33f6.0",txt)

end
