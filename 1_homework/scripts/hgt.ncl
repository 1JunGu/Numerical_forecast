begin
	url = "http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/"
	filename = "hgt.mon.mean.nc"
	f = addfile(url+filename,"r")
	
	time = cd_calendar(f->time,-1)
	time_idx = ind(time .eq. 199807)
	hgt = f->hgt(time_idx, {lev|850}, {lat|5:60}, {lon|60:140})
	printVarSummary(hgt)


end
