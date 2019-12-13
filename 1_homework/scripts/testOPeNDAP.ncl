;---Path to file served by OPeNDAP server at NOAA
begin
  url      = "http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/pressure/"
  filename = "hgt.mon.mean.nc"

  exists = isfilepresent(url+filename)
  if(.not.exists) then
    print("OPeNDAP test unsuccessful.")
    print("Either the file doesn't exist, or NCL does")
    print("not have OPeNDAP cabilities on this system.")
  else
    f = addfile(url + filename,"r")
    printVarSummary(f)
  end if
end


