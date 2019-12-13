        program interpolation
                implicit none
                real data(5,177),lon_sta(177),lat_sta(177),hgt_sta(177)
                real lon_grid(121), lat_grid(81),hgt(121,81)
                real w,fw,sum_w,sum_fw
                integer i,j,k 
                real google_distance

                open(10,file="../input/ex4.txt")
                read(10,*)((data(i,j), i=1, 5), j=1, 177)
                lon_sta = data(1,:)
                lat_sta = data(2,:)
                hgt_sta = data(3,:)
                close(10)
                do j=1, 121
                        do i=1, 81
                                lon_grid(j) = 75+(j-1)*0.5
                                lat_grid(i) = 15+(i-1)*0.5
                                sum_fw = 0
                                sum_w = 0
                                do k=1, 177
                                        if(hgt_sta(k).ne.-999) then
                w = 1./(0.001+(lon_sta(k)-lon_grid(j))**2+&
                    (lat_sta(k)-lat_grid(i))**2)
                                                fw = hgt_sta(k)*w 
                                                sum_fw = sum_fw + fw
                                                sum_w = sum_w + w
                                        end if
                                end do
                                if(sum_w.ne.0) then
                                        hgt(j,i) = sum_fw/sum_w
                                else
                                        hgt(j,i) = -999
                                end if
                        end do
                end do
                call write_nc(hgt)

        end program interpolation
        real function google_distance(lon1,lat1,lon2,lat2)
                implicit none
                real,parameter :: RE=6378.137
                real,parameter:: Radian=1.745329E-2
                real lon1,lat1,lon2,lat2
                real a1,b1,a2,b2,cs,a,b
!               radian=4*atan(1.0)/180.0
                a1=lon1*radian
                a2=lon2*radian
                b1=lat1*radian
                b2=lat2*radian
                a =b1-b2
                b =a1-a2
                cs =2*asin(sqrt(sin(a/2)**2+cos(b1)*cos(b2)*sin(b/2)**2))
                google_distance =(RE*cs)**2! square
                return 
        end function google_distance
        
        subroutine write_nc(hgt)
                implicit none
                include '/usr/include/netcdf.inc'
                character(len=256) ::newfile
                integer,parameter :: nlat = 81,nlon = 121
                integer ::statu,ncid,nlatid,nlonid,hgtid,latid,lonid,&
                          i,j
                integer ::hgtdims(2),latdims(1),lonims(1)
                real    ::hgt(nlon,nlat),lat(nlat),lon(nlon)
                do i=1,81
                       lat(i)=15+(i-1)*0.5
                end do
                do j=1,121
                       lon(j)=75+(j-1)*0.5
                end do 
                newfile = "interpolated_hgt.nc"
                write(*,*) "Creat New NETcdf", trim(newfile)
                statu=nf_create(trim(newfile),nf_noclobber,ncid)
                if(statu /= nf_noerr) call handle_err(statu)

                write(*,*) "Define Dimensions"
                statu=nf_def_dim(ncid,'lon',nlon,nlonid)
                if(statu /= nf_noerr) call handle_err(statu)
                statu=nf_def_dim(ncid,'lat',nlat,nlatid)
                if(statu /= nf_noerr) call handle_err(statu)

                write(*,*)"Define Variables"
                hgtdims(1)=nlonid
                hgtdims(2)=nlatid
                statu=nf_def_var(ncid,"hgt",nf_float,2,hgtdims,hgtid)
                if(statu /= nf_noerr) call handle_err(statu)
                statu=nf_def_var(ncid,"lon",nf_float,1,nlonid,lonid)
                if(statu /= nf_noerr) call handle_err(statu)
                statu=nf_def_var(ncid,"lat",nf_float,1,nlatid,latid)
                if(statu /= nf_noerr) call handle_err(statu)

                write(*,*)"Put Attributes"
                if(statu /= nf_noerr) call handle_err(statu)
                statu=nf_put_att_text(ncid,hgtid,"description",19,&
                "Geopotential Height")
                if(statu /= nf_noerr) call handle_err(statu)
                statu=nf_put_att_text(ncid,hgtid,"units",1,"m")
                if(statu /= nf_noerr) call handle_err(statu)
                statu=nf_put_att_text(ncid,hgtid,"coordinates",11,&
                "XLONG XLAT")

                statu=nf_put_att_text(ncid,latid,"units",12,"degree_north")
                if(statu /= nf_noerr) call handle_err(statu)
                statu=nf_put_att_text(ncid,latid,"description",27,&
                "LATITUDE, SOUTH IS NEGATIVE")
                if(statu /= nf_noerr) call handle_err(statu)
                statu=nf_put_att_text(ncid,latid,"axis",1,"X")
                if(statu /= nf_noerr) call handle_err(statu)

                statu=nf_put_att_text(ncid,lonid,"units",11,"degree_east")
                if(statu /= nf_noerr) call handle_err(statu)
                statu=nf_put_att_text(ncid,lonid,"description",27,&
                "LONGITUDE, WEST IS NEGATIVE")
                if(statu /= nf_noerr) call handle_err(statu)
                statu=nf_put_att_text(ncid,lonid,"axis",1,"y")
                if(statu /= nf_noerr) call handle_err(statu)


                write(*,*)"End Defination"
                statu=nf_enddef(ncid)
                if(statu /= nf_noerr) call handle_err(statu)

                write(*,*)"Put Variable"
                statu=nf_put_var_real(ncid,hgtid,hgt)
                if(statu /= nf_noerr) call handle_err(statu)
                statu=nf_put_var_real(ncid,lonid,lon)
                if(statu /= nf_noerr) call handle_err(statu)
                statu=nf_put_var_real(ncid,latid,lat)
                if(statu /= nf_noerr) call handle_err(statu)

                write(*,*)"Close nc"
                statu=nf_close(ncid)
                if(statu /= nf_noerr) call handle_err(statu)
        end

        subroutine handle_err(statu)
                implicit none
                integer statu
                write(*,*)"error!"
                stop
        end
