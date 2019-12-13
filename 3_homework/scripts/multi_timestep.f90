        program time_integration
         implicit none
         real z0(65,41), u0(65,41), v0(65,41)
         real z1(65,41), u1(65,41), v1(65,41)
         real z2(65,41), u2(65,41), v2(65,41)
         integer i, j, t, dt, sn
         real, parameter ::g=9.8, omega=7.29e-5, H=8000, Re=6371000,&
                           d2r=3.14159/180,c=300!wave velocity
         real :: lat, f, dx, dy,p
         real time_step
         character(len=256) :: serial_number, cmd
       
         !cmd="rm -r ./ncfiles/zuv_*"
         !call system(cmd)
         
         dy = Re*1.25*d2r
         write(*,*)"resolution :",dy,"kilos"
        
         ! step 1 determine delatT
         dy = Re*1.25*d2r
                !(c*dt)/dx .le. 1
         dx = dy*cos(60*d2r)
         dt=dx/c
         write(*,*)"max timestep is ",dt
         dt = 120               !deltaT = 3mins
         sn = 24*3600/dt        !determine calculate steps
         !print*, sn
        
         ! step 2:read z0,u0,v0
         open(10,file="../input/ex3.txt")
         read(10,*) ((z0(i,j), i=1, 65), j=1, 41)
         read(10,*) ((u0(i,j), i=1, 65), j=1, 41)
         read(10,*) ((v0(i,j), i=1, 65), j=1, 41)
         close(10)
         z1=z0; z2=z0;
         u1=u0; u2=u0;
         v1=v0; v2=v0;
         
         call write_nc(z0,u0,v0,dt,0)
         do j=2,40
                lat = (10+(j-1)*1.25)*d2r
                dx = dy*cos(lat)
                f = 2*omega*sin(lat)
                do i=2, 64
           u2(i,j)=u0(i,j)+dt*(-u0(i,j)*(u0(i+1,j)-u0(i-1,j))/(2*dx)&
                   -v0(i,j)*(u0(i,j+1)-u0(i,j-1))/(2*dy)&
                   -g*(z0(i+1,j)-z0(i-1,j))/(2*dx)+f*v0(i,j))
           v2(i,j)=v0(i,j)+dt*(-u0(i,j)*(v0(i+1,j)-v0(i-1,j))/(2*dx)&
                   -v0(i,j)*(v0(i,j+1)-v0(i,j-1))/(2*dy)&
                   -g*(z0(i,j+1)-z0(i,j-1))/(2*dy)-f*u0(i,j))
           z2(i,j)=z0(i,j)+dt*(-u0(i,j)*(z0(i+1,j)-z0(i-1,j))/(2*dx)&
                   -v0(i,j)*(z0(i,j+1)-z0(i,j-1))/(2*dy)&
                   -H*((u0(i+1,j)-u0(i-1,j))/(2*dx)&
                   +(v0(i,j+1)-v0(i,j-1))/(2*dy)))
                end do
         end do
         
         do j=2,40
                lat = (10+(j-1)*1.25)*d2r
                dx = dy*cos(lat)
                f = 2*omega*sin(lat)
                do i=2, 64
          u1(i,j)=u0(i,j)+dt*(-u2(i,j)*(u2(i+1,j)-u2(i-1,j))/(2*dx)&
                  -v2(i,j)*(u2(i,j+1)-u2(i,j-1))/(2*dy)&
                  -g*(z2(i+1,j)-z2(i-1,j))/(2*dx)+f*v2(i,j))
          v1(i,j)=v0(i,j)+dt*(-u2(i,j)*(v2(i+1,j)-v2(i-1,j))/(2*dx)&
                  -v2(i,j)*(v2(i,j+1)-v2(i,j-1))/(2*dy)&
                  -g*(z2(i,j+1)-z2(i,j-1))/(2*dy)-f*u2(i,j))
          z1(i,j)=z0(i,j)+dt*(-u2(i,j)*(z2(i+1,j)-z2(i-1,j))/(2*dx)&
                  -v2(i,j)*(z2(i,j+1)-z2(i,j-1))/(2*dy)&
                  -H*((u2(i+1,j)-u2(i-1,j))/(2*dx)&
                  +(v2(i,j+1)-v2(i,j-1))/(2*dy)))
                end do
        end do
        call write_nc(z1,u1,v1,dt,1)
        do t=2,sn 
           do j=2,40
                lat = (10+(j-1)*1.25)*d2r
                dx = dy*cos(lat)
                f = 2*omega*sin(lat)
                do i=2, 64
          u2(i,j)=u0(i,j)+dt*(-u1(i,j)*(u1(i+1,j)-u1(i-1,j))/(1*dx)&
                  -v1(i,j)*(u1(i,j+1)-u1(i,j-1))/(1*dy)&
                  -g*(z1(i+1,j)-z1(i-1,j))/(1*dx)+f*v1(i,j))
          v2(i,j)=v0(i,j)+dt*(-u1(i,j)*(v1(i+1,j)-v1(i-1,j))/(1*dx)&
                  -v1(i,j)*(v1(i,j+1)-v1(i,j-1))/(1*dy)&
                  -g*(z1(i,j+1)-z1(i,j-1))/(1*dy)-f*u1(i,j))
          z2(i,j)=z0(i,j)+dt*(-u1(i,j)*(z1(i+1,j)-z1(i-1,j))/(1*dx)&
                  -v1(i,j)*(z1(i,j+1)-z1(i,j-1))/(1*dy)&
                  -H*((u1(i+1,j)-u1(i-1,j))/(1*dx)&
                  +(v1(i,j+1)-v1(i,j-1))/(1*dy)))
                end do
           end do
          
           !write(serial_number,'(i0)')t 
           !write(*,*),"run "//AdjustL(trim(serial_number))//"circle"
           call smooth5(z2,0.5)
           call smooth5(u2,0.5)
           call smooth5(u2,0.5)
           
           call system("sleep 0.1")     
           call write_nc(z2,u2,v2,dt,t)
           z0=z1; u0=u1; v0=v1
           z1=z2; u1=u2; v1=v2
        end do
        end program


        subroutine write_nc(z,u,v,time_step,sn)
         implicit none
         include '/usr/include/netcdf.inc'
         character(len=256) ::newfile,serial_number,time_step_str,&
                              output_dir
         integer,parameter :: nlat = 41,nlon = 65
         integer ::statu,ncid,nlatid,nlonid,zid,uid,vid,latid,lonid,&
                   i,j
         integer ::time_step,sn
         integer ::zdims(2),udims(2),vdims(2),latdims(1),lonims(1)
         real    ::z1(nlat,nlon),u1(nlat,nlon),v1(nlat,nlon)
         real    ::z(nlon,nlat),u(nlon,nlat),v(nlon,nlat),&
                   lat(nlat),lon(nlon)
         do i=1,41
                lat(i)=10+(i-1)*1.25
         end do
         do j=1,65
                lon(j)=60+(j-1)*1.25
         end do 
         
         write(time_step_str,'(i0)')time_step
         output_dir="../output/"//"time_step="&
                    //AdjustL(trim(time_step_str))//"/"
         call system("mkdir "//AdjustL(trim(output_dir)))
         
         write(serial_number,'(i0)')sn
         newfile = AdjustL(trim(output_dir))//"zuv_"//&
                   AdjustL(trim(serial_number))//".nc"
         write(*,*),trim(newfile)
         write(*,*) "Creat New NETcdf", trim(newfile)
         statu=nf_create(trim(newfile),nf_noclobber,ncid)
         if(statu /= nf_noerr) call handle_err(statu)
         
         write(*,*) "Define Dimensions"
         statu=nf_def_dim(ncid,'lon',nlon,nlonid)
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_def_dim(ncid,'lat',nlat,nlatid)
         if(statu /= nf_noerr) call handle_err(statu)

         write(*,*)"Define Variables"
         zdims(1)=nlonid
         zdims(2)=nlatid
         statu=nf_def_var(ncid,"z",nf_float,2,zdims,zid)
         if(statu /= nf_noerr) call handle_err(statu)
         udims(1)=nlonid
         udims(2)=nlatid
         statu=nf_def_var(ncid,"u",nf_float,2,udims,uid)
         if(statu /= nf_noerr) call handle_err(statu)
         vdims(1)=nlonid
         vdims(2)=nlatid
         statu=nf_def_var(ncid,"v",nf_float,2,vdims,vid)
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_def_var(ncid,"lon",nf_float,1,nlonid,lonid)
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_def_var(ncid,"lat",nf_float,1,nlatid,latid)
         if(statu /= nf_noerr) call handle_err(statu)
         
         write(*,*)"Put Attributes"
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_put_att_text(ncid,zid,"description",19,&
         "Geopotential Height")
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_put_att_text(ncid,zid,"units",1,"m")
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_put_att_text(ncid,zid,"coordinates",11,&
         "XLONG XLAT")
         
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_put_att_text(ncid,uid,"description",18,&
         "Geostrophic wind_u")
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_put_att_text(ncid,uid,"units",3,"m/s")
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_put_att_text(ncid,uid,"coordinates",11,&
         "XLONG XLAT")
        
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_put_att_text(ncid,vid,"description",18,&
         "Geostrophic wind_v")
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_put_att_text(ncid,vid,"units",3,"m/s")
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_put_att_text(ncid,vid,"coordinates",11,&
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
         statu=nf_put_var_real(ncid,uid,u)
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_put_var_real(ncid,vid,v)
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_put_var_real(ncid,zid,z)
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
         subroutine smooth5(x,p)
          implicit none 
          real x(65,41)
          real p
          integer i,j
          do j=2,40     
                do i=2,64
                x(i,j)=x(i,j)+p/4*(x(i+1,j)+x(i-1,j)+x(i,j+1)+x(i,j-1)&
                        -4*x(i,j))
                end do
          end do
         end 
