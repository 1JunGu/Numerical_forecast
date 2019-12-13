
        program Difference
         implicit none
         include '/usr/include/netcdf.inc'
         character(len=256) ::filename=&
         "/lustre/data/jungu/homework_zrj/homework_2/output/& 
         hgt.850hPa.199807.nc"
         integer i,j,ncid,statu,latid,lonid,hgtid,Ugid,Vgid,&
                 dimlon,dimlat
         real, parameter :: g=9.8,omega=7.292e-5,Re=6.371e6,Pi=3.14159
         real, parameter :: radian = Pi/180
        
         real hgt(23,33),Ug(23,33),Vg(23,33),lon(33),lat(23),lat_r,&
                dx,dy,f
         integer Ugdims(2)
         !read nc  
         write(*,*) "open ",trim(filename)
         statu=nf_open(trim(filename),nf_write,ncid)
         if(statu /= nf_noerr) call handle_err(statu)
         
         !read hgt
         write(*,*) "Read hgt"
         statu=nf_inq_varid(ncid,"hgt",hgtid)
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_get_var_real(ncid,hgtid,hgt)
         if(statu /= nf_noerr) call handle_err(statu)
                
         write(*,*) "Read lon"
         statu=nf_inq_varid(ncid,"lon",lonid)
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_get_var_real(ncid,lonid,lon)
         if(statu /= nf_noerr) call handle_err(statu)
         
         write(*,*) "Read lat"
         statu=nf_inq_varid(ncid,"lat",latid)
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_get_var_real(ncid,latid,lat)
         if(statu /= nf_noerr) call handle_err(statu) 
        
         write(*,*)"Inq Dimensions"
         statu=nf_inq_dimid(ncid,"lon",dimlon)
         if(statu /= nf_noerr) call handle_err(statu)
         statu=nf_inq_dimid(ncid,"lat",dimlat)
         if(statu /= nf_noerr) call handle_err(statu)
         
         dy = Re*2.5*radian
         Ug=-999.
         Vg=-999.
         
         do i=2,22
                lat_r = lat(i)*radian
                f = 2*omega*sin(lat_r)
                dx = dy*cos(lat_r)
                do j=2,32
                        Ug(i,j)= -g*(hgt(i+1,j)-hgt(i-1,j))/(f*2*dy) 
                        Vg(i,j)= g*(hgt(i,j+1)-hgt(i,j-1))/(f*2*dx)
                enddo
        enddo
        !outer difference
        i=1
        lat_r = lat(i)*radian
        f = 2*omega*sin(lat_r)
        dx = dy*cos(lat_r)
        do j=2,32
                Ug(i,j)=-g*(hgt(i+1,j)-hgt(i,j))/(f*dy)
                Vg(i,j)= g*(hgt(i,j+1)-hgt(i,j-1))/(f*2*dx)
        enddo
        
        i=23
        lat_r = lat(i)*radian
        f = 2*omega*sin(lat_r)
        dx = dy*cos(lat_r)
        do j=2,32
                Ug(i,j)=-g*(hgt(i,j)-hgt(i-1,j))/(f*dy)
                Vg(i,j)= g*(hgt(i,j+1)-hgt(i,j-1))/(f*2*dx)
        enddo
        
        j=1
        do i=2,22
                lat_r = lat(i)*radian
                f =2*omega*sin(lat_r)
                dx = dy*cos(lat_r)
                Ug(i,j)=-g*(hgt(i+1,j)-hgt(i-1,j))/(f*2*dy)
                Vg(i,j)= g*(hgt(i,j+1)-hgt(i,j))/(f*dx)
        enddo        
         
        j=33
        do i=2,22
                lat_r = lat(i)*radian
                f =2*omega*sin(lat_r)
                dx = dy*cos(lat_r)
                Ug(i,j)=-g*(hgt(i+1,j)-hgt(i-1,j))/(f*2*dy)
                Vg(i,j)= g*(hgt(i,j)-hgt(i,j-1))/(f*dx)
        enddo        
        !four corners
        i=1
        lat_r=lat(i)*radian
        f =2*omega*sin(lat_r)
        dx = dy*cos(lat_r)
        Ug(1,1)=-g*(hgt(2,1)-hgt(1,1))/(f*dy)
        Vg(1,1)= g*(hgt(1,2)-hgt(1,1))/(f*dx)
        Ug(1,33)=-g*(hgt(2,33)-hgt(1,33))/(f*dy)
        Vg(1,33)=g*(hgt(1,33)-hgt(1,32))/(f*dx)
        
        i=23
        lat_r=lat(i)*radian
        f =2*omega*sin(lat_r)
        dx = dy*cos(lat_r)
        Ug(23,1)=-g*(hgt(23,1)-hgt(22,1))/(f*dy)
        Vg(23,1)= g*(hgt(23,2)-hgt(23,1))/(f*dx)
        Ug(23,33)=-g*(hgt(23,33)-hgt(22,33))/(f*dy)
        Vg(23,33)=g*(hgt(23,33)-hgt(23,32))/(f*dx)
       
        
        !add two variables named Ug&Vg 
        
        
        write(*,*) "Add geo_wind to",trim(filename)
        statu=nf_redef(ncid)
        if(statu /= nf_noerr) call handle_err(statu)
        
        !write(*,*)"Define Dimensions"
        !statu=nf_def_dim(ncid,'Lat',23,dimlat)
        !if(statu /= nf_noerr) call handle_err(statu)
        !statu=nf_def_dim(ncid,'Lon',33,dimLon)
        !if(statu /= nf_noerr) call handle_err(statu)
        
        write(*,*) "Define Variables"
        Ugdims(2)=dimlat
        Ugdims(1)=dimlon
        statu = nf_def_var(ncid,"Ug",nf_float,2,Ugdims,Ugid)
        if(statu /= nf_noerr) call handle_err(statu)
        statu = nf_def_var(ncid,"Vg",nf_float,2,Ugdims,Vgid)
        if(statu /= nf_noerr) call handle_err(statu)
        
        write(*,*)"Put Attributes"
        statu=nf_put_att_text(ncid,Ugid,"units",3,"m/s")
        if(statu /= nf_noerr) call handle_err(statu)
        statu=nf_put_att_text(ncid,Vgid,"units",3,"m/s")
        if(statu /= nf_noerr) call handle_err(statu)
        statu=nf_put_att_text(ncid,Ugid,"long_name",13,"geostrophic_u")
        if(statu /= nf_noerr) call handle_err(statu)
        statu=nf_put_att_text(ncid,Vgid,"long_name",13,"geostrophic_v")
        if(statu /= nf_noerr) call handle_err(statu)
        statu=nf_put_att_real(ncid,Ugid,"Level",nf_float,1,850.)
        if(statu /= nf_noerr) call handle_err(statu)
        statu=nf_put_att_real(ncid,Vgid,"Level",nf_float,1,850.)
        if(statu /= nf_noerr) call handle_err(statu)
       
        write(*,*)"End Defination"
        statu=nf_enddef(ncid)
        if(statu /= nf_noerr) call handle_err(statu)
        
        write(*,*)"Put Variable"
        statu=nf_put_var_real(ncid,Ugid,Ug)
        if(statu /= nf_noerr) call handle_err(statu)
        statu=nf_put_var_real(ncid,Vgid,Vg)
        if(statu /= nf_noerr) call handle_err(statu)
        
        write(*,*)"Close nc"
        statu=nf_close(ncid)
        if(statu /= nf_noerr) call handle_err(statu)
        
        end program
        subroutine handle_err(statu)
         implicit none
         integer statu
         write(*,*)"Error"
         stop
        end
