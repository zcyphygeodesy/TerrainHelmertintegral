!  TerrainHelmertintgrl.f90 
!
!  FUNCTIONS:
!  TerrainHelmertintgrl - Entry point of console application.
!
!****************************************************************************

      program TerrainHelmertintgrl
      !Numerical integral of terrain Helmert condensation effects on various field elements
      !dr - the integral radius (m)
      implicit none
 	character*800::line,line0
      integer row,nk,sn,len,astat(8),i,j,nlon,nlat,knd(6),kk
	real*8::dr,hd(6),hd1(6),hd2(6),rec(800),pi,RAD,mr,ksi,gra,rga,vm1,vm2,grr,tp,zch
	real*8::GRS(6),BLH(3),rln(3),ter(7),cmpn(7),rr,rst(5)
	real*8,allocatable::mu(:,:),dtm(:,:),sfh(:,:)
	integer::status=0
!---------------------------------------------------------------------
      GRS(1)= 3.986004415d14; GRS(2)=6378136.3d0; GRS(3)=1.0826359d-3
      GRS(4) = 7.292115d-5; GRS(5) = 1.d0/298.2564619427d0
   	pi=datan(1.d0)*4.d0;RAD=pi/180.d0;tp=2.67d3;mr=36.d2/RAD
	dr=90.d3!dr - 积分半径 the integral radius (m)
      open(unit=8,file="landtm1m.dat",status="old",iostat=status)
      if(status/=0)goto 902
      read(8,'(a)') line
      call PickRecord(line,len,rec,sn)
      hd(1:6)=rec(1:6)
	nlat=nint((hd(4)-hd(3))/hd(6))
	nlon=nint((hd(2)-hd(1))/hd(5))
	hd(5)=(hd(2)-hd(1))/real(nlon)
	hd(6)=(hd(4)-hd(3))/real(nlat)
 	allocate(dtm(nlat,nlon), stat=astat(1))!the ground digital elevation model
 	allocate(sfh(nlat,nlon), stat=astat(2))!the ground ellipsoidal height grid
 	allocate(mu(nlat,nlon), stat=astat(3))
	if (sum(astat(1:3)) /= 0) then
          close(8);goto 902
      endif
 	do i=1,nlat
	   read(8,*,end=903)(dtm(i,j),j=1,nlon)
         do j=1,nlon
           if(dtm(i,j)<0)dtm(i,j)=0.d0 !set zero in sea area
         enddo
      enddo
903   close(8)
      open(unit=10,file="landbmsurfhgt.dat",status="old",iostat=status)
      if(status/=0)goto 904
      read(10,'(a)') line
      call PickRecord(line,len,rec,sn)
      hd1(1:6)=rec(1:6)
      if(sum(hd1-hd)>1.d-5)then  !格网规格不同 The grid specifications are different
         close(10);goto 904
      endif
 	do i=1,nlat
	   read(10,*,end=905)(sfh(i,j),j=1,nlon)
      enddo
905   close(10)
 	do i=1,nlat
        BLH(1)=hd(3)+(real(i)-0.5d0)*hd(6)
        do j=1,nlon
	    BLH(2)=hd(1)+(real(j)-0.5d0)*hd(5);BLH(3)=sfh(i,j)
          call BLH_RLAT(GRS,BLH,rln);rr=rln(1);zch=dtm(i,j)
	    mu(i,j)=tp*zch*(1.d0+zch/rr+(zch/rr)**2/3.d0)
        enddo
      enddo
      open(unit=8,file="surfhgt.txt",status="old",iostat=status)
      if(status/=0)goto 904
      open(unit=10,file="reslt.txt",status="replace")
      read(8,'(a)') line  !读取头文件 read the file header
      write(10,101)trim(line);kk=0
      do while(.not.eof(8))  
         read(8,'(a)') line0
         call PickRecord(line0,len,rec,sn)
         if(sn<4.or.sn<row)goto 906
         BLH(2)=rec(2);BLH(1)=rec(3);BLH(3)=rec(4)!大地高ellipsoidal height
         call LTerAllBLH(BLH,dtm,sfh,nlat,nlon,hd,dr,GRS,ter)
         call TercmpnBLH(BLH,mu,sfh,nlat,nlon,hd,dr,GRS,cmpn)
         rst(1)=ter(1)-cmpn(1);rst(2)=(ter(2)-cmpn(2))*1.d5
         rst(3)=(ter(3)-cmpn(3))*mr;rst(4)=(ter(4)-cmpn(4))*mr
         rst(5)=(ter(5)-cmpn(5))*1.d9;kk=kk+1
         write(10,101)trim(line0),(rst(i),i=1,5)
         if(kk/25*25==kk)write(*, '(a,i9)'), '    Calculated point number: ',kk
      enddo
906   close(8)
      close(10)
904   deallocate(dtm,sfh,mu)
902   continue
101   format(a,40F12.4)
      write (*,*)'  Complete the computation! The results are saved in the file reslt.txt.'
      pause
      end
