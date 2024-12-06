      subroutine TercmpnBLH(BLH,mu,sfh,nlat,nlon,hd,dr,GRS,ter)
      !dr-积分半径m
!-------------------------------------------------------------
      implicit none
	integer::i,j,nlat,nlon,i0,j0,ni,nj,astat(5)
	real*8::dr,mu(nlat,nlon),sfh(nlat,nlon),gr,rln(3),NFD(5)
	real*8::hd(6),gp,pi,RAD,ds,mdr,tt,rr,r0,r1,r2,rst(7),ter(7)
	real*8::BLH(3),XYZ(3),BLH0(3),XYZ0(3),BLH1(3),XYZ1(3),rln1(3)
	real*8::CGrdPntD2,mGrdPntD2,up,u1,L0,L1,L2
	real*8::GRS(6),ge
	real*8 rlat,rlon,rlat1,rlon1,sin2f,cos2f,sina,cosa,sinf
!-----------------------------------------------------------------
      ge=6.67428d-11; pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      BLH0=BLH;BLH0(3)=CGrdPntD2(BLH(2),BLH(1),sfh,nlat,nlon,hd)!计算点正下方地面点位置
      call BLH_XYZ(GRS,BLH,XYZ);call BLH_XYZ(GRS,BLH0,XYZ0)
      call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      r0=dsqrt(XYZ0(1)**2+XYZ0(2)**2+XYZ0(3)**2)
      call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
      up=mGrdPntD2(BLH(2),BLH(1),mu,nlat,nlon,hd)
      ter=0.d0;mdr=r0*hd(5)*RAD*dcos(rln(2)*RAD)/2.d0 !奇异点判断
      ni=nint(dr/r0/RAD/hd(6)+1.d0) !积分半径dr对应的地面格网数
      nj=nint(dr/r0/RAD/hd(5)/dcos(rln(2)*RAD)+1.d0)
	i0=nint((BLH(1)-hd(3))/hd(6)+0.5d0)
	j0=nint((BLH(2)-hd(1))/hd(5)+0.5d0)!计算点所在的地面格网i0,j0
      rlat=rln(2)*RAD;rlon=rln(3)*RAD
	do i=i0-ni,i0+ni
	  if(i<1.or.i>nlat)goto 9100
        BLH1(1)=hd(3)+(real(i)-0.5d0)*hd(6)
	  do j=j0-nj,j0+nj
	    if(j<1.or.j>nlon)goto 9101
	    BLH1(2)=hd(1)+(real(j)-0.5d0)*hd(5)
          BLH1(3)=sfh(i,j);call BLH_XYZ(GRS,BLH1,XYZ1)
          L0=dsqrt((XYZ1(1)-XYZ0(1))**2+(XYZ1(2)-XYZ0(2))**2+(XYZ1(3)-XYZ0(3))**2)
          if(L0>dr)goto 9101
          L1=dsqrt((XYZ1(1)-XYZ(1))**2+(XYZ1(2)-XYZ(2))**2+(XYZ1(3)-XYZ(3))**2)
          if(L1<mdr)then!计算奇异积分
             call Tercmpnsgn(BLH,mu,sfh,nlat,nlon,hd,i,j,4,GRS,rst)
             ter=ter+rst; goto 9101 
          endif
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          u1=mu(i,j);tt=1.d0-2.d0*(L0/r1/2.d0)**2
          rlat1=rln1(2)*RAD;rlon1=rln1(3)*RAD
          ds=hd(5)*hd(6)*RAD**2*dcos(rln1(2)*RAD)*r1**2
          sin2f=dsqrt((1.d0-tt)/2.d0);sinf=2.d0*sin2f**dsqrt((1.d0+tt)/2.d0)
	    cosa=(dcos(rlat)*dsin(rlat1)-dsin(rlat)*dcos(rlat1)*dcos(rlon1-rlon))/sinf
	    sina=dcos(rlat1)*dsin(rlon1-rlon)/sinf
          ter(1)=ter(1)+(u1-up)*ge/L1*ds/gr
          ter(2)=ter(2)+(u1-up)*ge*(rr-r1*tt)/L1**3*ds
          if(dabs(sin2f)>1.d-12)then
            ter(3)=ter(3)+(u1-up)*ge*r1*sinf/L1**3/gr*ds*cosa
            ter(4)=ter(4)+(u1-up)*ge*r1*sinf/L1**3/gr*ds*sina
          endif
          ter(5)=ter(5)+(u1-up)*ge*(3.d0*(rr-r1*tt)/L1**5-1.d0/L1**3)*ds
9101      continue
	  enddo
9100    continue
	enddo
9002	return
      end
!--------------------------------------------------------------------------------
      subroutine Tercmpnsgn(BLH,mu,sfh,nlat,nlon,hd,i0,j0,m,GRS,rst)
      !m-核函数细化为m*m
      !i0,j0-奇异点格网位置
!-------------------------------------------------------------
      implicit none
	integer::m,i,j,nlat,nlon,i0,j0,ni,nj
	real*8::dr,mu(nlat,nlon),sfh(nlat,nlon),gr,rln(3),NFD(5)
	real*8::hd(6),gp,pi,RAD,ds,mdr,tt,rr,r0,r1,r2,rst(7),rv,lon,lat
	real*8::BLH(3),XYZ(3),BLH0(3),XYZ0(3),BLH1(3),XYZ1(3),rln1(3)
	real*8 CGrdPntD2,mGrdPntD2,L0,L1,up,u1
	real*8::GRS(6),ge
	real*8 rlat,rlon,rlat1,rlon1,sin2f,cos2f,sina,cosa,sinf
!-----------------------------------------------------------------
      ge=6.67428d-11; rv=hd(5)/dble(m);pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      lat=hd(3)+real(i0-1)*hd(6);lon=hd(1)+real(j0-1)*hd(5)!格网左下角经纬度
      BLH0=BLH;BLH0(3)=CGrdPntD2(BLH(2),BLH(1),sfh,nlat,nlon,hd)!计算点正下方地面点位置
      call BLH_XYZ(GRS,BLH,XYZ);call BLH_XYZ(GRS,BLH0,XYZ0)
      call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      r0=dsqrt(XYZ0(1)**2+XYZ0(2)**2+XYZ0(3)**2)
      call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
      up=mGrdPntD2(BLH(2),BLH(1),mu,nlat,nlon,hd)
      rst=0.d0;mdr=r0*rv*RAD*dcos(rln(2)*RAD)/dble(m)/2.d0  !奇异点判断
      rlat=rln(2)*RAD;rlon=rln(3)*RAD
	do i=1,m
        BLH1(1)=lat+(real(i)-0.5d0)*rv
	  do j=1,m
	    BLH1(2)=lon+(real(j)-0.5d0)*rv
          BLH1(3)=CGrdPntD2(BLH1(2),BLH1(1),sfh,nlat,nlon,hd)
          call BLH_XYZ(GRS,BLH1,XYZ1)
          L0=dsqrt((XYZ1(1)-XYZ0(1))**2+(XYZ1(2)-XYZ0(2))**2+(XYZ1(3)-XYZ0(3))**2)
          L1=dsqrt((XYZ1(1)-XYZ(1))**2+(XYZ1(2)-XYZ(2))**2+(XYZ1(3)-XYZ(3))**2)
          if(L1<mdr)L1=mdr
          u1=mGrdPntD2(BLH1(2),BLH1(1),mu,nlat,nlon,hd)
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          ds=rv**2*RAD**2*dcos(rln1(2)*RAD)*r1**2
          tt=1.d0-2.d0*(L0/r1/2.d0)**2
          sin2f=dsqrt((1.d0-tt)/2.d0);sinf=2.d0*sin2f**dsqrt((1.d0+tt)/2.d0)
          if(sinf<1.d-10)sinf=1.d-10
	    cosa=(dcos(rlat)*dsin(rlat1)-dsin(rlat)*dcos(rlat1)*dcos(rlon1-rlon))/sinf
	    sina=dcos(rlat1)*dsin(rlon1-rlon)/sinf
          rst(1)=rst(1)+(u1-up)*ge/L1*ds/gr
          rst(2)=rst(2)+(u1-up)*ge*(rr-r1*tt)/L1**3*ds
          rst(5)=rst(5)+(u1-up)*ge*(3.d0*(rr-r1*tt)/L1**5-1.d0/L1**3)*ds
	  enddo
	enddo
9002	return
      end
