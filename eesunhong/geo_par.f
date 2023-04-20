CC   Original program from Andy Gould
CC   Calling sequence of original:
CC	      tau = (t-t0)/te
CC	      call geta(qn,qe,t)
CC	      qw = qe*cos(0.2) + qn*sin(0.2)
CC	      dtau =  piex*qn + piey*qe
CC	      dbeta= -piex*qe + piey*qn
CC	      dbeta = -dbeta
CC	      taup = tau + dtau
CC	      betap = beta + dbeta
CC	      x2 = betap**2 + taup**2
CC	      x = sqrt(x2)
CC              z = x/rho
CC
CC      Info for OGLE 2003-BLG-235/MOA 2003 2003-BLG-53
CCc       alpha, delta = 18:05:16.35  -28:53:42.0 (J2000)
CCc                    = 18:05:18.52  -28:53:41.74  2000.57
CC	alpha = (18 +  12./60 + 2.24/3600.)*15
CC	delta = -(29 + 01./60 + 1.3/3600.)
CC        tfix = 2848.064
c       
c
	subroutine geo_par(qn,qe,hjd,alpha,delta,tfix)
	implicit real*8 (a-h,o-z)
	real*8 sun(3),xpos(3),ypos(3),rad(3),north(3),east(3)
	real*8 spring(3),summer(3)
	data spring/1.,0.,0./
	data summer/0.,0.9174,0.3971/
	data pi/3.1415926535897932385d0/
	ecc = 0.0167
c	vernal = 2719.0 + 1000.
c	vernal = 2719.0
	vernal = 2719.55
	offset = 75
c	offset = 0
	peri   = vernal - offset
	phi = (1 - offset/365.25)*2*pi
c	phi = (1 - 91/365.25)*2*pi
c	phi = 0
	call getpsi(psi,phi,ecc)
c	write(6,*)psi,phi
	costh = (cos(psi) - ecc)/(1-ecc*cos(psi))
	sinth = -sqrt(1-costh**2)
	do 3 i = 1,3
c	   xpos(i) = spring(i)*cos(psi) + summer(i)*sin(psi)
c	   ypos(i) =-spring(i)*sin(psi) + summer(i)*cos(psi)
	   xpos(i) = spring(i)*costh + summer(i)*sinth
	   ypos(i) =-spring(i)*sinth + summer(i)*costh
 3	continue
c	write(6,4)xpos
c	write(6,4)ypos
c	read(5,*)xyz
 4	format(3f10.4)
	north(1) = 0
	north(2) = 0
	north(3) = 1
	radian = 180/pi
ccc	t0 = 2880.6953
	t0 = tfix
c	t0 = t0-1
	rad(1) = cos(alpha/radian)*cos(delta/radian)
	rad(2) = sin(alpha/radian)*cos(delta/radian)
	rad(3) = sin(delta/radian)
	call cross(east,north,rad)
	call dot(e2,east,east)
	do 5 i=1,3
	   east(i) = east(i)/sqrt(e2)
 5	continue
	call cross(north,rad,east)
 6	format(3f7.3)
c	theta = (t0+1 - peri)/365.25*360.
	phi   = (t0+1 - peri)/365.25*2*pi
	call getpsi(psi,phi,ecc)
	qn2 = 0
	qe2 = 0
	do 10 i=1,3
c	   sun(i) = xpos(i)*cos(theta/radian)+ypos(i)*sin(theta/radian)
	   sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
	   qn2 = qn2 + sun(i)*north(i)
	   qe2 = qe2 + sun(i)*east(i)
 10	continue
c	theta = (t0-1 - peri)/365.25*360.
	phi   = (t0-1 - peri)/365.25*2*pi
	call getpsi(psi,phi,ecc)
	qn1 = 0
	qe1 = 0
	do 20 i=1,3
c	   sun(i) = xpos(i)*cos(theta/radian)+ypos(i)*sin(theta/radian)
	   sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
	   qn1 = qn1 + sun(i)*north(i)
	   qe1 = qe1 + sun(i)*east(i)
 20	continue
c	theta = (t0 - peri)/365.25*360.
	phi   = (t0 - peri)/365.25*2*pi
	call getpsi(psi,phi,ecc)
	qn0 = 0
	qe0 = 0
	do 30 i=1,3
c	   sun(i) = xpos(i)*cos(theta/radian)+ypos(i)*sin(theta/radian)
	   sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
	   qn0 = qn0 + sun(i)*north(i)
	   qe0 = qe0 + sun(i)*east(i)
 30	continue
	vn0 = (qn2-qn1)/2
	ve0 = (qe2-qe1)/2
	factor = 365.25*4.74
c	write(6,*)qn0,qe0,vn0*factor,ve0*factor
c	read(5,*)xyz
	t = hjd
c	theta = (t - peri)/365.25*360.
	phi   = (t - peri)/365.25*2*pi
	call getpsi(psi,phi,ecc)
	qn = -qn0 - vn0*(t-t0)
	qe = -qe0 - ve0*(t-t0)
	do 40 i=1,3
c	   sun(i) = xpos(i)*cos(theta/radian)+ypos(i)*sin(theta/radian)
	   sun(i) = xpos(i)*(cos(psi)-ecc) +
     *              ypos(i)*sin(psi)*sqrt(1-ecc**2)
	   qn = qn + sun(i)*north(i)
	   qe = qe + sun(i)*east(i)
 40	continue
c	write(31,11)njd,qn,qe
 11	format(i6,2f9.5)
 100	continue
	return
	end

        subroutine geo_tpar(qtn,qte,t,alpha,delta,flon_obs,flat_obs,
     &                      icheck)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	data pi/3.1415926535897932385d0/
        data t_equinox/4180.50486d0/
ccc        data t_equinox/5276.23056d0/  bug!
        data RearthAU/4.26352325d-5/
        data day_sid/0.997269566319d0/

        radian = 180.d0/pi
        delta_rad = delta/radian

        dphase = (t-t_equinox)/day_sid + 0.5d0
        if(dphase.lt.0.d0) then
          dphase = dphase + int(-dphase) + 1.d0
        else
          dphase = dphase - int(dphase)
        endif
        dalpha = (alpha - 360.d0*dphase)
        cos_delta = cos(delta_rad)
        xsrc = 0.d0
        ysrc = cos_delta
        zsrc = sin(delta_rad)
        cos_latobs = cos(flat_obs/radian)
        dlon_alpha = flon_obs - dalpha
        xobs = sin(dlon_alpha/radian)*cos_latobs
        yobs = cos(dlon_alpha/radian)*cos_latobs
        zobs = sin(flat_obs/radian)

        dot = ysrc*yobs + zsrc*zobs
        if(icheck.eq.1.and.dot.lt.0.d0) then
          write(6,*) 'dot = ',dot
          stop
        endif
c       rotate about the x-axis
        x = xobs
        y = cos(delta_rad)*yobs + sin(delta_rad)*zobs
        z = cos(delta_rad)*zobs - sin(delta_rad)*yobs

ccc        qtn = -z*RearthAU
ccc        qte = x*RearthAU
ccc        qtn = z*RearthAU
ccc        qte = -x*RearthAU
        qtn = -z*RearthAU
        qte = -x*RearthAU

	return
	end

        subroutine sat_par(qsn,qse,alpha,delta,flon_obs,flat_obs,Rsat)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	data pi/3.1415926535897932385d0/

        radian = 180.d0/pi
        delta_rad = delta/radian

ccc        dphase = (t-t_equinox)/day_sid + 0.5d0
        dphase = 0.5d0
        dalpha = (alpha - 360.d0*dphase)
        cos_delta = cos(delta_rad)
        xsrc = 0.d0
        ysrc = cos_delta
        zsrc = sin(delta_rad)
        cos_latobs = cos(flat_obs/radian)
        dlon_alpha = flon_obs - dalpha
        xobs = sin(dlon_alpha/radian)*cos_latobs
        yobs = cos(dlon_alpha/radian)*cos_latobs
        zobs = sin(flat_obs/radian)

        dot = ysrc*yobs + zsrc*zobs
c       rotate about the x-axis
        x = xobs
        y = cos(delta_rad)*yobs + sin(delta_rad)*zobs
        z = cos(delta_rad)*zobs - sin(delta_rad)*yobs

ccc        qsn = -z*Rsat
ccc        qse = -x*Rsat
ccc        qsn = -z*Rsat
ccc        qse = x*Rsat
        qsn = z*Rsat
        qse = x*Rsat

	return
	end

	subroutine getpsi(psi,phi,ecc)
	implicit real*8 (a-h,o-z)
	data pi/3.1415926535897932385d0/
	psi= phi
	do 10 i=1,4
	   fun = psi - ecc*sin(psi)
	   dif = phi - fun
	   der = 1 - ecc*cos(psi)
	   psi = psi + dif/der
 10	continue
	return
	end

c       
	subroutine cross(c,a,b)
	implicit real*8 (a-h,o-z)
	dimension a(3),b(3),c(3)
	c(1) = a(2)*b(3) - b(2)*a(3)
	c(2) = a(3)*b(1) - b(3)*a(1)
	c(3) = a(1)*b(2) - b(1)*a(2)
	return
	end
c
      subroutine dot(c,a,b)
      implicit real*8 (a-h,o-z)
      dimension a(3),b(3)
      c = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end
c
