c==============================================================================

       subroutine trilens(eps2_in,eps3_in,sep_in,sep2_in,angle)

c==============================================================================
c
c  Routines in this file written by David Bennett and Sun Hong Rhie.
c
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       double precision amp(10)
       double complex zz(10),z(10),cff(11),aa,bb,cc,dd
       double complex zx(10),work(100)
       double complex zk,ztos,czkaa,czkbb,cjac
       double complex zta,zta1,zta2,zta0,cx1,cx2
       double complex zb,zb0,zb1,zb2,ctempb,czer0
       double complex hh39,hh38,hh37,hh36,hh35,hh34,hh33,hh32,hh31 
       double complex hh30,hh28,hh27,hh26,hh25,hh24,hh23,hh22,hh21 
       double complex hh20,hh17,hh16,hh15,hh14,hh13,hh12,hh11,hh10 
       double complex hh06,hh05,hh04,hh03,hh02,hh01,hh00 
       double complex ww,ww1,ww2,ww3,wwbar,ww1bar,ww2bar,ww3bar
       double complex wwaa,wwbb,wwcc,wwdd,cc1,cc2,cc3,cc4
       double complex czkcc1,czkcc2,czkcc3,ctemp,ctemp1,ctemp2
       double complex ctemp3
       common/lenscoords/xx1,xx2,xx3,yy1,yy2,yy3
       common/q_count/n_q,n_noq,noquad
       integer iimage(11)
       logical polish

       save

c  yes, polish the roots
c  ---------------------
       polish = .true.

ccccc  CHANGED FROM HERE ON   ccccccccccccc 
ccccc  lens parameters  ccccccccccccccccccc
       eps2=eps2_in
       eps3=eps3_in
       sep=sep_in
       sep2 = sep2_in
       ang =angle

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  CENTER OF MASS SYSTEM           ccccccccccccccccc
cccccc  epsilon_1 is on the real axis at xx1    ccccccccc 
cccccc  Make sure eps4 does NOT  VANISH        cccccccccc
cccccc  because I am going to divide by it        ccccccc
cccccc  eps4 = 0 corresponds to  a single lens    ccccccc
cccccc  where  eps1 = 1.         cccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       eps4 = eps2+eps3

       if(eps4.le.0.d0) return
       

       czer0 = dcmplx(0.d0,0.d0)
       eps1 = 1.-eps4
       xx1 = -eps4*sep
       xx4 =  eps1*sep
       xx2 = xx4 +eps3/eps4*sep2* cos(ang)
       yy2 = eps3/eps4*sep2* sin(ang)  
       xx3 = xx4 - eps2/eps4*sep2* cos(ang) 
       yy3 = - eps2/eps4*sep2* sin(ang)

       cc1 = dcmplx(xx1,0.d0)
       cc2 = dcmplx(xx2,yy2)
       cc3 = dcmplx(xx3,yy3)
       cc4 = dcmplx(xx4,0.d0)
      
       aa = -(cc1+cc2+cc3) 
       bb = cc1*cc2 + cc1*cc3 + cc2*cc3 
       cc = -cc1*cc2*cc3
       dd = eps1*cc2*cc3 + eps2*cc1*cc3 + eps3*cc1*cc2

       hh39 = dcmplx(1.d0,0.d0)
       hh38 = 3.d0*aa
       hh37 = 3.d0*bb + 3.d0*aa*aa 
       hh36 = 3.d0*cc + 6.d0*aa*bb + aa*aa*aa
       hh35 = 6.d0*aa*cc + 3.d0*bb*bb + 3.d0*aa*aa*bb
       hh34 = 6.d0*bb*cc + 3.d0*aa*aa*cc + 3.d0*aa*bb*bb
       hh33 = 3.d0*cc*cc + 6.d0*aa*bb*cc + bb*bb*bb
       hh32 = 3.d0*aa*cc*cc + 3.d0*bb*bb*cc
       hh31 = 3.d0*bb*cc*cc
       hh30 = cc*cc*cc

       hh28 = dcmplx(1.d0,0.d0)
       hh27 = 3.d0*aa
       hh26 = dd + 2.d0*bb + 3.d0*aa*aa
       hh25 = 2.d0*aa*dd + 4.d0*aa*bb + aa*aa*aa + 2.d0*cc
       hh24 = 2.d0*dd*bb + dd*aa*aa + 4.d0*aa*cc 
     &       +2.d0*aa*aa*bb + bb*bb
       hh23 = 2.d0*dd*cc + 2.d0*dd*aa*bb + 2.d0*aa*aa*cc 
     &       +aa*bb*bb + 2.d0*bb*cc
       hh22 = 2.d0*cc*aa*dd + dd*bb*bb + 2.d0*aa*bb*cc + cc*cc
       hh21 = 2.d0*bb*cc*dd + aa*cc*cc
       hh20 = cc*cc*dd

       hh17 = dcmplx(1.d0,0.d0)
       hh16 = 3.d0*aa
       hh15 = 2.d0*dd + 3.d0*aa*aa + bb
       hh14 = 4.d0*aa*dd + aa*aa*aa + 2.d0*aa*bb + cc
       hh13 = dd*dd + 2.d0*aa*aa*dd + 2.d0*bb*dd 
     &       +bb*aa*aa + 2.d0*aa*cc
       hh12 = aa*dd*dd + 2.d0*aa*bb*dd + 2.d0*cc*dd + cc*aa*aa
       hh11 = bb*dd*dd + 2.d0*aa*cc*dd
       hh10 = cc*dd*dd

       hh06 = dcmplx(1.d0,0.d0)
       hh05 = 3.d0*aa
       hh04 = 3.d0*dd + 3.d0*aa*aa
       hh03 = 6.d0*aa*dd + aa*aa*aa
       hh02 = 3.d0*dd*dd + 3.d0*aa*aa*dd
       hh01 = 3.d0*aa*dd*dd
       hh00 = dd*dd*dd

cccccccccccccccccccccccccccccccccccccccccccc
ccccccc  MORE CHANGES MADE (5/04/99)  cccccc
ccccccc  COEFFICIENTS BELOW  ccccccccccccccc
cccccckccccccccccccccccccccccccccccccccccccc

       tol = 1.d-4
       nimage0=3

       return

c    entry trilens_im calculates the amplification of each image
c    -----------------------------------------------------------
      entry trilens_im(sx,sy,nimage,iimage,zz,amp)

      if(eps4.le.0.) then
c       single lens case
c       ----------------
        nimage = 2
        iimage(1) = 1
        u2 = sx**2+sy**2
        u = sqrt(u)
        a = (u2+2.d0)/(u*sqrt(u2+4.d0))
        amp(1) = 0.5d0*(a+1.d0)
        amp(2) = 0.5d0*(a-1.d0)
        r1fac = 0.5d0*(1.d0+sqrt(u2+4.d0)/u)
        r2fac = 0.5d0*(1.d0-sqrt(u2+4.d0)/u)
        zr1 = sx*r1fac
        zi1 = sy*r1fac
        zr2 = sx*r2fac
        zi2 = sy*r2fac
        zz(1) = dcmplx(zr1,zi1)
        zz(2) = dcmplx(zr2,zi2)
        return
      endif

c     special case (avoid singular solution)
c     --------------------------------------
      if(sep.le.1.1d-4.and.sx.eq.0..and.sy.eq.0.) sy=1.d-5

cccccccccccccccccccccccccccccccccccccccccccc
cccccccccc  MORE LINES from here ccccccccccc
cccccccccccc  COEFFICIENTS ccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccc

       ww  = dcmplx(sx,sy)
       ww1 = dcmplx(sx-xx1,sy)
       ww2 = dcmplx(sx-xx2,sy-yy2)
       ww3 = dcmplx(sx-xx3,sy-yy3)
 
       wwbar = conjg(ww)
       ww1bar = conjg(ww1)
       ww2bar = conjg(ww2)
       ww3bar = conjg(ww3)

       wwaa = ww1bar+ww2bar+ww3bar
       wwbb = ww1bar*ww2bar + ww2bar*ww3bar 
     &        + ww1bar*ww3bar
       wwcc = ww1bar*ww2bar*ww3bar
       wwdd = eps1*ww2bar*ww3bar + eps2*ww1bar*ww3bar 
     &        + eps3*ww1bar*ww2bar

ccccccc  so far, definitions cccccccccccccc
ccccccc  Finally the coefficients ccccccccc

       cff(11) = hh39*wwcc
       cff(10) = hh38*wwcc + hh28*wwbb - (ww*wwcc+wwdd)*hh39
       cff(9)  = hh37*wwcc + hh27*wwbb + hh17*wwaa
     &     - (ww*wwcc + wwdd)*hh38 - (ww*wwbb + wwaa - wwbar)*hh28
       cff(8)  = hh36*wwcc + hh26*wwbb + hh16*wwaa + hh06 
     &     - (ww*wwcc + wwdd)*hh37 - (ww*wwbb + wwaa-wwbar)*hh27  
     &     - (ww*wwaa + 1)*hh17 
       cff(7)  = hh35*wwcc + hh25*wwbb + hh15*wwaa + hh05
     &     - (ww*wwcc + wwdd)*hh36 - (ww*wwbb + wwaa-wwbar)*hh26
     &     - (ww*wwaa + 1)*hh16  - ww*hh06 
       cff(6)  = hh34*wwcc + hh24*wwbb + hh14*wwaa + hh04
     &     - (ww*wwcc + wwdd)*hh35 - (ww*wwbb + wwaa-wwbar)*hh25
     &     - (ww*wwaa + 1)*hh15  - ww*hh05
       cff(5)  = hh33*wwcc + hh23*wwbb + hh13*wwaa + hh03
     &     - (ww*wwcc + wwdd)*hh34 - (ww*wwbb + wwaa-wwbar)*hh24
     &     - (ww*wwaa + 1)*hh14  - ww*hh04
       cff(4)  = hh32*wwcc + hh22*wwbb + hh12*wwaa + hh02
     &     - (ww*wwcc + wwdd)*hh33 - (ww*wwbb + wwaa-wwbar)*hh23
     &     - (ww*wwaa + 1)*hh13  - ww*hh03
       cff(3)  = hh31*wwcc + hh21*wwbb + hh11*wwaa + hh01
     &     - (ww*wwcc + wwdd)*hh32 - (ww*wwbb + wwaa-wwbar)*hh22
     &     - (ww*wwaa + 1)*hh12  - ww*hh02
       cff(2)  = hh30*wwcc + hh20*wwbb + hh10*wwaa + hh00
     &     - (ww*wwcc + wwdd)*hh31 - (ww*wwbb + wwaa-wwbar)*hh21
     &     - (ww*wwaa + 1)*hh11  - ww*hh01
       cff(1)  = 
     &     - (ww*wwcc + wwdd)*hh30 - (ww*wwbb + wwaa-wwbar)*hh20
     &     - (ww*wwaa + 1)*hh10  - ww*hh00



       if(cff(11).ne.0.) then
          m=10
       elseif(cff(10).ne.0.) then
          m=9
       elseif(cff(9).ne.0.) then
          m=8
       elseif(cff(8).ne.0.) then
          m=7
       elseif(cff(7).ne.0.) then
          m=6
       elseif(cff(6).ne.0.) then
          m=5
       elseif(cff(5).ne.0.) then
          m=4
       elseif(cff(4).ne.0.) then
          m=3
       elseif(cff(3).ne.0.) then
          m=2
       elseif(cff(2).ne.0.) then
          m=1
       else
          m=0
       endif


cccccccccc  COMPLETE !!!!          cccccccccccc
ccccccccccc  I think            ccccccccccccccc
cccccc   If it doesn't work, let me know  ccccc
ccccccccccccccccccccccccccccccccccccccccccccccc



c   find the roots
c   --------------
       call zroots(cff,m,z,polish)
ccc       call czero(zx,cff,m,work)

c   find the amplification
c   ----------------------
       ttol_last = -1.
       ttol=tol
       nimage=0
       ampsum=0.d0
 44    continue

ccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  MAY 12 1999  ---  NEW                    cccccccccccccccccc
cccccc  There are up to  TEN images                     ccccccccccc
cc  Lens positions cc2 and cc3 are genunie  complex numbers ccc 
cc  That is why I do not use cx1 any more  
cc  cx1 and cx2 give you an impression that they are actually real
ccccccccccccccccccccccccccccccccccccccccccccccccc

       do 50 k = 1, 10 
          zk = z(k)
c         the following is the lens equation
c         ztos = zk - conjg((zk-cc)/((zk-cc1)*(zk-cc2)))  
c              !! correct only for binary  not correct for ternary.

             czkcc1 = conjg(zk-cc1)
             czkcc2 = conjg(zk-cc2)
             czkcc3 = conjg(zk-cc3)

          ctemp= czkcc1 * czkcc2 * czkcc3
          if(ctemp.eq.czer0) then
            ztos = zk - dcmplx(1.d30,0.d0)
          else
            ztos = zk - eps1/czkcc1 - eps2/czkcc2 
     &                - eps3/czkcc3
          endif
cc           ww=complx(sx,sy) is the input source position
c         check if that (ww - ztos) = 0
ccc          zerr = abs(imag(ww-ztos)) + abs(real(ww-ztos))     
          zerr = abs(dimag(ww-ztos)) + abs(dreal(ww-ztos))     
c         possibly calculate the amplitude
c         --------------------------------
          if(zerr.le.ttol.and.zerr.gt.ttol_last) then
             czkcc1 = conjg(zk-cc1)
             czkcc2 = conjg(zk-cc2)
             czkcc3 = conjg(zk-cc3) 
c             cjac = eps1/(czkcc1*czkcc1)+eps2/(czkcc2*czkcc2)
c                    + eps3/(czkcc3*czkcc3)
             ctemp1=(czkcc1*czkcc1)
             ctemp2=(czkcc2*czkcc2)
             ctemp3=(czkcc3*czkcc3)
             if(ctemp1.eq.czer0) ctemp1= dcmplx(1.d-30,0.d0)
             if(ctemp2.eq.czer0) ctemp2= dcmplx(1.d-30,0.d0)
             if(ctemp3.eq.czer0) ctemp3= dcmplx(1.d-30,0.d0)
             cjac = eps1/ctemp1+eps2/ctemp2+eps3/ctemp3
ccc             fjac = 1. - ((real(cjac))**2 + (imag(cjac))**2)
             fjac = 1. - ((dreal(cjac))**2 + (dimag(cjac))**2)
             nimage=nimage+1
             iimage(nimage)=k
             amp(nimage) = abs(1./(fjac+1.d-20))
             ampsum = ampsum + amp(nimage)
             zz(nimage) = z(k)
          endif
 50    continue

c      if we've missed an image, try again
cccccc   minimum number of images for ternary is 4  cccccccccc
cccccc   it can be 4, 6, 8, or 10                ccccccccccccc
c      -------------------------------------------------------
ccc       if(nimage.lt.nimage0.and.nimage.ne.4) then
cc       if(nimage.lt.nimage0) then
       if(nimage.lt.nimage0.or.ampsum.lt.0.99999d0.or.
     &   (ampsum.gt.1.003d0.and.nimage.lt.4) )         then
         if(ttol.lt.3.*tol) then
           ttol_last = ttol
           ttol=ttol*2.
ccc           write(6,*) 'WARNING:',nimage,
ccc     &                ' images; trying again w/ tol =',ttol
           go to 44
         elseif(ampsum.lt.0.99999d0.or.ampsum.gt.1.003d0.or.
     &                                ampsum0.gt.1.003d0) then
           call trilensq_im(sx,sy,nimage,iimage,zz,amp)
ccc           write(6,*) 'trilensq called at sx sy =',sx,sy
           ampsum=0.d0
           do im = 1,nimage
             ampsum = ampsum + amp(im)
           enddo
           n_q = n_q + 1
           n_noq = n_noq - 1
         endif
       endif
       nimage0=nimage
       ampsum0 = ampsum
       n_noq = n_noq + 1

       return
       end

c==============================================================================

       subroutine trilenso(eps2_in,eps3_in,sep_in,sep2_in,angle)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       double precision amp(10)
       double complex zz(10),z(10),cff(11),aa,bb,cc,dd
       double complex zx(10),work(100)
       double complex zk,ztos,czkaa,czkbb,cjac
       double complex zta,zta1,zta2,zta0,cx1,cx2
       double complex zb,zb0,zb1,zb2,ctempb,czer0
       double complex hh39,hh38,hh37,hh36,hh35,hh34,hh33,hh32,hh31 
       double complex hh30,hh28,hh27,hh26,hh25,hh24,hh23,hh22,hh21 
       double complex hh20,hh17,hh16,hh15,hh14,hh13,hh12,hh11,hh10 
       double complex hh06,hh05,hh04,hh03,hh02,hh01,hh00 
       double complex ww,ww1,ww2,ww3,wwbar,ww1bar,ww2bar,ww3bar
       double complex wwaa,wwbb,wwcc,wwdd,cc1,cc2,cc3,cc4
       double complex czkcc1,czkcc2,czkcc3,ctemp,ctemp1,ctemp2
       double complex ctemp3
       common/lenscoords/xx1,xx2,xx3,yy1,yy2,yy3
       common/q_count/n_q,n_noq,noquad
       integer iimage(11)
       logical polish

       save

c  yes, polish the roots
c  ---------------------
       polish = .true.

ccccc  CHANGED FROM HERE ON   ccccccccccccc 
ccccc  lens parameters  ccccccccccccccccccc
       eps2=eps2_in
       eps3=eps3_in
       sep=sep_in
       sep2 = sep2_in
       ang =angle

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  CENTER OF MASS SYSTEM           ccccccccccccccccc
cccccc  epsilon_1 is on the real axis at xx1    ccccccccc 
cccccc  Make sure eps4 does NOT  VANISH        cccccccccc
cccccc  because I am going to divide by it        ccccccc
cccccc  eps4 = 0 corresponds to  a single lens    ccccccc
cccccc  where  eps1 = 1.         cccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       eps4 = eps2+eps3

       if(eps4.le.0.d0) return
       

       czer0 = dcmplx(0.d0,0.d0)
       eps1 = 1.-eps4
       xx1 = -eps4*sep
       xx4 =  eps1*sep
       xx2 = xx4 +eps3/eps4*sep2* cos(ang)
       yy2 = eps3/eps4*sep2* sin(ang)  
       xx3 = xx4 - eps2/eps4*sep2* cos(ang) 
       yy3 = - eps2/eps4*sep2* sin(ang)

       cc1 = dcmplx(xx1,0.d0)
       cc2 = dcmplx(xx2,yy2)
       cc3 = dcmplx(xx3,yy3)
       cc4 = dcmplx(xx4,0.d0)
      
       aa = -(cc1+cc2+cc3) 
       bb = cc1*cc2 + cc1*cc3 + cc2*cc3 
       cc = -cc1*cc2*cc3
       dd = eps1*cc2*cc3 + eps2*cc1*cc3 + eps3*cc1*cc2

       hh39 = dcmplx(1.d0,0.d0)
       hh38 = 3.d0*aa
       hh37 = 3.d0*bb + 3.d0*aa*aa 
       hh36 = 3.d0*cc + 6.d0*aa*bb + aa*aa*aa
       hh35 = 6.d0*aa*cc + 3.d0*bb*bb + 3.d0*aa*aa*bb
       hh34 = 6.d0*bb*cc + 3.d0*aa*aa*cc + 3.d0*aa*bb*bb
       hh33 = 3.d0*cc*cc + 6.d0*aa*bb*cc + bb*bb*bb
       hh32 = 3.d0*aa*cc*cc + 3.d0*bb*bb*cc
       hh31 = 3.d0*bb*cc*cc
       hh30 = cc*cc*cc

       hh28 = dcmplx(1.d0,0.d0)
       hh27 = 3.d0*aa
       hh26 = dd + 2.d0*bb + 3.d0*aa*aa
       hh25 = 2.d0*aa*dd + 4.d0*aa*bb + aa*aa*aa + 2.d0*cc
       hh24 = 2.d0*dd*bb + dd*aa*aa + 4.d0*aa*cc 
     &       +2.d0*aa*aa*bb + bb*bb
       hh23 = 2.d0*dd*cc + 2.d0*dd*aa*bb + 2.d0*aa*aa*cc 
     &       +aa*bb*bb + 2.d0*bb*cc
       hh22 = 2.d0*cc*aa*dd + dd*bb*bb + 2.d0*aa*bb*cc + cc*cc
       hh21 = 2.d0*bb*cc*dd + aa*cc*cc
       hh20 = cc*cc*dd

       hh17 = dcmplx(1.d0,0.d0)
       hh16 = 3.d0*aa
       hh15 = 2.d0*dd + 3.d0*aa*aa + bb
       hh14 = 4.d0*aa*dd + aa*aa*aa + 2.d0*aa*bb + cc
       hh13 = dd*dd + 2.d0*aa*aa*dd + 2.d0*bb*dd 
     &       +bb*aa*aa + 2.d0*aa*cc
       hh12 = aa*dd*dd + 2.d0*aa*bb*dd + 2.d0*cc*dd + cc*aa*aa
       hh11 = bb*dd*dd + 2.d0*aa*cc*dd
       hh10 = cc*dd*dd

       hh06 = dcmplx(1.d0,0.d0)
       hh05 = 3.d0*aa
       hh04 = 3.d0*dd + 3.d0*aa*aa
       hh03 = 6.d0*aa*dd + aa*aa*aa
       hh02 = 3.d0*dd*dd + 3.d0*aa*aa*dd
       hh01 = 3.d0*aa*dd*dd
       hh00 = dd*dd*dd

cccccccccccccccccccccccccccccccccccccccccccc
ccccccc  MORE CHANGES MADE (5/04/99)  cccccc
ccccccc  COEFFICIENTS BELOW  ccccccccccccccc
cccccckccccccccccccccccccccccccccccccccccccc

       tol = 1.d-4
       nimage0=3

       return

c    entry trilenso_im calculates the amplification of each image
c    ------------------------------------------------------------
      entry trilenso_im(sx,sy,nimage,iimage,zz,amp)

      if(eps4.le.0.) then
c       single lens case
c       ----------------
        nimage = 2
        iimage(1) = 1
        u2 = sx**2+sy**2
        u = sqrt(u)
        a = (u2+2.d0)/(u*sqrt(u2+4.d0))
        amp(1) = 0.5d0*(a+1.d0)
        amp(2) = 0.5d0*(a-1.d0)
        r1fac = 0.5d0*(1.d0+sqrt(u2+4.d0)/u)
        r2fac = 0.5d0*(1.d0-sqrt(u2+4.d0)/u)
        zr1 = sx*r1fac
        zi1 = sy*r1fac
        zr2 = sx*r2fac
        zi2 = sy*r2fac
        zz(1) = dcmplx(zr1,zi1)
        zz(2) = dcmplx(zr2,zi2)
        return
      endif

c     special case (avoid singular solution)
c     --------------------------------------
      if(sep.le.1.1d-4.and.sx.eq.0..and.sy.eq.0.) sy=1.d-5

cccccccccccccccccccccccccccccccccccccccccccc
cccccccccc  MORE LINES from here ccccccccccc
cccccccccccc  COEFFICIENTS ccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccc

       ww  = dcmplx(sx,sy)
       ww1 = dcmplx(sx-xx1,sy)
       ww2 = dcmplx(sx-xx2,sy-yy2)
       ww3 = dcmplx(sx-xx3,sy-yy3)
 
       wwbar = conjg(ww)
       ww1bar = conjg(ww1)
       ww2bar = conjg(ww2)
       ww3bar = conjg(ww3)

       wwaa = ww1bar+ww2bar+ww3bar
       wwbb = ww1bar*ww2bar + ww2bar*ww3bar 
     &        + ww1bar*ww3bar
       wwcc = ww1bar*ww2bar*ww3bar
       wwdd = eps1*ww2bar*ww3bar + eps2*ww1bar*ww3bar 
     &        + eps3*ww1bar*ww2bar

ccccccc  so far, definitions cccccccccccccc
ccccccc  Finally the coefficients ccccccccc

       cff(11) = hh39*wwcc
       cff(10) = hh38*wwcc + hh28*wwbb - (ww*wwcc+wwdd)*hh39
       cff(9)  = hh37*wwcc + hh27*wwbb + hh17*wwaa
     &     - (ww*wwcc + wwdd)*hh38 - (ww*wwbb + wwaa - wwbar)*hh28
       cff(8)  = hh36*wwcc + hh26*wwbb + hh16*wwaa + hh06 
     &     - (ww*wwcc + wwdd)*hh37 - (ww*wwbb + wwaa-wwbar)*hh27  
     &     - (ww*wwaa + 1)*hh17 
       cff(7)  = hh35*wwcc + hh25*wwbb + hh15*wwaa + hh05
     &     - (ww*wwcc + wwdd)*hh36 - (ww*wwbb + wwaa-wwbar)*hh26
     &     - (ww*wwaa + 1)*hh16  - ww*hh06 
       cff(6)  = hh34*wwcc + hh24*wwbb + hh14*wwaa + hh04
     &     - (ww*wwcc + wwdd)*hh35 - (ww*wwbb + wwaa-wwbar)*hh25
     &     - (ww*wwaa + 1)*hh15  - ww*hh05
       cff(5)  = hh33*wwcc + hh23*wwbb + hh13*wwaa + hh03
     &     - (ww*wwcc + wwdd)*hh34 - (ww*wwbb + wwaa-wwbar)*hh24
     &     - (ww*wwaa + 1)*hh14  - ww*hh04
       cff(4)  = hh32*wwcc + hh22*wwbb + hh12*wwaa + hh02
     &     - (ww*wwcc + wwdd)*hh33 - (ww*wwbb + wwaa-wwbar)*hh23
     &     - (ww*wwaa + 1)*hh13  - ww*hh03
       cff(3)  = hh31*wwcc + hh21*wwbb + hh11*wwaa + hh01
     &     - (ww*wwcc + wwdd)*hh32 - (ww*wwbb + wwaa-wwbar)*hh22
     &     - (ww*wwaa + 1)*hh12  - ww*hh02
       cff(2)  = hh30*wwcc + hh20*wwbb + hh10*wwaa + hh00
     &     - (ww*wwcc + wwdd)*hh31 - (ww*wwbb + wwaa-wwbar)*hh21
     &     - (ww*wwaa + 1)*hh11  - ww*hh01
       cff(1)  = 
     &     - (ww*wwcc + wwdd)*hh30 - (ww*wwbb + wwaa-wwbar)*hh20
     &     - (ww*wwaa + 1)*hh10  - ww*hh00



       if(cff(11).ne.0.) then
          m=10
       elseif(cff(10).ne.0.) then
          m=9
       elseif(cff(9).ne.0.) then
          m=8
       elseif(cff(8).ne.0.) then
          m=7
       elseif(cff(7).ne.0.) then
          m=6
       elseif(cff(6).ne.0.) then
          m=5
       elseif(cff(5).ne.0.) then
          m=4
       elseif(cff(4).ne.0.) then
          m=3
       elseif(cff(3).ne.0.) then
          m=2
       elseif(cff(2).ne.0.) then
          m=1
       else
          m=0
       endif


cccccccccc  COMPLETE !!!!          cccccccccccc
ccccccccccc  I think            ccccccccccccccc
cccccc   If it doesn't work, let me know  ccccc
ccccccccccccccccccccccccccccccccccccccccccccccc



c   find the roots
c   --------------
       call zroots(cff,m,z,polish)
ccc       call czero(zx,cff,m,work)

c   find the amplification
c   ----------------------
       ttol_last = -1.
       ttol=tol
       nimage=0
       ampsum=0.d0
 44    continue

ccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  MAY 12 1999  ---  NEW                    cccccccccccccccccc
cccccc  There are up to  TEN images                     ccccccccccc
cc  Lens positions cc2 and cc3 are genunie  complex numbers ccc 
cc  That is why I do not use cx1 any more  
cc  cx1 and cx2 give you an impression that they are actually real
ccccccccccccccccccccccccccccccccccccccccccccccccc

       do 50 k = 1, 10 
          zk = z(k)
c         the following is the lens equation
c         ztos = zk - conjg((zk-cc)/((zk-cc1)*(zk-cc2)))  
c              !! correct only for binary  not correct for ternary.

             czkcc1 = conjg(zk-cc1)
             czkcc2 = conjg(zk-cc2)
             czkcc3 = conjg(zk-cc3)

          ctemp= czkcc1 * czkcc2 * czkcc3
          if(ctemp.eq.czer0) then
            ztos = zk - dcmplx(1.d30,0.d0)
          else
            ztos = zk - eps1/czkcc1 - eps2/czkcc2 
     &                - eps3/czkcc3
          endif
cc           ww=complx(sx,sy) is the input source position
c         check if that (ww - ztos) = 0
ccc          zerr = abs(imag(ww-ztos)) + abs(real(ww-ztos))     
          zerr = abs(dimag(ww-ztos)) + abs(dreal(ww-ztos))     
c         possibly calculate the amplitude
c         --------------------------------
          if(zerr.le.ttol.and.zerr.gt.ttol_last) then
             czkcc1 = conjg(zk-cc1)
             czkcc2 = conjg(zk-cc2)
             czkcc3 = conjg(zk-cc3) 
c             cjac = eps1/(czkcc1*czkcc1)+eps2/(czkcc2*czkcc2)
c                    + eps3/(czkcc3*czkcc3)
             ctemp1=(czkcc1*czkcc1)
             ctemp2=(czkcc2*czkcc2)
             ctemp3=(czkcc3*czkcc3)
             if(ctemp1.eq.czer0) ctemp1= dcmplx(1.d-30,0.d0)
             if(ctemp2.eq.czer0) ctemp2= dcmplx(1.d-30,0.d0)
             if(ctemp3.eq.czer0) ctemp3= dcmplx(1.d-30,0.d0)
             cjac = eps1/ctemp1+eps2/ctemp2+eps3/ctemp3
ccc             fjac = 1. - ((real(cjac))**2 + (imag(cjac))**2)
             fjac = 1. - ((dreal(cjac))**2 + (dimag(cjac))**2)
             nimage=nimage+1
             iimage(nimage)=k
             amp(nimage) = abs(1./(fjac+1.d-20))
             ampsum = ampsum + amp(nimage)
             zz(nimage) = z(k)
          endif
 50    continue

c      if we've missed an image, try again
cccccc   minimum number of images for ternary is 4  cccccccccc
cccccc   it can be 4, 6, 8, or 10                ccccccccccccc
c      -------------------------------------------------------
ccc       if(nimage.lt.nimage0.and.nimage.ne.4) then
cc       if(nimage.lt.nimage0) then
       if(nimage.lt.nimage0.or.ampsum.lt.0.99999d0.or.
     &   (ampsum.gt.1.06d0.and.nimage.lt.4) )         then
         if(ttol.lt.3.*tol) then
           ttol_last = ttol
           ttol=ttol*2.
ccc           write(6,*) 'WARNING:',nimage,
ccc     &                ' images; trying again w/ tol =',ttol
           go to 44
         elseif((ampsum.lt.0.99999d0.or.ampsum.gt.1.06d0.or.
     &           ampsum0.gt.1.06d0).and.noquad.ne.1) then
           call trilensq(eps2,eps3,sep,sep2,ang)
           call trilensq_im(sx,sy,nimage,iimage,zz,amp)
ccc           write(6,*) 'trilensq called at sx sy =',sx,sy
           ampsum=0.d0
           do im = 1,nimage
             ampsum = ampsum + amp(im)
           enddo
           n_q = n_q + 1
           n_noq = n_noq - 1
         endif
       endif
       nimage0=nimage
       ampsum0 = ampsum
       n_noq = n_noq + 1

       return
       end

c==============================================================================

       subroutine trilensqp(eps2_in,eps3_in,sep_in,sep2_in,angle)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       parameter(diffmin2=1.e-25)

       double precision amp(10)
       complex*32 zq(10),cffq(11)
       complex*32 zqp(10)
       complex*32 ctempq,czkqcc1,czkqcc2,czkqcc3
       complex*32 czer0q,zqtos
       complex*32 ctemp1q,ctemp2q,ctemp3q,cjacq
       complex*32 aaq,bbq,ccq,ddq
       complex*32 hh39q,hh38q,hh37q,hh36q,hh35q,hh34q,hh33q,hh32q,hh31q
       complex*32 hh30q,hh28q,hh27q,hh26q,hh25q,hh24q,hh23q,hh22q,hh21q
       complex*32 hh20q,hh17q,hh16q,hh15q,hh14q,hh13q,hh12q,hh11q,hh10q
       complex*32 hh06q,hh05q,hh04q,hh03q,hh02q,hh01q,hh00q
       complex*32 wwq,ww1q,ww2q,ww3q,wwbarq,ww1barq,ww2barq,ww3barq
       complex*32 wwaaq,wwbbq,wwccq,wwddq,cc1q,cc2q,cc3q,cc4q
       complex*32 czkcc1q,czkcc2q,czkcc3q
       complex*32 czkqpcc1q,czkqpcc2q,czkqpcc3q
       complex*32 ctemp1qpq,ctemp2qpq,ctemp3qpq,cjacqpq
       real*16 eps1q,eps2q,eps3q
       double complex zz(10),z(10),cff(11),aa,bb,cc,dd
       double complex zx(10),work(100)
       double complex znp(10)
       double complex zk,ztos,czkaa,czkbb,cjac
       double complex zkqp,zqptos
       double complex zkq,zdiff
       double complex zta,zta1,zta2,zta0,cx1,cx2
       double complex zb,zb0,zb1,zb2,ctempb,czer0
       double complex hh39,hh38,hh37,hh36,hh35,hh34,hh33,hh32,hh31 
       double complex hh30,hh28,hh27,hh26,hh25,hh24,hh23,hh22,hh21 
       double complex hh20,hh17,hh16,hh15,hh14,hh13,hh12,hh11,hh10 
       double complex hh06,hh05,hh04,hh03,hh02,hh01,hh00 
       double complex ww,ww1,ww2,ww3,wwbar,ww1bar,ww2bar,ww3bar
       double complex wwaa,wwbb,wwcc,wwdd,cc1,cc2,cc3,cc4
       double complex czkcc1,czkcc2,czkcc3,ctemp,ctemp1,ctemp2
       double complex czkqpcc1,czkqpcc2,czkqpcc3
       double complex ctemp3
       double complex ctemp1qp,ctemp2qp,ctemp3qp,cjacqp
       common/lenscoords/xx1,xx2,xx3,yy1,yy2,yy3
       common/q_count/n_q,n_noq,noquad
       integer iimage(11)
       logical polish,nopolish

       save

c  yes, polish the roots
c  ---------------------
       polish = .true.
       nopolish = .false.

ccccc  CHANGED FROM HERE ON   ccccccccccccc 
ccccc  lens parameters  ccccccccccccccccccc
       eps2=eps2_in
       eps3=eps3_in
       sep=sep_in
       sep2 = sep2_in
       ang =angle

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  CENTER OF MASS SYSTEM           ccccccccccccccccc
cccccc  epsilon_1 is on the real axis at xx1    ccccccccc 
cccccc  Make sure eps4 does NOT  VANISH        cccccccccc
cccccc  because I am going to divide by it        ccccccc
cccccc  eps4 = 0 corresponds to  a single lens    ccccccc
cccccc  where  eps1 = 1.         cccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       eps4 = eps2+eps3

       if(eps4.le.0.d0) return

       czer0 = dcmplx(0.d0,0.d0)
       eps1 = 1.-eps4
       xx1 = -eps4*sep
       xx4 =  eps1*sep
       xx2 = xx4 +eps3/eps4*sep2* cos(ang)
       yy2 = eps3/eps4*sep2* sin(ang)  
       xx3 = xx4 - eps2/eps4*sep2* cos(ang) 
       yy3 = - eps2/eps4*sep2* sin(ang)

       cc1 = dcmplx(xx1,0.d0)
       cc2 = dcmplx(xx2,yy2)
       cc3 = dcmplx(xx3,yy3)
       cc4 = dcmplx(xx4,0.d0)
      
       aa = -(cc1+cc2+cc3) 
       bb = cc1*cc2 + cc1*cc3 + cc2*cc3 
       cc = -cc1*cc2*cc3
       dd = eps1*cc2*cc3 + eps2*cc1*cc3 + eps3*cc1*cc2

       hh39 = dcmplx(1.d0,0.d0)
       hh38 = 3.d0*aa
       hh37 = 3.d0*bb + 3.d0*aa*aa 
       hh36 = 3.d0*cc + 6.d0*aa*bb + aa*aa*aa
       hh35 = 6.d0*aa*cc + 3.d0*bb*bb + 3.d0*aa*aa*bb
       hh34 = 6.d0*bb*cc + 3.d0*aa*aa*cc + 3.d0*aa*bb*bb
       hh33 = 3.d0*cc*cc + 6.d0*aa*bb*cc + bb*bb*bb
       hh32 = 3.d0*aa*cc*cc + 3.d0*bb*bb*cc
       hh31 = 3.d0*bb*cc*cc
       hh30 = cc*cc*cc

       hh28 = dcmplx(1.d0,0.d0)
       hh27 = 3.d0*aa
       hh26 = dd + 2.d0*bb + 3.d0*aa*aa
       hh25 = 2.d0*aa*dd + 4.d0*aa*bb + aa*aa*aa + 2.d0*cc
       hh24 = 2.d0*dd*bb + dd*aa*aa + 4.d0*aa*cc 
     &       +2.d0*aa*aa*bb + bb*bb
       hh23 = 2.d0*dd*cc + 2.d0*dd*aa*bb + 2.d0*aa*aa*cc 
     &       +aa*bb*bb + 2.d0*bb*cc
       hh22 = 2.d0*cc*aa*dd + dd*bb*bb + 2.d0*aa*bb*cc + cc*cc
       hh21 = 2.d0*bb*cc*dd + aa*cc*cc
       hh20 = cc*cc*dd

       hh17 = dcmplx(1.d0,0.d0)
       hh16 = 3.d0*aa
       hh15 = 2.d0*dd + 3.d0*aa*aa + bb
       hh14 = 4.d0*aa*dd + aa*aa*aa + 2.d0*aa*bb + cc
       hh13 = dd*dd + 2.d0*aa*aa*dd + 2.d0*bb*dd 
     &       +bb*aa*aa + 2.d0*aa*cc
       hh12 = aa*dd*dd + 2.d0*aa*bb*dd + 2.d0*cc*dd + cc*aa*aa
       hh11 = bb*dd*dd + 2.d0*aa*cc*dd
       hh10 = cc*dd*dd

       hh06 = dcmplx(1.d0,0.d0)
       hh05 = 3.d0*aa
       hh04 = 3.d0*dd + 3.d0*aa*aa
       hh03 = 6.d0*aa*dd + aa*aa*aa
       hh02 = 3.d0*dd*dd + 3.d0*aa*aa*dd
       hh01 = 3.d0*aa*dd*dd
       hh00 = dd*dd*dd

c      quad precision version
c      ----------------------

       cc1q = cmplx(xx1,0.d0, kind=16)
       cc2q = cmplx(xx2,yy2, kind=16)
       cc3q = cmplx(xx3,yy3, kind=16)
       cc4q = cmplx(xx4,0.d0, kind=16)
      
       aaq = -(cc1q+cc2q+cc3q) 
       bbq = cc1q*cc2q + cc1q*cc3q + cc2q*cc3q 
       ccq = -cc1q*cc2q*cc3q
       ddq = eps1*cc2q*cc3q + eps2*cc1q*cc3q + eps3*cc1q*cc2q

       hh39q = dcmplx(1.d0,0.d0)
       hh38q = 3.d0*aaq
       hh37q = 3.d0*bbq + 3.d0*aaq*aaq 
       hh36q = 3.d0*ccq + 6.d0*aaq*bbq + aaq*aaq*aaq
       hh35q = 6.d0*aaq*ccq + 3.d0*bbq*bbq + 3.d0*aaq*aaq*bbq
       hh34q = 6.d0*bbq*ccq + 3.d0*aaq*aaq*ccq + 3.d0*aaq*bbq*bbq
       hh33q = 3.d0*ccq*ccq + 6.d0*aaq*bbq*ccq + bbq*bbq*bbq
       hh32q = 3.d0*aaq*ccq*ccq + 3.d0*bbq*bbq*ccq
       hh31q = 3.d0*bbq*ccq*ccq
       hh30q = ccq*ccq*ccq

       hh28q = dcmplx(1.d0,0.d0)
       hh27q = 3.d0*aaq
       hh26q = ddq + 2.d0*bbq + 3.d0*aaq*aaq
       hh25q = 2.d0*aaq*ddq + 4.d0*aaq*bbq + aaq*aaq*aaq + 2.d0*ccq
       hh24q = 2.d0*ddq*bbq + ddq*aaq*aaq + 4.d0*aaq*ccq 
     &       +2.d0*aaq*aaq*bbq + bbq*bbq
       hh23q = 2.d0*ddq*ccq + 2.d0*ddq*aaq*bbq + 2.d0*aaq*aaq*ccq 
     &       +aaq*bbq*bbq + 2.d0*bbq*ccq
       hh22q = 2.d0*ccq*aaq*ddq + ddq*bbq*bbq+2.d0*aaq*bbq*ccq+ccq*ccq
       hh21q = 2.d0*bbq*ccq*ddq + aaq*ccq*ccq
       hh20q = ccq*ccq*ddq

       hh17q = dcmplx(1.d0,0.d0)
       hh16q = 3.d0*aaq
       hh15q = 2.d0*ddq + 3.d0*aaq*aaq + bbq
       hh14q = 4.d0*aaq*ddq + aaq*aaq*aaq + 2.d0*aaq*bbq + ccq
       hh13q = ddq*ddq + 2.d0*aaq*aaq*ddq + 2.d0*bbq*ddq 
     &       +bbq*aaq*aaq + 2.d0*aaq*ccq
       hh12q = aaq*ddq*ddq + 2.d0*aaq*bbq*ddq+2.d0*ccq*ddq+ccq*aaq*aaq
       hh11q = bbq*ddq*ddq + 2.d0*aaq*ccq*ddq
       hh10q = ccq*ddq*ddq

       hh06q = dcmplx(1.d0,0.d0)
       hh05q = 3.d0*aaq
       hh04q = 3.d0*ddq + 3.d0*aaq*aaq
       hh03q = 6.d0*aaq*ddq + aaq*aaq*aaq
       hh02q = 3.d0*ddq*ddq + 3.d0*aaq*aaq*ddq
       hh01q = 3.d0*aaq*ddq*ddq
       hh00q = ddq*ddq*ddq

cccccccccccccccccccccccccccccccccccccccccccc
ccccccc  MORE CHANGES MADE (5/04/99)  cccccc
ccccccc  COEFFICIENTS BELOW  ccccccccccccccc
cccccckccccccccccccccccccccccccccccccccccccc

       tol = 1.d-4
       nimage0=3

       return

c    entry trilenso_im calculates the amplification of each image
c    ------------------------------------------------------------
      entry trilensqp_im(sx,sy,nimage,iimage,zz,amp)

      if(eps4.le.0.) then
c       single lens case
c       ----------------
        nimage = 2
        iimage(1) = 1
        u2 = sx**2+sy**2
        u = sqrt(u)
        a = (u2+2.d0)/(u*sqrt(u2+4.d0))
        amp(1) = 0.5d0*(a+1.d0)
        amp(2) = 0.5d0*(a-1.d0)
        r1fac = 0.5d0*(1.d0+sqrt(u2+4.d0)/u)
        r2fac = 0.5d0*(1.d0-sqrt(u2+4.d0)/u)
        zr1 = sx*r1fac
        zi1 = sy*r1fac
        zr2 = sx*r2fac
        zi2 = sy*r2fac
        zz(1) = dcmplx(zr1,zi1)
        zz(2) = dcmplx(zr2,zi2)
        return
      endif

c     special case (avoid singular solution)
c     --------------------------------------
      if(sep.le.1.1d-4.and.sx.eq.0..and.sy.eq.0.) sy=1.d-5

cccccccccccccccccccccccccccccccccccccccccccc
cccccccccc  MORE LINES from here ccccccccccc
cccccccccccc  COEFFICIENTS ccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccc

       ww  = dcmplx(sx,sy)
       ww1 = dcmplx(sx-xx1,sy)
       ww2 = dcmplx(sx-xx2,sy-yy2)
       ww3 = dcmplx(sx-xx3,sy-yy3)
 
       wwbar = conjg(ww)
       ww1bar = conjg(ww1)
       ww2bar = conjg(ww2)
       ww3bar = conjg(ww3)

       wwaa = ww1bar+ww2bar+ww3bar
       wwbb = ww1bar*ww2bar + ww2bar*ww3bar 
     &        + ww1bar*ww3bar
       wwcc = ww1bar*ww2bar*ww3bar
       wwdd = eps1*ww2bar*ww3bar + eps2*ww1bar*ww3bar 
     &        + eps3*ww1bar*ww2bar

c      quad precision version
       wwq  = cmplx(sx,sy, kind=16)
       ww1q = cmplx(sx-xx1,sy, kind=16)
       ww2q = cmplx(sx-xx2,sy-yy2, kind=16)
       ww3q = cmplx(sx-xx3,sy-yy3, kind=16)

       wwbarq = conjg(wwq)
       ww1barq = conjg(ww1q)
       ww2barq = conjg(ww2q)
       ww3barq = conjg(ww3q)

       wwaaq = ww1barq+ww2barq+ww3barq
       wwbbq = ww1barq*ww2barq + ww2barq*ww3barq
     &        + ww1barq*ww3barq
       wwccq = ww1barq*ww2barq*ww3barq
       wwddq = eps1*ww2barq*ww3barq + eps2*ww1barq*ww3barq
     &        + eps3*ww1barq*ww2barq

ccccccc  so far, definitions cccccccccccccc
ccccccc  Finally the coefficients ccccccccc

       cff(11) = hh39*wwcc
       cff(10) = hh38*wwcc + hh28*wwbb - (ww*wwcc+wwdd)*hh39
       cff(9)  = hh37*wwcc + hh27*wwbb + hh17*wwaa
     &     - (ww*wwcc + wwdd)*hh38 - (ww*wwbb + wwaa - wwbar)*hh28
       cff(8)  = hh36*wwcc + hh26*wwbb + hh16*wwaa + hh06 
     &     - (ww*wwcc + wwdd)*hh37 - (ww*wwbb + wwaa-wwbar)*hh27  
     &     - (ww*wwaa + 1)*hh17 
       cff(7)  = hh35*wwcc + hh25*wwbb + hh15*wwaa + hh05
     &     - (ww*wwcc + wwdd)*hh36 - (ww*wwbb + wwaa-wwbar)*hh26
     &     - (ww*wwaa + 1)*hh16  - ww*hh06 
       cff(6)  = hh34*wwcc + hh24*wwbb + hh14*wwaa + hh04
     &     - (ww*wwcc + wwdd)*hh35 - (ww*wwbb + wwaa-wwbar)*hh25
     &     - (ww*wwaa + 1)*hh15  - ww*hh05
       cff(5)  = hh33*wwcc + hh23*wwbb + hh13*wwaa + hh03
     &     - (ww*wwcc + wwdd)*hh34 - (ww*wwbb + wwaa-wwbar)*hh24
     &     - (ww*wwaa + 1)*hh14  - ww*hh04
       cff(4)  = hh32*wwcc + hh22*wwbb + hh12*wwaa + hh02
     &     - (ww*wwcc + wwdd)*hh33 - (ww*wwbb + wwaa-wwbar)*hh23
     &     - (ww*wwaa + 1)*hh13  - ww*hh03
       cff(3)  = hh31*wwcc + hh21*wwbb + hh11*wwaa + hh01
     &     - (ww*wwcc + wwdd)*hh32 - (ww*wwbb + wwaa-wwbar)*hh22
     &     - (ww*wwaa + 1)*hh12  - ww*hh02
       cff(2)  = hh30*wwcc + hh20*wwbb + hh10*wwaa + hh00
     &     - (ww*wwcc + wwdd)*hh31 - (ww*wwbb + wwaa-wwbar)*hh21
     &     - (ww*wwaa + 1)*hh11  - ww*hh01
       cff(1)  = 
     &     - (ww*wwcc + wwdd)*hh30 - (ww*wwbb + wwaa-wwbar)*hh20
     &     - (ww*wwaa + 1)*hh10  - ww*hh00

ccccccc  quad precision coefficients ccccccccc

       cffq(11) = hh39q*wwccq
       cffq(10) = hh38q*wwccq + hh28q*wwbbq - (ww*wwccq+wwddq)*hh39q
       cffq(9)  = hh37q*wwccq + hh27q*wwbbq + hh17q*wwaaq
     &     - (wwq*wwccq + wwddq)*hh38q - (wwq*wwbbq+wwaaq-wwbarq)*hh28q
       cffq(8)  = hh36q*wwccq + hh26q*wwbbq + hh16q*wwaaq + hh06q
     &     - (wwq*wwccq + wwddq)*hh37q - (wwq*wwbbq+wwaaq-wwbarq)*hh27q
     &     - (wwq*wwaaq + 1)*hh17q
       cffq(7)  = hh35q*wwccq + hh25q*wwbbq + hh15q*wwaaq + hh05q
     &     - (wwq*wwccq + wwddq)*hh36q - (wwq*wwbbq+wwaaq-wwbarq)*hh26q
     &     - (wwq*wwaaq + 1)*hh16q  - wwq*hh06q
       cffq(6)  = hh34q*wwccq + hh24q*wwbbq + hh14q*wwaaq + hh04q
     &     - (wwq*wwccq + wwddq)*hh35q - (wwq*wwbbq+wwaaq-wwbarq)*hh25q
     &     - (wwq*wwaaq + 1)*hh15q  - wwq*hh05q
       cffq(5)  = hh33q*wwccq + hh23q*wwbbq + hh13q*wwaaq + hh03q
     &     - (wwq*wwccq + wwddq)*hh34q - (wwq*wwbbq+wwaaq-wwbarq)*hh24q
     &     - (wwq*wwaaq + 1)*hh14q  - wwq*hh04q
       cffq(4)  = hh32q*wwccq + hh22q*wwbbq + hh12q*wwaaq + hh02q
     &     - (wwq*wwccq + wwddq)*hh33q - (wwq*wwbbq+wwaaq-wwbarq)*hh23q
     &     - (wwq*wwaaq + 1)*hh13q  - wwq*hh03q
       cffq(3)  = hh31q*wwccq + hh21q*wwbbq + hh11q*wwaaq + hh01q
     &     - (wwq*wwccq + wwddq)*hh32q - (wwq*wwbbq+wwaaq-wwbarq)*hh22q
     &     - (wwq*wwaaq + 1)*hh12q  - wwq*hh02q
       cffq(2)  = hh30q*wwccq + hh20q*wwbbq + hh10q*wwaaq + hh00q
     &     - (wwq*wwccq + wwddq)*hh31q - (wwq*wwbbq+wwaaq-wwbarq)*hh21q
     &     - (wwq*wwaaq + 1)*hh11q  - wwq*hh01q
       cffq(1)  = 
     &     - (wwq*wwccq + wwddq)*hh30q - (wwq*wwbbq + wwaaq-wwbar)*hh20q
     &     - (wwq*wwaaq + 1)*hh10q  - wwq*hh00q

       if(cff(11).ne.0.) then
          m=10
       elseif(cff(10).ne.0.) then
          m=9
       elseif(cff(9).ne.0.) then
          m=8
       elseif(cff(8).ne.0.) then
          m=7
       elseif(cff(7).ne.0.) then
          m=6
       elseif(cff(6).ne.0.) then
          m=5
       elseif(cff(5).ne.0.) then
          m=4
       elseif(cff(4).ne.0.) then
          m=3
       elseif(cff(3).ne.0.) then
          m=2
       elseif(cff(2).ne.0.) then
          m=1
       else
          m=0
       endif

cccccccccc  COMPLETE !!!!          cccccccccccc
ccccccccccc  I think            ccccccccccccccc
cccccc   If it doesn't work, let me know  ccccc
ccccccccccccccccccccccccccccccccccccccccccccccc

c   find the roots
c   --------------
       call zroots(cff,m,znp,nopolish)
       zqp = znp
       call zrootsqp(cffq,m,zqp,polish)
ccc       call czero(zx,cff,m,work)

c   now check for double counted roots
c   ----------------------------------
       iq_root = 0
       do ir = 1,9
         do jr = ir+1,10
           zdiff = zqp(ir)-zqp(jr)
           abszdiff = dconjg(zdiff)*zdiff
           if(abszdiff.lt.diffmin2) then
c            double counted roots!
             iq_root = 1
             go to 33
           endif
         enddo
       enddo

 33    continue
       if(iq_root.gt.0) then
         call zrootsq(cffq,m,zq,polish)
         zqp = zq
       endif

c   find the amplification
c   ----------------------
       ttol_last = -1.
       ttol=tol
       nimage=0
       ampsum=0.d0
 44    continue

ccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  MAY 12 1999  ---  NEW                    cccccccccccccccccc
cccccc  There are up to  TEN images                     ccccccccccc
cc  Lens positions cc2 and cc3 are genunie  complex numbers ccc 
cc  That is why I do not use cx1 any more  
cc  cx1 and cx2 give you an impression that they are actually real
ccccccccccccccccccccccccccccccccccccccccccccccccc

       do 50 k = 1, 10 
c         the following is the lens equation

c         quad-precision polish version
c         -----------------------------
          zkqp = zqp(k)
          czkqpcc1 = conjg(zkqp-cc1)
          czkqpcc2 = conjg(zkqp-cc2)
          czkqpcc3 = conjg(zkqp-cc3)
          ctempqp = czkqpcc1 * czkqpcc2 * czkqpcc3

          if(ctempqp.eq.czer0) then
            zqptos = zkqp - dcmplx(1.d30,0.d0)
          else
            zqptos = zkqp - eps1/czkqpcc1 - eps2/czkqpcc2
     &                    - eps3/czkqpcc3
          endif
cc           ww=complx(sx,sy) is the input source position
c         check if that (ww - zqptos) = 0
          zqperr = abs(dimag(ww-zqptos)) + abs(dreal(ww-zqptos))

c         possibly calculate the amplitude
c         --------------------------------
          if(zqperr.le.ttol.and.zqperr.gt.ttol_last) then
c             cjac = eps1/(czkcc1*czkcc1)+eps2/(czkcc2*czkcc2)
c                    + eps3/(czkcc3*czkcc3)
            ctemp1qp=(czkqpcc1*czkqpcc1)
            ctemp2qp=(czkqpcc2*czkqpcc2)
            ctemp3qp=(czkqpcc3*czkqpcc3)
            if(ctemp1qp.eq.czer0) ctemp1qp= dcmplx(1.d-30,0.d0)
            if(ctemp2qp.eq.czer0) ctemp2qp= dcmplx(1.d-30,0.d0)
            if(ctemp3qp.eq.czer0) ctemp3qp= dcmplx(1.d-30,0.d0)
            cjacqp = eps1/ctemp1qp+eps2/ctemp2qp+eps3/ctemp3qp
            fjacqp = 1.d0 - ((dreal(cjacqp))**2 + (dimag(cjacqp))**2)
            nimage=nimage+1
            iimage(nimage)=k
            amp(nimage) = abs(1.d0/(fjacqp+1.d-20))
            ampsum = ampsum + amp(nimage)
            zz(nimage) = zqp(k)
          endif
 50    continue

c      if we've missed an image, try again
cccccc   minimum number of images for ternary is 4  cccccccccc
cccccc   it can be 4, 6, 8, or 10                ccccccccccccc
c      -------------------------------------------------------
ccc       if(nimage.lt.nimage0.and.nimage.ne.4) then
cc       if(nimage.lt.nimage0) then
       if(nimage.lt.nimage0.or.ampsum.lt.0.99999d0.or.
     &   (ampsum.gt.1.06d0.and.nimage.lt.4) )         then
         if(ttol.lt.3.*tol) then
           ttol_last = ttol
           ttol=ttol*2.
ccc           write(6,*) 'WARNING:',nimage,
ccc     &                ' images; trying again w/ tol =',ttol
           go to 44
         endif
       endif
       nimage0=nimage
       ampsum0 = ampsum
       n_noq = n_noq + 1

       return
       end

c==============================================================================

       subroutine trilensdqp(eps2_in,eps3_in,sep_in,sep2_in,angle)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       parameter(diffmin2=1.e-25)

       double precision amp(10)
       complex*32 zq(10),cffq(11)
       complex*32 zqp(10)
       complex*32 ctempq,czkqcc1,czkqcc2,czkqcc3
       complex*32 czer0q,zqtos
       complex*32 ctemp1q,ctemp2q,ctemp3q,cjacq
       complex*32 aaq,bbq,ccq,ddq
       complex*32 hh39q,hh38q,hh37q,hh36q,hh35q,hh34q,hh33q,hh32q,hh31q
       complex*32 hh30q,hh28q,hh27q,hh26q,hh25q,hh24q,hh23q,hh22q,hh21q
       complex*32 hh20q,hh17q,hh16q,hh15q,hh14q,hh13q,hh12q,hh11q,hh10q
       complex*32 hh06q,hh05q,hh04q,hh03q,hh02q,hh01q,hh00q
       complex*32 wwq,ww1q,ww2q,ww3q,wwbarq,ww1barq,ww2barq,ww3barq
       complex*32 wwaaq,wwbbq,wwccq,wwddq,cc1q,cc2q,cc3q,cc4q
       complex*32 czkcc1q,czkcc2q,czkcc3q
       complex*32 czkqpcc1q,czkqpcc2q,czkqpcc3q
       complex*32 ctemp1qpq,ctemp2qpq,ctemp3qpq,cjacqpq
       real*16 eps1q,eps2q,eps3q
       double complex zz(10),z(10),cff(11),aa,bb,cc,dd
       double complex zx(10),work(100)
       double complex znp(10),zdp(10)
       double complex zk,ztos,czkaa,czkbb,cjac
       double complex zkdp,zdptos
       double complex zkqp,zqptos
       double complex zkq,zdiff
       double complex zta,zta1,zta2,zta0,cx1,cx2
       double complex zb,zb0,zb1,zb2,ctempb,czer0
       double complex hh39,hh38,hh37,hh36,hh35,hh34,hh33,hh32,hh31 
       double complex hh30,hh28,hh27,hh26,hh25,hh24,hh23,hh22,hh21 
       double complex hh20,hh17,hh16,hh15,hh14,hh13,hh12,hh11,hh10 
       double complex hh06,hh05,hh04,hh03,hh02,hh01,hh00 
       double complex ww,ww1,ww2,ww3,wwbar,ww1bar,ww2bar,ww3bar
       double complex wwaa,wwbb,wwcc,wwdd,cc1,cc2,cc3,cc4
       double complex czkcc1,czkcc2,czkcc3,ctemp,ctemp1,ctemp2
       double complex czkqpcc1,czkqpcc2,czkqpcc3
       double complex czkdpcc1,czkdpcc2,czkdpcc3
       double complex ctemp3,ctempdp
       double complex ctemp1qp,ctemp2qp,ctemp3qp,cjacqp
       double complex ctemp1dp,ctemp2dp,ctemp3dp,cjacdp
       common/lenscoords/xx1,xx2,xx3,yy1,yy2,yy3
       common/q_count/n_q,n_noq,noquad
       integer iimage(11)
       logical polish,nopolish

       save

c  yes, polish the roots
c  ---------------------
       polish = .true.
       nopolish = .false.

ccccc  CHANGED FROM HERE ON   ccccccccccccc 
ccccc  lens parameters  ccccccccccccccccccc
       eps2=eps2_in
       eps3=eps3_in
       sep=sep_in
       sep2 = sep2_in
       ang =angle

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  CENTER OF MASS SYSTEM           ccccccccccccccccc
cccccc  epsilon_1 is on the real axis at xx1    ccccccccc 
cccccc  Make sure eps4 does NOT  VANISH        cccccccccc
cccccc  because I am going to divide by it        ccccccc
cccccc  eps4 = 0 corresponds to  a single lens    ccccccc
cccccc  where  eps1 = 1.         cccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       eps4 = eps2+eps3

       if(eps4.le.0.d0) return

       czer0 = dcmplx(0.d0,0.d0)
       eps1 = 1.-eps4
       xx1 = -eps4*sep
       xx4 =  eps1*sep
       xx2 = xx4 +eps3/eps4*sep2* cos(ang)
       yy2 = eps3/eps4*sep2* sin(ang)  
       xx3 = xx4 - eps2/eps4*sep2* cos(ang) 
       yy3 = - eps2/eps4*sep2* sin(ang)

       cc1 = dcmplx(xx1,0.d0)
       cc2 = dcmplx(xx2,yy2)
       cc3 = dcmplx(xx3,yy3)
       cc4 = dcmplx(xx4,0.d0)
      
       aa = -(cc1+cc2+cc3) 
       bb = cc1*cc2 + cc1*cc3 + cc2*cc3 
       cc = -cc1*cc2*cc3
       dd = eps1*cc2*cc3 + eps2*cc1*cc3 + eps3*cc1*cc2

       hh39 = dcmplx(1.d0,0.d0)
       hh38 = 3.d0*aa
       hh37 = 3.d0*bb + 3.d0*aa*aa 
       hh36 = 3.d0*cc + 6.d0*aa*bb + aa*aa*aa
       hh35 = 6.d0*aa*cc + 3.d0*bb*bb + 3.d0*aa*aa*bb
       hh34 = 6.d0*bb*cc + 3.d0*aa*aa*cc + 3.d0*aa*bb*bb
       hh33 = 3.d0*cc*cc + 6.d0*aa*bb*cc + bb*bb*bb
       hh32 = 3.d0*aa*cc*cc + 3.d0*bb*bb*cc
       hh31 = 3.d0*bb*cc*cc
       hh30 = cc*cc*cc

       hh28 = dcmplx(1.d0,0.d0)
       hh27 = 3.d0*aa
       hh26 = dd + 2.d0*bb + 3.d0*aa*aa
       hh25 = 2.d0*aa*dd + 4.d0*aa*bb + aa*aa*aa + 2.d0*cc
       hh24 = 2.d0*dd*bb + dd*aa*aa + 4.d0*aa*cc 
     &       +2.d0*aa*aa*bb + bb*bb
       hh23 = 2.d0*dd*cc + 2.d0*dd*aa*bb + 2.d0*aa*aa*cc 
     &       +aa*bb*bb + 2.d0*bb*cc
       hh22 = 2.d0*cc*aa*dd + dd*bb*bb + 2.d0*aa*bb*cc + cc*cc
       hh21 = 2.d0*bb*cc*dd + aa*cc*cc
       hh20 = cc*cc*dd

       hh17 = dcmplx(1.d0,0.d0)
       hh16 = 3.d0*aa
       hh15 = 2.d0*dd + 3.d0*aa*aa + bb
       hh14 = 4.d0*aa*dd + aa*aa*aa + 2.d0*aa*bb + cc
       hh13 = dd*dd + 2.d0*aa*aa*dd + 2.d0*bb*dd 
     &       +bb*aa*aa + 2.d0*aa*cc
       hh12 = aa*dd*dd + 2.d0*aa*bb*dd + 2.d0*cc*dd + cc*aa*aa
       hh11 = bb*dd*dd + 2.d0*aa*cc*dd
       hh10 = cc*dd*dd

       hh06 = dcmplx(1.d0,0.d0)
       hh05 = 3.d0*aa
       hh04 = 3.d0*dd + 3.d0*aa*aa
       hh03 = 6.d0*aa*dd + aa*aa*aa
       hh02 = 3.d0*dd*dd + 3.d0*aa*aa*dd
       hh01 = 3.d0*aa*dd*dd
       hh00 = dd*dd*dd

c      quad precision version
c      ----------------------

       cc1q = cmplx(xx1,0.d0, kind=16)
       cc2q = cmplx(xx2,yy2, kind=16)
       cc3q = cmplx(xx3,yy3, kind=16)
       cc4q = cmplx(xx4,0.d0, kind=16)
      
       aaq = -(cc1q+cc2q+cc3q) 
       bbq = cc1q*cc2q + cc1q*cc3q + cc2q*cc3q 
       ccq = -cc1q*cc2q*cc3q
       ddq = eps1*cc2q*cc3q + eps2*cc1q*cc3q + eps3*cc1q*cc2q

       hh39q = dcmplx(1.d0,0.d0)
       hh38q = 3.d0*aaq
       hh37q = 3.d0*bbq + 3.d0*aaq*aaq 
       hh36q = 3.d0*ccq + 6.d0*aaq*bbq + aaq*aaq*aaq
       hh35q = 6.d0*aaq*ccq + 3.d0*bbq*bbq + 3.d0*aaq*aaq*bbq
       hh34q = 6.d0*bbq*ccq + 3.d0*aaq*aaq*ccq + 3.d0*aaq*bbq*bbq
       hh33q = 3.d0*ccq*ccq + 6.d0*aaq*bbq*ccq + bbq*bbq*bbq
       hh32q = 3.d0*aaq*ccq*ccq + 3.d0*bbq*bbq*ccq
       hh31q = 3.d0*bbq*ccq*ccq
       hh30q = ccq*ccq*ccq

       hh28q = dcmplx(1.d0,0.d0)
       hh27q = 3.d0*aaq
       hh26q = ddq + 2.d0*bbq + 3.d0*aaq*aaq
       hh25q = 2.d0*aaq*ddq + 4.d0*aaq*bbq + aaq*aaq*aaq + 2.d0*ccq
       hh24q = 2.d0*ddq*bbq + ddq*aaq*aaq + 4.d0*aaq*ccq 
     &       +2.d0*aaq*aaq*bbq + bbq*bbq
       hh23q = 2.d0*ddq*ccq + 2.d0*ddq*aaq*bbq + 2.d0*aaq*aaq*ccq 
     &       +aaq*bbq*bbq + 2.d0*bbq*ccq
       hh22q = 2.d0*ccq*aaq*ddq + ddq*bbq*bbq+2.d0*aaq*bbq*ccq+ccq*ccq
       hh21q = 2.d0*bbq*ccq*ddq + aaq*ccq*ccq
       hh20q = ccq*ccq*ddq

       hh17q = dcmplx(1.d0,0.d0)
       hh16q = 3.d0*aaq
       hh15q = 2.d0*ddq + 3.d0*aaq*aaq + bbq
       hh14q = 4.d0*aaq*ddq + aaq*aaq*aaq + 2.d0*aaq*bbq + ccq
       hh13q = ddq*ddq + 2.d0*aaq*aaq*ddq + 2.d0*bbq*ddq 
     &       +bbq*aaq*aaq + 2.d0*aaq*ccq
       hh12q = aaq*ddq*ddq + 2.d0*aaq*bbq*ddq+2.d0*ccq*ddq+ccq*aaq*aaq
       hh11q = bbq*ddq*ddq + 2.d0*aaq*ccq*ddq
       hh10q = ccq*ddq*ddq

       hh06q = dcmplx(1.d0,0.d0)
       hh05q = 3.d0*aaq
       hh04q = 3.d0*ddq + 3.d0*aaq*aaq
       hh03q = 6.d0*aaq*ddq + aaq*aaq*aaq
       hh02q = 3.d0*ddq*ddq + 3.d0*aaq*aaq*ddq
       hh01q = 3.d0*aaq*ddq*ddq
       hh00q = ddq*ddq*ddq

cccccccccccccccccccccccccccccccccccccccccccc
ccccccc  MORE CHANGES MADE (5/04/99)  cccccc
ccccccc  COEFFICIENTS BELOW  ccccccccccccccc
cccccckccccccccccccccccccccccccccccccccccccc

       tol = 1.d-4
       nimage0=3

       return

c    entry trilenso_im calculates the amplification of each image
c    ------------------------------------------------------------
      entry trilensdqp_im(sx,sy,nimage,iimage,zz,amp)

      if(eps4.le.0.) then
c       single lens case
c       ----------------
        nimage = 2
        iimage(1) = 1
        u2 = sx**2+sy**2
        u = sqrt(u)
        a = (u2+2.d0)/(u*sqrt(u2+4.d0))
        amp(1) = 0.5d0*(a+1.d0)
        amp(2) = 0.5d0*(a-1.d0)
        r1fac = 0.5d0*(1.d0+sqrt(u2+4.d0)/u)
        r2fac = 0.5d0*(1.d0-sqrt(u2+4.d0)/u)
        zr1 = sx*r1fac
        zi1 = sy*r1fac
        zr2 = sx*r2fac
        zi2 = sy*r2fac
        zz(1) = dcmplx(zr1,zi1)
        zz(2) = dcmplx(zr2,zi2)
        return
      endif

c     special case (avoid singular solution)
c     --------------------------------------
      if(sep.le.1.1d-4.and.sx.eq.0..and.sy.eq.0.) sy=1.d-5

cccccccccccccccccccccccccccccccccccccccccccc
cccccccccc  MORE LINES from here ccccccccccc
cccccccccccc  COEFFICIENTS ccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccc

       ww  = dcmplx(sx,sy)
       ww1 = dcmplx(sx-xx1,sy)
       ww2 = dcmplx(sx-xx2,sy-yy2)
       ww3 = dcmplx(sx-xx3,sy-yy3)
 
       wwbar = conjg(ww)
       ww1bar = conjg(ww1)
       ww2bar = conjg(ww2)
       ww3bar = conjg(ww3)

       wwaa = ww1bar+ww2bar+ww3bar
       wwbb = ww1bar*ww2bar + ww2bar*ww3bar 
     &        + ww1bar*ww3bar
       wwcc = ww1bar*ww2bar*ww3bar
       wwdd = eps1*ww2bar*ww3bar + eps2*ww1bar*ww3bar 
     &        + eps3*ww1bar*ww2bar

ccccccc  so far, definitions cccccccccccccc
ccccccc  Finally the coefficients ccccccccc

       cff(11) = hh39*wwcc
       cff(10) = hh38*wwcc + hh28*wwbb - (ww*wwcc+wwdd)*hh39
       cff(9)  = hh37*wwcc + hh27*wwbb + hh17*wwaa
     &     - (ww*wwcc + wwdd)*hh38 - (ww*wwbb + wwaa - wwbar)*hh28
       cff(8)  = hh36*wwcc + hh26*wwbb + hh16*wwaa + hh06 
     &     - (ww*wwcc + wwdd)*hh37 - (ww*wwbb + wwaa-wwbar)*hh27  
     &     - (ww*wwaa + 1)*hh17 
       cff(7)  = hh35*wwcc + hh25*wwbb + hh15*wwaa + hh05
     &     - (ww*wwcc + wwdd)*hh36 - (ww*wwbb + wwaa-wwbar)*hh26
     &     - (ww*wwaa + 1)*hh16  - ww*hh06 
       cff(6)  = hh34*wwcc + hh24*wwbb + hh14*wwaa + hh04
     &     - (ww*wwcc + wwdd)*hh35 - (ww*wwbb + wwaa-wwbar)*hh25
     &     - (ww*wwaa + 1)*hh15  - ww*hh05
       cff(5)  = hh33*wwcc + hh23*wwbb + hh13*wwaa + hh03
     &     - (ww*wwcc + wwdd)*hh34 - (ww*wwbb + wwaa-wwbar)*hh24
     &     - (ww*wwaa + 1)*hh14  - ww*hh04
       cff(4)  = hh32*wwcc + hh22*wwbb + hh12*wwaa + hh02
     &     - (ww*wwcc + wwdd)*hh33 - (ww*wwbb + wwaa-wwbar)*hh23
     &     - (ww*wwaa + 1)*hh13  - ww*hh03
       cff(3)  = hh31*wwcc + hh21*wwbb + hh11*wwaa + hh01
     &     - (ww*wwcc + wwdd)*hh32 - (ww*wwbb + wwaa-wwbar)*hh22
     &     - (ww*wwaa + 1)*hh12  - ww*hh02
       cff(2)  = hh30*wwcc + hh20*wwbb + hh10*wwaa + hh00
     &     - (ww*wwcc + wwdd)*hh31 - (ww*wwbb + wwaa-wwbar)*hh21
     &     - (ww*wwaa + 1)*hh11  - ww*hh01
       cff(1)  = 
     &     - (ww*wwcc + wwdd)*hh30 - (ww*wwbb + wwaa-wwbar)*hh20
     &     - (ww*wwaa + 1)*hh10  - ww*hh00


       if(cff(11).ne.0.) then
          m=10
       elseif(cff(10).ne.0.) then
          m=9
       elseif(cff(9).ne.0.) then
          m=8
       elseif(cff(8).ne.0.) then
          m=7
       elseif(cff(7).ne.0.) then
          m=6
       elseif(cff(6).ne.0.) then
          m=5
       elseif(cff(5).ne.0.) then
          m=4
       elseif(cff(4).ne.0.) then
          m=3
       elseif(cff(3).ne.0.) then
          m=2
       elseif(cff(2).ne.0.) then
          m=1
       else
          m=0
       endif


cccccccccc  COMPLETE !!!!          cccccccccccc
ccccccccccc  I think            ccccccccccccccc
cccccc   If it doesn't work, let me know  ccccc
ccccccccccccccccccccccccccccccccccccccccccccccc

c   find the roots
c   --------------
       call zroots(cff,m,znp,nopolish)
       zdp = znp
       call zrootsdp(cff,m,zdp,polish)
ccc       call zrootsqp(cffq,m,zqp,polish)
ccc       call czero(zx,cff,m,work)

c   find the amplification
c   ----------------------
       ttol_last = -1.
       ttol=tol
       nimage=0
       ampsum=0.d0
 44    continue

ccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  MAY 12 1999  ---  NEW                    cccccccccccccccccc
cccccc  There are up to  TEN images                     ccccccccccc
cc  Lens positions cc2 and cc3 are genunie  complex numbers ccc 
cc  That is why I do not use cx1 any more  
cc  cx1 and cx2 give you an impression that they are actually real
ccccccccccccccccccccccccccccccccccccccccccccccccc

       do 50 k = 1, 10 
c         the following is the lens equation

c         double-precision polish version
c         -------------------------------
          zkdp = zdp(k)
          czkdpcc1 = conjg(zkdp-cc1)
          czkdpcc2 = conjg(zkdp-cc2)
          czkdpcc3 = conjg(zkdp-cc3)
          ctempdp = czkdpcc1 * czkdpcc2 * czkdpcc3

          if(ctempdp.eq.czer0) then
            zdptos = zkdp - dcmplx(1.d30,0.d0)
          else
            zdptos = zkdp - eps1/czkdpcc1 - eps2/czkdpcc2
     &                    - eps3/czkdpcc3
          endif
cc           ww=complx(sx,sy) is the input source position
c         check if that (ww - zdptos) = 0
          zdperr = abs(dimag(ww-zdptos)) + abs(dreal(ww-zdptos))

c         possibly calculate the amplitude
c         --------------------------------
          if(zdperr.le.ttol.and.zdperr.gt.ttol_last) then
c             cjac = eps1/(czkcc1*czkcc1)+eps2/(czkcc2*czkcc2)
c                    + eps3/(czkcc3*czkcc3)
            ctemp1dp=(czkdpcc1*czkdpcc1)
            ctemp2dp=(czkdpcc2*czkdpcc2)
            ctemp3dp=(czkdpcc3*czkdpcc3)
            if(ctemp1dp.eq.czer0) ctemp1dp= dcmplx(1.d-30,0.d0)
            if(ctemp2dp.eq.czer0) ctemp2dp= dcmplx(1.d-30,0.d0)
            if(ctemp3dp.eq.czer0) ctemp3dp= dcmplx(1.d-30,0.d0)
            cjacdp = eps1/ctemp1dp+eps2/ctemp2dp+eps3/ctemp3dp
            fjacdp = 1.d0 - ((dreal(cjacdp))**2 + (dimag(cjacdp))**2)
            nimage=nimage+1
            iimage(nimage)=k
            amp(nimage) = abs(1.d0/(fjacdp+1.d-20))
            ampsum = ampsum + amp(nimage)
            zz(nimage) = zdp(k)
          endif
 50    continue

c      if we've missed an image, try again
cccccc   minimum number of images for ternary is 4  cccccccccc
cccccc   it can be 4, 6, 8, or 10                ccccccccccccc
c      -------------------------------------------------------
ccc       if(nimage.lt.nimage0.and.nimage.ne.4) then
cc       if(nimage.lt.nimage0) then
       if(nimage.lt.nimage0.or.ampsum.lt.0.99999d0.or.
     &   (ampsum.gt.1.06d0.and.nimage.lt.4) )         then
         if(ttol.lt.3.*tol) then
           ttol_last = ttol
           ttol=ttol*2.
ccc           write(6,*) 'WARNING:',nimage,
ccc     &                ' images; trying again w/ tol =',ttol
           go to 44
         elseif((ampsum.lt.0.99999d0.or.ampsum.gt.1.06d0.or.
     &           ampsum0.gt.1.06d0).and.noquad.ne.1) then
c          ccccccccccc  COEFFICIENTS ccccccccccccccccc

c          quad precision version
           wwq  = cmplx(sx,sy, kind=16)
           ww1q = cmplx(sx-xx1,sy, kind=16)
           ww2q = cmplx(sx-xx2,sy-yy2, kind=16)
           ww3q = cmplx(sx-xx3,sy-yy3, kind=16)

           wwbarq = conjg(wwq)
           ww1barq = conjg(ww1q)
           ww2barq = conjg(ww2q)
           ww3barq = conjg(ww3q)

           wwaaq = ww1barq+ww2barq+ww3barq
           wwbbq = ww1barq*ww2barq + ww2barq*ww3barq
     &           + ww1barq*ww3barq
           wwccq = ww1barq*ww2barq*ww3barq
           wwddq = eps1*ww2barq*ww3barq + eps2*ww1barq*ww3barq
     &           + eps3*ww1barq*ww2barq

c          cccccc  quad precision coefficients ccccccccc

           cffq(11) = hh39q*wwccq
           cffq(10) = hh38q*wwccq + hh28q*wwbbq - (ww*wwccq+wwddq)*hh39q
           cffq(9)  = hh37q*wwccq + hh27q*wwbbq + hh17q*wwaaq
     &     - (wwq*wwccq + wwddq)*hh38q - (wwq*wwbbq+wwaaq-wwbarq)*hh28q
           cffq(8)  = hh36q*wwccq + hh26q*wwbbq + hh16q*wwaaq + hh06q
     &     - (wwq*wwccq + wwddq)*hh37q - (wwq*wwbbq+wwaaq-wwbarq)*hh27q
     &     - (wwq*wwaaq + 1)*hh17q
           cffq(7)  = hh35q*wwccq + hh25q*wwbbq + hh15q*wwaaq + hh05q
     &     - (wwq*wwccq + wwddq)*hh36q - (wwq*wwbbq+wwaaq-wwbarq)*hh26q
     &     - (wwq*wwaaq + 1)*hh16q  - wwq*hh06q
           cffq(6)  = hh34q*wwccq + hh24q*wwbbq + hh14q*wwaaq + hh04q
     &     - (wwq*wwccq + wwddq)*hh35q - (wwq*wwbbq+wwaaq-wwbarq)*hh25q
     &     - (wwq*wwaaq + 1)*hh15q  - wwq*hh05q
           cffq(5)  = hh33q*wwccq + hh23q*wwbbq + hh13q*wwaaq + hh03q
     &     - (wwq*wwccq + wwddq)*hh34q - (wwq*wwbbq+wwaaq-wwbarq)*hh24q
     &     - (wwq*wwaaq + 1)*hh14q  - wwq*hh04q
           cffq(4)  = hh32q*wwccq + hh22q*wwbbq + hh12q*wwaaq + hh02q
     &     - (wwq*wwccq + wwddq)*hh33q - (wwq*wwbbq+wwaaq-wwbarq)*hh23q
     &     - (wwq*wwaaq + 1)*hh13q  - wwq*hh03q
           cffq(3)  = hh31q*wwccq + hh21q*wwbbq + hh11q*wwaaq + hh01q
     &     - (wwq*wwccq + wwddq)*hh32q - (wwq*wwbbq+wwaaq-wwbarq)*hh22q
     &     - (wwq*wwaaq + 1)*hh12q  - wwq*hh02q
           cffq(2)  = hh30q*wwccq + hh20q*wwbbq + hh10q*wwaaq + hh00q
     &     - (wwq*wwccq + wwddq)*hh31q - (wwq*wwbbq+wwaaq-wwbarq)*hh21q
     &     - (wwq*wwaaq + 1)*hh11q  - wwq*hh01q
           cffq(1)  = 
     &     - (wwq*wwccq + wwddq)*hh30q - (wwq*wwbbq + wwaaq-wwbar)*hh20q
     &     - (wwq*wwaaq + 1)*hh10q  - wwq*hh00q

c          find the roots
c          --------------
           zqp = znp
           call zrootsqp(cffq,m,zqp,polish)

c          now check for double counted roots
c          ----------------------------------
           iq_root = 0
           do ir = 1,9
             do jr = ir+1,10
               zdiff = zqp(ir)-zqp(jr)
               abszdiff = dconjg(zdiff)*zdiff
               if(abszdiff.lt.diffmin2) then
c                double counted roots!
                 iq_root = 1
                 go to 33
               endif
             enddo
           enddo

 33        continue
           if(iq_root.gt.0) then
             call zrootsq(cffq,m,zq,polish)
             zqp = zq
           endif

c          find the amplification
c          ----------------------
           ttol_last = -1.
           ttol=tol
           nimage=0
           ampsum=0.d0

           do 150 k = 1, 10 
c             the following is the lens equation

c             quad-precision polish version
c             -----------------------------
              zkqp = zqp(k)
              czkqpcc1 = conjg(zkqp-cc1)
              czkqpcc2 = conjg(zkqp-cc2)
              czkqpcc3 = conjg(zkqp-cc3)
              ctempqp = czkqpcc1 * czkqpcc2 * czkqpcc3

              if(ctempqp.eq.czer0) then
                zqptos = zkqp - dcmplx(1.d30,0.d0)
              else
                zqptos = zkqp - eps1/czkqpcc1 - eps2/czkqpcc2
     &                        - eps3/czkqpcc3
              endif
cc           ww=complx(sx,sy) is the input source position
c            check if that (ww - zqptos) = 0
             zqperr = abs(dimag(ww-zqptos)) + abs(dreal(ww-zqptos))

c            possibly calculate the amplitude
c            --------------------------------
             if(zqperr.le.ttol.and.zqperr.gt.ttol_last) then
c                cjac = eps1/(czkcc1*czkcc1)+eps2/(czkcc2*czkcc2)
c                       + eps3/(czkcc3*czkcc3)
               ctemp1qp=(czkqpcc1*czkqpcc1)
               ctemp2qp=(czkqpcc2*czkqpcc2)
               ctemp3qp=(czkqpcc3*czkqpcc3)
               if(ctemp1qp.eq.czer0) ctemp1qp= dcmplx(1.d-30,0.d0)
               if(ctemp2qp.eq.czer0) ctemp2qp= dcmplx(1.d-30,0.d0)
               if(ctemp3qp.eq.czer0) ctemp3qp= dcmplx(1.d-30,0.d0)
               cjacqp = eps1/ctemp1qp+eps2/ctemp2qp+eps3/ctemp3qp
               fjacqp = 1.d0 - ((dreal(cjacqp))**2 + (dimag(cjacqp))**2)
               nimage=nimage+1
               iimage(nimage)=k
               amp(nimage) = abs(1.d0/(fjacqp+1.d-20))
               ampsum = ampsum + amp(nimage)
               zz(nimage) = zqp(k)
             endif
 150       continue

         endif
       endif
       nimage0=nimage
       ampsum0 = ampsum
       n_noq = n_noq + 1

       return
       end

c==============================================================================

       subroutine trilensq(eps2_in,eps3_in,sep_in,sep2_in,angle)

c==============================================================================
       IMPLICIT REAL*16 (A-H,O-Z)

       double precision amp(10),scp(10)
       double precision eps2_in,eps3_in,sep_in,sep2_in,angle,sx_in,sy_in
       double complex zz(10)
       complex*16 zd(10),cffd(11)
       complex*32 z(10),cff(11),aa,bb,cc,dd
       complex*32 q1,aaaa,aaaaaa,aa3,aaaa3,bbbb,bb3,bbbb3,cccc,dddd
       complex*32 aabb,aacc,aadd,bbcc,bbdd,ccdd,aabbcc,aabbdd,aaccdd
       complex*32 zx(10),work(1000),zcp(10),zq(20),zc(10)
       complex*32 cffcp(20),cffq(20),cffc(11)
       complex*32 zk,ztos,czkaa,czkbb,cjac
       complex*32 zkcp,zcptos,zkq,zqtos,zkc,zctos
       complex*32 zta,zta1,zta2,zta0,cx1,cx2
       complex*32 zb,zb0,zb1,zb2,ctempb,czer0
       complex*32 hh39,hh38,hh37,hh36,hh35,hh34,hh33,hh32,hh31 
       complex*32 hh30,hh28,hh27,hh26,hh25,hh24,hh23,hh22,hh21 
       complex*32 hh20,hh17,hh16,hh15,hh14,hh13,hh12,hh11,hh10 
       complex*32 hh06,hh05,hh04,hh03,hh02,hh01,hh00 
       complex*32 ww,ww1,ww2,ww3,wwbar,ww1bar,ww2bar,ww3bar
       complex*32 wwaa,wwbb,wwcc,wwdd,cc1,cc2,cc3,cc4
       complex*32 wwcwd,wwbwawbar,wwa1
       complex*32 czkcc1,czkcc2,czkcc3,ctemp,ctemp1,ctemp2
       complex*32 czkcpcc1,czkcpcc2,czkcpcc3,ccptemp
       complex*32 czkqcc1,czkqcc2,czkqcc3,cqtemp
       complex*32 czkccc1,czkccc2,czkccc3,cctemp
       complex*32 ctemp3
ccc       common/lenscoords/xx1,xx2,xx3,yy1,yy2,yy3
       integer iimage(11)
       logical polish

       save

c  yes, polish the roots
c  ---------------------
       polish = .true.

ccccc  CHANGED FROM HERE ON   ccccccccccccc 
ccccc  lens parameters  ccccccccccccccccccc
       eps2=eps2_in
       eps3=eps3_in
       sep=sep_in
       sep2 = sep2_in
       ang =angle

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  CENTER OF MASS SYSTEM           ccccccccccccccccc
cccccc  epsilon_1 is on the real axis at xx1    ccccccccc 
cccccc  Make sure eps4 does NOT  VANISH        cccccccccc
cccccc  because I am going to divide by it        ccccccc
cccccc  eps4 = 0 corresponds to  a single lens    ccccccc
cccccc  where  eps1 = 1.         cccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       eps4 = eps2+eps3

       if(eps4.le.0.q0) return
       

       czer0 = cmplx(0.q0,0.q0, kind=16)
       eps1 = 1.q0-eps4
       xx1 = -eps4*sep
       xx4 =  eps1*sep
       xx2 = xx4 +eps3/eps4*sep2* cos(ang)
       yy2 = eps3/eps4*sep2* sin(ang)  
       xx3 = xx4 - eps2/eps4*sep2* cos(ang) 
       yy3 = - eps2/eps4*sep2* sin(ang)

       cc1 = cmplx(xx1,0.q0, kind=16)
       cc2 = cmplx(xx2,yy2, kind=16)
       cc3 = cmplx(xx3,yy3, kind=16)
       cc4 = cmplx(xx4,0.q0, kind=16)
      
       aa = -(cc1+cc2+cc3) 
       bb = cc1*cc2 + cc1*cc3 + cc2*cc3 
       cc = -cc1*cc2*cc3
       dd = eps1*cc2*cc3 + eps2*cc1*cc3 + eps3*cc1*cc2

cc       q1 = cmplx(1.q0,0.q0, kind=16)
cc       aaaa = aa*aa
cc       aaaaaa = aaaa*aa
cc       aa3 = 3.q0*aa
cc       aaaa3 = 3.q0*aaaa
cc       bbbb = bb*bb
cc       bb3 = 3.q0*bb
cc       bbbb3 = 3.q0*bbbb
cc       cccc = cc*cc
cc       dddd = dd*dd
cc       aabb = aa*bb
cc       aacc = aa*cc
cc       aadd = aa*dd
cc       bbcc = bb*cc
cc       bbdd = bb*dd
cc       ccdd = cc*dd
cc       aabbcc = aabb*cc
cc       aabbdd = aabb*dd
cc       aaccdd = aacc*dd

cc       hh39 = q1
cc       hh38 = aa3
cc       hh37 = bb3 + aaaa3 
cc       hh36 = 3.q0*cc + 6.q0*aabb + aaaaaa
cc       hh35 = 6.q0*aacc + bbbb3 + aaaa3*bb
cc       hh34 = 6.q0*bbcc + aaaa3*cc + aa3*bbbb
cc       hh33 = 3.q0*cccc + 6.q0*aabbcc + bbbb*bb
cc       hh32 = aa3*cccc + bbbb3*cc
cc       hh31 = bb3*cccc
cc       hh30 = cccc*cc

       hh39 = cmplx(1.q0,0.q0, kind=16)
       hh38 = 3.q0*aa
       hh37 = 3.q0*bb + 3.q0*aa*aa 
       hh36 = 3.q0*cc + 6.q0*aa*bb + aa*aa*aa
       hh35 = 6.q0*aa*cc + 3.q0*bb*bb + 3.q0*aa*aa*bb
       hh34 = 6.q0*bb*cc + 3.q0*aa*aa*cc + 3.q0*aa*bb*bb
       hh33 = 3.q0*cc*cc + 6.q0*aa*bb*cc + bb*bb*bb
       hh32 = 3.q0*aa*cc*cc + 3.q0*bb*bb*cc
       hh31 = 3.q0*bb*cc*cc
       hh30 = cc*cc*cc

cc       hh28 = q1
cc       hh27 = aa3
cc       hh26 = dd + 2.q0*bb + aaaa3
cc       hh25 = 2.q0*aadd + 4.q0*aabb + aaaaaa + 2.q0*cc
cc       hh24 = 2.q0*bbdd + dd*aaaa + 4.q0*aacc 
cc     &       +2.q0*aaaa*bb + bbbb
cc       hh23 = 2.q0*ccdd + 2.q0*aabbdd + 2.q0*aaaa*cc 
cc     &       +aa*bbbb + 2.q0*bbcc
cc       hh22 = 2.q0*aaccdd + dd*bbbb + 2.q0*aabbcc + cccc
cc       hh21 = 2.q0*bbcc*dd + aa*cccc
cc       hh20 = cccc*dd

       hh28 = cmplx(1.q0,0.q0, kind=16)
       hh27 = 3.q0*aa
       hh26 = dd + 2.q0*bb + 3.q0*aa*aa
       hh25 = 2.q0*aa*dd + 4.q0*aa*bb + aa*aa*aa + 2.q0*cc
       hh24 = 2.q0*dd*bb + dd*aa*aa + 4.q0*aa*cc 
     &       +2.q0*aa*aa*bb + bb*bb
       hh23 = 2.q0*dd*cc + 2.q0*dd*aa*bb + 2.q0*aa*aa*cc 
     &       +aa*bb*bb + 2.q0*bb*cc
       hh22 = 2.q0*cc*aa*dd + dd*bb*bb + 2.q0*aa*bb*cc + cc*cc
       hh21 = 2.q0*bb*cc*dd + aa*cc*cc
       hh20 = cc*cc*dd

cc       hh17 = q1
cc       hh16 = aa3
cc       hh15 = 2.q0*dd + aaaa3 + bb
cc       hh14 = 4.q0*aadd + aaaaaa + 2.q0*aabb + cc
cc       hh13 = dddd + 2.q0*aaaa*dd + 2.q0*bbdd 
cc     &       +bb*aaaa + 2.q0*aacc
cc       hh12 = aa*dddd + 2.q0*aabbdd + 2.q0*ccdd + cc*aaaa
cc       hh11 = bb*dddd + 2.q0*aaccdd
cc       hh10 = cc*dddd

       hh17 = cmplx(1.q0,0.q0, kind=16)
       hh16 = 3.q0*aa
       hh15 = 2.q0*dd + 3.q0*aa*aa + bb
       hh14 = 4.q0*aa*dd + aa*aa*aa + 2.q0*aa*bb + cc
       hh13 = dd*dd + 2.q0*aa*aa*dd + 2.q0*bb*dd 
     &       +bb*aa*aa + 2.q0*aa*cc
       hh12 = aa*dd*dd + 2.q0*aa*bb*dd + 2.q0*cc*dd + cc*aa*aa
       hh11 = bb*dd*dd + 2.q0*aa*cc*dd
       hh10 = cc*dd*dd

cc       hh06 = q1
cc       hh05 = aa3
cc       hh04 = 3.q0*dd + aaaa3
cc       hh03 = 6.q0*aadd + aaaaaa
cc       hh02 = 3.q0*dddd + aaaa3*dd
cc       hh01 = aa3*dddd
cc       hh00 = dddd*dd

       hh06 = cmplx(1.q0,0.q0, kind=16)
       hh05 = 3.q0*aa
       hh04 = 3.q0*dd + 3.q0*aa*aa
       hh03 = 6.q0*aa*dd + aa*aa*aa
       hh02 = 3.q0*dd*dd + 3.q0*aa*aa*dd
       hh01 = 3.q0*aa*dd*dd
       hh00 = dd*dd*dd

cccccccccccccccccccccccccccccccccccccccccccc
ccccccc  MORE CHANGES MADE (5/04/99)  cccccc
ccccccc  COEFFICIENTS BELOW  ccccccccccccccc
cccccckccccccccccccccccccccccccccccccccccccc

       tol = 1.q-4
       nimage0=3

       return

c    entry trilensq_im calculates the amplification of each image
c    -----------------------------------------------------------
      entry trilensq_im(sx_in,sy_in,nimage,iimage,zz,amp)

      sx=sx_in
      sy=sy_in

      if(eps4.le.0.) then
c       single lens case
c       ----------------
        nimage = 2
        iimage(1) = 1
        u2 = sx**2+sy**2
        u = sqrt(u2)
        a = (u2+2.q0)/(u*sqrt(u2+4.q0))
        amp(1) = 0.5q0*(a+1.q0)
        amp(2) = 0.5q0*(a-1.q0)
        r1fac = 0.5q0*(1.q0+sqrt(u2+4.q0)/u)
        r2fac = 0.5q0*(1.q0-sqrt(u2+4.q0)/u)
        zr1 = sx*r1fac
        zi1 = sy*r1fac
        zr2 = sx*r2fac
        zi2 = sy*r2fac
        zz(1) = dcmplx(zr1,zi1)
        zz(2) = dcmplx(zr2,zi2)
        return
      endif

c     special case (avoid singular solution)
c     --------------------------------------
      if(sep.le.1.1q-4.and.sx.eq.0..and.sy.eq.0.) sy=1.q-5

cccccccccccccccccccccccccccccccccccccccccccc
cccccccccc  MORE LINES from here ccccccccccc
cccccccccccc  COEFFICIENTS ccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccc

       ww  = cmplx(sx,sy, kind=16)
       ww1 = cmplx(sx-xx1,sy, kind=16)
       ww2 = cmplx(sx-xx2,sy-yy2, kind=16)
       ww3 = cmplx(sx-xx3,sy-yy3, kind=16)
 
       wwbar = conjg(ww)
       ww1bar = conjg(ww1)
       ww2bar = conjg(ww2)
       ww3bar = conjg(ww3)

       wwaa = ww1bar+ww2bar+ww3bar
       wwbb = ww1bar*ww2bar + ww2bar*ww3bar 
     &        + ww1bar*ww3bar
       wwcc = ww1bar*ww2bar*ww3bar
       wwdd = eps1*ww2bar*ww3bar + eps2*ww1bar*ww3bar 
     &        + eps3*ww1bar*ww2bar

       wwcwd = (ww*wwcc + wwdd)
       wwbwawbar = (ww*wwbb + wwaa-wwbar)
       wwa1 = (ww*wwaa + 1.q0)

ccccccc  so far, definitions cccccccccccccc
ccccccc  Finally the coefficients ccccccccc

       cff(11) = hh39*wwcc
       cff(10) = hh38*wwcc + hh28*wwbb - (ww*wwcc+wwdd)*hh39
       cff(9)  = hh37*wwcc + hh27*wwbb + hh17*wwaa
     &     - wwcwd*hh38 - wwbwawbar*hh28
       cff(8)  = hh36*wwcc + hh26*wwbb + hh16*wwaa + hh06 
     &     - wwcwd*hh37 - wwbwawbar*hh27  
     &     - wwa1*hh17 
       cff(7)  = hh35*wwcc + hh25*wwbb + hh15*wwaa + hh05
     &     - wwcwd*hh36 - wwbwawbar*hh26
     &     - wwa1*hh16  - ww*hh06 
       cff(6)  = hh34*wwcc + hh24*wwbb + hh14*wwaa + hh04
     &     - wwcwd*hh35 - wwbwawbar*hh25
     &     - wwa1*hh15  - ww*hh05
       cff(5)  = hh33*wwcc + hh23*wwbb + hh13*wwaa + hh03
     &     - wwcwd*hh34 - wwbwawbar*hh24
     &     - wwa1*hh14  - ww*hh04
       cff(4)  = hh32*wwcc + hh22*wwbb + hh12*wwaa + hh02
     &     - wwcwd*hh33 - wwbwawbar*hh23
     &     - wwa1*hh13  - ww*hh03
       cff(3)  = hh31*wwcc + hh21*wwbb + hh11*wwaa + hh01
     &     - wwcwd*hh32 - wwbwawbar*hh22
     &     - wwa1*hh12  - ww*hh02
       cff(2)  = hh30*wwcc + hh20*wwbb + hh10*wwaa + hh00
     &     - wwcwd*hh31 - wwbwawbar*hh21
     &     - wwa1*hh11  - ww*hh01
       cff(1)  = 
     &     - wwcwd*hh30 - wwbwawbar*hh20
     &     - wwa1*hh10  - ww*hh00


       if(cff(11).ne.0.) then
          m=10
       elseif(cff(10).ne.0.) then
          m=9
       elseif(cff(9).ne.0.) then
          m=8
       elseif(cff(8).ne.0.) then
          m=7
       elseif(cff(7).ne.0.) then
          m=6
       elseif(cff(6).ne.0.) then
          m=5
       elseif(cff(5).ne.0.) then
          m=4
       elseif(cff(4).ne.0.) then
          m=3
       elseif(cff(3).ne.0.) then
          m=2
       elseif(cff(2).ne.0.) then
          m=1
       else
          m=0
       endif


cccccccccc  COMPLETE !!!!          cccccccccccc
ccccccccccc  I think            ccccccccccccccc
cccccc   If it doesn't work, let me know  ccccc
ccccccccccccccccccccccccccccccccccccccccccccccc



c   find the roots
c   --------------
ccc       call coeff(m+1,cff,cffcp,cffc)
       call zrootsq(cff,m,z,polish)
ccc       call cpzero(m,icpflag,cffcp,zcp,scp)
ccc       call czero(zc,cffc,m,work)
ccc       call cpqr79(m,cffcp,zq,iqerr)

c   find the amplification
c   ----------------------
       ttol=tol
       ttol_last = -1.q0
       nimage=0
 44    continue

ccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  MAY 12 1999  ---  NEW                    cccccccccccccccccc
cccccc  There are up to  TEN images                     ccccccccccc
cc  Lens positions cc2 and cc3 are genunie  complex numbers ccc 
cc  That is why I do not use cx1 any more  
cc  cx1 and cx2 give you an impression that they are actually real
ccccccccccccccccccccccccccccccccccccccccccccccccc

       do 50 k = 1, 10 
          zk = z(k)
          zkcp = zcp(k)
          zkc = zc(k)
          zkq = zq(k)
c         the following is the lens equation
c         ztos = zk - conjg((zk-cc)/((zk-cc1)*(zk-cc2)))  
c              !! correct only for binary  not correct for ternary.

             czkcc1 = conjg(zk-cc1)
             czkcc2 = conjg(zk-cc2)
             czkcc3 = conjg(zk-cc3)

ccc             czkcpcc1 = conjg(zkcp-cc1)
ccc             czkcpcc2 = conjg(zkcp-cc2)
ccc             czkcpcc3 = conjg(zkcp-cc3)

ccc             czkqcc1 = conjg(zkq-cc1)
ccc             czkqcc2 = conjg(zkq-cc2)
ccc             czkqcc3 = conjg(zkq-cc3)

ccc             czkccc1 = conjg(zkc-cc1)
ccc             czkccc2 = conjg(zkc-cc2)
ccc             czkccc3 = conjg(zkc-cc3)

          ctemp= czkcc1 * czkcc2 * czkcc3
ccc          ccptemp= czkcpcc1 * czkcpcc2 * czkcpcc3
ccc          cqtemp= czkqcc1 * czkqcc2 * czkqcc3
ccc          cctemp= czkccc1 * czkccc2 * czkccc3
          if(ctemp.eq.czer0) then
            ztos = zk - cmplx(1.q30,0.q0, kind=16)
          else
            ztos = zk - eps1/czkcc1 - eps2/czkcc2 
     &                - eps3/czkcc3
          endif
ccc          if(ccptemp.eq.czer0) then
ccc            zcptos = zkcp - cmplx(1.q30,0.q0, kind=16)
ccc          else
ccc            zcptos = zkcp - eps1/czkcpcc1 - eps2/czkcpcc2 
ccc     &                - eps3/czkcpcc3
ccc          endif
ccc          if(cqtemp.eq.czer0) then
ccc            zqtos = zkq - cmplx(1.q30,0.q0, kind=16)
ccc          else
ccc            zqtos = zkq - eps1/czkqcc1 - eps2/czkqcc2 
ccc     &                - eps3/czkqcc3
ccc          endif
ccc          if(cctemp.eq.czer0) then
ccc            zctos = zkc - cmplx(1.q30,0.q0, kind=16)
ccc          else
ccc            zctos = zkc - eps1/czkccc1 - eps2/czkccc2 
ccc     &                - eps3/czkccc3
ccc          endif
cc           ww=complx(sx,sy) is the input source position
c         check if that (ww - ztos) = 0
ccc          zerr = abs(imag(ww-ztos)) + abs(real(ww-ztos))     
          zerr = abs(imag(ww-ztos)) + abs(real(ww-ztos, kind=16))
ccc          zcperr = abs(imag(ww-zcptos)) + abs(real(ww-zcptos, kind=16))
ccc          zqerr = abs(imag(ww-zqtos)) + abs(real(ww-zqtos, kind=16))
ccc          zcerr = abs(imag(ww-zctos)) + abs(real(ww-zctos, kind=16))
c         possibly calculate the amplitude
c         --------------------------------
          if(zerr.le.ttol.and.zerr.gt.ttol_last) then
             czkcc1 = conjg(zk-cc1)
             czkcc2 = conjg(zk-cc2)
             czkcc3 = conjg(zk-cc3) 
c             cjac = eps1/(czkcc1*czkcc1)+eps2/(czkcc2*czkcc2)
c                    + eps3/(czkcc3*czkcc3)
             ctemp1=(czkcc1*czkcc1)
             ctemp2=(czkcc2*czkcc2)
             ctemp3=(czkcc3*czkcc3)
             if(ctemp1.eq.czer0) ctemp1= cmplx(1.q-30,0.q0, kind=16)
             if(ctemp2.eq.czer0) ctemp2= cmplx(1.q-30,0.q0, kind=16)
             if(ctemp3.eq.czer0) ctemp3= cmplx(1.q-30,0.q0, kind=16)
             cjac = eps1/ctemp1+eps2/ctemp2+eps3/ctemp3
ccc             fjac = 1.q0 - ((real(cjac))**2 + (imag(cjac))**2)
             fjac = 1.q0 - ((real(cjac, kind=16))**2 + (aimag(cjac))**2)
             nimage=nimage+1
             iimage(nimage)=k
             amp(nimage) = abs(1.d0/(fjac*1.d0+1.d-20))
             zz(nimage) = z(k)
          endif
 50    continue

c      if we've missed an image, try again
cccccc   minimum number of images for ternary is 4  cccccccccc
cccccc   it can be 4, 6, 8, or 10                ccccccccccccc
c      -------------------------------------------------------
ccc       if(nimage.lt.nimage0.and.nimage.ne.4) then
       if(nimage.lt.nimage0) then
         if(ttol.lt.3.q0*tol) then
           ttol_last = ttol
           ttol=ttol*2.q0
ccc           write(6,*) 'WARNING:',nimage,
ccc     &                ' images; trying again w/ tol =',ttol
           go to 44
         endif
       endif
       nimage0=nimage

       return
       end

ccc       subroutine coeff(mp1,cff,cffcp,cffc)
ccc
ccc       complex*32 cff(mp1),cffcp(mp1),cffc(mp1)
ccc
ccc       do 10 i=1,mp1
ccc         cffc(i)=cff(i)/cff(mp1)
ccc         cffcp(i)=cff(mp1+1-i)
ccc 10    continue
ccc
ccc       return
ccc       end
