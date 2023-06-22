      module eesunhong_bilens
          use stdlib_kinds, only : dp, int32
          use polyroots_cmplx_roots_gen, only: cmplx_roots_gen_laguerre
          implicit none

          integer(int32), parameter :: allowed_image_number = 5
          double precision amp(allowed_image_number)
          double complex zz(allowed_image_number)
          double complex z(allowed_image_number)
          double complex cff(allowed_image_number + 1)
          double complex aa,bb,cc,alph,beta,gamm
          double complex zk,ztos,czkaa,czkbb,cjac
          double complex zta,zta1,zta2,zta0,cx1,cx2
          double complex zb,zb0,zb1,zb2,ctemp,ctempb,czero
          double complex aa2, bb2, cc2, aapcc, aaaa, bbbb, cccc
          double complex aaaapbb2
          double complex aabb2, aaccpbb, aamcc, bbcc, cff2term
          real(dp) :: xx0, xx1, xx2
          real(dp) :: eps1, eps2
          real(dp) :: sep, tol
          integer(int32) :: nimage0
          integer iimage(allowed_image_number)
          logical polish
          save

      contains
c==============================================================================

       subroutine bilens(eps1_in,sep_in)

c==============================================================================
c
c   Authors: David P. Bennett and Sun Hong Rhie
c
       use polyroots_cmplx_roots_gen, only: cmplx_roots_gen_laguerre
          implicit none
          real(dp), intent(in) :: eps1_in
          real(dp), intent(in) :: sep_in

c  yes, polish the roots
c  ---------------------
       polish = .true.

       eps1=eps1_in
       sep=sep_in

       czero = dcmplx(0.,0.)
       eps2=1.-eps1
       xx1 = -eps2*sep
       xx2 =  eps1*sep
c      xx0 = "anti" CM
       xx0 = eps1*xx2 + eps2*xx1
       cx1 = dcmplx(xx1,0.)
       cx2 = dcmplx(xx2,0.)
       aa = dcmplx(xx1+xx2,0.)
       bb = dcmplx(xx1*xx2,0.)
       cc = dcmplx(xx0,0.)
       aa2 = 2.*aa
       bb2 = 2.*bb
       cc2 = 2.*cc
       aapcc = aa+cc
       aaaa = aa*aa
       bbbb = bb*bb
       cccc = cc*cc
       aaaapbb2 = aaaa+bb2
       aabb2 = aa*bb2
       aaccpbb = aa*cc+bb
       aamcc = aa-cc
       bbcc = bb*cc
       cff2term = cc*cc-aa*cc-bb
       tol = 3.d-4
       nimage0=3

       return
       end subroutine bilens

c    entry bilens_im calculates the amplification of each image
c    ----------------------------------------------------------
      subroutine bilens_im(sx,sy,nimage,iimage_in,zz_in,amp_in)
          use stdlib_kinds, only : dp, int32
          implicit none
          real(dp), intent(in) :: sx
          real(dp), intent(inout) :: sy
          integer(int32), intent(inout) :: nimage
          integer(int32), intent(inout) ::
     &        iimage_in(allowed_image_number)
          complex(dp), intent(inout) :: zz_in(allowed_image_number)
          real(dp), intent(inout) :: amp_in(allowed_image_number)

          integer(int32) :: k, m
          real(dp) :: fjac, ttol, zerr

          iimage = iimage_in
          zz = zz_in
          amp = amp_in


c     special case (avoid singular solution)
c     --------------------------------------
      if(sep.le.1.1d-4.and.sx.eq.0..and.sy.eq.0.) sy=1.d-5

c     analytic polynomial coefficients
c     --------------------------------
       zta  = dcmplx(sx,sy) 
       zta1 = dcmplx(sx-xx1,sy)
       zta2 = dcmplx(sx-xx2,sy)
       zta0 = dcmplx(sx-xx0,sy) 
       alph = conjg(zta1) * conjg(zta2)
       beta = conjg(zta1) + conjg(zta2) 
       gamm = conjg(zta0)
            
       cff(6) = alph
       cff(5) = -alph*(aa2+zta) + beta - gamm
       cff(4) = alph*(aaaapbb2 + aa2*zta)
     &        - beta*(aapcc+zta) + gamm*aa2
       cff(3) = -alph*((aaaapbb2)*zta + aabb2)
     &                 + beta*(zta*cc+zta*aa+aaccpbb)
     &                 - gamm*(aaaapbb2) + aa-cc-zta
       cff(2) = alph*(bbbb+aabb2*zta) 
     &                 - beta*(aaccpbb*zta+bbcc)
     &                 + gamm*aabb2 + cff2term + zta*cc2
       cff(1) = -alph*bbbb*zta + beta*bbcc*zta
     &                 -gamm*bbbb + bbcc-zta*cccc

ccc       cff(5) = -alph*(2.*aa+zta) + beta - gamm
ccc       cff(4) = alph*(aa*aa + 2.*bb + 2.*aa*zta)
ccc     &        - beta*(aa+zta+cc) + gamm*2.*aa
ccc       cff(3) = -alph*((aa*aa+2.*bb)*zta + 2.*aa*bb)
ccc     &                 + beta*(bb+zta*cc+zta*aa+aa*cc)
ccc     &                 - gamm*(aa*aa+2.*bb) + aa-cc-zta
ccc       cff(2) = alph*(bb*bb+2.*aa*bb*zta) 
ccc     &                 - beta*(aa*zta*cc+bb*zta+bb*cc)
ccc     &                 + gamm*2.*aa*bb 
ccc     &                 + cc*cc-aa*cc-bb+2.*zta*cc
ccc       cff(1) = -alph*bb*bb*zta + beta*bb*zta*cc
ccc     &                 -gamm*bb*bb + cc*(bb-zta*cc) 
       if(cff(6).ne.0.) then
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
          m=4
       endif

c   find the roots
c   --------------
        call cmplx_roots_gen_laguerre(m, cff, z, polish)

c   find the amplification
c   ----------------------
       ttol=tol
 44    nimage=0
       do 50 k = 1, 5
          zk = z(k)
c         the following is the lens equation
c         ztos = zk - conjg((zk-cc)/((zk-cx1)*(zk-cx2)))
          ctemp=(zk-cx1)*(zk-cx2)
          if(ctemp.eq.czero) then
            ztos = zk - dcmplx(1.d30,0.)
          else
            ztos = zk - conjg((zk-cc)/ctemp)
          endif
c         check if that (zta - ztos) = 0
ccc          zerr = abs(imag(zta-ztos)) + abs(real(zta-ztos))     
          zerr = abs(dimag(zta-ztos)) + abs(dreal(zta-ztos))     
c         possibly calculate the amplitude
c         --------------------------------
          if(zerr.le.ttol) then
             czkaa = conjg(zk-cx1)
             czkbb = conjg(zk-cx2)
c             cjac = eps1/(czkaa*czkaa)+eps2/(czkbb*czkbb)
             ctemp=(czkaa*czkaa)
             ctempb=(czkbb*czkbb)
             if(ctemp.eq.czero) ctemp= dcmplx(1.d30,0.d0)
             if(ctempb.eq.czero) ctempb= dcmplx(1.d-30,0.d0)
             cjac = eps1/ctemp+eps2/ctempb
ccc             fjac = 1. - ((real(cjac))**2 + (imag(cjac))**2)
             fjac = 1. - ((dreal(cjac))**2 + (dimag(cjac))**2)
             nimage=nimage+1
             iimage(nimage)=k
             amp(nimage) = abs(1./(fjac+1.d-20))
             zz(nimage) = z(k)
          endif
 50    continue

c      if we've missed an image, try again
c      -----------------------------------
       if(nimage.lt.nimage0.and.nimage.ne.3) then
         if(ttol.lt.10.*tol) then
           ttol=ttol*2.
ccc           write(6,*) 'WARNING:',nimage,
ccc     &                ' images; trying again w/ tol =',ttol
           go to 44
         endif
       endif
       nimage0=nimage
       iimage_in = iimage
       zz_in = zz
       amp_in = amp

       return
       end subroutine bilens_im

      end module eesunhong_bilens
