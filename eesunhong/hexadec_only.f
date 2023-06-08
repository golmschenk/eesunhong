c==============================================================================

       subroutine hexadec(sx,sy,Ustar,nimage,z,ampsum,amp,igoodhex,
     &                    hexthresh,iclr)

c==============================================================================
c
c  Hexadecapole calculation for binary lenses - no caustic crossing checks
c  Routines in this file written by David Bennett.
c
       use stdlib_kinds, only : dp, int32
       use eesunhong_bilens, only: bilens, bilens_im
       use eesunhong_real_complex_conversion,
     &     only: from_2d_real_to_complex, from_complex_to_2d_real
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       double precision amp(5),ampq(5)
       double precision amppt(5,12),amptot(12)
       double precision sxpt(12),sypt(12)
       double precision z(2,5),zh(2,5,12)
       complex(dp) :: zh_complex(5,12)
       double precision zchk(2,5),ampchk(5)
       integer iim(11),nim(12)
       integer inew(5)
       integer im_mag(5),nim_magpt(12), im_magpt(5,12)
       parameter (mxclr=59,ndark=5000)
       common/hexlimb/qlimb(0:mxclr),hlimb(0:mxclr)
       data icall/0/

       save

       halfUstar = 0.5d0*Ustar
       hexpt_thresh = 0.125d0*sqrt(hexthresh)

       if(icall.eq.0) then
         fmag_thresh = 3.d-4
         rat_thresh = 1.d-4
         sqrthalf = sqrt(0.5d0)
         pi = acos(-1.d0)
         pis8 = 0.125d0*pi
         pis80 = 0.1d0*pis8
         cos_pis8 = cos(pis8)
         sin_pis8 = sin(pis8)
         icall = 1
       endif

       if(Ustar.ne.Ustar0) then
         Ustar2 = Ustar*Ustar
         Ustar4 = Ustar2*Ustar2
         hUstar = 0.5d0*Ustar
         rt2Ustar = sqrthalf*Ustar
         Ustar_11 = 1.15d0*Ustar
         cpis8_U11 = cos_pis8*Ustar_11
         spis8_U11 = sin_pis8*Ustar_11
         Ustar0 = Ustar
       endif

c      calculate the source positions
       if(sx0.ne.sx.and.sy0.ne.sy) then
         sxpt(1) = sx + hUstar        !  phi = 0
         sypt(1) = sy
         sxpt(2) = sx                 !  phi = pi/2
         sypt(2) = sy + hUstar
         sxpt(3) = sx - hUstar        !  phi = pi
         sypt(3) = sy
         sxpt(4) = sx                 !  phi = 3*pi/2
         sypt(4) = sy - hUstar
         sxpt(5) = sx + Ustar         !  phi = 0
         sypt(5) = sy
         sxpt(6) = sx                 !  phi = pi/2
         sypt(6) = sy + Ustar
         sxpt(7) = sx - Ustar         !  phi = pi
         sypt(7) = sy
         sxpt(8) = sx                 !  phi = 3*pi/2
         sypt(8) = sy - Ustar
         sxpt(9) = sx + rt2Ustar      !  phi = pi/4
         sypt(9) = sy + rt2Ustar
         sxpt(10)= sx - rt2Ustar      !  phi = 3*pi/4
         sypt(10)= sy + rt2Ustar
         sxpt(11)= sx + rt2Ustar      !  phi = 5*pi/4
         sypt(11)= sy - rt2Ustar
         sxpt(12)= sx - rt2Ustar      !  phi = 7*pi/4
         sypt(12)= sy - rt2Ustar
         sx0 = sx
         sy0 = sy
       endif

       nmx = 12

       do j = 1,nmx
         call from_2d_real_to_complex(zh(:, :, j), zh_complex(:, j))
         call bilens_im(sxpt(j),sypt(j),nim(j),iim,
     &                  zh_complex(:,j),amppt(1,j))
         call from_complex_to_2d_real(zh_complex(:, j), zh(:, :, j))
         amptot(j) = 0.d0
         do im = 1,nim(j)
           amptot(j) = amptot(j) + amppt(im,j)
         enddo
       enddo

c      do the quadrupole and hexadecapole sums
c      ---------------------------------------
       A_hplus = -ampsum
       A_plus = -ampsum
       A_cross = -ampsum
       do i = 1,4
         j = i + 4
         k = i + 8
         A_hplus = A_hplus + 0.25d0*amptot(i)
         A_plus = A_plus + 0.25d0*amptot(j)
         A_cross = A_cross + 0.25d0*amptot(k)
       enddo

       A_2 = (16.d0*A_hplus-A_plus)/3.d0
       A_quad = ampsum + qlimb(iclr)*A_2
       dA_hex = hlimb(iclr)*(0.5d0*(A_plus+A_cross) - A_2)
       A_hex = A_quad + dA_hex
       err_hex = abs(dA_hex/A_hex)
       err_hex_sng = abs((A_hex-ampsum)/A_hex)

c      if hexadecaopole doesn't work for all images, forget it
c      -------------------------------------------------------
       if(err_hex.gt.hexthresh.or.err_hex_sng.gt.hexpt_thresh) then
         igoodhex = 0
       else
         igoodhex = 1
       endif

c      possibly replace ampsum by hexadecapole approximation
       if(nimnew.gt.0) igoodhex = 0
       if(igoodhex.eq.1) ampsum = A_hex

       return
       end

c==============================================================================

       subroutine hexlimbfac()

c==============================================================================

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       parameter (mxclr=59,ndark=5000)
       common/limbdark/dark(ndark,0:mxclr),ald(0:mxclr),bld(0:mxclr),
     &                 darklimb(0:mxclr),hcut,bcut,iend,ildtab(0:mxclr)
       common/hexlimb/qlimb(0:mxclr),hlimb(0:mxclr)

       do i=0,mxclr
         if(ildtab(i).eq.0) then
           gam = 2.d0*ald(i)/(3.d0-ald(i))
           qlimb(i) = 0.5d0*(1.d0-0.2d0*gam)
           hlimb(i) = (1.d0-11.d0*gam/35.d0)/3.d0
         else
           qlimb(i) = 0.d0
           hlimb(i) = 0.d0
           do j = 1,ndark
             r2 = (j-0.5)/ndark
             qlimb(i) = qlimb(i) + r2*dark(j,i)
             hlimb(i) = hlimb(i) + r2*r2*dark(j,i)
           enddo
           qlimb(i) = qlimb(i)/ndark
           hlimb(i) = hlimb(i)/ndark
         endif
       enddo

       return
       end
