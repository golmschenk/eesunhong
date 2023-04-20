c==============================================================================

       subroutine critical(sep,eps1,npts4,bmat,smat)

c==============================================================================
c
c   Authors: David P. Bennett and Sun Hong Rhie
c
       use polyroots_cmplx_roots_gen, only: cmplx_roots_gen_laguerre
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       parameter (pi = 3.14159265358979)
 
       double precision bmat(2,npts4), smat(2,npts4)
       double complex w, ph, cff(5), zrt(4), xx1, xx2, xx0, u0 
       double complex zk, zbk1, zbk2, zbk0, zdenom
       logical polish

       save

c  yes, polish the roots
c  ---------------------
       polish = .true.

       npts=npts4/4
ccc       ctheta=cos(theta)
ccc       stheta=sin(theta)
ccc       w=1./thats2
ccc       t00=t0-3*thats2
ccc       vfac=w*(t00-t0)
ccc       s0x=ctheta*vfac-stheta*umin
ccc       s0y=stheta*vfac+ctheta*umin
ccc       t11=t0+3*thats2
ccc       vfac=w*(t11-t0)
ccc       s1x=ctheta*vfac-stheta*umin
ccc       s1y=stheta*vfac+ctheta*umin
ccc       write(6,*) t00,s0x,s0y,t11,s1x,s1y

       eps = eps1
       eps2 = 1. - eps1
ccc       delph = 2.*pi/float(npts)
       delph = 2.*pi/dble(npts)

c   c.s. = cm  system  
c   ------------------
       xx1 = -eps2*sep
       xx2 =  eps1*sep
       xx0 = eps1*xx2 + eps2*xx1

c   c.s. = anti-cm system
c   ---------------------
c       xx1 = -eps1*sep
c       xx2 =  eps2*sep 
c       xx0 = 0. 

       u0  = eps1*eps2*sep*sep 

       do j = 1, npts
ccc         phx = cos(delph*float(j-1)) 
ccc         phy = sin(delph*float(j-1))
         phx = cos(delph*dble(j-1)) 
         phy = sin(delph*dble(j-1))
         ph  = dcmplx(phx, phy)
           
         cff(5) = 1.
         cff(4) = -2.*(xx1+xx2)
         cff(3) = (xx1+xx2)**2 + 2.*xx1*xx2 - ph 
         cff(2) = -2.*(xx1+xx2)*xx1*xx2 + 2.*ph*xx0
         cff(1) = xx1*xx1*xx2*xx2 - ph*(xx0*xx0+u0)

c        find the roots
c        --------------
         call cmplx_roots_gen_laguerre(4, cff, zrt, polish)
         do k = 1, 4
             zk = zrt(k)
             zbk1 = conjg(zk - xx1)
             zbk2 = conjg(zk - xx2)
             zbk0 = conjg(zk - xx0) 
             zdenom = zbk1*zbk2
             if(zdenom.eq.dcmplx(0.d0,0.d0)) zdenom=dcmplx(1.d-30,0.d0)
             w = zk - zbk0/zdenom
ccc             bmat(1,j+(k-1)*npts)= real(zk)
ccc             bmat(2,j+(k-1)*npts)= imag(zk)
             bmat(1,j+(k-1)*npts)= dreal(zk)
             bmat(2,j+(k-1)*npts)= dimag(zk)
ccc             smat(1,j+(k-1)*npts)= real(w)
ccc             smat(2,j+(k-1)*npts)= imag(w)
             smat(1,j+(k-1)*npts)= dreal(w)
             smat(2,j+(k-1)*npts)= dimag(w)
           enddo 
        enddo    

        return
        end

