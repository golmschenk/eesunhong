c==============================================================================

       subroutine microcurve(t,a,yfit,iclr,nimage,alpha,delta,
     &                       Ein_R,xcc,flon_obs,flat_obs,icheck,
     &                       brgrid,bphigrid,ngr,nphimax,sxg,syg,sx,sy)

c
c   Author: David P. Bennett
c
c   routines in this file:
c       subroutine microcurve(t,a,yfit,nimage)
c       subroutine ray_shoot(zr,zphi,iclr,tol,
c     &                      eps1,sep,fampsum,famp_bndc,Ein_R,xcc,sx,sy,
c     &                      Ustar,brgrid,bphigrid,im,izr,izphi,nimage,
c     &                      included,ngr,nphimax,sxg,syg)
c       subroutine ray_shoot_row(id0,jds,id_cen,zr_cen,bphi,iclr,tol,
c     &                      eps1,sep,famp_row,fbndc_row,xcc,sx,sy,
c     &                      id_mbnd,id_pbnd,Ustar2,brgrid,nin,ngr,
c     &                      nphimax,sxg,syg)
c       subroutine check_bktrack(jds,bphi,jd0,id1,zr_cen,brgrid,
c     &                      id_cen,eps1,sep,xcc,sx,sy,Ustar2,nbpts,
c     &                      nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
c     &                      ngr,nphimax,sxg,syg)
c       subroutine limb_dark()
c       subroutine sourceloc(newcall,sep,eps1,bx,by,sx,sy)
c       function endpoint(famp0,h,br,bgridfine,jclr)
c       function endpoint2(famp1,famp15,famp2,h,br,bgridfine,jclr)
c       function endpoint3(famp1,famp2,h,br,bgridfine,jclr)
c       function crosslocr(bra,brb,cosbphi,sinbphi,
c     &                    tol,ssx,ssy,Ustar2,sep,eps1,bgridfine)
c       function crosslocphi(bphia,bphib,br,brdum,
c     &                    tol,ssx,ssy,Ustar2,sep,eps1,bgridfine)
c       function starbndr(br,cosbphi,sinbphi,ssx,ssy,Ustar2,sep,eps1)
c       function starbndphi(bphi,br,brdum,ssx,ssy,Ustar2,sep,eps1)
c      function lenc(s)
c
c       subroutine calls
c             call ray_shoot(zr(im),zphi(im),iclr,tol,eps1,sep,
c     &                      fampsum(im),famp_bndc(im),Ein_R,xcc,sx,sy,
c     &                      Ustar,brgrid,bphigrid,im,izr,izphi,nimage,
c     &                      included,ngr,nphimax,sxg,syg)
c         call bilens(eps1db,sepdb)
c         call bilens_im(sx,sy,nimage,iimage,z,ampim)
c         call critical(sep,eps1,ncpts,crmat,camat)
c         call critical(sepp,eps1,ncpts,crmats,camats)
c         call bilens_im(dsx,dsy,nimage,iimage,z,ampim)
c                 call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
c                   call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
c                   call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
c                 call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
c                 call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
c       call sourceloc(n,sep,eps1,bx,by,sx,sy)
c------------------------------------------------------------------------------
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
ccc       parameter(ncpts=16000,ndark=20000,nllmax=200,ngxmx=800,ngymx=200)
ccc       parameter(ncpts=16000,ndark=20000,nllmax=200)
       parameter(ncpts=108000,ndark=5000,nllmax=200)
ccc       parameter(ncpts=36000,ndark=5000,nllmax=200)
       parameter(ncaxgrid=700,ncaygrid=700)
       parameter(nbgridmax=6000000)
       parameter(ncrgrid=nllmax)
ccc       double precision bgridmax(2),bgridmin(2),dd(0:3)
       double precision dd(0:3)
       double precision crmat(2,ncpts), camat(2,ncpts)
       double precision crmats(2,ncpts),camats(2,ncpts),camatp(2,ncpts)
       common/critcaus/crmin(2,3),crmax(2,3),camin(2,3),camax(2,3),
     &      cacen(2,3),cadel(2),crdel(2),xcaust(3),critmin(2),
     &      critmax(2),causmin(2),causmax(2),critinmin(2),critinmax(2)
       common/integrate/gridUstar
       common/locate_src/newcall
       double precision sxg(-ngr:ngr,0:nphimax),syg(-ngr:ngr,0:nphimax)
       double precision fampsum(10),famp_bndc(10)
       double precision famp0(-nbgridmax:nbgridmax),
     &                  fampn(-nbgridmax:nbgridmax)
       double precision a(30)
       double precision z(2,10),ampim(10),zr(10),zphi(10)
       integer icrhead(nllmax,nllmax),nxtcr(ncpts)
       integer icahead(ncaxgrid,ncaygrid),nxtca(ncpts)
       double precision ca_zr(ncpts),ca_zphi(ncpts)
       integer ica_zr(ncpts),ica_zphi(ncpts),nxtcain(ncpts),
     &         ica_not_inc(ncpts),ipt_caus(ncpts)
       integer icapt_jds(-nbgridmax:nbgridmax)
       integer iimage(10),iminside(10),included(10)
       integer izr(10),izphi(10)
ccc       integer nbgrid(2),iimage(5),iminside(5)
ccc       character*80 fileb,critfile,gridfile
       parameter (max_sat=40000)
       common/satcoords/tsat(max_sat,9),rsat(max_sat,9),
     &        rasat(max_sat,9),decsat(max_sat,9),ndatsat(max_sat)

       logical polish
       external starbndr,starbndphi

       save

c  skip the calculation if t=t_old
c  -------------------------------
       if(t.eq.t_old.and.flon_obs_old.eq.flon_obs.and.
     &                   flat_obs_old.eq.flat_obs.and.
     &                   (a(9).le.0..or.famp.eq.ptamp0)) then
          nimage=nimage_old
          yfit=yfit_old
          go to 800
       endif
       t_old=t
       flon_obs_old = flon_obs
       flat_obs_old = flat_obs

c      trap for too large nphimax
c      --------------------------
       if(nphimax.gt.nbgridmax) then
         write(6,*) 'nphimax =',nphimax,' > nbgridmax =',nbgridmax
         yfit = 0.1d0
         return
       endif

c      hardwired parameters
c      --------------------
       ncptss4=ncpts/4
       pi=acos(-1.d0)
       twopi=2.d0*pi
       d_sep=0.02d0
       sep_thresh=0.08d0

       dt=0.01d0
       third=1.d0/3.d0

ccc       hexthresh = 1.d-4
       hexthresh = 2.d-4
ccc       gridUstar=0.05
ccc       dcritUmin=7.d0
ccc       dcausUmin=7.d0
       dcritUmin=20.d0
       dcausUmin=20.d0
       dgrUstar=20.d0
       dgrfnUstar=5.d0

       if(gridUstar.le.0.03d0) then
         ihires = 1
       else
         ihires = 0
       endif

c      polish the roots
c      ----------------
       polish = .true.


c    notation:
c
c    a(1) = inverse of Einstein radius crossing time = 2/t_hat.
c    a(2) = time of umin w.r.t. CM
c    a(3) = umin w.r.t. CM
c    a(4) = lens separation at time t_fix in units of Einstein radius
c    a(5) = angle between velocity vectory and separation vector at time t_fix
c    a(6) = epsilon_1 = mass fraction of mass 1
c    a(7) = 1/Tbin = inverse period of binary orbit
c    a(8) = v_sep = velocity in seperation direction.
c    a(9) = Rstar = stellar radius in days
c    a(10)= t_fix (not a fit parameter)
c    a(11)= piEr
c    a(12)= piEtheta
c    
       eps1 = a(6)
       eps2 = 1.d0-eps1
       w=twopi*a(7)
       v_sep=a(8)
       sep_in=a(4)
       theta = a(5)
       ctheta=cos(theta)
       stheta=sin(theta)
       tfix=a(10)
       piEr = a(11)
       piEtheta = a(12)
       piex = piEr*cos(piEtheta)
       piey = piEr*sin(piEtheta)
       if(v_sep.eq.0..and.tfix.eq.tfix0.and.sep_in.eq.sep_in0) then
         isep0=1
       else
         isep0=0
         sep_in0=sep_in
         tfix0=tfix
       endif

c      solve for orbital position as a function of time
c      ------------------------------------------------
       if(w.ne.0.or.v_sep.ne.0.) then
ccc       if(w.ne.0..or.gamma.ne.0.) then
         wtt = w*(t-tfix)
         sep = sep_in+v_sep*(t-tfix)
         thetaproj = wtt
cccc        don't recalculate caustic due to round-off errors!
ccc         if(isep0.eq.1) sep0=sep
       else
         sep=sep_in
         thetaproj=0.
       endif
       psi_0=psi
       sepnp_0=sepnp
       thetaorb_0=thetaorb

c      call Andy Gould's Parallax code
c      -------------------------------
       call geo_par(qn,qe,t,alpha,delta,tfix)

       if(flon_obs.eq.0.d0.and.flat_obs.eq.0.d0) then
         qtn = 0.d0
         qte = 0.d0
       elseif(flat_obs.ge.1000.d0) then
         isat = nint(flat_obs*0.001d0)
         do is = 2,ndatsat(isat)
           if(tsat(is-1,isat).le.t.and.tsat(is,isat).ge.t) then
             f0 = (tsat(is,isat)-t)/(tsat(is,isat)-tsat(is-1,isat))
             f1 = 1.d0 - f0
             rsat_t = f0*rsat(is-1,isat) + f1*rsat(is,isat)
             rasat_t = f0*rasat(is-1,isat) + f1*rasat(is,isat)
             decsat_t = f0*decsat(is-1,isat) + f1*decsat(is,isat)
             go to 5
           endif
         enddo
         if(t.le.tsat(1,isat)) then
           rsat_t = rsat(1,isat)
           rasat_t = rasat(1,isat)
           decsat_t = decsat(1,isat)
         elseif(t.ge.tsat(ndatsat(isat),isat)) then
           rsat_t = rsat(ndatsat(isat),isat)
           rasat_t = rasat(ndatsat(isat),isat)
           decsat_t = decsat(ndatsat(isat),isat)
         endif
 5       continue
         call sat_par(qtn,qte,alpha,delta,rasat_t,decsat_t,rsat_t)
       else
         call geo_tpar(qtn,qte,t,alpha,delta,flon_obs,flat_obs,icheck)
       endif

       vfac=a(1)*(t-a(2))
       dvfac = piex*(qn+qtn) + piey*(qe+qte)
       da3 =  -piex*(qe+qte) + piey*(qn+qtn)
       vfacp = vfac + dvfac
       a3p   = a(3) + da3

       sx_norot=ctheta*vfacp-stheta*a3p
       sy_norot=stheta*vfacp+ctheta*a3p

       cosfac=cos(thetaproj)
       sinfac=sin(thetaproj)
       sx=cosfac*sx_norot-sinfac*sy_norot
       sy=cosfac*sy_norot+sinfac*sx_norot

c  initialize the fit parameter variables
c  --------------------------------------
       newcall=1

c      yes, polish the roots
c      ---------------------
       polish = .true.

       eps1=min(0.9999999d0,max(1.d-7,a(6)))
       eps2=1.d0-eps1
       t0 = a(2)
       sepdb=sep
       eps1db=eps1
       eps2db=eps2

       Ustar=a(9)*a(1)
       bgriddel=gridUstar*Ustar
       that=2./a(1)
       thats2i=a(1)
       Ustar2 = Ustar**2
       delcaus=0.2*Ustar
       thats2=0.5*that

       if(sep.ne.sep0.or.eps1.ne.eps10) then
c        initialzation for lens parameters:
ccc         write(6,*) 'calling bilens for sep eps1=',sepdb,eps1db
         call bilens(eps1db,sepdb)
       endif

c      skip finite size calculation if Rstar = 0.
c      ------------------------------------------
       if(a(9).le.0.) then
c        find pt. source image locations and amplifications
c        --------------------------------------------------
         call bilens_im(sx,sy,nimage,iimage,z,ampim)

c        loop over images
c        ----------------
         nimages = nimage
         ptamp=0.
         nin=0
         do 20 im=1,nimage
           iminside(im)=0
           ptamp=ptamp+ampim(im)
 20      continue
ccc         write(4,444) t,nimage,ptamp,(ampim(im),
ccc     &                z(1,im),z(2,im), im=1,nimage)
         famp=ptamp
         ptamp0 = ptamp
         yfit=famp
         go to 800
       endif

 444   format(f11.4,i2,f9.4,10(3f9.4,1x))
    
c      determine the number of caustic regions
c      ---------------------------------------
c      c.s. = cm  system  
 21    if(abs(sep-sep0).gt.sep_thresh.or.eps1.ne.eps10.or.
     &          Ustar.ne.Ustar0.or.sep_in.ne.sep_in_0) then
         if(sep.gt.sep0) then
           sepp=sep+d_sep
         else
           sepp=sep+d_sep
         endif
         eps10=eps1
         Ustar0=Ustar
         sep0=sep
         sep_in_0=sep_in

         newcall=1

c        find the critical curves and caustics
c        -------------------------------------
c        ncpts MUST be divisible by 4
ccc         write(6,*) 'calling critical for sep eps1=',sep,eps1
         call critical(sep,eps1,ncpts,crmat,camat)
         call critical(sepp,eps1,ncpts,crmats,camats)
         dca_dsep_mx=0.
         do 25 j=1,ncptss4
           do 22 k=0,3
             camatp(1,j+k*ncptss4)=camats(1,j+k*ncptss4)
     &                             -camat(1,j+k*ncptss4)
             camatp(2,j+k*ncptss4)=camats(2,j+k*ncptss4)
     &                             -camat(2,j+k*ncptss4)
             dd(k)=camatp(1,j+k*ncptss4)**2+camatp(2,j+k*ncptss4)**2
 22        continue
           do 23 k=0,2
             do 23 kk=k+1,3
               ddcheck=(camats(1,j+kk*ncptss4)-camat(1,j+k*ncptss4))**2
     &                +(camats(2,j+kk*ncptss4)-camat(2,j+k*ncptss4))**2
               if(ddcheck.lt.dd(k)) then
c                swap solutions
                 dum1=camats(1,j+k*ncptss4)
                 dum2=camats(2,j+k*ncptss4)
                 camats(1,j+k*ncptss4)=camats(1,j+kk*ncptss4)
                 camats(2,j+k*ncptss4)=camats(2,j+kk*ncptss4)
                 camats(1,j+kk*ncptss4)=dum1
                 camats(2,j+kk*ncptss4)=dum2
                 camatp(1,j+k*ncptss4)=camats(1,j+k*ncptss4)
     &                                 -camat(1,j+k*ncptss4)
                 camatp(2,j+k*ncptss4)=camats(2,j+k*ncptss4)
     &                                 -camat(2,j+k*ncptss4)
                 dd(k)=camatp(1,j+k*ncptss4)**2+camatp(2,j+k*ncptss4)**2
                 camatp(1,j+kk*ncptss4)=camats(1,j+kk*ncptss4)
     &                                 -camat(1,j+kk*ncptss4)
                 camatp(2,j+kk*ncptss4)=camats(2,j+kk*ncptss4)
     &                                 -camat(2,j+kk*ncptss4)
                 dd(kk)=camatp(1,j+kk*ncptss4)**2
     &                 +camatp(2,j+kk*ncptss4)**2
               endif
 23        continue
           do 24 k=0,3
             camatp(1,j+k*ncptss4)=camatp(1,j+k*ncptss4)/d_sep
             camatp(2,j+k*ncptss4)=camatp(2,j+k*ncptss4)/d_sep
             dcheck=sqrt(camatp(1,j+k*ncptss4)**2
     &                  +camatp(2,j+k*ncptss4)**2)
             if(dcheck.gt.dca_dsep_mx) then
               dca_dsep_mx=dcheck
             endif
 24        continue
 25      continue

cccc        determine the number of caustic regions
cccc        ---------------------------------------
cccc        c.s. = cm  system  
ccc         xx1 = -eps2db*sep
ccc         xx2 =  eps1db*sep
ccc         xx0 = eps1db*xx2 + eps2db*xx1

c        for general case try just 1 caustic grid (even if there's > 1 caustic)
c        ------------------------------------------
cc         fls=(eps1db**third+eps2db**third)**1.5
cc         flm=1./sqrt(fls)
cc         if(sep.lt.fls) then
cc           bx=xx0
cc           by=sqrt(eps1*eps2)*sep
cc           call sourceloc(newcall,sep,eps1,bx,by,sxc,syc)
cc           xcaust(2)=sxc
cc           ncausgrid=2
cc         elseif(sep.gt.flm) then
cc           y2=0.5*(sqrt((sep**2-1.)**2+4.*eps1db*sep**2)+1.-sep**2)
cc           xcaust(2)=xx2+eps1db*(xx1-xx2)/y2
cc           ncausgrid=2
cc         else
           ncausgrid=1
cc         endif
cc         xcaust(1)=0.

c        find range of critical curve and caustic positions
c        --------------------------------------------------
cc         do 30 ic=1,2
           do 30 j=1,2
cc             crmax(j,ic)=-1.d9
cc             crmin(j,ic)=1.d9
cc             camax(j,ic)=-1.d9
cc             camin(j,ic)=1.d9
             critmax(j)=-1.d9
             critmin(j)=1.d9
             causmax(j)=-1.d9
             causmin(j)=1.d9
 30      continue
         do 50 i=1,ncpts
cc           if(abs(camat(1,i)-xcaust(1)).le.
cc       &      abs(camat(1,i)-xcaust(2))) then
             ic=1
cc           else
cc             ic=2
cc           endif
cc           if(camat(1,i).gt.camax(1,ic)) camax(1,ic)=camat(1,i)
cc           if(camat(1,i).lt.camin(1,ic)) camin(1,ic)=camat(1,i)
cc           if(abs(camat(2,i)).gt.camax(2,ic)) camax(2,ic)=abs(camat(2,i))
cc           if(abs(camat(2,i)).lt.camin(2,ic)) camin(2,ic)=abs(camat(2,i))
           if(camat(1,i).gt.causmax(1)) causmax(1)=camat(1,i)
           if(camat(1,i).lt.causmin(1)) causmin(1)=camat(1,i)
           if(abs(camat(2,i)).gt.causmax(2)) causmax(2)=abs(camat(2,i))
           if(abs(camat(2,i)).lt.causmin(2)) causmin(2)=abs(camat(2,i))
           do 50 j=1,2
cc             if(crmat(j,i).gt.crmax(j,ic)) crmax(j,ic)=crmat(j,i)
cc             if(crmat(j,i).lt.crmin(j,ic)) crmin(j,ic)=crmat(j,i)
             if(crmat(j,i).gt.critmax(j)) critmax(j)=crmat(j,i)
             if(crmat(j,i).lt.critmin(j)) critmin(j)=crmat(j,i)
 50      continue

         do 60 j=1,2
cc           critmax(j)=max(crmax(j,1),crmax(j,2))
cc           critmin(j)=min(crmin(j,1),crmin(j,2))
           crdel(j)=(critmax(j)-critmin(j))/ncrgrid
cc           causmax(j)=max(camax(j,1),camax(j,2))
cc           causmin(j)=min(camin(j,1),camin(j,2))
 60      continue
         cadel(1)=(causmax(1)-causmin(1))/ncaxgrid
         cadel(2)=(2.*causmax(2))/ncaygrid

c        determine the dimensions of the source grid(s)
c        ----------------------------------------------
         fpad=dgrUstar*Ustar
         ffpad=dgrfnUstar*Ustar
         if(ncausgrid.eq.1) then
         elseif(camax(1,2).lt.camin(1,2).or.
     &          camax(2,2).lt.camin(2,2)) then
           ncausgrid=1
         elseif(camin(1,1).gt.camax(1,2)+2.*fpad) then
           ncausgrid=2
         elseif(camin(1,2).gt.camax(1,1)+2.*fpad) then
           ncausgrid=2
         else
           ncausgrid=1
         endif
         if(ncausgrid.eq.1) then
           camin(1,1)=min(camin(1,1),camin(1,2))
           camax(1,1)=max(camax(1,1),camax(1,2))
           camin(2,1)=min(camin(2,1),camin(2,2))
           camax(2,1)=max(camax(2,1),camax(2,2))
         endif
         do 70 ic=1,ncausgrid
           do 68 j=1,2
             cacen(j,ic)=0.5*(camin(j,ic)+camax(j,ic))
             nnn=nint(cacen(j,ic)/delcaus)
             cacen(j,ic)=nnn*delcaus
             range=0.5*(camax(j,ic)-camin(j,ic))+fpad
             rangefn=0.5*(camax(j,ic)-camin(j,ic))+ffpad
ccc             nca(j,ic)=2*nint(range/(2.*delcaus))
ccc             fac=nca(j,ic)*delcaus
             camin(j,ic)=cacen(j,ic)-rangefn
             camax(j,ic)=cacen(j,ic)+rangefn
 68        continue
ccc           if(nca(1,ic).gt.ngxmx) then
ccc             write(6,*) 'NCA > ngxmx at NCA, ic =',nca(1,ic),ic
ccc             STOP
ccc           endif
ccc           if(nca(2,ic).gt.ngymx) then
ccc             write(6,*) 'NCA > ngymx at NCA, ic =',nca(2,ic),ic
ccc             STOP
ccc           endif
 70      continue

c        zero the icahead grid
c        ---------------------
         do 80 icay=1,ncaygrid
           do 80 icax=1,ncaxgrid
             icahead(icax,icay)=0
 80      continue

c        zero the icrhead grid
c        ---------------------
         do 90 icry=1,ncrgrid
           do 90 icrx=1,ncrgrid
             icrhead(icrx,icry)=0
 90      continue

c        setup the critical curve/caustic linked lists
c        ---------------------------------------------
         r2check=1.d-8*Ustar2
         do 100 i=1,ncpts
           icrx=min(ncrgrid,1+int((crmat(1,i)-critmin(1))/crdel(1)) )
           icry=min(ncrgrid,1+int((crmat(2,i)-critmin(2))/crdel(2)) )
           icax=min(ncaxgrid,1+int((camat(1,i)-causmin(1))/cadel(1)) )
           icay=min(ncaygrid,1+int((camat(2,i)+causmax(2))/cadel(2)) )
           nxtcr(i)=icrhead(icrx,icry)
           icrhead(icrx,icry)=i
           if(icahead(icax,icay).eq.0) then
             nxtca(i)=icahead(icax,icay)
             icahead(icax,icay)=i
             cagridx=camat(1,i)
             cagridy=camat(2,i)
           else
             r2=(camat(1,i)-cagridx)**2+(camat(2,i)-cagridy)**2
             if(r2.ge.r2check) then
               nxtca(i)=icahead(icax,icay)
               icahead(icax,icay)=i
               cagridx=camat(1,i)
               cagridy=camat(2,i)
             endif
           endif
 100     continue

       endif

c==============================================================================

ccc       t00=t0-3*thats2
ccc       vf=thats2i*(t00-t0)
ccc       s0x=ctheta*vf-stheta*umin
ccc       s0y=stheta*vf+ctheta*umin
ccc       t11=t0+3*thats2
ccc       vf=thats2i*(t11-t0)
ccc       s1x=ctheta*vf-stheta*umin
ccc       s1y=stheta*vf+ctheta*umin

ccc       open(unit=11,file=critfile,status='unknown')
ccc 20      format(4(f13.5))
ccc 21      format(10(f13.5))
ccc       write(11,21) crmat(1,1),crmat(2,1),camat(1,1),camat(2,1),
ccc   &                t00,s0x,s0y,t11,s1x,s1y
ccc       write(11,20)(crmat(1,j),crmat(2,j),camat(1,j),camat(2,j),
ccc   &                j=2,ncpts)
ccc       close(11)

ccc      write(6, *) 'eps1, sep, Ustar, delUstar, ncausgrid:'
ccc      write(6,*) eps1,sep,Ustar,delUstar,ncausgrid
      ngymin=0
      inone=1

c       now, compute the magnification
c       ------------------------------
         dsy=sy
         dsx=sx

c        find pt. source image locations and amplifications
c        --------------------------------------------------
ccc         write(6,*) 'calling bilens_im for dsx dsy=',dsx,dsy
         call bilens_im(dsx,dsy,nimage,iimage,z,ampim)
c        loop over images
c        ----------------
         ptamp=0.
         do 110 im=1,nimage
           ptamp=ptamp+ampim(im)
c          determine the image locations in offset polar coords
c          ----------------------------------------------------
           z1mxcc = z(1,im)-xcc
           zr(im) = sqrt(z1mxcc**2+z(2,im)**2)
           if(z1mxcc.gt.0.d0) then
             zphi(im) = atan(z(2,im)/z1mxcc)
           elseif(z1mxcc.lt.0.d0) then
             zphi(im) = atan(z(2,im)/z1mxcc) + pi
           elseif(z(2,im).ge.0.d0) then
             zphi(im) = 0.5d0*pi
           else
             zphi(im) = -0.5d0*pi
           endif
           if(zphi(im).gt.pi) zphi(im)=zphi(im)-twopi
           if(zphi(im).lt.-pi) zphi(im)=zphi(im)+twopi
           izphi(im) =nint(zphi(im)/bphigrid)
           if(zphi(im).gt.nphimax) then
             zphi(im)=zphi(im)-2*nphimax-1
           elseif(zphi(im).lt.-nphimax) then
             zphi(im)=zphi(im)+2*nphimax+1
           endif
           izr(im)=nint((zr(im)-Ein_R)/brgrid)
 110     continue

         ptamp0 = ptamp

c        determine if seperation from caustics is large
c        enough to skip finite source calculation
c        ----------------------------------------------
         none_inside=1
ccc         dcaus=dcausUmin*sqrt(ptamp)*Ustar
         dcaus=dcausUmin*ptamp*Ustar
         rchk=dcaus+dca_dsep_mx*abs(sep-sep0)
         drchk2=(dcaus+abs(sep-sep0))**2
         if(sx-rchk.le.causmax(1).and.sx+rchk.ge.causmin(1).and.
     &      sy-rchk.le.causmax(2).and.sy+rchk.ge.-causmax(2)) then
           icaxmin=1+int((sx-rchk-causmin(1))/cadel(1))
           icaxmax=1+int((sx+rchk-causmin(1))/cadel(1))
           icaymin=1+int((sy-rchk+causmax(2))/cadel(2))
           icaymax=1+int((sy+rchk+causmax(2))/cadel(2))

           icaxmin=max(1,icaxmin)
           icaxmax=min(ncaxgrid,icaxmax)
           icaymin=max(1,icaymin)
           icaymax=min(ncaygrid,icaymax)
           dsep=sep-sep0
           do 120 icay=icaymin,icaymax
             do 120 icax=icaxmin,icaxmax
               i=icahead(icax,icay)
 115           if(i.ne.0) then
                 cax=camat(1,i)+dsep*camatp(1,i)
                 cay=camat(2,i)+dsep*camatp(2,i)
                 r2=(cax-sx)**2+(cay-sy)**2
                 if(r2.le.drchk2) then
                   none_inside=0
                   go to 122
                 endif
                 i=nxtca(i)
                 go to 115
               endif
 120       continue
         endif
 122     if(none_inside.eq.1) then
           famp=ptamp
           ptamp0 = ptamp
           yfit=famp
           go to 800
         elseif(sep.ne.sep0) then
           sep0=sep-1.
           go to 21
         endif

c        loop over images
c        ----------------
         ptamp=0.
         nin=0
         do 180 im=1,nimage
           iminside(im)=0
           ptamp=ptamp+ampim(im)
           icrx=min(ncrgrid,1+int((z(1,im)-critmin(1))/crdel(1)) )
           icry=min(ncrgrid,1+int((z(2,im)-critmin(2))/crdel(2)) )
           if(ihires.eq.1) then
             crdthmin=Ustar*(max(ampim(im),1.d0))
           else
             crdthmin=Ustar*sqrt(ampim(im))
           endif
           crdthresh=dcritUmin*crdthmin
           crdthresh2=crdthresh**2
           CAdthresh=dcausUmin*crdthmin
           CAdthresh2=CAdthresh**2
           ndcr=(2.d0*crdthresh)/min(crdel(1),crdel(2))+1
           icy0=max(icry-ndcr,1)
           icy1=min(icry+ndcr,ncrgrid)
           icx0=max(icrx-ndcr,1)
           icx1=min(icrx+ndcr,ncrgrid)
           if(icy1.ge.1.and.icy0.le.ncrgrid.and.
     &        icx1.ge.1.and.icx0.le.ncrgrid.and.
     &                  crdthmin.gt.1.5*bgriddel) then
             do 140 icy=icy0,icy1
               do 140 icx=icx0,icx1
                 i=icrhead(icx,icy)
 135             if(i.gt.0) then
                   if(i.gt.ncpts) then
                     write(6,*) 'ERROR ',
     &                 'icrhead, nxtcr out of range at i,icx,icy=',
     &                 i,icx,icy
                     write(6,*) 'skipping to next grid cell'
                   else
                     d2=(crmat(1,i)-z(1,im))**2+(crmat(2,i)-z(2,im))**2
                     if(d2.le.crdthresh2) then
                       iminside(im)=1
                       go to 141
                     endif
                     dca2=(camat(1,i)-sx)**2+(camat(2,i)-sy)**2
                     if(dca2.le.CAdthresh2) then
                       iminside(im)=1
                       go to 141
                     endif
                     i=nxtcr(i)
                     go to 135
                   endif
                 endif
 140         continue
 141         continue
           endif
 146       continue
           
 180     continue

c        determine if there are any caustics inside the source
c        -----------------------------------------------------
         ncausin=0
         if(ncausin_save.gt.0) then
c          zero the linked list array
           do in = 1, ncausin_save
             nxtcain(in) = 0
           enddo
         endif
c        and the linked list header array
         do jds = -nphimax,nphimax
           icapt_jds(jds) = 0
         enddo
         critinmin(1)=1.d9
         critinmin(2)=1.d9
         critinmax(1)=-1.d9
         critinmax(2)=-1.d9
         if(nimage.lt.7) then
           if(sx-Ustar.le.causmax(1).and.sx+Ustar.ge.causmin(1).and.
     &        sy-Ustar.le.causmax(2).and.sy+Ustar.ge.-causmax(2)) then
             icaxmin=1+int((sx-Ustar-causmin(1))/cadel(1))
             icaxmax=1+int((sx+Ustar-causmin(1))/cadel(1))
             icaymin=1+int((sy-Ustar+causmax(2))/cadel(2))
             icaymax=1+int((sy+Ustar+causmax(2))/cadel(2))

             icaxmin=max(1,icaxmin)
             icaxmax=min(ncaxgrid,icaxmax)
             icaymin=max(1,icaymin)
             icaymax=min(ncaygrid,icaymax)
             Us2_m_r2_max = -1.d0
             do 200 icay=icaymin,icaymax
               do 200 icax=icaxmin,icaxmax
                 i=icahead(icax,icay)
 195             if(i.ne.0) then
                   r2=(camat(1,i)-sx)**2+(camat(2,i)-sy)**2
                   Us2_m_r2 = Ustar2 - r2
                   if(Us2_m_r2.ge.0.d0) then
c                    convert to offset polar coords!
c                    -------------------------------
                     crmat1mxcc= crmat(1,i) - xcc
                     crr = sqrt(crmat1mxcc**2+crmat(2,i)**2)
                     if(crmat1mxcc.gt.0.d0) then
                       crphi = atan(crmat(2,i)/crmat1mxcc)
                     elseif(crmat1mxcc.lt.0.d0) then
                       crphi = atan(crmat(2,i)/crmat1mxcc) + pi
                     elseif(crmat(2,i).ge.0.d0) then
                       crphi = 0.5d0*pi
                     else
                       crphi = -0.5d0*pi
                     endif
c                    critinmin and critinmax are in polar coords (1,2)=(r,phi)
c                    ---------------------------------------------------------
                     if(crr.lt.critinmin(1))
     &                  critinmin(1)=crr
                     if(crr.gt.critinmax(1))
     &                  critinmax(1)=crr
                     if(crphi.lt.critinmin(2))
     &                  critinmin(2)=crphi
                     if(crphi.gt.critinmax(2))
     &                  critinmax(2)=crphi
c                    determine the integer coordinates
c                    ---------------------------------
                     ica_zr_tmp = nint((crr-Ein_R)/brgrid)
                     ica_zphi_tmp = nint(crphi/bphigrid)
                     if(ica_zphi_tmp.gt.nphimax) then
                       ica_zphi_tmp=ica_zphi_tmp-2*nphimax-1
                     elseif(ica_zphi_tmp.lt.-nphimax) then
                       ica_zphi_tmp=ica_zphi_tmp+2*nphimax+1
                     endif
c                    check to see if this point is already included
c                    ----------------------------------------------
                     new = 1
                     call check_caust_list(new,ica_zr_tmp,ica_zphi_tmp,
     &                          ncausin,nxtcain,ica_zr,ica_zphi,
     &                          nbgridmax,icapt_jds)
                     if(new.eq.1) then
c                      check to see if the integer point is included
c                      ---------------------------------------------
                       ica_zr_tmp2 = ica_zr_tmp
                       ica_zphi_tmp2 = ica_zphi_tmp
                       call cen_in_image(crr,crphi,ica_zr_tmp,
     &                         ica_zphi_tmp,sx,sy,Ustar2,brgrid,
     &                         bphigrid,eps1,sep,Ein_R,xcc,nphimax,ipt)
                       if(ica_zphi_tmp2.ne.ica_zphi_tmp.or.
     &                    ica_zr_tmp2.ne.ica_zr_tmp)        then
                         call check_caust_list(new,ica_zr_tmp,
     &                          ica_zphi_tmp,ncausin,nxtcain,
     &                          ica_zr,ica_zphi,nbgridmax,icapt_jds)
                       endif
                       if(new.eq.1) then
                         ncausin=ncausin+1
                         ipt_caus(ncausin) = ipt
                         ca_zr(ncausin) = crr
                         ca_zphi(ncausin) = crphi
                         ica_zr(ncausin) = ica_zr_tmp
                         ica_zphi(ncausin) = ica_zphi_tmp
c                        save multiple points in the same row using linked list
c                        ------------------------------------------------------
                         nxtcain(ncausin) = icapt_jds(ica_zphi(ncausin))
                         icapt_jds(ica_zphi(ncausin)) = ncausin
                         ica_not_inc(ncausin) = 1
                       endif
                     endif
                   endif
                   i=nxtca(i)
                   go to 195
                 endif
 200         continue
           endif
         endif

         ncausin_save = ncausin
         ptamp0 = ptamp

c        calculate hexadecapole approximation if there is no caustic crossing
c        --------------------------------------------------------------------
         if(ncausin.eq.0) then
           call hexadec(sx,sy,Ustar,nimage,z,ptamp,ampim,igoodhex,
     &                  hexthresh,iclr)
           if(igoodhex.eq.1.and.nimcaus.lt.1) then
             famp=ptamp
             yfit=famp
             go to 800
           endif
         endif

c        ensure that center points are in the images
c        -------------------------------------------
         call cen_in_images(nimage,zr,zphi,izr,izphi,iminside,
     &             sx,sy,Ustar,brgrid,bphigrid,eps1,sep,
     &             Ein_R,xcc,nphimax)

c        sum over the included image points
c        ----------------------------------
         tol=0.1d0*bgriddel*gridUstar
         do im = 1,5
           included(im) = 0
         enddo
         do 600 im=1,nimage
           fampsum(im)=0.
           famp_bndc(im)=0.
           if(iminside(im).eq.1.and.included(im).eq.0) then
             imi=iimage(im)

c            calculate image brightness by ray shooting
c            ------------------------------------------
             call ray_shoot(zr(im),zphi(im),iclr,tol,eps1,sep,
     &                      fampsum(im),famp_bndc(im),Ein_R,xcc,sx,sy,
     &                      Ustar,brgrid,bphigrid,im,izr,izphi,nimage,
     &                      included,iminside,
     &                      ngr,nphimax,sxg,syg,ncausin,ica_zr,
     &                      ica_zphi,ica_not_inc,nxtcain,icapt_jds,0)

           endif
 600     continue

c        now add a possible caustic crossing image
c        -----------------------------------------
         ncausin2 = 0
         icausin = 0
         igood_pt = 0
         nimcaus = 0
         nimages = nimage
         ds2_min = Ustar2
         if(ncausin.gt.0) then
           do in = 1,ncausin
             if(ica_not_inc(in).gt.0) then
               ncausin2 = ncausin2 + 1
ccc               if(icausin.eq.0) then
ccc                 icausin = in
ccc               elseif(igood_pt.eq.0) then
                 zrc = (ica_zr(in)*brgrid) + Ein_R
                 zphic = bphigrid*ica_zphi(in)
                 bx = zrc*cos(zphic) + xcc
                 by = zrc*sin(zphic)
                 call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
                 ds2=((sx-sxe)**2+(sy-sye)**2)
                 if(ds2.lt.ds2_min) then
                   icausin = in
                   igood_pt = 1
                   ds2_min = ds2
                 endif
ccc               endif
             endif
           enddo
         endif
         if(ncausin2.eq.0) nimcaus = 0
ccc         iextra_caus = 0
ccc         do while (ncausin2.gt.0)
         do
           if(ncausin2.le.0) exit
ccc           iextra_caus = iextra_caus + 1
           nimcaus = nimcaus + 1
           im = nimage + nimcaus
           nimages = nimage + nimcaus
           if(nimages.gt.5) then
             write(6,*) nimages,' images at t =',t
             if(nimages.gt.10) then
               nimages = 10
               exit
             endif
           endif
           zr(im) = ca_zr(icausin)
           zphi(im) = ca_zphi(icausin)
           if(zphi(im).gt.pi) zphi(im)=zphi(im)-twopi
           if(zphi(im).lt.-pi) zphi(im)=zphi(im)+twopi
           z(1,im)=zr(im)*cos(zphi(im)) + xcc
           z(2,im)=zr(im)*sin(zphi(im))
           iminside(im)=1
           ampim(im)=1.d-10
           izr(im) = ica_zr(icausin)
           izphi(im) = ica_zphi(icausin)
           if(izphi(im).gt.nphimax) then
             izphi(im)=izphi(im)-2*nphimax-1
           elseif(zphi(im).lt.-nphimax) then
             izphi(im)=izphi(im)+2*nphimax+1
           endif
c          calculate image brightness by ray shooting
c          ------------------------------------------
           if(ipt_caus(icausin).eq.0) then
c            normal caustic
             call ray_shoot(zr(im),zphi(im),iclr,tol,eps1,sep,
     &                  fampsum(im),famp_bndc(im),Ein_R,xcc,sx,sy,
     &                  Ustar,brgrid,bphigrid,im,izr,izphi,nimages,
     &                  included,iminside,
     &                  ngr,nphimax,sxg,syg,ncausin,ica_zr,
     &                  ica_zphi,ica_not_inc,nxtcain,icapt_jds,0)
           else
c            single grid-point caustic
             call ray_shoot(zr(im),zphi(im),iclr,tol,eps1,sep,
     &                  fampsum(im),famp_bndc(im),Ein_R,xcc,sx,sy,
     &                  Ustar,brgrid,bphigrid,im,izr,izphi,nimages,
     &                  included,iminside,
     &                  0,nphimax,sxg,syg,ncausin,ica_zr,
     &                  ica_zphi,ica_not_inc,nxtcain,icapt_jds,1)
           endif

c          there should be no caustic points left now, but there could be?
c          ------------------------------------------
           ncausin2 = 0
           ncausin3 = 0
ccc           if(iextra_caus.gt.1) then
ccc             write(6,*)
ccc     &         'iextra_caus fampsum(im) =',iextra_caus,fampsum(im)
ccc           endif
           ds2_min = Ustar2
           igood_pt = 0
           do in = 1,ncausin
             if(ica_not_inc(in).gt.0) then
               ncausin3 = ncausin3 + 1
               zrc = (ica_zr(in)*brgrid) + Ein_R
               zphic = bphigrid*ica_zphi(in)
               bx = zrc*cos(zphic) + xcc
               by = zrc*sin(zphic)
               call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
               ds2=((sx-sxe)**2+(sy-sye)**2)
               if(ds2.lt.ds2_min) then
                 icausin = in
                 igood_pt = 1
                 ds2_min = ds2
               endif
             endif
           enddo
ccc           if(iextra_caus.gt.3.or.igood_pt.le.0) ncausin3 = 0
           if(igood_pt.le.0) ncausin3 = 0
           ncausin2 = ncausin3
         enddo

         famptot=0.
         do 635 im=1,nimages
           if(included(im).eq.0)
     &                famptot=famptot+fampsum(im)+famp_bndc(im)
 635     continue
         famp=famptot

c        add the "point" images
c        ----------------------
         do 650 iptim=1,nimage
           if(iminside(iptim).eq.0) famp=famp+ampim(iptim)
 650     continue

       yfit=famp

 800   continue
       yfit_old=yfit
       nimage_old=nimage

       return
       end

c==============================================================================

       subroutine microcurve_init(a,Ein_R,xcc,brgrid,bphigrid,
     &                            ngr,ngphi,sxg,syg)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       double precision sxg(-ngr:ngr,0:ngphi),syg(-ngr:ngr,0:ngphi)
       double precision a(30)
       common/locate_src/newcall

       save eps1_0,sep_0,ngr_0,ngphi_0

c    notation:
c
c    a(1) = inverse of Einstein radius crossing time = 2/t_hat.
c    a(2) = time of umin w.r.t. CM
c    a(3) = umin w.r.t. CM
c    a(4) = lens separation at time t_fix in units of Einstein radius
c    a(5) = angle between velocity vectory and separation vector at time t_fix
c    a(6) = epsilon_1 = mass fraction of mass 1
c    a(7) = 1/Tbin = inverse period of binary orbit
c    a(8) = v_sep = velocity in seperation direction.
c    a(9) = Rstar = stellar radius in days
c    a(10)= t_fix (not a fit parameter)
c
       eps1 = a(6)
       sep=a(4)
       if(eps1_0.eq.eps1.and.sep_0.eq.sep.and.ngr.eq.ngr_0.and.
     &                                        ngphi.eq.ngphi_0) return

       newcall=1

       do 100 jd=0,ngphi
         bphi=bphigrid*jd
         cos_bphi=cos(bphi)
         sin_bphi=sin(bphi)
         do 100 id=-ngr,ngr,1
           br=Ein_R+id*brgrid
           bx=br*cos_bphi + xcc
           by=br*sin_bphi
           call sourceloc(newcall,sep,eps1,bx,by,sxg(id,jd),syg(id,jd))
 100   continue

       eps1_0 = eps1
       sep_0 = sep
       ngr_0 = ngr
       ngphi_0 = ngphi

       return
       end

c==============================================================================

       subroutine ray_shoot(zr,zphi,iclr,tol,
     &                    eps1,sep,fampsum,famp_bndc,Ein_R,xcc,sx,sy,
     &                    Ustar,brgrid,bphigrid,im,izr,izphi,nimage,
     &                    included,iminside,
     &                    ngr,nphimax,sxg,syg,ncausin,ica_zr,
     &                    ica_zphi,ica_not_inc,nxtcain,icapt_jds,no_int)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       parameter(nbgridmax=6000000)
       double precision sxg(-ngr:ngr,0:nphimax),syg(-ngr:ngr,0:nphimax)
       double precision famp0(-nbgridmax:nbgridmax),
     &                  fampn(-nbgridmax:nbgridmax)
       integer izr(nimage),izphi(nimage),included(nimage),
     &         iminside(nimage)
       integer jds_btrack(100),idmn_btrack(100),idmx_btrack(100),
     &         ibtrack_dir(100),idmn_btrk_last(100),idmx_btrk_last(100)
       integer ica_zr(ncausin),ica_zphi(ncausin),ica_not_inc(ncausin)
       integer nxtcain(ncausin)
       integer icapt_jds(-nbgridmax:nbgridmax)
       parameter (mxclr=59,ndark=5000)
       common/limbdark/dark(ndark,0:mxclr),ald(0:mxclr),bld(0:mxclr),
     &                 darklimb(0:mxclr),hcut,bcut,iend,ildtab(0:mxclr)
       common/locate_src/newcall

       Ustar2 = Ustar*Ustar
       pi=acos(-1.d0)
       twopi=2.d0*pi
       newcall = 1

       brgridfine=brgrid
       bphigridfine=bphigrid

c      star area in grid units
c      -----------------------
       area=pi*Ustar2/(brgridfine*bphigridfine)
       area_inv=1./area

       jd_cen = izphi(im)
       id_cen = izr(im)
       if(no_int.eq.1) then
         zphi_cen = zphi
         zr_cen = zr
       else
         zphi_cen = bphigrid*jd_cen
         zr_cen = Ein_R+id_cen*brgrid
       endif

       fampsum=0.d0
       famp_bndc=0.d0
       nptin = 0
       nbtrack = 0

c      start at the center and move upward
c      -----------------------------------
       id0 = 0
       id_mbnd = 0
       id_pbnd = 0
       famp_last_row = 0.d0
       famp_last_row2 = 0.d0
       iEinR = 0
       do 50 jd = 0,2*nphimax
         famp_row=0.d0
         fbndc_row=0.d0

         jds = jd_cen + jd
         if(jds.gt.nphimax) then
           jds=jds-2*nphimax-1
         elseif(jds.lt.-nphimax) then
           jds=jds+2*nphimax+1
         endif
         jdsabs = abs(jds)
         signjds = sign(1,jds)
         bphi = zphi_cen+jd*bphigridfine
         bbx=cos(bphi)
         bby=sin(bphi)

         id_mbnd_last = id_mbnd
         id_pbnd_last = id_pbnd
         call ray_shoot_row(id0,jds,id_cen,zr_cen,bphi,iclr,tol,
     &                      eps1,sep,famp_row,fbndc_row,xcc,sx,sy,
     &                      id_mbnd,id_pbnd,Ustar2,brgrid,nin,ngr,
     &                      nphimax,sxg,syg)

c        save the first row values
c        -------------------------
         if(jd.eq.0) then
           id_mbnd0 = id_mbnd
           id_pbnd0 = id_pbnd
         endif

c        check if the center of any other image is included
c        --------------------------------------------------
         do img = 1,nimage
           if(img.ne.im.and.izphi(img).eq.jds) then
             if(id_cen+id_mbnd.lt.izr(img).and.
     &          id_cen+id_pbnd.gt.izr(img)) then
               if(img.gt.im) then
c                skip the subsequent image
                 included(img) = 1
               elseif(iminside(img).eq.0) then
c                previous image was included as point source
c                reject it
                 included(img) = 1
                 iminside(img) = 1
               else
c                skip this image
                 included(im) = 1
                 fampsum = 0.d0
                 famp_bndc = 0.d0
                 return
               endif
             endif
           endif
         enddo

c        remove any caustic points that have been included
c        -------------------------------------------------
         if(icapt_jds(jds).gt.0) then
           icausin = icapt_jds(jds)
           do
             if(id_cen+id_mbnd.le.ica_zr(icausin).and.
     &          id_cen+id_pbnd.ge.ica_zr(icausin)) then
               ica_not_inc(icausin) = 0
             endif
             icausin = nxtcain(icausin)
             if(icausin.le.0) exit
           enddo
         endif

         fampsum=fampsum+famp_row
         famp_bndc=famp_bndc+fbndc_row
         id0 = (id_mbnd+id_pbnd)/2

c        exit loop if the last row included no image points
c        --------------------------------------------------
         if(jd.eq.0) famp_first_row = famp_row+fbndc_row
         if(nin.le.0) go to 51
         famp_last_row2 = famp_last_row
         famp_last_row = famp_row+fbndc_row

c        check for a "backtrack" image
c        -----------------------------
         nbtrack0 = nbtrack
         if(jd.ne.0.and.id_mbnd.lt.id_mbnd_last) then
           jdsm1 = jds - 1
           if(jdsm1.lt.-nphimax) jdsm1=jdsm1+2*nphimax+1
           bphim = bphi - bphigridfine
           call check_bktrack(jdsm1,bphim,id_mbnd,id_mbnd_last,
     &               id_cen,zr_cen,brgrid,eps1,sep,xcc,sx,sy,Ustar2,
     &               nbpts,nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
     &               ngr,nphimax,sxg,syg)
         endif
         if(jd.ne.0.and.id_pbnd.gt.id_pbnd_last) then
           jdsm1 = jds - 1
           if(jdsm1.lt.-nphimax) jdsm1=jdsm1+2*nphimax+1
           bphim = bphi - bphigridfine
           call check_bktrack(jdsm1,bphim,id_pbnd_last,id_pbnd,
     &               id_cen,zr_cen,brgrid,eps1,sep,xcc,sx,sy,Ustar2,
     &               nbpts,nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
     &               ngr,nphimax,sxg,syg)
         endif
         if(nbtrack.gt.nbtrack0) then
           do ibk = nbtrack0+1,nbtrack
             ibtrack_dir(ibk) = -1
             idmn_btrk_last(ibk) = id_mbnd
             idmx_btrk_last(ibk) = id_pbnd
           enddo
         endif
         nbtrack0 = nbtrack
cccc        ray_shoot_row will not miss this side
ccc         if(jd.ne.0.and.id_mbnd.gt.id_mbnd_last) then
ccc           call check_bktrack(jds,bphi,id_mbnd_last,id_mbnd,
ccc     &               id_cen,zr_cen,brgrid,eps1,sep,xcc,sx,sy,Ustar2,
ccc     &               nbpts,nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
ccc     &               ngr,nphimax,sxg,syg)
ccc         endif
         if(jd.ne.0.and.id_pbnd.lt.id_pbnd_last) then
           call check_bktrack(jds,bphi,id_pbnd,id_pbnd_last,
     &               id_cen,zr_cen,brgrid,eps1,sep,xcc,sx,sy,Ustar2,
     &               nbpts,nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
     &               ngr,nphimax,sxg,syg)
         endif
         if(nbtrack.gt.nbtrack0) then
           do ibk = nbtrack0+1,nbtrack
             ibtrack_dir(ibk) = 2
             idmn_btrk_last(ibk) = id_mbnd_last
             idmx_btrk_last(ibk) = id_pbnd_last
           enddo
         endif
         id_mbnd_last = id_mbnd
         id_pbnd_last = id_pbnd
 50    continue
c      we've completed a circle: is it an Einstein Ring?
c      -------------------------------------------------
       if(id_pbnd.gt.id_mbnd0.and.id_pbnd0.gt.id_mbnd) then
c        Einstein ring
c        -------------
         iEinR = 1
         id_mbnd_EinR = id_mbnd
         id_pbnd_EinR = id_pbnd
         go to 81
       else
c        image covers all phi without completing the ring
c        check for backtrack
c        ------------------------------------------------
         jdsp1 = jds + 1
         if(jdsp1.gt.nphimax) jdsp1=jdsp1-2*nphimax-1
         bphip = bphi + bphigridfine
         nbtrack0 = nbtrack
         call check_bktrack(jdsp1,bphip,id_mbnd,id_pbnd,
     &               id_cen,zr_cen,brgrid,eps1,sep,xcc,sx,sy,Ustar2,
     &               nbpts,nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
     &               ngr,nphimax,sxg,syg)
         if(nbtrack.gt.nbtrack0) then
           do ibk = nbtrack0+1,nbtrack
             ibtrack_dir(ibk) = 1
             idmn_btrk_last(ibk) = id_mbnd
             idmx_btrk_last(ibk) = id_pbnd
           enddo
         endif
c        check for backtrack from the first row
c        --------------------------------------
         jdsm1 = jd_cen - 1
         if(jdsm1.lt.-nphimax) jdsm1=jdsm1+2*nphimax+1
         bphim = bphi - bphigridfine
         nbtrack0 = nbtrack
         call check_bktrack(jdsm1,bphim,id_mbnd0,id_pbnd0,
     &               id_cen,zr_cen,brgrid,eps1,sep,xcc,sx,sy,Ustar2,
     &               nbpts,nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
     &               ngr,nphimax,sxg,syg)
         if(nbtrack.gt.nbtrack0) then
           do ibk = nbtrack0+1,nbtrack
             ibtrack_dir(ibk) = -1
             idmn_btrk_last(ibk) = id_mbnd0
             idmx_btrk_last(ibk) = id_pbnd0
           enddo
         endif
         write(6,*)
     &    'WARNING: unusual case - phi wrapping without Einstein Ring'
       endif
 51    continue
       jd00 = jd
       jdmin = jd - 2*nphimax - 1
       id_mbnd_top = id_mbnd
       id_pbnd_top = id_pbnd

c      find the maximum phi-extent of the image
c      ----------------------------------------
       jd = jd - 1
       bphi = zphi_cen+jd*bphigridfine
       bphi1 = bphi+bphigridfine
       bbx=cos(bphi)
       bby=sin(bphi)
       hmax=0.d0
       jds = jd+jd_cen
       if(jds.gt.nphimax) then
         jds=jds-2*nphimax-1
       elseif(jds.lt.-nphimax) then
         jds=jds+2*nphimax+1
       endif
       jdsabs = abs(jds)
       signjds = sign(1,jds)
       do id = id_mbnd_last+1,id_pbnd_last-1
         ids = id+id_cen
         br = zr_cen+id*brgrid
         if(abs(ids).le.ngr.and.ngr.ne.0) then
c          recall the saved source locations
c          ---------------------------------
           sxe=sxg(ids,jdsabs)
           sye=signjds*syg(ids,jdsabs)
         else
c          grid value not saved - must calculate
c          -------------------------------------
           bx=br*bbx + xcc
           by=br*bby
           call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
         endif
         ds2=((sx-sxe)**2+(sy-sye)**2)
         if(ds2.le.Ustar2) then
           ilimb=min(int(ndark*ds2/Ustar2)+1,ndark)
           h=crosslocphi(bphi,bphi1,br,brdum,tol,
     &                   sx,sy,Ustar2,sep,eps1,bphigridfine)
           if(h.gt.hmax) hmax=h
         endif
       enddo

       if(iend.eq.0.and.jd00.ne.0) then
c        used closed formula assuming f(boundary) = 0 (old formula)
c        --------------------------------------------
         famp_bndc=famp_bndc-0.5d0*(1.d0-hmax)*famp_last_row
       elseif((famp_last_row2.eq.0.d0.or.iend.le.1.).and.jd00.ne.0) then
c        new single pt formula
c        ---------------------
         famp_bndc=famp_bndc+(hmax-0.5d0)*famp_last_row
       elseif(jd00.ne.0) then
c        new 2-pt formula
c        ----------------
         w = (9.d0/8.d0) - 0.5d0*hmax**2
         v = 1.5d0 + hmax - w
         famp_bndc=famp_bndc + (v-1.d0)*famp_last_row
     &                       + (w-1.d0)*famp_last_row2
       endif

c      start 1 row below the center and move downward
c      -----------------------------------------------
       id_mbnd = id_mbnd0
       id_pbnd = id_pbnd0
       id0 = 0
       famp_last_row = famp_first_row
       famp_last_row2 = 0.d0
       do 70 jd = -1,jdmin,-1
         famp_row=0.d0
         fbndc_row=0.d0

         jds = jd_cen + jd
         if(jds.gt.nphimax) then
           jds=jds-2*nphimax-1
         elseif(jds.lt.-nphimax) then
           jds=jds+2*nphimax+1
         endif
         jdsabs = abs(jds)
         signjds = sign(1,jds)
         bphi = zphi_cen+jd*bphigridfine
         bbx=cos(bphi)
         bby=sin(bphi)

         id_mbnd_last = id_mbnd
         id_pbnd_last = id_pbnd
         call ray_shoot_row(id0,jds,id_cen,zr_cen,bphi,iclr,tol,
     &                      eps1,sep,famp_row,fbndc_row,xcc,sx,sy,
     &                      id_mbnd,id_pbnd,Ustar2,brgrid,nin,ngr,
     &                      nphimax,sxg,syg)

c        check if the center of any other image is included
c        --------------------------------------------------
         do img = 1,nimage
           if(img.ne.im.and.izphi(img).eq.jds) then
             if(id_cen+id_mbnd.lt.izr(img).and.
     &          id_cen+id_pbnd.gt.izr(img)) then
               if(img.gt.im) then
c                skip the subsequent image
                 included(img) = 1
               elseif(iminside(img).eq.0) then
c                previous image was included as point source
c                reject it
                 included(img) = 1
                 iminside(img) = 1
               else
c                skip this image
                 included(im) = 1
                 fampsum = 0.d0
                 famp_bndc = 0.d0
                 return
               endif
             endif
           endif
         enddo

c        remove any caustic points that have been included
c        -------------------------------------------------
         if(icapt_jds(jds).gt.0) then
           icausin = icapt_jds(jds)
           do
             if(id_cen+id_mbnd.le.ica_zr(icausin).and.
     &          id_cen+id_pbnd.ge.ica_zr(icausin)) then
               ica_not_inc(icausin) = 0
             endif
             icausin = nxtcain(icausin)
             if(icausin.le.0) exit
           enddo
         endif

         fampsum=fampsum+famp_row
         famp_bndc=famp_bndc+fbndc_row
         id0 = (id_mbnd+id_pbnd)/2

c        exit loop if the last row included no image points
c        --------------------------------------------------
         if(nin.le.0) go to 71
         famp_last_row2 = famp_last_row
         famp_last_row = famp_row+fbndc_row

c        check for a "backtrack" image
c        -----------------------------
         nbtrack0 = nbtrack
         if(id_mbnd.lt.id_mbnd_last) then
           jdsp1 = jds + 1
           if(jdsp1.gt.nphimax) jdsp1=jdsp1-2*nphimax-1
           bphip = bphi + bphigridfine
           call check_bktrack(jdsp1,bphip,id_mbnd,id_mbnd_last,
     &               id_cen,zr_cen,brgrid,eps1,sep,xcc,sx,sy,Ustar2,
     &               nbpts,nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
     &               ngr,nphimax,sxg,syg)
         endif
         if(id_pbnd.gt.id_pbnd_last) then
           jdsp1 = jds + 1
           if(jdsp1.gt.nphimax) jdsp1=jdsp1-2*nphimax-1
           bphip = bphi + bphigridfine
           call check_bktrack(jdsp1,bphip,id_pbnd_last,id_pbnd,
     &               id_cen,zr_cen,brgrid,eps1,sep,xcc,sx,sy,Ustar2,
     &               nbpts,nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
     &               ngr,nphimax,sxg,syg)
         endif
         if(nbtrack.gt.nbtrack0) then
           do ibk = nbtrack0+1,nbtrack
             ibtrack_dir(ibk) = 1
             idmn_btrk_last(ibk) = id_mbnd
             idmx_btrk_last(ibk) = id_pbnd
           enddo
         endif
         nbtrack0 = nbtrack
cccc        ray_shoot_row will not miss this side
ccc         if(id_mbnd.gt.id_mbnd_last) then
ccc           call check_bktrack(jds,bphi,id_mbnd_last,id_mbnd,
ccc     &               id_cen,zr_cen,brgrid,eps1,sep,xcc,sx,sy,Ustar2,
ccc     &               nbpts,nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
ccc     &               ngr,nphimax,sxg,syg)
ccc         endif
         if(id_pbnd.lt.id_pbnd_last) then
           call check_bktrack(jds,bphi,id_pbnd,id_pbnd_last,
     &               id_cen,zr_cen,brgrid,eps1,sep,xcc,sx,sy,Ustar2,
     &               nbpts,nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
     &               ngr,nphimax,sxg,syg)
         endif
         if(nbtrack.gt.nbtrack0) then
           do ibk = nbtrack0+1,nbtrack
             ibtrack_dir(ibk) = -2
             idmn_btrk_last(ibk) = id_mbnd_last
             idmx_btrk_last(ibk) = id_pbnd_last
           enddo
         endif
         id_mbnd_last = id_mbnd
         id_pbnd_last = id_pbnd
 70    continue
c      we've completed a circle: is it an Einstein Ring?
c      -------------------------------------------------
       if(id_pbnd.gt.id_mbnd_top.and.id_pbnd_top.gt.id_mbnd) then
c        Einstein ring
c        -------------
         iEinR = -1
         id_mbnd_EinR = id_mbnd
         id_pbnd_EinR = id_pbnd
         go to 81
       elseif(jdmin.gt.-1) then
c        fully wrapped in do 50 loop
c        ---------------------------
         go to 81
       else
c        image covers all phi without completing the ring
c        check for backtrack
c        ------------------------------------------------
         jdsm1 = jds - 1
         if(jdsm1.gt.nphimax) jdsm1=jdsm1+2*nphimax+1
         bphim = bphi - bphigridfine
         nbtrack0 = nbtrack
         call check_bktrack(jdsm1,bphim,id_mbnd,id_pbnd,
     &               id_cen,zr_cen,brgrid,eps1,sep,xcc,sx,sy,Ustar2,
     &               nbpts,nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
     &               ngr,nphimax,sxg,syg)
         if(nbtrack.gt.nbtrack0) then
           do ibk = nbtrack0+1,nbtrack
             ibtrack_dir(ibk) = -1
             idmn_btrk_last(ibk) = id_mbnd
             idmx_btrk_last(ibk) = id_pbnd
           enddo
         endif
         write(6,*)
     &    'WARNING: unusual case2 - phi wrapping without Einstein Ring'
       endif

 71    continue
       jd00 = jd

c      find the maximum phi-extent of the image
c      ----------------------------------------
       jd = jd + 1
       bphi = zphi_cen+jd*bphigridfine
       bphi0 = bphi-bphigridfine
       bbx=cos(bphi)
       bby=sin(bphi)
       hmin=1.d0
       jds = jd+jd_cen
       if(jds.gt.nphimax) then
         jds=jds-2*nphimax-1
       elseif(jds.lt.-nphimax) then
         jds=jds+2*nphimax+1
       endif
       jdsabs = abs(jds)
       signjds = sign(1,jds)
       do id = id_mbnd_last+1,id_pbnd_last-1
         ids = id+id_cen
         br = zr_cen+id*brgrid
         if(abs(ids).le.ngr.and.ngr.ne.0) then
c          recall the saved source locations
c          ---------------------------------
           sxe=sxg(ids,jdsabs)
           sye=signjds*syg(ids,jdsabs)
         else
c          grid value not saved - must calculate
c          -------------------------------------
           bx=br*bbx + xcc
           by=br*bby
           call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
         endif
         ds2=((sx-sxe)**2+(sy-sye)**2)
         if(ds2.le.Ustar2) then
           ilimb=min(int(ndark*ds2/Ustar2)+1,ndark)
           h=crosslocphi(bphi0,bphi,br,brdum,tol,
     &                   sx,sy,Ustar2,sep,eps1,bphigridfine)
           if(h.lt.hmin) hmin=h
         endif
       enddo

       if(iend.eq.0.and.jd00.ne.0) then
c        used closed formula assuming f(boundary) = 0 (old formula)
c        --------------------------------------------
         famp_bndc=famp_bndc-0.5d0*hmin*famp_last_row
       elseif((famp_last_row2.eq.0.d0.or.iend.le.1.).and.jd00.ne.0) then
c        new single pt formula
c        ---------------------
         famp_bndc=famp_bndc+(0.5d0-hmin)*famp_last_row
       elseif(jd00.ne.0) then
c        new 2-pt formula
c        ----------------
         w = (9.d0/8.d0) - 0.5d0*(1.d0-hmin)**2
         v = 1.5d0 + (1.d0-hmin) - w
         famp_bndc=famp_bndc + (v-1.d0)*famp_last_row
     &                       + (w-1.d0)*famp_last_row2
       endif

c      81 is entry for Einstein ring case with no boundary
 81    continue

       if(nbtrack.gt.0) then
         if(nbtrack.gt.1) then
c          check for backtrack duplicates
c          ------------------------------
           do ib = 1,nbtrack-1
             do jb = ib+1,nbtrack
               if(jds_btrack(ib).eq.jds_btrack(jb).and.
     &            idmn_btrack(ib).eq.idmn_btrack(jb).and.
     &            idmx_btrack(ib).eq.idmx_btrack(jb).and.
     &            ibtrack_dir(ib)*ibtrack_dir(jb).gt.0)  then
c                duplicate!
                 ibtrack_dir(jb) = 0
               endif
             enddo
           enddo
         endif

c        now integrate the "backtrack" parts of the image
         ibtrack = 0
         do
 75        ibtrack = ibtrack + 1
           if(ibtrack.gt.nbtrack) exit
           if(ibtrack_dir(ibtrack).eq.0)  go to 75
           jd0 = jds_btrack(ibtrack) - jd_cen
           if(jd0.gt.nphimax) then
             jd0=jd0-2*nphimax-1
           elseif(jd0.lt.-nphimax) then
             jd0=jd0+2*nphimax+1
           endif
           id_mbnd = idmn_btrack(ibtrack)
           id_pbnd = idmx_btrack(ibtrack)
           id_mbnd_last = idmn_btrk_last(ibtrack)
           id_pbnd_last = idmx_btrk_last(ibtrack)
           if(abs(ibtrack_dir(ibtrack)).gt.1) then
             ibtrack_dir(ibtrack) = ibtrack_dir(ibtrack)/2
             nocheck = 2
           else
             nocheck = 1
           endif

c          start at the center and move upward
c          -----------------------------------
           id0 = (id_mbnd+id_pbnd)/2
           jd_lim = ibtrack_dir(ibtrack)*nphimax
           famp_last_row = 0.d0
           famp_last_row2 = 0.d0
           do 100 jd = jd0,jd_lim,ibtrack_dir(ibtrack)
             famp_row=0.d0
             fbndc_row=0.d0

             jds = jd_cen + jd
             if(jds.gt.nphimax) then
               jds=jds-2*nphimax-1
             elseif(jds.lt.-nphimax) then
               jds=jds+2*nphimax+1
             endif
             jdsabs = abs(jds)
             signjds = sign(1,jds)

c            escape from the loop if we hit the first row of an Einstein ring
c            ----------------------------------------------------------------
             if(jds.eq.jd_cen) then
               if(ibtrack_dir(ibtrack)*iEinR.eq.-1) then
                 if(id_mbnd.lt.id_pdnd0.and.id_mbnd0.lt.id_pbnd)
     &                         go to 101
               endif
             endif

c            escape from the loop if we hit the last row of an Einstein ring
c            ---------------------------------------------------------------
             if(jds.eq.jd_cen-1.or.jds.eq.jd_cen+2*nphimax.or.
     &          jds.eq.jd_cen-2*nphimax-2) then
               if(ibtrack_dir(ibtrack)*iEinR.eq.-1) then
                 if(id_mbnd.lt.id_pdnd_EinR.and.id_mbnd_EinR.lt.id_pbnd)
     &                         go to 101
               endif
             endif

             bphi = zphi_cen+jd*bphigridfine
             bbx=cos(bphi)
             bby=sin(bphi)

             if(jd.ne.jd0) then
               id_mbnd_last = id_mbnd
               id_pbnd_last = id_pbnd
             endif
             call ray_shoot_row(id0,jds,id_cen,zr_cen,bphi,iclr,tol,
     &                      eps1,sep,famp_row,fbndc_row,xcc,sx,sy,
     &                      id_mbnd,id_pbnd,Ustar2,brgrid,nin,ngr,
     &                      nphimax,sxg,syg)

c            check if the center of any other image is included
c            --------------------------------------------------
             do img = 1,nimage
               if(img.ne.im.and.izphi(img).eq.jds) then
                 if(id_cen+id_mbnd.lt.izr(img).and.
     &              id_cen+id_pbnd.gt.izr(img)) then
                   if(img.gt.im) then
c                    skip the subsequent image
                     included(img) = 1
                   elseif(iminside(img).eq.0) then
c                    previous image was included as point source
c                    reject it
                     included(img) = 1
                     iminside(img) = 1
                   else
c                    skip this image
                     included(im) = 1
                     fampsum = 0.d0
                     famp_bndc = 0.d0
                     return
                   endif
                 endif
               endif
             enddo

c            remove any caustic points that have been included
c            -------------------------------------------------
             if(icapt_jds(jds).gt.0) then
               icausin = icapt_jds(jds)
               do
                 if(id_cen+id_mbnd.le.ica_zr(icausin).and.
     &              id_cen+id_pbnd.ge.ica_zr(icausin)) then
                   ica_not_inc(icausin) = 0
                 endif
                 icausin = nxtcain(icausin)
                 if(icausin.le.0) exit
               enddo
             endif

             fampsum=fampsum+famp_row
             famp_bndc=famp_bndc+fbndc_row
             id0 = (id_mbnd+id_pbnd)/2

c            check to see if we have hit another "backtrack" image boundary
c            --------------------------------------------------------------
             if(ibtrack.lt.nbtrack) then
               do jbtrack = ibtrack+1,nbtrack
                 if(jds_btrack(jbtrack).eq.jds) then
                   if(idmn_btrack(jbtrack).lt.id_pbnd.and.
     &                idmx_btrack(jbtrack).gt.id_mbnd)  then
c                    back-track duplicate! - forget it
                     if(ibtrack_dir(ibtrack)*ibtrack_dir(jbtrack).gt.0)
     &                                     then
c                      running over another "backtrack" image in the same 
c                      direction probably means double counting
c                      --------------------------------------------------
                       write(6,*)'double counting backtrack images'
                       stop
                     else
                       ibtrack_dir(jbtrack) = 0
c                      exit the loop
c                      and skip the "last row" boundary correction
                       go to 121
                     endif
                   endif
                 endif
               enddo
             endif

c            exit loop if the last row included no image points
c            --------------------------------------------------
             if(nin.le.0) go to 101
             famp_last_row2 = famp_last_row
             famp_last_row = famp_row+fbndc_row

c            check for a "backtrack" image
c            -----------------------------
             nbtrack0 = nbtrack
cccc        ray_shoot_row will not miss this side
ccc             if(nocheck.eq.0.and.id_mbnd_last.lt.id_mbnd) then
ccc               call check_bktrack(jds,bphi,id_mbnd_last,id_mbnd,
ccc     &                 id_cen,zr_cen,brgrid,eps1,sep,xcc,sx,sy,Ustar2,
ccc     &                 nbpts,nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
ccc     &                 ngr,nphimax,sxg,syg)
ccc             endif
             if(nocheck.eq.0.and.id_pbnd_last.gt.id_pbnd) then
               call check_bktrack(jds,bphi,id_pbnd,id_pbnd_last,
     &                 id_cen,zr_cen,brgrid,eps1,sep,xcc,sx,sy,Ustar2,
     &                 nbpts,nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
     &                 ngr,nphimax,sxg,syg)
             endif
             if(nbtrack.gt.nbtrack0) then
               do ibk = nbtrack0+1,nbtrack
                 ibtrack_dir(ibk) = 2*ibtrack_dir(ibtrack)
                 idmn_btrk_last(ibk) = id_mbnd_last
                 idmx_btrk_last(ibk) = id_pbnd_last
                 if(ibk.gt.1) then
                   do jbk = 1,ibk-1
                     if(jds_btrack(ibk).eq.jds_btrack(jbk).and.
     &                  idmn_btrack(ibk).eq.idmn_btrack(jbk).and.
     &                  idmx_btrack(ibk).eq.idmx_btrack(jbk).and.
     &                  ibtrack_dir(ibk)*ibtrack_dir(jbk).gt.0)  then
c                      duplicate!
                       ibtrack_dir(ibk) = 0
                     endif
                   enddo
                 endif
               enddo
             endif
             nbtrack0 = nbtrack
             if(id_mbnd_last.gt.id_mbnd) then
               jds1 = jds - ibtrack_dir(ibtrack)
               if(jds1.lt.-nphimax) jds1=jds1+2*nphimax+1
               if(jds1.gt.nphimax) jds1=jds1-2*nphimax-1
               bphi1 = bphi - ibtrack_dir(ibtrack)*bphigridfine
               call check_bktrack(jds1,bphi1,id_mbnd,id_mbnd_last,
     &                 id_cen,zr_cen,brgrid,eps1,sep,xcc,sx,sy,Ustar2,
     &                 nbpts,nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
     &                 ngr,nphimax,sxg,syg)
             endif
             if(id_pbnd_last.lt.id_pbnd) then
               jds1 = jds - ibtrack_dir(ibtrack)
               if(jds1.lt.-nphimax) jds1=jds1+2*nphimax+1
               if(jds1.gt.nphimax) jds1=jds1-2*nphimax-1
               bphi1 = bphi - ibtrack_dir(ibtrack)*bphigridfine
               call check_bktrack(jds1,bphi1,id_pbnd_last,id_pbnd,
     &                 id_cen,zr_cen,brgrid,eps1,sep,xcc,sx,sy,Ustar2,
     &                 nbpts,nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
     &                 ngr,nphimax,sxg,syg)
             endif
             if(nbtrack.gt.nbtrack0) then
               do ibk = nbtrack0+1,nbtrack
                 ibtrack_dir(ibk) = -ibtrack_dir(ibtrack)
                 idmn_btrk_last(ibk) = id_mbnd
                 idmx_btrk_last(ibk) = id_pbnd
                 if(ibk.gt.1) then
                   do jbk = 1,ibk-1
                     if(jds_btrack(ibk).eq.jds_btrack(jbk).and.
     &                  idmn_btrack(ibk).eq.idmn_btrack(jbk).and.
     &                  idmx_btrack(ibk).eq.idmx_btrack(jbk).and.
     &                  ibtrack_dir(ibk)*ibtrack_dir(jbk).gt.0)  then
c                      duplicate!
                       ibtrack_dir(ibk) = 0
                     endif
                   enddo
                 endif
               enddo
             endif
             id_mbnd_last = id_mbnd
             id_pbnd_last = id_pbnd
             nocheck = 0
 100       continue
 101       continue

c          find the maximum phi-extent of the image
c          ----------------------------------------
           jd = jd - ibtrack_dir(ibtrack)
           bphi = zphi_cen+jd*bphigridfine
           if(ibtrack_dir(ibtrack).gt.0) then
             bphi0 = bphi
             bphi1 = bphi+bphigridfine
           else
             bphi1 = bphi
             bphi0 = bphi-bphigridfine
           endif
           bbx=cos(bphi)
           bby=sin(bphi)
           hmax=0.d0
           hmin=1.d0
           jds = jd+jd_cen
           if(jds.gt.nphimax) then
             jds=jds-2*nphimax-1
           elseif(jds.lt.-nphimax) then
             jds=jds+2*nphimax+1
           endif
           jdsabs = abs(jds)
           signjds = sign(1,jds)
           do id = id_mbnd_last+1,id_pbnd_last-1
             ids = id+id_cen
             br = zr_cen+id*brgrid
             if(abs(ids).le.ngr.and.ngr.ne.0) then
c              recall the saved source locations
c              ---------------------------------
               sxe=sxg(ids,jdsabs)
               sye=signjds*syg(ids,jdsabs)
             else
c              grid value not saved - must calculate
c              -------------------------------------
               bx=br*bbx + xcc
               by=br*bby
               call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
             endif
             ds2=((sx-sxe)**2+(sy-sye)**2)
             if(ds2.le.Ustar2) then
               ilimb=min(int(ndark*ds2/Ustar2)+1,ndark)
               h=crosslocphi(bphi0,bphi1,br,brdum,tol,
     &                       sx,sy,Ustar2,sep,eps1,bphigridfine)
               if(h.gt.hmax) hmax=h
               if(h.lt.hmin) hmin=h
             endif
           enddo

           if(iend.eq.0) then
c            used closed formula assuming f(boundary) = 0 (old formula)
c            --------------------------------------------
             if(ibtrack_dir(ibtrack).gt.0) then
               famp_bndc=famp_bndc-0.5d0*(1.d0-hmax)*famp_last_row
             else
               famp_bndc=famp_bndc+(0.5d0-hmin)*famp_last_row
             endif
           elseif(famp_last_row2.eq.0.d0.or.iend.le.1) then
c            new single pt formula
c            ---------------------
             if(ibtrack_dir(ibtrack).gt.0) then
               famp_bndc=famp_bndc+(hmax-0.5d0)*famp_last_row
             else
               famp_bndc=famp_bndc+(0.5d0-hmin)*famp_last_row
             endif
           else
c            new 2-pt formula
c            ----------------
             if(ibtrack_dir(ibtrack).gt.0) then
               w = (9.d0/8.d0) - 0.5d0*hmax**2
               v = 1.5d0 + hmax - w
             else
               w = (9.d0/8.d0) - 0.5d0*(1.d0-hmin)**2
               v = 1.5d0 + (1.d0-hmin) - w
             endif
             famp_bndc=famp_bndc + (v-1.d0)*famp_last_row
     &                           + (w-1.d0)*famp_last_row2
           endif
 121       continue
         enddo
       endif


       fampsum=fampsum*area_inv
       famp_bndc=famp_bndc*area_inv

       return
       end
c==============================================================================

       subroutine ray_shoot_row(id0,jds,id_cen,zr_cen,bphi,iclr,tol,
     &                      eps1,sep,famp_row,fbndc_row,xcc,sx,sy,
     &                      id_mbnd,id_pbnd,Ustar2,brgrid,nin,ngr,
     &                      nphimax,sxg,syg)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       parameter(nbgridmax=6000000)
       parameter(nrow_max=20000)
       double precision sxg(-ngr:ngr,0:nphimax),syg(-ngr:ngr,0:nphimax)
       double precision famp1(-nrow_max:nrow_max),
     &                 brsave(-nrow_max:nrow_max)
       parameter (mxclr=59,ndark=5000)
       common/limbdark/dark(ndark,0:mxclr),ald(0:mxclr),bld(0:mxclr),
     &                 darklimb(0:mxclr),hcut,bcut,iend,ildtab(0:mxclr)
       common/locate_src/newcall

       famp_row=0.d0
       fbndc_row=0.d0

       jdsabs = abs(jds)
       signjds = sign(1,jds)
       bbx=cos(bphi)
       bby=sin(bphi)

c      calculate the integrand at each grid point
c      ------------------------------------------
       id_mbnd_last = id_mbnd
       id_pbnd_last = id_pbnd
       id_min = id_mbnd
       id_max = id_pbnd
       nin = 0
       do 50 id = id_mbnd_last,9999999
         idd = id-id0
         br = zr_cen+id*brgrid
         brsave(idd) = br
         ids = id+id_cen
         if(abs(ids).le.ngr.and.ngr.ne.0) then
c          recall the saved source locations
c          ---------------------------------
           sxe=sxg(ids,jdsabs)
           sye=signjds*syg(ids,jdsabs)
         else
c          grid value not saved - must calculate
c          -------------------------------------
           bx=br*bbx + xcc
           by=br*bby
           call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
         endif
         ds2=((sx-sxe)**2+(sy-sye)**2)
         if(ds2.le.Ustar2) then
c          lensed point
           ilimb=min(int(ndark*ds2/Ustar2)+1,ndark)
           famp1(idd)=br*dark(ilimb,iclr)
           famp_row=famp_row+famp1(idd)
           nin = nin + 1
         else
           famp1(idd)=0.d0
ccc           if(id.ge.id_pbnd_last) then
           if(id.ge.id_pbnd_last.or.nin.ge.1) then
ccc             id_max = max(id_max,id)
             id_max = id
             go to 51
           endif
         endif
 50    continue
 51    continue

c      possibly move backwards from id_mbnd_last
c      -----------------------------------------
       if(famp1(id_mbnd_last-id0).gt.0.d0) then
         do 70 id = id_mbnd_last-1,-9999999,-1
           idd = id-id0
           br = zr_cen+id*brgrid
           brsave(idd) = br
           ids = id+id_cen
           if(abs(ids).le.ngr.and.ngr.ne.0) then
c            recall the saved source locations
c            ---------------------------------
             sxe=sxg(ids,jdsabs)
             sye=signjds*syg(ids,jdsabs)
           else
c            grid value not saved - must calculate
c            -------------------------------------
             bx=br*bbx + xcc
             by=br*bby
             call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
           endif
           ds2=((sx-sxe)**2+(sy-sye)**2)
           if(ds2.le.Ustar2) then
c            lensed point
             ilimb=min(int(ndark*ds2/Ustar2)+1,ndark)
             famp1(idd)=br*dark(ilimb,iclr)
             famp_row=famp_row+famp1(idd)
             nin = nin + 1
           else
             famp1(idd)=0.d0
ccc             if(id.le.id_mbnd_last) then 
             id_min = id
             go to 71
ccc             endif
           endif
 70      continue
 71      continue
       endif

       if(nin.eq.0) then
         id_mbnd = id0
         id_pbnd = id0
         return
       endif

c      find new grid boundaries
c      ------------------------
ccc       imin = 0
ccc       do id = id_min+1,id_max
ccc         idd = id-id0
ccc         if(imin.eq.0.and.famp1(idd-1).eq.0.d0.and.
ccc     &                    famp1(idd).gt.0.d0) then
ccc           id_mbnd = id-1
ccc           imin = 1
ccc         endif
ccc         if(famp1(idd-1).gt.0.d0.and.famp1(idd).eq.0.d0) id_pbnd = id
ccc       enddo
       do 80 id = id_min+1,id_max
         idd = id-id0
         if(famp1(idd-1).eq.0.d0.and.famp1(idd).gt.0.d0) then
           id_mbnd = id-1
           go to 81
         endif
 80    continue
 81    continue
       id_pbnd = id_max

c      do the boundary corrections
c      ---------------------------
       if(nin.gt.0) then
c        left boundary
         id = id_mbnd+1
         idd = id-id0
         br = brsave(idd)
         br_last = br - brgrid
         h=crosslocr(br_last,br,bbx,bby,tol,
     &               sx,sy,Ustar2,sep,eps1,brgrid)
         hbg=(h-1.d0)
         if(famp1(idd+1).eq.0.d0.or.iend.le.2) then
           fbndc_cor1=endpoint(famp1(idd),hbg,br,brgrid,iclr)
         elseif(iend.eq.3) then
           fbndc_cor1=endpoint3(famp1(idd),famp1(idd+1),hbg,
     &                          br,brgrid,iclr)
         else
           br15 = br+0.5d0*brgrid
           bx15 = br15*bbx + xcc
           by15 = br15*bby
           call sourceloc(newcall,sep,eps1,bx15,by15,sxe15,sye15)
           ds2=((sx-sxe15)**2+(sy-sye15)**2)
           ilimb=min(int(ndark*ds2/Ustar2)+1,ndark)
           famp15=br15*dark(ilimb,iclr)
           if(iend.eq.5) then
             fbndc_cor1=endpoint5(famp1(idd),famp15,hbg,br,brgrid,iclr)
           else
             fbndc_cor1=endpoint2(famp1(idd),famp15,famp1(idd+1),hbg,
     &                            br,brgrid,iclr)
           endif
         endif
         fbndc_row=fbndc_row + fbndc_cor1
c        right boundary
         id = id_pbnd-1
         idd = id-id0
         br = brsave(idd)
         br_next = br + brgrid
         h=crosslocr(br,br_next,bbx,bby,tol,
     &               sx,sy,Ustar2,sep,eps1,brgrid)
         hbg=h
         if(famp1(idd-1).eq.0.d0.or.iend.le.2) then
           fbndc_cor1=endpoint(famp1(idd),hbg,br,brgrid,iclr)
         elseif(iend.eq.3) then
           fbndc_cor1=endpoint3(famp1(idd),famp1(idd-1),hbg,
     &                          br,brgrid,iclr)
         else
           br15 = br-0.5d0*brgrid
           bx15 = br15*bbx + xcc
           by15 = br15*bby
           call sourceloc(newcall,sep,eps1,bx15,by15,sxe15,sye15)
           ds2=((sx-sxe15)**2+(sy-sye15)**2)
           ilimb=min(int(ndark*ds2/Ustar2)+1,ndark)
           famp15=br15*dark(ilimb,iclr)
           if(iend.eq.5) then
             fbndc_cor1=endpoint5(famp1(idd),famp15,hbg,br,brgrid,iclr)
           else
             fbndc_cor1=endpoint2(famp1(idd),famp15,famp1(idd-1),hbg,
     &                            br,brgrid,iclr)
           endif
         endif
         fbndc_row=fbndc_row + fbndc_cor1
       endif

       return
       end

c==============================================================================

       subroutine cen_in_images(nimages,zr,zphi,izr,izphi,iminside,
     &                          sx,sy,Ustar,brgrid,bphigrid,eps1,sep,
     &                          Ein_R,xcc,nphimax)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       double precision zr(nimages),zphi(nimages)
       integer izr(nimages),izphi(nimages),iminside(nimages)

       Ustar2 = Ustar*Ustar
       do 50 im = 1,nimages
         if(iminside(im).eq.1) then
           call cen_in_image(zr(im),zphi(im),izr(im),izphi(im),sx,sy,
     &                       Ustar2,brgrid,bphigrid,eps1,sep,
     &                       Ein_R,xcc,nphimax,ipt)
         endif
 50    continue

       return
       end

c==============================================================================

       subroutine cen_in_image(zr,zphi,izr,izphi,sx,sy,
     &                         Ustar2,brgrid,bphigrid,eps1,sep,
     &                         Ein_R,xcc,nphimax,ipt)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c      This routine checks if the grid points or its neighbors are inside
c      an image. If the grid point is not, but a neighbor is, then the
c      neighor grid coordinates are used as the integer coordinates. If none
c      are on the grid, a "point image" is declared.

       common/locate_src/newcall

       ipt = 0
       zr_int = Ein_R+izr*brgrid
       zphi_int = bphigrid*izphi
       bbx=cos(zphi_int)
       bby=sin(zphi_int)
       bx=zr_int*bbx + xcc
       by=zr_int*bby
       call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
       ds2=((sx-sxe)**2+(sy-sye)**2)
       if(ds2.le.Ustar2) then
c        lensed point
       else
c        central point is not on the grid - so find a neighbor
         do idr = -1,1,2
           br = zr_int +idr*brgrid
           bx=br*bbx + xcc
           by=br*bby
           call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
           ds2=((sx-sxe)**2+(sy-sye)**2)
           if(ds2.le.Ustar2) then
c            lensed point
             izr = izr + idr
             zr_int = zr_int + idr*brgrid
             go to 50
           endif
         enddo
         do idphi = -1,1,2
           bphi = zphi_int + idphi*bphigrid
           bbx=cos(bphi)
           bby=sin(bphi)
           bx=zr*bbx + xcc
           by=zr*bby
           call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
           ds2=((sx-sxe)**2+(sy-sye)**2)
           if(ds2.le.Ustar2) then
c            lensed point
             izphi = izphi + idphi
             if(izphi.gt.nphimax) then
               izphi=izphi-2*nphimax-1
             elseif(izphi.lt.-nphimax) then
               izphi=izphi+2*nphimax+1
             endif
             zphi_int = zphi_int + idphi*brgrid
             go to 50
           endif
         enddo
         do idphi = -1,1,2
           bphi = zphi_int + idphi*bphigrid
           bbx=cos(bphi)
           bby=sin(bphi)
           do idr = -1,1,2
             br = zr_int +idr*brgrid
             bx=br*bbx + xcc
             by=br*bby
             call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
             ds2=((sx-sxe)**2+(sy-sye)**2)
             if(ds2.le.Ustar2) then
c              lensed point
               izr = izr + idr
               zr_int = zr_int + idr*brgrid
               izphi = izphi + idphi
               if(izphi.gt.nphimax) then
                 izphi=izphi-2*nphimax-1
               elseif(izphi.lt.-nphimax) then
                 izphi=izphi+2*nphimax+1
               endif
               zphi_int = zphi_int + idphi*brgrid
               go to 50
             endif
           enddo
         enddo
ccc         write(6,*) 'WARNING: image falls between grip points'
         ipt = 1
       endif
 50    continue

       return
       end

c==============================================================================

       subroutine check_caust_list(new,ica_zr_tmp,ica_zphi_tmp,ncausin,
     &                  nxtcain,ica_zr,ica_zphi,nbgridmax,icapt_jds)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       integer ica_zr(ncausin),ica_zphi(ncausin),nxtcain(ncausin)
       integer icapt_jds(-nbgridmax:nbgridmax)

       icausin = icapt_jds(ica_zphi_tmp)
       do 
         if(icausin.eq.0) exit
         if(ica_zr(icausin).eq.ica_zr_tmp.and.
     &      ica_zphi(icausin).eq.ica_zphi_tmp) then
           new = 0
           icausin = 0
         else
           icausin = nxtcain(icausin)
         endif
       enddo

       return
       end

c==============================================================================

       subroutine check_bktrack(jds,bphi,id0,id1,id_cen,zr_cen,brgrid,
     &                      eps1,sep,xcc,sx,sy,Ustar2,nbpts,
     &                      nbtrack,jds_btrack,idmn_btrack,idmx_btrack,
     &                      ngr,nphimax,sxg,syg)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       parameter(nbgridmax=6000000)
       double precision sxg(-ngr:ngr,0:nphimax),syg(-ngr:ngr,0:nphimax)
       integer jds_btrack(100),idmn_btrack(100),idmx_btrack(100)

c      check for a "backtrack" image
c      -----------------------------
       nbtpts = 0
       jdsabs = abs(jds)
       if(jds.eq.0) then
         signjds = 1.d0
       else
         signjds = jds/jdsabs
       endif
       bbx=cos(bphi)
       bby=sin(bphi)
       do id = id0,id1
         br = zr_cen+id*brgrid
         ids = id+id_cen
         if(abs(ids).le.ngr.and.ngr.ne.0) then
c          recall the saved source locations
c          ---------------------------------
           sxe=sxg(ids,jdsabs)
           sye=signjds*syg(ids,jdsabs)
         else
c          grid value not saved - must calculate
c          -------------------------------------
           bx=br*bbx + xcc
           by=br*bby
           call sourceloc(newcall,sep,eps1,bx,by,sxe,sye)
         endif
         ds2=((sx-sxe)**2+(sy-sye)**2)
         if(ds2.le.Ustar2) then
c          lensed point
           if(nbtpts.eq.0) then
             nbtrack = nbtrack+1
             idmn_btrack(nbtrack) = id
             jds_btrack(nbtrack) = jds
           endif
           nbtpts = nbtpts + 1
           idmx_btrack(nbtrack) = id
         else
           nbtpts = 0
         endif
       enddo

       return
       end
c==============================================================================

       function endpoint(famp0,h,br,bgridfine,jclr)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       parameter (mxclr=59,ndark=5000)
       common/limbdark/dark(ndark,0:mxclr),ald(0:mxclr),bld(0:mxclr),
     &                 darklimb(0:mxclr),hcut,bcut,iend,ildtab(0:mxclr)

       habs=abs(h)
       brlimb=br+h*bgridfine
       dklimb=brlimb*darklimb(jclr)
       if(habs.ge.hcut) then
         habsp05=habs+0.5d0
         b=(2.d0/3.d0)*sqrt(habsp05/habs)
         endpoint=habsp05*(b*famp0+(1.d0-b)*dklimb)-famp0
       else
         b=2.d0/3.d0
         endpoint=habs*(b*famp0+(1.d0-b)*dklimb)-0.5d0*famp0
       endif

       return
       end

c==============================================================================

       function endpoint2(famp1,famp15,famp2,h,br,bgridfine,jclr)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       parameter (mxclr=59,ndark=5000)
       common/limbdark/dark(ndark,0:mxclr),ald(0:mxclr),bld(0:mxclr),
     &                 darklimb(0:mxclr),hcut,bcut,iend,ildtab(0:mxclr)

       habs=abs(h)
       if(habs.ge.hcut) then
         brlimb=br+h*bgridfine
         dklimb=brlimb*darklimb(jclr)
         habsp05=habs+0.5d0
         b=(2.d0/3.d0)*sqrt(habsp05/habs)
         endpoint2=habsp05*(b*famp1+(1.d0-b)*dklimb)-famp1
       else
         sqh = sqrt(habs)
         sqh05 = sqrt(habs+0.5d0)
         sqh1 = sqrt(habs+1.d0)
         sqh15 = sqrt(habs+1.5d0)
         b = (2.d0/3.d0)*sqh15 - sqh - (0.75d0-0.5d0*habs)*(sqh1-sqh)
         b = b/(sqh05 - 0.5d0*sqh1 - 0.5d0*sqh)
         c = 0.75d0 - 0.5d0*(habs + b)
         a = 1.d0 - b - c
         endpoint2 = (habs+1.5d0)*(a*famp1+b*famp15+c*famp2)
     &                         -(famp1+famp2)
       endif

       return
       end

c==============================================================================

       function endpoint3(famp1,famp2,h,br,bgridfine,jclr)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       parameter (mxclr=59,ndark=5000)
       common/limbdark/dark(ndark,0:mxclr),ald(0:mxclr),bld(0:mxclr),
     &                 darklimb(0:mxclr),hcut,bcut,iend,ildtab(0:mxclr)

       habs=abs(h)
       brlimb=br+h*bgridfine
       dklimb=brlimb*darklimb(jclr)
       habsp05=habs+0.5d0
       if(habs.ge.hcut) then
         b=(2.d0/3.d0)*sqrt(habsp05/habs)
         endpoint3=habsp05*(b*famp1+(1.d0-b)*dklimb)-famp1
       else
         c=((2.d0/3.d0)*sqrt(habsp05)-bcut*sqrt(habs))/sqrt(habs+1.d0)
         endpoint3=habsp05*(c*famp2+bcut*famp1+(1.d0-bcut-c)*dklimb)
     &              -famp1
       endif

       return
       end

c==============================================================================

       function endpoint5(famp1,famp15,h,br,bgridfine,jclr)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       parameter (mxclr=59,ndark=5000)
       common/limbdark/dark(ndark,0:mxclr),ald(0:mxclr),bld(0:mxclr),
     &                 darklimb(0:mxclr),hcut,bcut,iend,ildtab(0:mxclr)

       habs=abs(h)
       brlimb=br+h*bgridfine
       dklimb=brlimb*darklimb(jclr)
       habsp05=habs+0.5d0
       if(habs.ge.hcut) then
         b=(2.d0/3.d0)*sqrt(habsp05/habs)
         endpoint5=habsp05*(b*famp1+(1.d0-b)*dklimb)-famp1
       else
         c=(2.d0/3.d0)-bcut*sqrt(habs)/sqrt(habsp05)
         endpoint5=habsp05*(c*famp15+bcut*famp1+(1.d0-bcut-c)*dklimb)
     &              -famp1
       endif

       return
       end

c==============================================================================

       subroutine limb_dark(limb_tab,limbfile)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       parameter (mxclr=59,ndark=5000)
       common/limbdark/dark(ndark,0:mxclr),ald(0:mxclr),bld(0:mxclr),
     &                 darklimb(0:mxclr),hcut,bcut,iend,ildtab(0:mxclr)
       character*80 limbfile

       f=1.d0/ndark
       do jclr = 0,mxclr
         darklimb(jclr) = 0.d0
         ildtab(jclr) = 0
       enddo
       if(limb_tab.ne.0) then
c        read in limb darkening tables
c        -----------------------------
         open(unit=1,file=limbfile,status='old')
         do 20 j = 0,mxclr
           read(1,*,end=21,err=21) jclr
           do i = 1,ndark
             read(1,*) dark(i,jclr)
           enddo
           read(1,*) darklimb(jclr)
           write(6,*) 'read limb darkening table for color',jclr
           ildtab(jclr) = 1
 20      continue
 21      continue
       endif
       nonlin = 0
       do j = 0,mxclr
         if(bld(j).ne.0) nonlin = 1
       enddo
       if(nonlin.eq.1) write(6,*) 'using non-linear limb darkening'

       do 50 j=0,mxclr
         if(darklimb(j).eq.0.d0) then
           a=ald(j)
           if(nonlin.eq.1) then
             do i=1,ndark
               b=bld(j)
               dark(i,j)=(1.d0-a+a*sqrt(1.d0-r2)-b+b*(1.d0-r2)**0.25d0)
               dark(i,j)=dark(i,j)/(1.d0-a/3.d0-0.2d0*b)
             enddo
             darklimb(j)=(1.d0-a-b)/(1.d0-a/3.d0-0.2d0*b)
           else
             eps = a
             do i=1,ndark
               r2=(i-0.5)*f
c              the Eddington approximate solution to the Gray Amosphere
c              --------------------------------------------------------
               dark(i,j)=(1.d0-eps+eps*sqrt(1.d0-r2))
               dark(i,j)=dark(i,j)/(1.d0-eps/3.d0)
             enddo
             darklimb(j)=(1.d0-eps)/(1.d0-eps/3.d0)
           endif
         endif
 50    continue

       return
       end

c==============================================================================

       subroutine sourceloc(newcall,sep,eps1,bx,by,sx,sy)

c==============================================================================
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       double complex cc,cx1,cx2,zk,zta
       save

       if(newcall.gt.0) then
         eps2=1.d0-eps1
         xx1 = -eps2*sep
         xx2 =  eps1*sep
c        xx0 = "anti" CM
         xx0 = eps1*xx2 + eps2*xx1
         cx1 = dcmplx(xx1,0.d0)
         cx2 = dcmplx(xx2,0.d0)
         cc = dcmplx(xx0,0.d0)
         newcall=0
       endif

       zk=dcmplx(bx,by)

c      the following is the lens equation
c      ----------------------------------
       zta = zk - conjg((zk-cc)/((zk-cx1)*(zk-cx2)))
ccc       sx = real(zta)
ccc       sy = imag(zta)
       sx = dreal(zta)
       sy = dimag(zta)

       return
       end

c==============================================================================

       function crosslocr(bra,brb,cosbphi,sinbphi,
     &                    tol,ssx,ssy,Ustar2,sep,eps1,bgridfine)

c==============================================================================
       use eesunhong_recipes_replacements,
     &     only: brent_wrapper_with_additional_lens_arguments
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       external starbndr

c        find the boundary
c        -----------------
         crosslocr=brent_wrapper_with_additional_lens_arguments(
     &                 starbndr,bra,brb,cosbphi,sinbphi,tol,
     &                 ssx,ssy,Ustar2,sep,eps1) -bra
         crosslocr=crosslocr/bgridfine

       return
       end

c==============================================================================

       function crosslocphi(bphia,bphib,br,brdum,
     &                    tol,ssx,ssy,Ustar2,sep,eps1,bgridfine)

c==============================================================================
       use eesunhong_recipes_replacements,
     &     only: brent_wrapper_with_additional_lens_arguments
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       external starbndphi

c        find the boundary
c        -----------------
         crosslocphi=brent_wrapper_with_additional_lens_arguments(
     &                   starbndphi,bphia,bphib,br,brdum,tol,
     &                   ssx,ssy,Ustar2,sep,eps1) -bphia
         crosslocphi=crosslocphi/bgridfine

       return
       end

       function starbndr(br,cosbphi,sinbphi,ssx,ssy,Ustar2,sep,eps1)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       common/fudge_coords/ Ein_R,xcc

       bx=br*cosbphi + xcc
       by=br*sinbphi
       call sourceloc(n,sep,eps1,bx,by,sx,sy)
       starbndr=(ssx-sx)**2+(ssy-sy)**2-Ustar2

       return
       end

       function starbndphi(bphi,br,brdum,ssx,ssy,Ustar2,sep,eps1)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       common/fudge_coords/ Ein_R,xcc

       bx=br*cos(bphi) + xcc
       by=br*sin(bphi)
       call sourceloc(n,sep,eps1,bx,by,sx,sy)
       starbndphi=(ssx-sx)**2+(ssy-sy)**2-Ustar2

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      function lenc(s)
c
c  Get the length of a character variable.  Only necessary because of the
c  Unix lobotomy of Fortran.
c
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      character*(*) s
      lenc = len(s)
      do 900 i = len(s),1,-1
          if(s(i:i).ne.' ') then
              lenc = i
              return
          end if
900   continue
      lenc = 0
      return
      end
