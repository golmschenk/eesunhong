c==============================================================================

       subroutine fcn(nfit,grad,chi2,a,iflag)

c==============================================================================
c  This program fits light curves for microlensing by binary lenses.
c
c  Routines in this file written by David Bennett.
c
c------------------------------------------------------------------------------
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       parameter (mmax=30,maxdata=80000,ngmax=12000000)
cc       parameter (mmax=30,maxdata=4000,ngmax=400000)
       double precision x(maxdata),y(maxdata),sig(maxdata)
       integer iclr(maxdata),nimage(maxdata),iwork(maxdata),
     &         jwork(maxdata)
       double precision a(mmax),grad(mmax),dtcaustic(maxdata),
     &        work(maxdata)
       parameter(mxclr=59,ndark=5000)
       double precision yfits(0:mxclr),yyfits(0:mxclr),chi2clr(0:mxclr)
       double precision A2(0:mxclr),A0(0:mxclr),A2sig(0:mxclr),
     &                  A0sig(0:mxclr)
       double precision yfit(maxdata,0:mxclr),xclr(maxdata,0:mxclr),
     &                  yclr(maxdata,0:mxclr),eclr(maxdata,0:mxclr)
       double precision fmin(0:mxclr),fmax(0:mxclr),scaus(20),tcaus(20)
       double precision fudge(0:mxclr),errmin(0:mxclr),dayoff(0:mxclr)
       double precision sxg(ngmax),syg(ngmax)
       double precision fmass(2),extvvv(0:20)
       double precision lon_obs(0:mxclr),lat_obs(0:mxclr)
       double complex z(10)
       double precision ampim(10),flux_meas(10),flux_sig(10)
       double precision fMag_off(10),A_last_obs(10)
       double precision fMag_meas(10),sig_Mag_lens(10),A_lens(10)
       integer iimage(10),icband(10),Iup(10)
       character*6 sfx(0:mxclr)
       integer nclr(0:mxclr),iclrind(maxdata),lsfx(0:mxclr),
     &         jclrinc(0:mxclr)
       integer idcaustic(12)
       double complex zalph,zgamm,box
       double precision Date(maxdata),rMag(maxdata),rErr(maxdata),
     &                    bMag(maxdata),bErr(maxdata),
     &                    rFlux(maxdata),rFerr(maxdata),
     &                    bFlux(maxdata),bFerr(maxdata)
       character*1 c(0:9)
       character*19 char_mcmc(64)
       character*200 inline
       character*120 line
       character*80 infile,ctinfile,newinffile,lcoutfile,tmpfile,parfile
       character*80 oldinffile,residfile,tmpfile2,outline,mcmc_file
       character*80 satfile(0:mxclr),limbfile
       character*30 starname,fitname
       character*20 blank20
       character*10 parname(30)
       character*1 cband(10)
       common/mcmc_output/nchar_mcmc,char_mcmc
       common/seek_opt2/mcmc_file
       common/integrate/gridUstar
       common/limbdark/dark(ndark,0:mxclr),ald(0:mxclr),bld(0:mxclr),
     &                 darklimb(0:mxclr),hcut
       parameter (max_sat=40000)
       common/satcoords/tsat(max_sat,9),rsat(max_sat,9),
     &        rasat(max_sat,9),decsat(max_sat,9),ndatsat(max_sat)
       data c /'0','1','2','3','4','5','6','7','8','9'/
       data fudge /60*1.0d0/
       data errmin /2*0.014d0,28*0.003d0,10*0.d0,10*0.003d0,10*0.d0/
       data dayoff /60*0.d0/
       data fmin /60*0.0001d0/
       data fmax /60*10000.d0/
       data ald /60*0.6d0/
       data bld /60*0.d0/
       data blank20 /'                    '/

       data hcut / 0.15d0 /
       data grid_rat / 4.d0 /
       data icc / 3 /

       save

       if(iflag.eq.1) then
         pi=acos(-1.d0)
         twobypi=2.d0/pi
         twopi=2.d0*pi
         wfm=1.d0

c        hardwired data cut parameters
c        -----------------------------
         seemax=8.0d0
         seemax74=8.0d0
         itypemax=6
         errormax=3.0d0
         idskymax=200
         ipsfc2max=140
         psfcmagmin=-6.5d0
         ipsfc2max74=160
         psfcmagmin74=-9.1
         icrdmax=170
         icosmax=1
         imispmax=50

         fudgectio=1.5d0

         minorm=200
         ndata=0
         nread=0
         icolor=2

c      hardwired parameters
c      --------------------
       newcall=1
       dt=0.01d0
       dtfine=0.001d0
       third=1.d0/3.d0

       gridUstar=0.01d0
       dcritUmin=7.d0
       dgrUstar=20.d0
       dgrfnUstar=5.d0
       gridUstar=0.03d0

c ---
         write(6,*) 'enter star name'
         read(5,'(a)') starname
         lsn=lenc(starname)
         write(6,*) 'enter fit name'
         read(5,'(a)') fitname
         lfn=lenc(fitname)
         read(5,'(a)') limbfile
         if(limbfile(1:4).eq.'limb') then
           write(6,*) 'limb darkening fits not yet implemented'
           STOP
         elseif(limbfile(1:2).eq.'no') then
           write(6,*) 'limb darkening parameters assumed - not fit'
           limb_tab = 0
         else
           limb_tab = 1
         endif
         write(6,*)
     &    'Enter RA & DEC precessed to time of maximum magnification'
         write(6,*)
     &    'format: hh mm ss.s  deg mm ss.ss'
         read(5,*) rah,ram,ras,decd,decm,decs
         alpha = (rah +  ram/60. + ras/3600.)*15
         if(decd.lt.0.d0.or.decm.lt.0.d0.or.decs.lt.0.d0) then
           delta = -(abs(decd) + abs(decm)/60. + abs(decs)/3600.)
         else
           delta = (abs(decd) + abs(decm)/60. + abs(decs)/3600.)
         endif

         write(6,*) 'enter theta_suasC, Flux0, iflux0'
         read(5,*) theta_suasC, Flux0, iflux0

         write(6,*) 'enter Dsource, Ds_sig, ncband:'
         read(5,*) Dso,Ds_sig,ncband
         if(ncband.gt.10) then
           stop 'ncband too large'
         elseif(ncband.gt.0) then
           write(6,*) 'enter Band, Band #, Mag_off',
     &           ' Lens-Mag, Sigma, A_lens, Iup values:'
           do ic = 1,ncband
             read(5,*) cband(ic),icband(ic),fMag_off(ic),
     &             fMag_meas(ic),sig_Mag_lens(ic),A_lens(ic),Iup(ic)
           enddo
         endif
         write(6,*) 'enter 0, b_gal, h_dust for exponential dust or'
         write(6,*) 'enter 1 and 21 extinction values for VVV 3-d model'
         read(5,'(a)') inline
         read(inline,*) iextinct
         if(iextinct.eq.1) then
           read(inline,*) iex,(extvvv(ie), ie = 0,20)
           ie0 = int(2.d0*(Dso-0.25d0))
           ie1 = ie0 + 1
           if(ie1.gt.20.or.ie0.lt.0) stop 'Dso out of range'
           f0 = 2.d0*(0.5d0*ie1+0.25d0-Dso)
           f1 = 1.d0-f0
           ext_Dso = f0*extvvv(ie0)+f1*(extvvv(ie1))
         else
           read(inline,*) iex,b_gal,h_dust
           scale_dust = h_dust/sind(abs(b_gal))
         endif
         write(6,*) 'enter 1 for integration grid, 0 to skip;',
     &      ' and name of optional MCMC output file'
         read(5,'(a)') line
         read(line(1:2),*) make_grid_in
         len_filep2 = lenc(line)
         if(len_filep2.gt.2) then
           read(line(3:len_filep2),'(a)') mcmc_file
           open(unit=10,file=mcmc_file,position='append')
           nchar_mcmc=24
         endif
         write(6,*) 'Source reference size & Distance =',theta_suasC,Dso
         if(ncband.gt.0) then
           do ic = 1,ncband
             flux_meas(ic) = 10.d0**(0.4d0*(21.d0-fMag_meas(ic)))
             flux_sig(ic) = flux_meas(ic)*sig_Mag_lens(ic)/1.0857d0
             write(6,*) cband(ic),'-band Mag =',fMag_meas(ic),' +-',
     &                  sig_Mag_lens(ic)
             write(6,*) cband(ic)//'-band extinction A_lens =',
     &                  A_lens(ic)
             if(Iup(ic).eq.1) then
               write(6,*) 'Upper-limit Magnitude constraint'
             elseif(Iup(ic).eq.2) then
               write(6,*) 'Full Magnitude constraint using total flux'
             elseif(Iup(ic).eq.3) then
               write(6,*)
     &           'Upper-limit Magnitude constraint using total flux'
             else
               write(6,*) 'Full Magnitude constraint'
             endif
           enddo
         endif

c        open and read parameter file
c        ----------------------------
         do 3 iop=1,3
           if(iop.eq.1) then
             parfile=fitname(1:lfn)//'.par'
           elseif(iop.eq.2) then
             parfile='par'//starname(1:lsn)
           elseif(iop.eq.3) then
             parfile='par.default'
           endif
           open(unit=3,file=parfile,status='old',err=2)
           write(6,*) parfile,' opened'
           go to 4
 2         write(6,*) parfile,' not found'
 3       continue
         write(6,*) 'no parameter file found'
         write(6,*) 'default parameters used'
 4       continue

         read(3,*)
         read(3,'(a)') line
         read(line,*,err=71,end=71) daycausmin,daycausmax,deltcaus,
     &                delfine,gridUstar,hcut,grid_rat,icc
         go to 77
 71      read(line,*,err=74,end=74) daycausmin,daycausmax,deltcaus,
     &                delfine,gridUstar,hcut,grid_rat
         go to 77
 74      read(line,*,err=76,end=76) daycausmin,daycausmax,deltcaus,
     &                delfine,gridUstar,hcut
         go to 77
 76      read(line,*) daycausmin,daycausmax,deltcaus,delfine,gridUstar
 77      continue
         write(6,*) 'hcut =',hcut,' grid_rat =',grid_rat,' icc =',icc
         read(3,*)
         do 5 jjclr=0,mxclr
           read(3,'(a)',err=6,end=6) line
           lenline = lenc(line)
           read(line,*) jclr
           if(lenline.gt.80) then
             read(line,*) idum,fudge(jclr),errmin(jclr),fmin(jclr),
     &          fmax(jclr),ald(jclr),bld(jclr),dayoff(jclr),sfx(jclr),
     &          lon_obs(jclr),lat_obs(jclr)
             write(6,*) sfx(jclr),' long., lat. =',
     &                  lon_obs(jclr),lat_obs(jclr)
             if(lat_obs(jclr).ge.1000.d0) then
               isat = nint(lat_obs(jclr)*0.001d0)
               satfile(jclr) = 'sat_position_'//c(isat)//'.dat'
             endif
           else
             read(line,*) idum,fudge(jclr),errmin(jclr),fmin(jclr),
     &          fmax(jclr),ald(jclr),bld(jclr),dayoff(jclr),sfx(jclr)
           endif
           lsfx(jclr)=lenc(sfx(jclr))
 5       continue
 111       format('jclr  fudge  errmin     fmin     fmax',
     &        '      ald      bld        dayoff       sfx')
 112       format(i3,f8.3,f8.4,2g11.4,2f8.4,f15.4,a7)
 6       continue
         nnclr = min(jjclr-1,mxclr)
         close(3)
         if(ncband.gt.0) then
           write(6,*) 'Mag constraints:'
           write(6,*) ('using Iclr =',icband(ic),sfx(icband(ic)),
     &                 ic = 1,ncband)
         endif
         write(6,111)
         do jclr=0,nnclr
           write(6,112) jclr,fudge(jclr),errmin(jclr),fmin(jclr),
     &        fmax(jclr),ald(jclr),bld(jclr),dayoff(jclr),sfx(jclr)
           if(lat_obs(jclr).ge.1000.d0) then
             open(unit=3,file=satfile(jclr),status='old')
             isat = nint(lat_obs(jclr)*0.001d0)
             do is = 1,max_sat
ccc               read(3,*,end=7,err=7) tsat(is,isat),
ccc     &               rasat(is,isat),decsat(is,isat),rsat(is,isat)
               read(3,*,end=7,err=7) tsat0,rasat0,decsat0,rsat0
               tsat(is,isat) = tsat0
               rasat(is,isat) = rasat0
               decsat(is,isat) = decsat0
               rsat(is,isat) = rsat0
             enddo
 7           continue
             ndatsat(isat) = is - 1
             write(6,*) 'read in ',ndatsat(isat),
     &                  ' satellite position points'
             close(3)
           endif
         enddo

c        assign input data file names
c        ----------------------------
         infile='lc'//starname(1:lsn)//'.'//sfx(0)(1:lsfx(0))

c        find the iteration number
c        -------------------------
         do 10 it=1,999
           if(it.lt.10) then
             newinffile=fitname(1:lfn)//c(it)//'.in'
             tmpfile='fit.lc_'//fitname(1:lfn)//c(it)
             tmpfile2='resid.'//fitname(1:lfn)//c(it)
           elseif(it.lt.100) then
             i1=mod(it,10)
             i2=it/10
             newinffile=fitname(1:lfn)//c(i2)//c(i1)//'.in'
             tmpfile='fit.lc_'//fitname(1:lfn)//c(i2)//c(i1)
             tmpfile2='resid.'//fitname(1:lfn)//c(i2)//c(i1)
           else
             i1=mod(it,10)
             i2=mod(it,100)/10
             i3=mod(it,1000)/100
             newinffile=fitname(1:lfn)//c(i3)//c(i2)//c(i1)//'.in'
             tmpfile=
     &         'fit.lc_'//fitname(1:lfn)//c(i3)//c(i2)//c(i1)
             tmpfile2=
     &         'resid.'//fitname(1:lfn)//c(i3)//c(i2)//c(i1)
           endif
           open(unit=8,file=newinffile,status='old',err=11)
           oldinffile=newinffile
           lcoutfile=tmpfile
           residfile=tmpfile2
           close(8)
 10      continue
 11      continue

c        determine min and max days
c        --------------------------
         daymin=1.d10
         daymax=-1.d10

         open(unit=1,file=infile,status='old',err=101)
         read(1,'(a)') cstar
         read(1,*)
         read(1,'(a)') head2
         read(1,*)
         read(1,*)
         read(1,*)
         read(1,*)
         jclr=0
         do 100 j=1,9999
            read(1,*,end=101,err=101) Date(j),ObsId,Lot,OK,Airmass,
     &         dxC35,dyC35,SodDate,SodVers,TmplVers,rMag(j),rErr(j),
     &         irType,rDSky,rChi2,rCrwd,rCosm,rMiss,rSeeing,rSky,rTmpl,
     &         bMag(j),bErr(j),ibType,bDSky,bChi2,bCrwd,bCosm,bMiss,
     &         bSeeing,bSky,bTmpl
            Date(j)=Date(j)+dayoff(jclr)
            if(Date(j).lt.daymin) daymin=Date(j)
            if(Date(j).gt.daymax) daymax=Date(j)
            irt=mod(irtype,10)
            if(rSeeing.le.seemax.and.irt.ne.7.and.irt.ne.8.and.
     &         rErr(j).le.errormax.and.abs(rDSky).le.idskymax.and.
     &         (rChi2.le.ipsfc2max.or.rMag(j).le.psfcmagmin).and.
     &         rCrwd.le.icrdmax.and.rCosm.le.icosmax.and.
     &         rMiss.le.imispmax.and.rErr(j).gt.0.d0) then

               ifield=Lot/1000
               ndata=ndata+1
               nclr(0)=nclr(0)+1
               x(ndata)=Date(j)
               y(ndata)=10**(-0.4*rMag(j))
               sig(ndata)=sqrt((rErr(j)/1.086)**2+errmin(0)**2)
               sig(ndata)=fudge(0)*sig(ndata)*y(ndata)
               xclr(nclr(0),0)=x(ndata)
               yclr(nclr(0),0)=y(ndata)
               eclr(nclr(0),0)=sig(ndata)
               iclrind(ndata)=nclr(0)
               iclr(ndata)=0
               rFlux(j)=y(ndata)
               rFerr(j)=sig(ndata)
            else
               rFlux(j)=0.
               rFerr(j)=9999.999
            endif
            ibt=mod(ibtype,10)
            if(bSeeing.le.seemax.and.ibt.ne.7.and.ibt.ne.8.and.
     &         bErr(j).le.errormax.and.abs(bDSky).le.idskymax.and.
     &         (bChi2.le.ipsfc2max.or.bMag(j).le.psfcmagmin).and.
     &         bCrwd.le.icrdmax.and.bCosm.le.icosmax.and.
     &         bMiss.le.imispmax.and.bErr(j).gt.0.) then

               ifield=Lot/1000
               ndata=ndata+1
               if(Date(j).lt.daymin) daymin=Date(j)
               if(Date(j).gt.daymax) daymax=Date(j)
               nclr(1)=nclr(1)+1
               x(ndata)=Date(j)
               y(ndata)=10**(-0.4*bMag(j))
               sig(ndata)=sqrt((bErr(j)/1.086)**2+errmin(1)**2)
               sig(ndata)=fudge(1)*sig(ndata)*y(ndata)
               xclr(nclr(1),1)=x(ndata)
               yclr(nclr(1),1)=y(ndata)
               eclr(nclr(1),1)=sig(ndata)
               iclrind(ndata)=nclr(1)
               iclr(ndata)=1
               bFlux(j)=y(ndata)
               bFerr(j)=sig(ndata)
            else
               bFlux(j)=0.
               bFerr(j)=9999.999
            endif
            if(bErr(j).lt.-99.999.or.bErr(j).gt.10.) bErr(j)=-99.999
            if(rErr(j).lt.-99.999.or.rErr(j).gt.10.) rErr(j)=-99.999
 100     continue
 101     continue
          close(1)
          if(j.gt.0) then
            jnext=j
          else
            jnext=1
          endif
          write(6,*) 'read in ',ndata,' measurements from '
     &         ,infile(1:20)

c        possibly read in the GMAN data
c        ------------------------------
         do 210 jclr=2,7
           ctinfile='lc'//starname(1:lsn)//'.'//sfx(jclr)(1:lsfx(jclr))
           open(unit=1,file=ctinfile,status='old',err=201)
           read(1,*)
           read(1,*)
           read(1,*)
           do 200 j=jnext,99999
             read(1,*,end=201,err=201) Date(j),rMag(j),rErr(j),sky,xpix,
     &           ypix,amass,xFWHM,yFWHM,chi2,sharp,dMag,dErr,dChi2,
     &           obsid,norms,expsr
               Date(j)=Date(j)+dayoff(jclr)
               if(dChi2.le.4.0.and.sharp.ge.-0.5.and.sharp.le.1.0.and.
     &            abs((xFWHM-yFWHM)/(xFWHM+yFWHM)).le.0.4.and.
     &          norms.ge.3) then
                 ndata=ndata+1
                 if(Date(j).lt.daymin) daymin=Date(j)
                 if(Date(j).gt.daymax) daymax=Date(j)
                 x(ndata)=Date(j)
                 y(ndata)=10**(-0.4*(rMag(j)-21.))
                 sig(ndata)=sqrt((rErr(j)/1.086)**2+(dErr/1.086)**2
     &               +errmin(jclr)**2)
                 sig(ndata)=sig(ndata)*fudge(jclr)*y(ndata)
                 iclr(ndata)=jclr
                 nclr(jclr)=nclr(jclr)+1
                 xclr(nclr(jclr),jclr)=x(ndata)
                 yclr(nclr(jclr),jclr)=y(ndata)
                 eclr(nclr(jclr),jclr)=sig(ndata)
                 iclrind(ndata)=nclr(jclr)
                 if(sig(ndata).gt.3*y(ndata)) then
                   ndata=ndata-1
                   nclr(jclr)=nclr(jclr)-1
                 endif
               endif
 200       continue
 201       continue
           close(1)
           write(6,*) 'read in GMAN data from ctinfile:'
           write(6,*) ctinfile(1:50)
           write(6,*) 'read in ',ndata,' total measurements'
 210     continue

           jclr=8
           ctinfile='lc'//starname(1:lsn)//'.'//sfx(jclr)(1:lsfx(jclr))
           open(unit=1,file=ctinfile,status='old',err=221)
           read(1,*)
           read(1,*)
           read(1,*)
           do 220 j=jnext,99999
             read(1,*,end=221,err=221) Date(j),rMag(j),rErr(j),sky,
     &           amass,xFWHM,yFWHM,rDSky,irType,rCrwd,rChi2,rMiss,
     &           rCosm,dMag,dMagC,dErr,dChi2,norms
             rSeeing=sqrt(xFWHM*yFWHM)
             Date(j)=Date(j)+dayoff(jclr)
             if(Date(j).lt.daymin) daymin=Date(j)
             if(Date(j).gt.daymax) daymax=Date(j)
             irt=mod(irtype,10)
             if(rSeeing.le.seemax74.and.irt.ne.7.and.irt.ne.8.and.
     &         rErr(j).le.errormax.and.abs(rDSky).le.idskymax.and.
     &         (rChi2.le.ipsfc2max74.or.rMag(j).le.psfcmagmin74).and.
     &         rCrwd.le.icrdmax.and.rCosm.le.icosmax.and.
     &         rMiss.le.imispmax.and.rErr(j).gt.0.d0.and.
     &         norms.gt.minorm) then

               ndata=ndata+1
               nclr(jclr)=nclr(jclr)+1
               x(ndata)=Date(j)
               y(ndata)=10**(-0.4*rMag(j))
               sig(ndata)=sqrt((rErr(j)/1.086)**2+errmin(jclr)**2)
               sig(ndata)=sig(ndata)*fudge(jclr)*y(ndata)
               xclr(nclr(jclr),jclr)=x(ndata)
               yclr(nclr(jclr),jclr)=y(ndata)
               eclr(nclr(jclr),jclr)=sig(ndata)
               iclrind(ndata)=nclr(jclr)
               iclr(ndata)=jclr
               rFlux(j)=y(ndata)
               rFerr(j)=sig(ndata)
             else
               rFlux(j)=0.
               rFerr(j)=9999.999
             endif
 220       continue
 221       continue
           close(1)
           write(6,*) 'read in GMAN data from ctinfile:'
           write(6,*) ctinfile(1:50)
           write(6,*) 'read in ',ndata,' total measurements'

         do 240 jclr=9,14
           ctinfile='lc'//starname(1:lsn)//'.'//sfx(jclr)(1:lsfx(jclr))
           open(unit=1,file=ctinfile,status='old',err=231)
           do 230 j=jnext,99999
             read(1,*,end=231,err=231) Date(j),ia,ib,itype,xxx,yyy,
     &           rMag(j),dm,rErr(j),chi2,dsky,c1,crd,fmiss,cr
               if(cr.lt.1.e-4.and.fmiss.lt.0.1.and.crd.lt.2..and.
     &                     rErr(j).ge.0..and.
     &                     (chi2.lt.1500..or.rMag(j).le.-6.5)) then
                 ndata=ndata+1
                 Date(j)=Date(j)+dayoff(jclr)
                 if(Date(j).lt.daymin) daymin=Date(j)
                 if(Date(j).gt.daymax) daymax=Date(j)
                 x(ndata)=Date(j)
                 y(ndata)=10**(-0.4*rMag(j))
                 sig(ndata)=sqrt((rErr(j)/1.086)**2+errmin(jclr)**2)
                 sig(ndata)=sig(ndata)*fudge(jclr)*y(ndata)
                 iclr(ndata)=jclr
                 nclr(jclr)=nclr(jclr)+1
                 xclr(nclr(jclr),jclr)=x(ndata)
                 yclr(nclr(jclr),jclr)=y(ndata)
                 eclr(nclr(jclr),jclr)=sig(ndata)
                 iclrind(ndata)=nclr(jclr)
                 if(sig(ndata).gt.3*y(ndata)) then
                   ndata=ndata-1
                   nclr(jclr)=nclr(jclr)-1
                 endif
               endif
 230       continue
 231       continue
           close(1)
           write(6,*) 'read in GMAN data from ctinfile:'
           write(6,*) ctinfile(1:50)
           write(6,*) 'read in ',ndata,' total measurements'
 240     continue

         do 255 jclr=15,29
           ctinfile='lc'//starname(1:lsn)//'.'//sfx(jclr)(1:lsfx(jclr))
           open(unit=1,file=ctinfile,status='old',err=251)
           do 250 j=jnext,99999
             read(1,*,end=251,err=251) Date(j),rMag(j),rErr(j)
                 Date(j)=Date(j)+dayoff(jclr)
                 ndata=ndata+1
                 if(Date(j).lt.daymin) daymin=Date(j)
                 if(Date(j).gt.daymax) daymax=Date(j)
                 x(ndata)=Date(j)
                 y(ndata)=10**(-0.4*(rMag(j)-21.))
                 sig(ndata)=sqrt((rErr(j)/1.086)**2+errmin(jclr)**2)
                 sig(ndata)=sig(ndata)*fudge(jclr)*y(ndata)
                 iclr(ndata)=jclr
                 nclr(jclr)=nclr(jclr)+1
                 xclr(nclr(jclr),jclr)=x(ndata)
                 yclr(nclr(jclr),jclr)=y(ndata)
                 eclr(nclr(jclr),jclr)=sig(ndata)
                 iclrind(ndata)=nclr(jclr)
 250       continue
 251       continue
           close(1)
           write(6,*) 'read in GMAN data from ctinfile:'
           write(6,*) ctinfile(1:50)
           write(6,*) 'read in ',ndata,' total measurements'
 255     continue

         do 265 jclr=30,39
           ctinfile='lc'//starname(1:lsn)//'.'//sfx(jclr)(1:lsfx(jclr))
           open(unit=1,file=ctinfile,status='old',err=261)
           datamin=1.d20
           datamax=-1.d20
           ndata0=ndata
           do 260 j=jnext,99999
             read(1,*,end=261,err=261) Date(j),rMag(j),rErr(j)
                 ndata=ndata+1
                 Date(j)=Date(j)+dayoff(jclr)
                 if(Date(j).lt.daymin) daymin=Date(j)
                 if(Date(j).gt.daymax) daymax=Date(j)
                 x(ndata)=Date(j)
                 y(ndata)=rMag(j)
                 datamin=min(datamin,y(ndata))
                 datamax=max(datamax,y(ndata))
                 sig(ndata)=sqrt(rErr(j)**2+(errmin(jclr)*rMag(j))**2)
                 sig(ndata)=sig(ndata)*fudge(jclr)
                 iclr(ndata)=jclr
                 nclr(jclr)=nclr(jclr)+1
                 xclr(nclr(jclr),jclr)=x(ndata)
                 yclr(nclr(jclr),jclr)=y(ndata)
                 eclr(nclr(jclr),jclr)=sig(ndata)
                 iclrind(ndata)=nclr(jclr)
 260       continue
 261       continue
           close(1)
           write(6,*) 'read in GMAN data from ctinfile:'
           write(6,*) ctinfile(1:50)
           write(6,*) 'read in ',ndata,' total measurements'
           write(6,*) 'datamin =',datamin,'  datamax =',datamax
           if(datamin.lt.0.) then
             write(6,*) 'adding',-datamin,' to each flux value'
             do idat=ndata0+1,ndata
               y(idat)=y(idat)-datamin
             enddo
           endif
 265     continue

         if(nnclr.lt.40) go to 286
         do 275 jclr=40,49
           ctinfile='lc'//starname(1:lsn)//'.'//sfx(jclr)(1:lsfx(jclr))
           open(unit=1,file=ctinfile,status='old',err=271)
           do 270 j=jnext,99999
             read(1,*,end=271,err=271) Date(j),rMag(j),rErr(j)
                 Date(j)=Date(j)+dayoff(jclr)
                 ndata=ndata+1
                 if(Date(j).lt.daymin) daymin=Date(j)
                 if(Date(j).gt.daymax) daymax=Date(j)
                 x(ndata)=Date(j)
                 y(ndata)=10**(-0.4*(rMag(j)-21.))
                 sig(ndata)=sqrt((rErr(j)/1.0857362)**2+errmin(jclr)**2)
                 sig(ndata)=sig(ndata)*fudge(jclr)*y(ndata)
                 iclr(ndata)=jclr
                 nclr(jclr)=nclr(jclr)+1
                 xclr(nclr(jclr),jclr)=x(ndata)
                 yclr(nclr(jclr),jclr)=y(ndata)
                 eclr(nclr(jclr),jclr)=sig(ndata)
                 iclrind(ndata)=nclr(jclr)
 270       continue
 271       continue
           close(1)
           write(6,*) 'read in GMAN data from ctinfile:'
           write(6,*) ctinfile(1:50)
           write(6,*) 'read in ',ndata,' total measurements'
 275     continue

         if(nnclr.lt.50) go to 286
         do 285 jclr=50,59
           ctinfile='lc'//starname(1:lsn)//'.'//sfx(jclr)(1:lsfx(jclr))
           open(unit=1,file=ctinfile,status='old',err=281)
           datamin=1.d20
           datamax=-1.d20
           ndata0=ndata
           do 280 j=jnext,99999
             read(1,*,end=281,err=281) Date(j),rMag(j),rErr(j)
                 ndata=ndata+1
                 Date(j)=Date(j)+dayoff(jclr)
                 if(Date(j).lt.daymin) daymin=Date(j)
                 if(Date(j).gt.daymax) daymax=Date(j)
                 x(ndata)=Date(j)
                 y(ndata)=rMag(j)
                 datamin=min(datamin,y(ndata))
                 datamax=max(datamax,y(ndata))
                 sig(ndata)=sqrt(rErr(j)**2+(errmin(jclr)*rMag(j))**2)
                 sig(ndata)=sig(ndata)*fudge(jclr)
                 iclr(ndata)=jclr
                 nclr(jclr)=nclr(jclr)+1
                 xclr(nclr(jclr),jclr)=x(ndata)
                 yclr(nclr(jclr),jclr)=y(ndata)
                 eclr(nclr(jclr),jclr)=sig(ndata)
                 iclrind(ndata)=nclr(jclr)
 280       continue
 281       continue
           close(1)
           write(6,*) 'read in GMAN data from ctinfile:'
           write(6,*) ctinfile(1:50)
           write(6,*) 'read in ',ndata,' total measurements'
           write(6,*) 'datamin =',datamin,'  datamax =',datamax
           if(datamin.lt.0.) then
             write(6,*) 'adding',-datamin,' to each flux value'
             do idat=ndata0+1,ndata
               y(idat)=y(idat)-datamin
             enddo
           endif
 285     continue
 286     continue
cccc        sort the data by time
cccc        ---------------------
ccc         call sort4(ndata,x,y,sig,iclr,work,iwork,jwork)

c        sort the data by time
c        ---------------------
         call sort5(ndata,x,y,sig,iclr,iclrind,work,iwork,jwork)

c        make the limb darkening table
c        -----------------------------
         call limb_dark()
         call hexlimbfac()

       endif

c      specify triple lens fit parameters
c      ----------------------------------
       t_E = 1.d0/a(1)
       sep1cm = a(4)
       sep23_0 = a(5)
       theta1cm = a(6)
       theta = theta1cm
       eps1 = a(7)
       eps2 = a(8)
       ang23_0 = a(9)
       ds23xdt = a(10)
       ds23ydt = a(11)
       Rstar = a(12)
       Tstar = Rstar
       t_fix = a(13)
       tfix = t_fix
       T_bin_inv = a(14)
       dt_orb = a(15)
       pier = a(22)
       pietheta = a(23)
       ang23 = ang23_0
       sep23 = sep23_0

       eps1 = min(max(eps1,0.d0),1.d0)
       eps2 = min(max(eps2,0.d0),1.d0)
       eps3 = 1.d0-eps2-eps1
       eps4 = eps2+eps3
       cosang23_0 = cos(ang23)
       sinang23_0 = sin(ang23)
       sep23_x0 = cosang23_0*sep23
       sep23_y0 = sinang23_0*sep23
       fm2 = eps2/eps4
       fm3 = eps3/eps4

       w = twopi*T_bin_inv
       if(T_bin_inv.ne.0.d0) then
         T_bin = 1.d0/T_bin_inv
         sep23_xmx = sqrt(sep23_x0**2+(ds23xdt/w)**2)
         sep23_ymx = sqrt(sep23_y0**2+(ds23ydt/w)**2)
         tan_wtx = -ds23xdt/(w*sep23_x0)
         wtx = atan(tan_wtx)
         if(sep23_x0.lt.0.d0) wtx = wtx + pi
         dt23x = wtx/w
         tan_wty = -ds23ydt/(w*sep23_y0)
         wty = atan(tan_wty)
         if(sep23_y0.lt.0.d0) wty = wty + pi
         dt23y = wty/w
       else
         T_bin = 1.d20
         sep23_xmx = sep23_x0
         sep23_ymx = sep23_y0
       endif

       write(6,*) 'FCN call with a =',(a(k), k=1,23)

       chi2=0.
       r0=0.
       maximages=0
       do 290 i=0,nnclr
         chi2clr(i)=0.
 290   continue

c  initialize source grid storage
c  ------------------------------
       Ustar=a(12)*a(1)
       bgriddel0=gridUstar*Ustar
       ngfac=nint(2.d0*log(1.d5*bgriddel0)/log(2.d0))
       bgriddel=1.d-5*sqrt(2.d0**ngfac)
       nphimax=int(pi/(grid_rat*bgriddel))
       bphigrid=twopi/(2*nphimax+1)
       brgrid=bgriddel
       ngrtot=ngmax/(2*nphimax+1)
       ngr=(ngrtot-1)/2

         Ein_R = 1.d0
         xcc = 0.d0

       ang = ang23
       sep2 = sep23
       sep = sep1cm
       if(make_grid.gt.0.and.T_bin_inv.eq.0.d0.and.dt_orb.gt.0.d0) then
         call microcurve_init(a,Ein_R,xcc,brgrid,bphigrid,
     &           sep1cm,sep23_0,ang23_0,theta,ngr,nphimax,sxg,syg)
         ngrchk=ngr
       else
         ngrchk = -1
       endif

c  sum over data points
c  --------------------
       Del_t_orb0 = -1.e20
       do 300 i=1,ndata
          t=x(i)
          ii=iclrind(i)
          if(dt_orb.gt.0.d0) then
            Del_t_orb = dt_orb*nint((t-t_fix)/dt_orb)
          else
            Del_t_orb=t-t_fix
          endif
c         move the masses only at dt_orb intervals
          if(Del_t_orb.ne.Del_t_orb0.and.T_bin_inv.ne.0.d0) then
            Del_t_orb0 = Del_t_orb
            sep23_x = sep23_xmx*cos(w*(Del_t_orb+dt23x))
            sep23_y = sep23_ymx*cos(w*(Del_t_orb+dt23y))
            sep23 = sqrt(sep23_x**2+sep23_y**2)
            if(sep23_x.lt.0.d0) then
              ang23 = atan(sep23_y/sep23_x) + pi
            elseif(sep23_x.gt.0.d0) then
              ang23 = atan(sep23_y/sep23_x)
            else
              if(sep23_y.ge.0.d0) then
                ang23 = 0.5d0*pi
              else
                ang23 = -0.5d0*pi
              endif
            endif
            ang = ang23
            sep2 = sep23
            sep = sep1cm
            ngrchk=-1
          endif

          if(make_grid.eq.1.and.ngrchk.le.0.and.dt_orb.gt.0.d0) then
            call microcurve_init(a,Ein_R,xcc,brgrid,bphigrid,
     &                      sep,sep2,ang,theta,ngr,nphimax,sxg,syg)
            ngrchk=ngr
          endif

          jc = iclr(i)
          call microcurve(t,a,yfit(ii,jc),jc,nimage(i),
     &                    sep,sep2,ang,theta,alpha,delta,Ein_R,xcc,
ccc     &                    lon_obs(jc),lat_obs(jc),1,
     &                    lon_obs(jc),lat_obs(jc),0,
     &                    brgrid,bphigrid,ngr,nphimax,sxg,syg,
     &                    ngrchk,make_gridm,icc)
          make_grid=make_gridm*make_grid_in
          maximages=max(maximages,nimage(i))
          if(nclr(iclr(i)).eq.ii) then
            do ic = 1,ncband
              if(icband(ic).eq.iclr(i)) A_last_obs(ic) = yfit(ii,jc)
            enddo
          endif
 300   continue
       if(nchar_mcmc.gt.0) nchar_mcmc = 24
       do 350 j=0,nnclr
          if(nclr(j).gt.0) then
             call normfit(yfit(1,j),yclr(1,j),nclr(j),eclr(1,j),wfm,
     &             A2(j),A0(j),A2sig(j),A0sig(j),chi,fmin(j),fmax(j),
     &             j,iflag,xclr(1,j))
             chi2clr(j)=chi
             chi2=chi2+chi
             if(nchar_mcmc.gt.0) then
               nchar_mcmc = nchar_mcmc + 1
               write(char_mcmc(nchar_mcmc),998) A0(j)
               nchar_mcmc = nchar_mcmc + 1
               write(char_mcmc(nchar_mcmc),998) A2(j)
             endif
          else
             chi2clr(j)=0.
             A0(j)=0.
             A2(j)=0.
          endif
 350   continue

c      correct theta_suas for source brightness
c      ----------------------------------------
       theta_suas = theta_suasC*sqrt(A0(iflux0)/flux0)

       call orb_mass(cc,a,fmass,Dso,Ds_sig,Ds_calc,Dl_calc,
     &      theta_suas,fMmax,fMsig)

c      compute mass and magnitude from parallax and t_*
c      ------------------------------------------------
       thetaEmas = 0.001d0*theta_suas*(t_E/Tstar)
       fmass_par = thetaEmas/(8.1438533d0*abs(pier))
       write(6,*) 'thetaE =',thetaEmas,' mas'
       write(6,*) 'Mass_par =',fmass_par,' M_solar'
       fmass_par1 = eps1*fmass_par
       fmass_par2 = eps2*fmass_par
       fmass_par3 = eps3*fmass_par
       DM = 10.d0 + 5.d0*dlog10(Dl_calc)

c      estimate the extinction in the lens foreground
c      ----------------------------------------------
       if(iextinct.eq.1) then
         ie0 = int(2.d0*(Dl_calc-0.25d0))
         ie1 = ie0 + 1
         if(ie1.le.0) then
           ext_Dl = extvvv(0)
         elseif(ie0.ge.20) then
           ext_Dl = extvvv(20)
         else
           f0 = 2.d0*(0.5d0*ie1+0.25d0-Dl_calc)
           f1 = 1.d0-f0
           ext_Dl = f0*extvvv(ie0)+f1*(extvvv(ie1))
         endif
         ext_fac = ext_Dl/ext_Dso
       else
         ext_fac = (1.d0-exp(-Dl_calc/scale_dust))/
     &             (1.d0-exp(-Dso/scale_dust))
       endif
       chi2mag = 0.d0
       if(ncband.gt.0) then
         do ic = 1,ncband
           if(fmass_par1.gt.0.06d0) then
             absMag_lens1 = absmag(cband(ic),fmass_par1)
           else
             absMag_lens1 = 30.d0
           endif
           if(fmass_par2.gt.0.06d0) then
             absMag_lens2 = absmag(cband(ic),fmass_par2)
           else
             absMag_lens2 = 30.d0
           endif
           if(fmass_par3.gt.0.06d0) then
             absMag_lens3 = absmag(cband(ic),fmass_par3)
           else
             absMag_lens3 = 30.d0
           endif
           aflux1 = 10**(0.4d0*(21.d0-absMag_lens1))
           aflux2 = 10**(0.4d0*(21.d0-absMag_lens2))
           aflux3 = 10**(0.4d0*(21.d0-absMag_lens3))
           aflux = aflux1 + aflux2 + aflux3
           absMag_lens = 21.d0 - 2.5d0*dlog10(aflux)

           fMag_lens = DM + absMag_lens + ext_fac*A_lens(ic)
           flux_lens = 10.d0**(0.4d0*(21.d0-fMag_lens))

           write(6,777) cband(ic),fMag_lens

           if(Iup(ic).gt.1) then
             flux_src = A0(icband(ic))*10**(-0.4d0*fMag_off(ic))
             fI_src0 = 21.d0-2.5d0*dlog10(flux_src)
             flux_src = flux_src*A_last_obs(ic)
             fI_src = 21.d0-2.5d0*dlog10(flux_src)
             flux_fit = flux_lens + flux_src
             fI_tot = 21.d0-2.5d0*dlog10(flux_fit)
             write(6,776) cband(ic)//'_src0 =',fI_src0,
     &           cband(ic)//'_src =',fI_src,cband(ic)//'_tot =',fI_tot
           elseif(Iup(ic).gt.0) then
             flux_fit = A2(icband(ic))
             if(flux_fit.le.0.d0) flux_fit=0.001d0
           else
             flux_fit = flux_lens
           endif
           if(flux_fit.gt.flux_meas(ic).or.Iup(ic).le.0.or.Iup(ic).eq.2)
     &            then
              chi2m = ((flux_fit-flux_meas(ic))/flux_sig(ic))**2
              chi2mag = chi2mag + chi2m
              write(6,*) 'chi2_'//cband(ic)//' =',chi2m
           endif
         enddo
       endif
 776   format(a8,f9.5,2(3x,a7,f9.5))
 777   format(a1,'-band lens Mag =',f10.6)

       chi2lc = chi2

c      include the D_source constraint
c      -------------------------------
       chi2_Ds = ((Ds_calc-Dso)/Ds_sig)**2
       chi2 = chi2lc + chi2_Ds

c      include the Lens Brightness Constraint
c      --------------------------------------
       chi2 = chi2 + chi2mag

       write(6,*) ' chi2mag =',chi2mag,'  chi2_Ds =',chi2_Ds
       write(6,*) ' chi2 =',chi2,'  chi2lc =',chi2lc
       if(nchar_mcmc.gt.0) then
         write(char_mcmc(1),998) chi2
         do ipar=1,23
           write(char_mcmc(ipar+1),998) a(ipar)
         enddo
       endif
       call flush(6)
       call flush(10)

       if(iflag.eq.3) then
         open(unit=9,file=residfile)
         write(9,555) chi2,a(10),a(11),1./a(1),a(2),a(3),a(4),
     &         a(5),a(6),a(7),a(8),1.-a(7)-a(8),a(9),a(12),a(15),
     &         a(14),T_bin,a(16),a(17),a(18),a(19),a(20),a(21),a(13),
     &         a(22),a(23),alpha,delta,theta_suas,Dso,Ds_sig,Ds_calc
         if(nnclr.lt.40) then
           write(9,558) (sfx(i),sfx(i),sfx(i+1),sfx(i+1),
     &                 sfx(i+2),sfx(i+2),sfx(i+3),sfx(i+3),
     &                 A0(i),A2(i),A0(i+1),A2(i+1),
     &                 A0(i+2),A2(i+2),A0(i+3),A2(i+3), i=3,35,4)
         elseif(nnclr.lt.50) then
           write(9,558) (sfx(i),sfx(i),sfx(i+1),sfx(i+1),
     &                 sfx(i+2),sfx(i+2),sfx(i+3),sfx(i+3),
     &                 A0(i),A2(i),A0(i+1),A2(i+1),
     &                 A0(i+2),A2(i+2),A0(i+3),A2(i+3), i=3,45,4)
         else
           write(9,558) (sfx(i),sfx(i),sfx(i+1),sfx(i+1),
     &                 sfx(i+2),sfx(i+2),sfx(i+3),sfx(i+3),
     &                 A0(i),A2(i),A0(i+1),A2(i+1),
     &                 A0(i+2),A2(i+2),A0(i+3),A2(i+3), i=3,55,4)
         endif
         nclrinc=0
         do j=0,nnclr
           if(A0(j).ne.0.) then
             nclrinc=nclrinc+1
             jclrinc(nclrinc)=j
           endif
         enddo

         write(9,565)
         do 360 i=1,ndata
           t=x(i)
           ii=iclrind(i)
           ampdata=(yclr(ii,iclr(i))-A2(iclr(i)))/A0(iclr(i))
           resid = (yfit(ii,iclr(i))-ampdata)
           resid_sig = eclr(ii,iclr(i))/A0(iclr(i))
           chi2resid = (resid/resid_sig)**2
           write(9,566) t,yfit(ii,iclr(i)),resid,resid_sig,
     &                  chi2resid,iclr(i),sfx(iclr(i))
 360     continue
 565     format('       t          Fit-Ampl     Residual    Resid-Sig',
     &          '      chi^2    Iclr')
 566     format(f13.6,g15.6,2f13.7,g13.5,i4,1x,a6)

         do 410 i=0,nnclr
           if(nclr(i).gt.0) write(6,553) i,sfx(i),chi2clr(i),nclr(i)
           if(nclr(i).gt.0) write(9,553) i,sfx(i),chi2clr(i),nclr(i)
 410     continue
 553     format('color #',i3,1x,a6,1x,' chi2 =',g12.6,' for',
     &          i6,' points')
         close(9)

         write(6,*)
         write(6,559)
         if(A0(0)+A2(0).ne.0.)
     &        write(6,560) A0(0),A0sig(0),A2(0),A2sig(0)
         if(A0(1)+A2(1).ne.0.)
     &        write(6,561) A0(1),A0sig(1),A2(1),A2sig(1)
         do 415 i=2,nnclr
           if(A0(i)+A2(i).ne.0.) 
     &       write(6,562) sfx(i),A0(i),A0sig(i),sfx(i),A2(i),A2sig(i)      
 415     continue

 559     format('Normalization parameters with pseudo-errors follow:')
 560     format('A0rm     =',f12.4,' +/-',f11.4,
     &         '  A2rm     =',f12.4,' +/-',f11.4)
 561     format('A0bm     =',f12.4,' +/-',f11.4,
     &         '  A2bm     =',f12.4,' +/-',f11.4)
 562     format('A0',a6,' =',f12.4,' +/-',f11.4,
     &         '  A2',a6,' =',f12.4,' +/-',f11.4)

cccc        determine caustic crossing times
cccc        --------------------------------
         umin=a(3)
         ctheta=cos(theta)
         stheta=sin(theta)
ccc         nim=3
ccc         nim00=3
ccc         ncaus=0
         ds=0.001
         nstep=4000
ccc         dst2=ds*2.
ccc         dst22=dst2**2
ccc         dsx=ds*ctheta
ccc         dsy=ds*stheta
         cosfac=cos(thetaproj)
         sinfac=sin(thetaproj)
ccc         call trilens(eps2,eps3,sep,sep2,ang)
ccc         do 400 i=-nstep,nstep,1
ccc           vfac=ds*i
ccc           sx_norot=ctheta*vfac-stheta*a(3)
ccc           sy_norot=stheta*vfac+ctheta*a(3)
ccc           sx=cosfac*sx_norot-sinfac*sy_norot
ccc           sy=cosfac*sy_norot+sinfac*sx_norot
ccc
ccc 393       nim0=nim
ccc           call trilens_im(sx,sy,nim,iimage,z,ampim)
ccc          if(nim.ne.nim0) then
ccc             if(nim0.eq.3.or.nim0.eq.5) s0=ds*(i-1)
ccc             if(nim.eq.3.or.nim.eq.5) then
ccc               s1=vfac
ccc               do 395 j=1,20
ccc                 s=0.5*(s0+s1)
ccc                 sx_norot=ctheta*s-stheta*a(3)
ccc                 sy_norot=stheta*s+ctheta*a(3)
ccc                 sx=cosfac*sx_norot-sinfac*sy_norot
ccc                 sy=cosfac*sy_norot+sinfac*sx_norot
ccc                 call trilens_im(sx,sy,nim_mid,iimage,z,ampim)
ccc                 if(nim_mid.eq.nim) then
ccc                   s1=s
ccc                 elseif(nim_mid.eq.nim00) then
ccc                   s0=s
ccc                 else
ccc                   go to 396
ccc                 endif
ccc 395           continue
ccc 396           continue
ccc               ncaus=ncaus+1
ccc               scaus(ncaus)=s
ccc               tcaus(ncaus)=a(2)+scaus(ncaus)/a(1)
ccc               nim00=nim
ccc             endif
ccc           endif
ccc 400     continue
ccc
ccc         if(ncaus.gt.1) then
ccc           write(6,*) 'Caustic crossings found at',ncaus,' times:'
ccc           write(6,552) (tcaus(i),i=1,ncaus)
ccc         endif
ccc 552     format(6f12.5)
         do 400 i=-nstep,nstep,2*nstep
           vfac=ds*i
           sx_norot=ctheta*vfac-stheta*a(3)
           sy_norot=stheta*vfac+ctheta*a(3)
           sx=cosfac*sx_norot-sinfac*sy_norot
           sy=cosfac*sy_norot+sinfac*sx_norot
           if(i.eq.-nstep) then
             sx0=sx
             sy0=sy
           elseif(i.eq.nstep) then
             sx1=sx
             sy1=sy
           endif
 400     continue

         write(6,555) chi2,a(10),a(11),1./a(1),a(2),a(3),a(4),
     &         a(5),a(6),a(7),a(8),1.-a(7)-a(8),a(9),a(12),a(15),
     &         a(14),T_bin,a(16),a(17),a(18),a(19),a(20),a(21),a(13),
     &         a(22),a(23),alpha,delta,theta_suas,Dso,Ds_sig,Ds_calc
         if(nnclr.lt.40) then
           write(2,558) (sfx(i),sfx(i),sfx(i+1),sfx(i+1),
     &                 sfx(i+2),sfx(i+2),sfx(i+3),sfx(i+3),
     &                 A0(i),A2(i),A0(i+1),A2(i+1),
     &                 A0(i+2),A2(i+2),A0(i+3),A2(i+3), i=3,35,4)
         elseif(nnclr.lt.50) then
           write(2,558) (sfx(i),sfx(i),sfx(i+1),sfx(i+1),
     &                 sfx(i+2),sfx(i+2),sfx(i+3),sfx(i+3),
     &                 A0(i),A2(i),A0(i+1),A2(i+1),
     &                 A0(i+2),A2(i+2),A0(i+3),A2(i+3), i=3,45,4)
         else
           write(2,558) (sfx(i),sfx(i),sfx(i+1),sfx(i+1),
     &                 sfx(i+2),sfx(i+2),sfx(i+3),sfx(i+3),
     &                 A0(i),A2(i),A0(i+1),A2(i+1),
     &                 A0(i+2),A2(i+2),A0(i+3),A2(i+3), i=3,55,4)
         endif

 555     format('  chisq   ds23xdt   ds23ydt    t_E       t0     umin',
     &    '     sep1cm  sep23 theta1cm',/,f8.2,2f9.5,f9.3,f10.4,f10.6,
     &    3f8.4,/,'    eps1       eps2       eps3    ang23    Tstar   ',
     &    ' dt_orb   1/T_bin   T_bin',/,3(1pe11.4),0pf9.5,f9.6,f8.4,
     &     f10.6,f9.2,/,'tcch1_min tcch1_max tcch2_min tcch2_max ',
     &       'tcch3_min tcch3_max  t_fix',/,7f10.3,/,
     &     '    pi_E  theta_piE  alpha    delta  theta_suas  Dso     ',
     &     'Ds_sig  Ds_calc',/,f10.6,f9.5,2f9.3,f9.5,3f9.4)

 558     format(4('    A0',a6,'    A2',a6),/,8f12.3)

c        write out fit lightcurve
c        ------------------------
         dt=abs(1./a(1))
         t0=a(2)
         umin=a(3)
         r0=0.
         delc=1.
         if(delfine.le.0.) delfine=0.1
         dayminc1=daymin-50.
         daymaxc1=t0-1.5*max(dt,r0)
         dayminc2=t0+1.5*max(dt,r0)
         daymaxc2=daymax+200.
         open(unit=2,file=lcoutfile,status='new')
         write(2,555) chi2,a(10),a(11),1./a(1),a(2),a(3),a(4),
     &         a(5),a(6),a(7),a(8),1.-a(7)-a(8),a(9),a(12),a(15),
     &         a(14),T_bin,a(16),a(17),a(18),a(19),a(20),a(21),a(13),
     &         a(22),a(23),alpha,delta,theta_suas,Dso,Ds_sig,Ds_calc
         if(nnclr.lt.40) then
           write(2,558) (sfx(i),sfx(i),sfx(i+1),sfx(i+1),
     &                 sfx(i+2),sfx(i+2),sfx(i+3),sfx(i+3),
     &                 A0(i),A2(i),A0(i+1),A2(i+1),
     &                 A0(i+2),A2(i+2),A0(i+3),A2(i+3), i=3,35,4)
         elseif(nnclr.lt.50) then
           write(2,558) (sfx(i),sfx(i),sfx(i+1),sfx(i+1),
     &                 sfx(i+2),sfx(i+2),sfx(i+3),sfx(i+3),
     &                 A0(i),A2(i),A0(i+1),A2(i+1),
     &                 A0(i+2),A2(i+2),A0(i+3),A2(i+3), i=3,45,4)
         else
           write(2,558) (sfx(i),sfx(i),sfx(i+1),sfx(i+1),
     &                 sfx(i+2),sfx(i+2),sfx(i+3),sfx(i+3),
     &                 A0(i),A2(i),A0(i+1),A2(i+1),
     &                 A0(i+2),A2(i+2),A0(i+3),A2(i+3), i=3,55,4)
         endif
         nclrinc=0
         do 420 j=0,nnclr
           if(A0(j).ne.0.) then
             nclrinc=nclrinc+1
             jclrinc(nclrinc)=j
           endif
 420     continue

         t=dayminc1
         Del_t_orb0 = -1.e20
         ngrchk = -1
         do 460 i=1,9999999
           if(i.gt.1) then
             if(deltcaus.gt.0..and.t.ge.daycausmin
     &                        .and.t.le.daycausmax) then
               delt=deltcaus
             elseif(t.ge.daymaxc1.and.t.le.dayminc2) then
               delt=delfine
             else
               delt=delc
             endif
             t=t+delt
             if(t.gt.daymaxc2) go to 461
           endif

           if(dt_orb.gt.0.d0) then
             Del_t_orb = dt_orb*nint((t-t_fix)/dt_orb)
           else
             Del_t_orb=t-t_fix
           endif
c          move the masses only at dt_orb intervals
           if(Del_t_orb.ne.Del_t_orb0.and.T_bin_inv.ne.0.d0) then
             Del_t_orb0 = Del_t_orb
             sep23_x = sep23_xmx*cos(w*(Del_t_orb+dt23x))
             sep23_y = sep23_ymx*cos(w*(Del_t_orb+dt23y))
             sep23 = sqrt(sep23_x**2+sep23_y**2)
             if(sep23_x.lt.0.d0) then
               ang23 = atan(sep23_y/sep23_x) + pi
             elseif(sep23_x.gt.0.d0) then
               ang23 = atan(sep23_y/sep23_x)
             else
               if(sep23_y.ge.0.d0) then
                 ang23 = 0.5d0*pi
               else
                 ang23 = -0.5d0*pi
               endif
             endif
             ang = ang23
             sep2 = sep23
             sep = sep1cm
             ngrchk=-1
           endif
           if(make_grid.eq.1.and.ngrchk.le.0.and.
     &                             dt_orb.gt.0.d0) then
             call microcurve_init(a,Ein_R,xcc,brgrid,bphigrid,
     &                      sep,sep2,ang,theta,ngr,nphimax,sxg,syg)
             ngrchk=ngr
           endif
           do 450 j=0,nnclr
             if(A0(j).ne.0.) then
               if(j.gt.0) then
                 do 440 jj=0,j-1
                   if(ald(j).eq.ald(jj).and.bld(j).eq.bld(jj).and.
     &                A0(jj).gt.0.) then
                     yfits(j)=yfits(jj)
                     yyfits(j)=A0(j)*yfits(jj)+A2(j)
                     go to 450
                   endif
 440             continue
               endif
               call microcurve(t,a,yfits(j),j,nn,sep,sep2,ang,theta,
     &                  alpha,delta,Ein_R,xcc,lon_obs(j),lat_obs(j),0,
     &                  brgrid,bphigrid,ngr,nphimax,sxg,syg,
     &                  ngrchk,make_gridm,icc)
               yyfits(j)=A0(j)*yfits(j)+A2(j)
               make_grid=make_gridm*make_grid_in
             else
               yyfits(j)=A2(j)
             endif
 450       continue
           if(A0(32).ne.0.d0) then
             ilc = 32
           elseif(A0(16).ne.0.d0) then
             ilc = 16
           elseif(A0(22).ne.0.d0) then
             ilc = 22
           else
             do j=1,nclrinc
               if(A0(jclrinc(j)).ne.0.d0) then
                 ilc = jclrinc(j)
                 go to 457
               endif
             enddo
           endif
 457       continue
           write(2,556) t,yfits(ilc),(yyfits(jclrinc(j)),j=1,nclrinc)
 460     continue
 461     continue
 556     format(f9.4,f11.4,39f14.4)
         close(2)

         open(unit=11,file=oldinffile,status='old')
         open(unit=4,file=newinffile,status='new')
         read(11,'(a)') line
         len=lenc(line)
         write(4,'(a)') line(1:len)
         do i=1,23
           read(11,'(a)') line
           len=lenc(line)
           outline = blank20//blank20//blank20//blank20
           k1 = index(line,'''')
           if(k1.eq.0) then
             go to 480
           else
             k2 = k1 + index(line(k1+1:),'''')
           endif
           read(line(1:k1-1),*) iii
           istart = 0
           do j=k2+1,80
             if(line(j:j).ne.' ') then
               istart = 1
             elseif(istart.eq.1) then
               go to 475
             endif
           enddo
 475       continue
           write(outline(1:20),'(a)') line(1:k2+1)
           write(outline(21:45),999) a(i)
           write(outline(46:80),'(a)') line(j:len)
           leno = lenc(outline)
           write(4,'(a)') outline(1:leno)
         enddo
 480     continue
         read(11,*)
         write(4,*)
         do i=1,999
           read(11,'(a)',end=485,err=485) line
           len=lenc(line)
           write(4,'(a)') line(1:len)
         enddo
 485     continue
         close(4)
         close(11)
       endif
 998   format(1x,g18.11)
 999   format(g18.11)

       return
       end

      SUBROUTINE normFIT(X,Y,NDATA,SIG,w,A,B,SIGA,SIGB,CHI2,fmin,fmax,
     &           jclr,iflag,xclr)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision X(NDATA),Y(NDATA),SIG(NDATA),xclr(ndata)

      dmn=(1./fmax-1.)
      dmx=(1./fmin-1.)
      wdmn=w*dmn
      wdmx=w*dmx
      SX=0.
      SY=0.
      ST2=0.
      B=0.
      SS=0.
      sxpdmn=0.
      sxpdmx=0.
      DO 11 I=1,NDATA
        WT=1./(SIG(I)**2)
        SS=SS+WT
        xwt=X(I)*WT
        SX=SX+xwt
        SY=SY+Y(I)*WT
        sxpdmn=sxpdmn+xwt-wdmn*wt
        sxpdmx=sxpdmx+xwt-wdmx*wt
11    CONTINUE
      SXOSS=SX/SS
      DO 13 I=1,NDATA
        T=(X(I)-SXOSS)/SIG(I)
        ST2=ST2+T*T
        B=B+T*Y(I)/SIG(I)
13    CONTINUE
      b0=b
      if(st2.eq.0.) then
        b=0.
        a=sy/ss
        st2=1.e-20
      else
        B=B/ST2
        A=(SY-SX*B)/SS
        f=b/(b+a)
        if(f.lt.fmin.or.f.gt.fmax) then
          ssp=(1+w)*ss
          if(f.lt.fmin) then
            dm=dmx
            sxp=sxpdmx
          else
            dm=dmn
            sxp=sxpdmn
          endif
          sxpossp=sxp/ssp
          b=0.
          do 14 i=1,ndata
            t=(x(i)-sxpossp)/sig(i)
            b=b+t*y(i)/sig(i)
14        continue
          st2p=st2+(w/ssp)*(sx+dm*ss)**2
          b=b/st2p
          a=(sy-sxp*b)/ssp
        endif
      endif
      SIGA=SQRT((1.+SX*SX/(SS*ST2))/SS)
      SIGB=SQRT(1./ST2)
      CHI2=0.
      if(iflag.eq.3) then
        write(6,21)
      endif
      DO 16 I=1,NDATA
        CHI=((Y(I)-A-B*X(I))/SIG(I))**2
        CHI2=CHI2+chi
        if(iflag.eq.3) then
          write(6,22) xclr(i),chi,y(i),sig(i),jclr
        endif
16    CONTINUE
21    format(' date        chi2      amp        error    color')
22    format(f10.4,f9.3,2f11.3,i4)
      RETURN
      END

       subroutine orb_mass(c,a,fmass,Dso,Ds_sig,Ds_calc,Dl_calc,
     &                     theta_suas,fMmax,fMsig)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       parameter (mmax=30)

       double precision a(mmax),xrt(3),fmass(2),Dsrt(3)

c      calculates mass and other lens parameters from orbital motion
c      triple lens fit parameters.

c      units:
c      sep23_mx, s3 : dimensionless
c      T_bin : days
c      thetaE : radians
c      theta_st : radians
c      theta_suas : micro arc sec
c      Dso : kpc
c      Ds_calc : kpc

       pi = acos(-1.d0)
       twopi = 2.d0*pi

c        specify triple lens fit parameters
c        ----------------------------------
         thats2 = 1.d0/a(1)
         t_E = thats2
         t0 = a(2)
         umin = a(3)
         sep1cm = a(4)
         sep23_0 = a(5)
         theta1cm = a(6)
         theta = theta1cm
         eps1 = a(7)
         eps2 = a(8)
         ang23_0 = a(9)
         ds23xdt = a(10)
         ds23ydt = a(11)
         Rstar = a(12)
         t_fix = a(13)
         tfix = t_fix
         T_bin_inv = a(14)
         pier = a(22)
         pietheta = a(23)
         piex = pier*cos(pietheta)
         piey = pier*sin(pietheta)

         dt_orb = a(15)
         ang23 = ang23_0
         sep23 = sep23_0

       theta_st = theta_suas*pi/(180.d0*1.d6*3600.d0)
       tstar = Rstar

       thetaE = t_E*theta_st/tstar
       thetaEmas = 0.001d0*t_E*theta_suas/tstar

       eps3 = 1.d0-eps1-eps2
       epsorb = eps2+eps3
       eps4 = eps2+eps3
       cosang23_0 = cos(ang23)
       sinang23_0 = sin(ang23)
       sep23_x0 = cosang23_0*sep23
       sep23_y0 = sinang23_0*sep23
       fm2 = eps2/eps4
       fm3 = eps3/eps4

       w = twopi*T_bin_inv
       if(T_bin_inv.ne.0.d0) then
         T_bin = 1.d0/T_bin_inv
         sep23_xmx = sqrt(sep23_x0**2+(ds23xdt/w)**2)
         sep23_ymx = sqrt(sep23_y0**2+(ds23ydt/w)**2)
         tan_wtx = -ds23xdt/(w*sep23_x0)
         wtx = atan(tan_wtx)
         if(sep23_x0.lt.0.d0) wtx = wtx + pi
         dt23x = wtx/w
         tan_wty = -ds23ydt/(w*sep23_y0)
         wty = atan(tan_wty)
         if(sep23_y0.lt.0.d0) wty = wty + pi
         dt23y = wty/w
       else
         T_bin = 1.d20
         sep23_xmx = sep23_x0
         sep23_ymx = sep23_y0
       endif

cccXXX       sep23_mx = sqrt(sep23_xmx**2+sep23_ymx**2)

c      determine radius of 3d circular orbit
c      -------------------------------------
       delphi = w*(dt23y-dt23x)
       delphi2 = 2.d0*delphi
       denom = sep23_xmx**2+cos(delphi2)*sep23_ymx**2
       if(denom.eq.0.d0) then
         tan2alpha = 1.d30
       else
         tan2alpha = -sep23_ymx**2*sin(delphi2)/denom
       endif
       alpha = 0.5d0*atan(tan2alpha)
       alpha_pi2 = alpha + 0.5d0*pi
       sep23_x = sep23_xmx*cos(alpha)
       sep23_y = sep23_ymx*cos(alpha+delphi)
       sep23_mx = sqrt(sep23_x**2+sep23_y**2)
       sep23_x_pi2 = sep23_xmx*cos(alpha_pi2)
       sep23_y_pi2 = sep23_ymx*cos(alpha_pi2+delphi)
       sep23_pi2 = sqrt(sep23_x_pi2**2+sep23_y_pi2**2)
       sep23_mx=max(sep23_mx,sep23_pi2)

       sep14_mx = sep_rat*sep23_mx
       s14 = sep14_mx

       s = sep23_mx
       s3 = s**3

       al = thetaEmas*abs(pier)
       beta = 4.462202d-15*epsorb*T_bin**2/(s3*thetaE)
       c = beta/Dso**2
       write(6,*) 'sep23_mx =',sep23_mx,'   c =',c
       if(c.gt.4.d0/27.d0) then
         write(6,*) 'c > 4/27 : no solutions!'
         fmass(1)=-1.d0
         fmass(2)=c
         Ds_calc = 0.d0
         return
       endif

c      solve the cubic for x following sec 5.6 of Num Recipes
c      ------------------------------------------------------
       rnr = 0.5d0*c - 1.d0/27.d0
       qnr = 1.d0/9.d0
       sqrtq = sqrt(qnr)
       thetaNR = acos(rnr/sqrt(qnr**3))
       xrt(1) = -2.d0*sqrtq*cos(thetaNR/3.d0) + 1.d0/3.d0
       xrt(2) = -2.d0*sqrtq*cos((thetaNR+twopi)/3.d0) + 1.d0/3.d0
       xrt(3) = -2.d0*sqrtq*cos((thetaNR-twopi)/3.d0) + 1.d0/3.d0

c      solve the cubic equation for Ds_calc following sec 5.6 of Num Recipes
c      ---------------------------------------------------------------------
       xnr = 1.d0/(al**2-1.d0/beta)
       anr = 3.d0*al*xnr
       bnr = 3.d0*xnr
       cnr = xnr/al
       qnrd = (anr**2-3.d0*bnr)/9.d0
       rnrd = (2.d0*anr**3-9.d0*anr*bnr+27.d0*cnr)/54.d0
       q3 = qnrd**3
       r2 = rnrd**2
       if(q3.gt.r2) then
         sqrtqd = sqrt(qnrd)
         thetaNRd = acos(rnrd/sqrt(q3))
         Dsrt(1) = -2.d0*sqrtqd*cos(thetaNRd/3.d0)+1.d0/3.d0
         Dsrt(2) = -2.d0*sqrtqd*cos((thetaNRd+twopi)/3.d0)+1.d0/3.d0
         Dsrt(3) = -2.d0*sqrtqd*cos((thetaNRd-twopi)/3.d0)+1.d0/3.d0
         nDrt = 3
       else
         absr = abs(rnrd)
         signR = rnrd/absr
         AA = -signR*(absr+sqrt(r2-q3))**(1.d0/3.d0)
         if(AA.eq.0.d0) then
           BB = 0.d0
         else
           BB = qnrd/AA
         endif
         Dsrt(1) = (AA+BB)-anr/3.d0
         nDrt = 1
       endif

       write(6,98) theta_suas,Dso
       im=0
       do irt=1,3
         if(xrt(irt).ge.0.d0.and.xrt(irt).le.1.d0) then
           x=xrt(irt)
           im=im+1
           fmass(im)=thetaE**2*Dso*(3.085677581d16/5.90650016d0)
     &              *(x/(1.d0-x))
           ReAU = thetaE*x*Dso*(3.085677581d16/1.4959787966d8)
           vp = ReAU*1.4959787966d8/(t_E*24.d0*3600.d0)
           write(6,99) x,x*Dso,fmass(im),vp,ReAU,s*ReAU,sep1cm*ReAU
         endif
       enddo

       chi2_Ds = 1.d30
       do jrt = 1,nDrt
         chi2 = ((Dsrt(jrt)-Dso)/Ds_sig)**2
         x = 1.d0/(1.d0+al*Dsrt(jrt))
         fm = 
     &   thetaE**2*Dsrt(jrt)*(3.085677581d16/5.90650016d0)*(x/(1.d0-x))
         if(chi2.lt.chi2_Ds) then
           chi2_Ds = chi2
           Ds_calc = Dsrt(jrt)
         endif
       enddo
       x = 1.d0/(1.d0+al*Ds_calc)
       fmassD = 
     &   thetaE**2*Ds_calc*(3.085677581d16/5.90650016d0)*(x/(1.d0-x))
       ReAUd = thetaE*x*Ds_calc*(3.085677581d16/1.4959787966d8)
       vpd = ReAUd*1.4959787966d8/(t_E*24.d0*3600.d0)
       Dl_calc = x*Ds_calc
       write(6,97) theta_suas,Ds_calc,x,Dl_calc,fmassD,vpd,ReAUd,
     &             s*ReAUd,sep1cm*ReAUd

 97    format('Parallax Solution: theta_src =',f8.4,' uas',
     &        ' D_s =',f8.4,' kpc',/,
     &  '     x       D_l      M_tot    v_per    R_e(AU)   s_orb',
     &  '  s_plan',/,f9.6,2f10.5,f9.2,3f9.5)
 98    format('Solution assuming theta_src =',f8.4,' uas',
     &        ' D_s =',f8.4,' kpc',/,
     &  '     x       D_l      M_tot    v_per    R_e(AU)   s_orb',
     &  '  s_plan')
 99    format(f9.6,2f10.5,f9.2,3f9.5)

       return
       end

