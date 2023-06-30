c==============================================================================

       subroutine fcn(nfit,grad,chi2,a,iflag)

c==============================================================================
c  This program fits light curves for microlensing by binary lenses.
c  It is called as a subroutine by Minuit, and it is optimised for
c  high magnification events.
c
c   Author: David P. Bennett
c
c------------------------------------------------------------------------------
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       parameter (mmax=30,maxdata=50000,ngmax=12000000)
       parameter (nmcmc=100)
       double precision x(maxdata),y(maxdata),sig(maxdata)
       integer iclr(maxdata),nimage(maxdata),iwork(maxdata),
     &         jwork(maxdata)
       double precision a(mmax),grad(mmax),dtcaustic(maxdata),
     &        work(maxdata)
       parameter (mxclr=59,ndark=5000)
       double precision yfits(0:mxclr),yyfits(0:mxclr),
     &         chi2clr(0:mxclr)
       double precision A2(0:mxclr),A0(0:mxclr),A2sig(0:mxclr),
     &                  A0sig(0:mxclr),A0fix(0:mxclr)
       double precision yfit(maxdata,0:mxclr),xclr(maxdata,0:mxclr),
     &                  yclr(maxdata,0:mxclr),eclr(maxdata,0:mxclr)
       double precision fmin(0:mxclr),fmax(0:mxclr),scaus(20),tcaus(20)
       double precision fudge(0:mxclr),errmin(0:mxclr),dayoff(0:mxclr)
       double precision lon_obs(0:mxclr),lat_obs(0:mxclr)
       character*6 sfx(0:mxclr)
       integer nclr(0:mxclr),iclrind(maxdata),lsfx(0:mxclr),
     &         jclrinc(0:mxclr)
       integer idcaustic(12)
       double complex zalph,zgamm,box,z(5)
       double precision Date(maxdata),rMag(maxdata),rErr(maxdata),
     &                    bMag(maxdata),bErr(maxdata),
     &                    rFlux(maxdata),rFerr(maxdata),
     &                    bFlux(maxdata),bFerr(maxdata)
       double precision sxg(ngmax),syg(ngmax)
       character*1 c(0:9)
       character*19 char_mcmc(nmcmc)
       character*120 line,outline
       character*80 infile,ctinfile,newinffile,lcoutfile,tmpfile,parfile
       character*80 oldinffile,residfile,tmpfile2,mcmc_file
       character*80 satfile(0:mxclr),limbfile
       character*30 starname,fitname
       character*20 blank20
       character*10 parname(30)
       common/mcmc_output/nchar_mcmc,nrej_mcmc,char_mcmc
       common/seek_opt2/mcmc_file
       common/integrate/gridUstar
       common/limbdark/dark(ndark,0:mxclr),ald(0:mxclr),bld(0:mxclr),
     &                 darklimb(0:mxclr),hcut,bcut,iend,ildtab(0:mxclr)
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
       data iend /2/

       save

       if(iflag.eq.1) then
         pi=acos(-1.)
         twopi=2.d0*pi
         twobypi=2.d0/pi
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
     &    'Enter RA & DEC precessed to time of maximum magnification',
     &    'format: hh mm ss.s  deg mm ss.ss'
         read(5,*) rah,ram,ras,decd,decm,decs
         alpha = (rah +  ram/60. + ras/3600.)*15
         if(decd.lt.0.d0.or.decm.lt.0.d0.or.decs.lt.0.d0) then
           delta = -(abs(decd) + abs(decm)/60. + abs(decs)/3600.)
         else
           delta = (abs(decd) + abs(decm)/60. + abs(decs)/3600.)
         endif

         write(6,*) 'enter 1 for integration grid, 0 to skip;',
     &      ' and name of optional MCMC output file'
         read(5,'(a)') line
         read(line(1:2),*) make_grid_in
         make_grid = make_grid_in
         len_filep2 = lenc(line)
         if(len_filep2.gt.2) then
           read(line(3:len_filep2),'(a)') mcmc_file
           open(unit=10,file=mcmc_file,position='append')
           nchar_mcmc=13
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

         read(3,'(a)') line
         lline=lenc(line)
         do while (lline.gt.3.and.line(1:3).eq.'fix')
           read(line(4:lline),*) jclrfix,A0fix_in
           A0fix(jclrfix) = A0fix_in
           read(3,'(a)') line
           lline=lenc(line)
         enddo
         read(3,'(a)') line
         read(line,*,err=72,end=72) daycausmin,daycausmax,deltcaus,
     &                delfine,gridUstar,hcut,iend,grid_rat
         go to 77
 72      read(line,*,err=74,end=74) daycausmin,daycausmax,deltcaus,
     &                delfine,gridUstar,hcut,iend
         go to 77
 74      read(line,*,err=76,end=76) daycausmin,daycausmax,deltcaus,
     &                delfine,gridUstar,hcut
         go to 77
 76      read(line,*) daycausmin,daycausmax,deltcaus,delfine,gridUstar
 77      continue
         bcut = (2.d0/3.d0)*sqrt((hcut+0.5d0)/hcut)
         write(6,*) 'daycausmin =',daycausmin,' daycausmax =',daycausmax
         write(6,*) 'deltcaus =',deltcaus,' delfine =',delfine,
     &              ' gridUstar =',gridUstar
         write(6,*) 'hcut =',hcut,' grid_rat =',grid_rat,' iend =',iend
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
               sig(ndata)=sqrt((rErr(j)/1.0857362)**2+errmin(0)**2)
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
               sig(ndata)=sqrt((bErr(j)/1.0857362)**2+errmin(1)**2)
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
          jnext=max(j,1)
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
                 sig(ndata)=sqrt((rErr(j)/1.0857362)**2
     &               +(dErr/1.0857362)**2+errmin(jclr)**2)
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
           write(6,*) 'read in GMAN data from ctinfile:',ctinfile(1:50)
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
               sig(ndata)=sqrt((rErr(j)/1.0857362)**2+errmin(jclr)**2)
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
           write(6,*) 'read in GMAN data from ctinfile:',ctinfile(1:50)
           write(6,*) 'read in ',ndata,' total measurements'

         do 240 jclr=9,14
           ctinfile='lc'//starname(1:lsn)//'.'//sfx(jclr)(1:lsfx(jclr))
           open(unit=1,file=ctinfile,status='old',err=231)
           do 230 j=jnext,99999
             read(1,*,end=231,err=231) Date(j),ia,ib,itype,xxx,yyy,
     &           rMag(j),dm,rErr(j),chi2,dsky,c1,crd,fmiss,cr
               if(cr.lt.1.e-4.and.fmiss.lt.0.1.and.crd.lt.2..and.
     &                     rErr(j).ge.0..and.
     &                     (chi2.lt.500..or.rMag(j).le.-7.0)) then
                 ndata=ndata+1
                 Date(j)=Date(j)+dayoff(jclr)
                 if(Date(j).lt.daymin) daymin=Date(j)
                 if(Date(j).gt.daymax) daymax=Date(j)
                 x(ndata)=Date(j)
                 y(ndata)=10**(-0.4*(rMag(j)))
                 sig(ndata)=sqrt((rErr(j)/1.0857362)**2+errmin(jclr)**2)
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
           write(6,*) 'read in GMAN data from ctinfile:',ctinfile(1:50)
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
                 sig(ndata)=sqrt((rErr(j)/1.0857362)**2+errmin(jclr)**2)
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
           write(6,*) 'read in GMAN data from ctinfile:',ctinfile(1:50)
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
           write(6,*) 'read in GMAN data from ctinfile:',ctinfile(1:50)
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
           write(6,*) 'read in GMAN data from ctinfile:',ctinfile(1:50)
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
           write(6,*) 'read in GMAN data from ctinfile:',ctinfile(1:50)
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
         call sort_light_curve_data_by_time(ndata,x,y,sig,iclr,iclrind)

c        make the limb darkening table
c        -----------------------------
         call limb_dark(limb_tab,limbfile)
         call hexlimbfac()

       endif

       write(6,*) 'FCN call with a =',(a(k), k=1,12)
       call flush(6)

       chi2=0.
       r0=0.
       maximages=0
       do 290 i=0,nnclr
         chi2clr(i)=0.
 290   continue

c  initialize source grid storage
c  ------------------------------
       Ustar=a(9)*a(1)
       eps1 = a(6)
       eps2 = 1.d0 - eps1
       if(a(9).gt.0.) then
         bgriddel0=gridUstar*Ustar
         ngfac=nint(2.d0*log(1.d5*bgriddel0)/log(2.d0))
         bgriddel=1.d-5*sqrt(2.d0**ngfac)
c        bphigrid ~ grid_rat*bgriddel
         nphimax=int(pi/(grid_rat*bgriddel))
         bphigrid=twopi/(2*nphimax+1)
         ngrtot=ngmax/(nphimax+1)
         ngr=(ngrtot-1)/2
         brgrid=bgriddel

c        different Einstein Radius and grid origin for sep > 1
c        -----------------------------------------------------
         if(sep.gt.1.d0) then
           Ein_R = sqrt(eps2)
           brgrid = Ein_R*brgrid
           xcc = eps1*(sep-1.d0/sep)
         else
           Ein_R = 1.d0
           xcc = 0.d0
         endif
         if(make_grid.ne.0) then
           call microcurve_init(a,Ein_R,xcc,brgrid,bphigrid,
     &                      ngr,nphimax,sxg,syg)
         else
           ngr = 0
         endif
       endif

c  sum over data points
c  --------------------
       do 300 i=1,ndata
          t=x(i)
          ii=iclrind(i)
          jc = iclr(i)
          call microcurve(t,a,yfit(ii,jc),jc,nimage(i),
     &                 alpha,delta,Ein_R,xcc,lon_obs(jc),lat_obs(jc),0,
     &                 brgrid,bphigrid,ngr,nphimax,sxg,syg,ssx,ssy)
          maximages=max(maximages,nimage(i))
 300   continue
c      nchar_mcmc = # of fit parameters + 1 for chi2
       if(nchar_mcmc.gt.0) nchar_mcmc = 13
       do 350 j=0,nnclr
          if(nclr(j).gt.0.and.A0fix(j).le.0.d0) then
             if(iflag.eq.3) write(6,*) 'chi2 values for color:',j
             call normfit(yfit(1,j),yclr(1,j),nclr(j),eclr(1,j),wfm,
     &          A2(j),A0(j),A2sig(j),A0sig(j),chi,fmin(j),fmax(j),iflag)
             chi2clr(j)=chi
             chi2=chi2+chi
             if(nchar_mcmc.gt.0) then
               nchar_mcmc = nchar_mcmc + 1
               write(char_mcmc(nchar_mcmc),998) A0(j)
               nchar_mcmc = nchar_mcmc + 1
               write(char_mcmc(nchar_mcmc),998) A2(j)
             endif
          elseif(nclr(j).gt.0) then
             A0(j) = A0fix(j)
             if(iflag.eq.3) write(6,*) 'chi2 values for color:',j
             call normfit_fix1(yfit(1,j),yclr(1,j),nclr(j),eclr(1,j),
     &          A2(j),A0(j),A2sig(j),A0sig(j),chi,iflag)
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

       write(6,*) ' chi2 =',chi2
       if(nchar_mcmc.gt.0) then
         write(char_mcmc(1),998) chi2
         do ipar=1,12
           write(char_mcmc(ipar+1),998) a(ipar)
         enddo
       endif
       call flush(6)

       if(iflag.eq.3) then
         open(unit=9,file=residfile)
         if(abs(a(7)).lt.1.e-20) then
          write(9,557) chi2,A0(0),A0(1),1./a(1),a(2),a(3),a(4),a(5),
     &                a(6),1.-a(6),a(9),a(7),a(8),a(11),a(12),
     &                A0(2),A2(0),A2(1),A2(2),a(10)
         else
          write(9,555) chi2,A0(0),A0(1),1./a(1),a(2),a(3),a(4),a(5),
     &                a(6),1.-a(6),a(9),1./a(7),a(7),a(8),a(11),a(12),
     &         A0(2),A2(0),A2(1),A2(2),a(10)
         endif
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
         close(9)
 565     format('       t          Fit-Ampl     Residual    Resid-Sig',
     &          '      chi^2    Iclr')
 566     format(f13.6,g15.6,2f13.7,g13.5,i4,1x,a6)

         do 410 i=0,nnclr
           if(nclr(i).gt.0) write(6,553) i,chi2clr(i),nclr(i)
 410     continue
 553     format('color #',i3,' chi2 =',g12.6,' for',i6,' points')

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
 560     format('A0rm     =',f12.3,' +/-',f11.3,
     &         '  A2rm     =',f12.3,' +/-',f11.3)
 561     format('A0bm     =',f12.3,' +/-',f11.3,
     &         '  A2bm     =',f12.3,' +/-',f11.3)
 562     format('A0',a6,' =',f12.3,' +/-',f11.3,
     &         '  A2',a6,' =',f12.3,' +/-',f11.3)

c        determine caustic crossing times
c        --------------------------------
         umin=a(3)
         theta = a(5)
         ctheta=cos(theta)
         stheta=sin(theta)
         nim=3
         nim00=3
         ncaus=0
         ds=0.001
         nstep=4000
         dst2=ds*2.
         dst22=dst2**2
         dsx=ds*ctheta
         dsy=ds*stheta
         cosfac=cos(thetaproj)
         sinfac=sin(thetaproj)
         do 400 i=-nstep,nstep,1
           vfac=ds*i
           sx_norot=ctheta*vfac-stheta*a(3)
           sy_norot=stheta*vfac+ctheta*a(3)
           sx=cosfac*sx_norot-sinfac*sy_norot
           sy=cosfac*sy_norot+sinfac*sx_norot

 393       nim0=nim
           call bilens_im(sx,sy,nim,iimage,z,ampim)
          if(nim.ne.nim0) then
             if(nim0.eq.3.or.nim0.eq.5) s0=ds*(i-1)
             if(nim.eq.3.or.nim.eq.5) then
               s1=vfac
               do 395 j=1,20
                 s=0.5*(s0+s1)
                 sx_norot=ctheta*s-stheta*a(3)
                 sy_norot=stheta*s+ctheta*a(3)
                 sx=cosfac*sx_norot-sinfac*sy_norot
                 sy=cosfac*sy_norot+sinfac*sx_norot
                 call bilens_im(sx,sy,nim_mid,iimage,z,ampim)
                 if(nim_mid.eq.nim) then
                   s1=s
                 elseif(nim_mid.eq.nim00) then
                   s0=s
                 else
                   go to 396
                 endif
 395           continue
 396           continue
               ncaus=ncaus+1
               scaus(ncaus)=s
               tcaus(ncaus)=a(2)+scaus(ncaus)/a(1)
               nim00=nim
             endif
           endif
 400     continue

         if(ncaus.gt.1) then
           write(6,*) 'Caustic crossings found at',ncaus,' times:'
           write(6,552) (tcaus(i),i=1,ncaus)
         endif
 552     format(6f12.5)

         if(abs(a(7)).lt.1.e-20) then
          write(6,557) chi2,A0(0),A0(1),1./a(1),a(2),a(3),a(4),a(5),
     &                a(6),1.-a(6),a(9),a(7),a(8),a(11),a(12),
     &                A0(2),A2(0),A2(1),A2(2),a(10)
         else
          write(6,555) chi2,A0(0),A0(1),1./a(1),a(2),a(3),a(4),a(5),
     &                a(6),1.-a(6),a(9),1./a(7),a(7),a(8),a(11),a(12),
     &         A0(2),A2(0),A2(1),A2(2),a(10)
         endif
         if(nnclr.lt.40) then
           write(6,558) (sfx(i),sfx(i),sfx(i+1),sfx(i+1),
     &                 sfx(i+2),sfx(i+2),sfx(i+3),sfx(i+3),
     &                 A0(i),A2(i),A0(i+1),A2(i+1),
     &                 A0(i+2),A2(i+2),A0(i+3),A2(i+3), i=3,35,4)
         elseif(nnclr.lt.50) then
           write(6,558) (sfx(i),sfx(i),sfx(i+1),sfx(i+1),
     &                 sfx(i+2),sfx(i+2),sfx(i+3),sfx(i+3),
     &                 A0(i),A2(i),A0(i+1),A2(i+1),
     &                 A0(i+2),A2(i+2),A0(i+3),A2(i+3), i=3,45,4)
         else
           write(6,558) (sfx(i),sfx(i),sfx(i+1),sfx(i+1),
     &                 sfx(i+2),sfx(i+2),sfx(i+3),sfx(i+3),
     &                 A0(i),A2(i),A0(i+1),A2(i+1),
     &                 A0(i+2),A2(i+2),A0(i+3),A2(i+3), i=3,55,4)
         endif

 555     format('    chisq       A0ogle    A0saao   t_E         t0',
     &        '       umin    sep      theta',/,
     &        f12.2,2f9.2,f9.3,f10.4,f11.7,2f9.5,/,
     &        '    eps1      eps2    Tstar     Tbin    1/Tbin  ',
     &        '  v_sep    piEr   piEtheta',/,
     &        2(1pe11.4),0pf8.5,f9.3,4f9.5,/,
     &        '    A0ms74     A2ms74     A2ogle    A2saao  t_fix',/,
     &        4f11.3,f11.3)

 557     format('    chisq       A0ogle    A0saao   t_E         t0',
     &        '       umin    sep      theta',/,
     &        f12.2,2f9.2,f9.3,f10.4,f11.7,2f9.5,/,
     &        '    eps1      eps2    Tstar     Tbin    1/Tbin  ',
     &        '  v_sep    piEr   piEtheta',/,
     &        2(1pe11.4),0pf8.5,'   1.e20 ',4f9.5,/,
     &        '    A0ms74     A2ms74     A2ogle    A2saao  t_fix',/,
     &        4f11.3,f11.3)

 558     format(4('    A0',a6,'    A2',a6),/,8f12.3)
ccc 554     format('   A0',a6,'   A2',a6,/,2f12.3)

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
         if(abs(a(7)).lt.1.e-20) then
          write(2,557) chi2,A0(0),A0(1),1./a(1),a(2),a(3),a(4),a(5),
     &                a(6),1.-a(6),a(9),a(7),a(8),a(11),a(12),
     &                A0(2),A2(0),A2(1),A2(2),a(10)
         else
          write(2,555) chi2,A0(0),A0(1),1./a(1),a(2),a(3),a(4),a(5),
     &                a(6),1.-a(6),a(9),1./a(7),a(7),a(8),a(11),a(12),
     &         A0(2),A2(0),A2(1),A2(2),a(10)
         endif
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
           yfit0=0.
           do 450 j=0,nnclr
             if(A0(j).ne.0.) then
               call microcurve(t,a,yfits(j),j,nn,
     &                 alpha,delta,Ein_R,xcc,lon_obs(j),lat_obs(j),0,
     &                 brgrid,bphigrid,ngr,nphimax,sxg,syg,ssx,ssy)
ccc               call microcurve(t,a,yfits(j),cfits(j),j,nn,alpha,delta,-9.,
ccc     &                       zalph,zgamm)
               yyfits(j)=A0(j)*yfits(j)+A2(j)
               if(yfit0.eq.0.) yfit0=yfits(j)
             else
               yyfits(j)=A2(j)
             endif
 450       continue
           write(2,556) t,yfit0,(yyfits(jclrinc(j)),j=1,nclrinc),ssx,ssy
 460     continue
 461     continue
 556     format(f9.4,f10.4,41f15.5)
         close(2)

         open(unit=11,file=oldinffile,status='old')
         open(unit=4,file=newinffile,status='new')
         read(11,'(a)') line
         len=lenc(line)
         write(4,'(a)') line(1:len)
         do i=1,12
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
     &           iflag)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision X(NDATA),Y(NDATA),SIG(NDATA)

      dmn=(1./fmax-1.)
      dmx=(1./fmin-1.)
      wdmn=w*dmn
      wdmx=w*dmx
      SX=0.d0
      SY=0.d0
      ST2=0.d0
      B=0.d0
      SS=0.d0
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
        b=0.d0
        a=sy/ss
        st2=1.e-20
      else
        B=B/ST2
        A=(SY-SX*B)/SS
        f=b/(b+a)
        if((f.lt.fmin.or.f.gt.fmax).and.abs(fmin).lt.1.d8) then
          ssp=(1+w)*ss
          if(f.lt.fmin) then
            dm=dmx
            sxp=sxpdmx
          else
            dm=dmn
            sxp=sxpdmn
          endif
          sxpossp=sxp/ssp
          b=0.d0
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
      DO 16 I=1,NDATA
        chi21=((Y(I)-A-B*X(I))/SIG(I))**2
        CHI2=CHI2+chi21
        if(iflag.eq.3) write(6,*) X(I),y(i),chi21,CHI2
16    CONTINUE
      RETURN
      END

      SUBROUTINE normFIT_fix1(X,Y,NDATA,SIG,A,B,SIGA,SIGB,CHI2,iflag)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision X(NDATA),Y(NDATA),SIG(NDATA)

      SX=0.d0
      SY=0.d0
      ST2=0.d0
      SS=0.d0
      DO 11 I=1,NDATA
        WT=1./(SIG(I)**2)
        SS=SS+WT
        xwt=X(I)*WT
        SX=SX+xwt
        SY=SY+Y(I)*WT
11    CONTINUE
      A=(SY-SX*B)/SS
      CHI2=0.
      DO 16 I=1,NDATA
        chi21=((Y(I)-A-B*X(I))/SIG(I))**2
        CHI2=CHI2+chi21
        if(iflag.eq.3) write(6,*) X(I),chi21,CHI2
16    CONTINUE
      RETURN
      END
