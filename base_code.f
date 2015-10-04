c^^^^^^^^^^^^^^^^^^^^^^^^^^^
      Program Channel                                                     
c^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
c  A Three-Dimensional Unsteady Navier-Stokes Solver
c  developed by Daniel Chiu-Leung Chan, 1993
c
c  homogeneous in x and y direction
c  inhomogeneous in z direction
c
c  baseline code: r99.f
c  modifications: (1) subroutine update, making sure w=0 at walls
c                 (2) full skew symmetric form
c                 (3) subroutine convec, 4-stage R-K scheme
c                 (4) modify restart file format
c                 (5) source term in main program
c
c  uses CGS solver instead of direct solvers, save in memory
c
      include 'dim_2d.h'
c
      common  /wave/  xw(nxhp),xsq(nxhp),alpha,u00
      common  /leng/  zpts(0:nzm),wg(0:nzm),d(0:nzm,0:nzm)                                 
      common  /time5/  t,dt,nsteps                                                 
      common  /flow/  u(ntot),w(ntot),temp(ntot),                 
     1                un(ntot),wn(ntot),tn(ntot)
      common  /fft/   trigsx(nxf),work(ntot),ifaxx(13)
      common /matrix/ su(ntot),sw(ntot),st(ntot)
      common /coeff/  diff(nz,nz), amass(nz)
      common /param/  re
      common /bc  /   ubc(nxpp),wbc(nxpp)
      common /rotat/  omz
      common /scratch/ ox(ntot),oz(ntot)
      common /width/   ybar  
      common /source/ body(ntot)
      common /solver/ cgstol
      common /exforce/ wavlen
c
c
      dimension uw1(nxpp),ww1(nxpp)
      dimension uw2(nxpp),ww2(nxpp)
      character*2 f11
      character*12 fil
c
c-------------------------------------------------------------------
c
c  irstart: set to 1 to restart from a previous run
c  nsteps:  number of time steps to run
c  dt:      size of time step
c  alpha:   wavelength in the x-direction
c  beta:    wavelength in the y-direction
c  re:      Reynolds number (based on channel half width and utau
c  cgstol was 1e-4 
c-------------------------------------------------------------------
      namelist/inputs/ istart,dt,nsteps,nwrt,iform,iles
      namelist/params/ xlen,ylen,re,ta,ybar,cgstol,cs,u00,wavlen
      data istart,dt,nsteps,alpha,beta,re,ta
     &    /0,0.01,100,.0,1.0,180.,0./  
      data nwrt,iform/10,0/
      data ybar/1.0/
      data cgstol/1.0e-4/
      data iles/0/
      data cs/0.1/
      data wavlen/1.0/
      data ylen/1.0/
c
c  set up the working arrays for the FFT calls
c
      call fftfax(nx,ifaxx,trigsx)
c
      pi2 = 8.*atan(1.0)
c
      do i=1,ntot
       su(i) = 0.
       sw(i) = 0.
      enddo
c
      do i=1,nxpp
       uw1(i) = 0.
       ww1(i) = 0.
       uw2(i) = 0.
       ww2(i) = 0.
      enddo
c
      open(7,file='input.dat',status='old')
      open(10,file='time.dat',status='unknown')
c
      read(7,inputs)                                                  
      read(7,params)                                                  
c                                                                       
      alpha = pi2/xlen
      retau = re
c
      print *,'alpha',alpha
      print *,'u00/vel ratio',u00
      print *,'Wavelength :',wavlen
c
      if(istart.eq.0) then 
        call setup                                                      
        call init(u,w,temp)                                             
        t = 0.                                                             
        fac1 = 1.0                                                          
        fac2 = 0.                                                          
      else                                         
       open(8,file='start.dat',form='unformatted',status='unknown')
       call restart(8,t,uw1,ww1,uw2,ww2)
       call setup                                                      
       fac1 = 1.5                                                      
       fac2 = -.5                                                      
      end if                                                            
c                                                                       
      if(istart.eq.1) then
      rewind(7)
      read(7,inputs)                                                  
      read(7,params)                                                  
      endif
c
      write(6,inputs)
      write(6,params)
c
      fact = 2./dt
      facvel = fact*re
c
      print *,(xw(i),i=1,nxhp)
c
c
      ntavg = 0
c
      do 1000 istep = 1,nsteps
c
c                                                                       
       call convec(fac1,fac2)
c

       fac1 = 1.5                                                      
       fac2 = -.5                                                      
c
       do ijk=1,ntot
        su(ijk) = (u(ijk)+su(ijk))*re/dt    
        sw(ijk) = (w(ijk)+sw(ijk))*re/dt 
       enddo 
c

      call helmv(u,su,uw1,uw2,
     &           w,sw,ww1,ww2,
     &           temp,st,xw,
     &           xsq,facvel,fact,diff,amass,
     &           dt,nz)
c
      call edge(u,w)
c
      t = t + dt
      write(6,*) 'at time ',t
c
      do m=1,ntot
       su(m) = 0.
      enddo
c
      call div(u,w,xw,su)
      call fbig(su)
c
      if(mod(istep,nwrt).eq.0) then
      open(2,file='run.dat',form='unformatted',status='unknown')
      call output(2,t,uw1,ww1,uw2,ww2)
      close(2)
c
      f11 = 'p_'
      write(fil,15)f11,t
15    format(a2,e10.4)
      if (iform.eq.0) then
      open(3,file=fil,form='unformatted',status='unknown')
      else
      open(3,file=fil,status='unknown')
      endif
      call plotf(3,iform,t,u,w,temp)
      close(3)
c
      endif
c
1000  continue
c
      open(2,file='run.dat',form='unformatted',status='unknown')
      call output(2,t,uw1,ww1,uw2,ww2)
      if (iform.eq.0) then
      open(3,file='plot.dat',form='unformatted',status='unknown')
      else
      open(3,file='plot.dat',status='unknown')
      endif
      call plotf(3,iform,t,u,w,temp)
c
62     format(5e14.6)
       close(1)
c
      stop
      end               
c
c*********************************************************************  
      subroutine fbig(a)
c*********************************************************************  
      include 'dim_2d.h'
      dimension a(*)
c
c  transform to physical space
c
       call horfft(a,1,nz)
c           
      do i=1,ntot
       a(i) = abs(a(i))
      enddo
c
      big = 0
      imax = 0
      jmax = 0
      kmax = 0
c
      do k=2,nzm
        kk = (k-1)*nxpp
          do i=1,nx
           ik = i + kk
           if(a(ik).gt.big) then 
             big = a(ik)
             imax = i
             kmax = k
           endif
          enddo
      enddo 
c
      write(6,*) 'largest divergence defect ',big,' at ',imax,jmax,kmax
c
      return
      end
c
c*********************************************************************  
      subroutine exbc(su,uw1,ww1,uw2,ww2)
c*********************************************************************  
c
c   extrapolate the velocity boundary condition
c   assuming the wall velocity is zero
c
      include 'dim_2d.h'
      common  /leng/  zpts(0:nzm),wg(0:nzm),d(0:nzm,0:nzm)                                 
      common  /time5/  t,dt,nsteps                                                 
      common  /wave/  xw(nxhp),xsq(nxhp),alpha,u00
      dimension  ox(ntot),oz(ntot)
      complex uw1(*),ww1(*),
     1        uw2(*),ww2(*),
     1        o1(nxhp,nz),o2(nxhp,nz),unit
      common /width/ ybar
      dimension su(*)
      equivalence (ox,o1),(oz,o2)
c
      call copy(su,oy,ntot)
      call dpdz(su,ox,d,nzm,ybar)
c
      unit = (0.,1.)
       do i=1,nxhp
        uw1(i)  = 0.5*unit*xw(i)*o2(i,1)*dt
        uw2(i)  = 0.5*unit*xw(i)*o2(i,nz)*dt
        ww1(i)  = 0.5*dt*o1(i,1)
        ww2(i)  = 0.5*dt*o1(i,nz)
       enddo
c
       return
       end
c
c*********************************************************************  
      subroutine update(u,w,ub,wb,xw,su,dt)
c*********************************************************************  
c
c   update the velocity field
c
      include 'dim_2d.h'
c
      dimension su(*),xw(*)     
      common  /leng/  zpts(0:nzm),wg(0:nzm),d(0:nzm,0:nzm)                                 
      dimension ox(ntot),oy(ntot)
      common /width/   ybar
      complex o1(nxhp,nz),o2(nxhp,nz)
      complex u(nxhp,nz),w(nxhp,nz),
     1        ub(nxhp,nz),wb(nxhp,nz)
      equivalence (ox,o1), (oy,o2)
c
      call copy (su,oy,ntot)
      call dpdz (su,ox,d,nzm,ybar)
c
       do k=2,nz-1
         do i=1,nxhp
           u(i,k) = ub(i,k) - dt*(0.,1.)*xw(i)*o2(i,k) 
           w(i,k) = wb(i,k) - dt*o1(i,k) 
         enddo
        enddo
c
        do i=1,nxhp
         w(i,1)  = cmplx(0.0,0.0)
         w(i,nz) = cmplx(0.0,0.0)
        enddo
c
       return
       end
c
c*********************************************************************  
      subroutine div(u,w,xw,su)
c*********************************************************************  
c
c   calculate the divergence of a velocity field
c
      include 'dim_2d.h'
c
      dimension w(ntot),xw(*)      
      common  /leng/  zpts(0:nzm),wg(0:nzm),d(0:nzm,0:nzm)                                 
      common /scratch/ ox(ntot),oz(ntot)
      common /width/   ybar
      complex o3(nxhp,nz)
      complex su(nxhp,nz),u(nxhp,nz)
      equivalence (oz,o3)
c
      call dpdz (w,oz,d,nzm,ybar)
c
       do k=1,nz
         do i=1,nxhp
           su(i,k) = (0.,1.)*xw(i)*u(i,k) +
     1                  o3(i,k) 
         enddo
        enddo
       return
       end 
c
c*********************************************************************  
       subroutine zero(a,n)
c*********************************************************************  
       dimension a(*)
       do i=1,n
        a(i) = 0.
       enddo
       return
       end
c
c*********************************************************************  
      subroutine convec(fac1,fac2)
c*********************************************************************  
c
      include 'dim_2d.h'
c
c
      common  /profile/ ubar(ntot)
      common  /flow/  u(ntot),w(ntot),temp(ntot),                 
     1                un(ntot),wn(ntot),tn(ntot)                                  
      common /matrix/ su(ntot),sw(ntot),st(ntot)
      common  /time5/  t,dt,nsteps                                                 
      common  /wave/  xw(nxhp),xsq(nxhp),alpha,u00
      common  /leng/  zpts(0:nzm),wg(0:nzm),d(0:nzm,0:nzm)                                 
      common  /param/ re
      common  /coeff/ diff(nz,nz), amass(nz)
      common  /fft/   trigsx(nxf),work(ntot),ifaxx(13)
      common /rotat/ omz
      common /scratch/ ox(ntot),oz(ntot)
      common /width/   ybar
      common /bc/      ubc(nxpp),wbc(nxpp)
      dimension dudx(ntot),dudz(ntot),
     1          dwdx(ntot),dwdz(ntot),
     1          uu(ntot),uw(ntot),ww(ntot)
c
      call copy(u,un,ntot)
      call copy(w,wn,ntot)
c
      do irk=4,1,-1
       alf = 1.0/float(irk)
       call dpdz(u,dudz,d,nzm,ybar)
       call dpdz(w,dwdz,d,nzm,ybar)
c
       call cdpdx(u,dudx,xw)
       call cdpdx(w,dwdx,xw)
c
c  spectral-to-physcial
c
      call  horfft (u,1,nz)
      call  horfft (w,1,nz)
      call  horfft (dudx,1,nz)
      call  horfft (dudz,1,nz)
      call  horfft (dwdx,1,nz)
      call  horfft (dwdz,1,nz)
c
      if (irk.eq.4) then
        call stavg(u,w,t,dudz)
      endif

c
c  udx + wdz in physical space                            
c
      call zero(su,ntot)
      call zero(sw,ntot)
c
      do i=1,ntot                                                  
       su(i) = -u(i)*dudx(i)  - w(i)*dudz(i) - 
c
c  it is assumed that dU/dz=1
c
     &         ubar(i)*dudx(i) - w(i)
c
       sw(i) = -u(i)*dwdx(i) - w(i)*dwdz(i)
     &         -ubar(i)*dwdx(i)

       uu(i) = u(i)*u(i)
       uw(i) = u(i)*w(i)
       ww(i) = w(i)*w(i)
      enddo 
c
c  physical-to-spectral
c
      call  horfft (su,-1,nz)
      call  horfft (sw,-1,nz)
      call  horfft (u,-1,nz)
      call  horfft (w,-1,nz)
      call  horfft (uu,-1,nz)
      call  horfft (uw,-1,nz)
      call  horfft (ww,-1,nz)
c
      call cdpdx(uu,dudx,xw)
      call dpdz(uw,dudz,d,nzm,ybar)
c
      do i=1,ntot
       su(i) = 0.5*( su(i) - dudx(i)  - dudz(i) ) 
      enddo
c
      call cdpdx(uw,dudx,xw)
      call dpdz(ww,dudz,d,nzm,ybar)
c
      do i=1,ntot
       sw(i) = 0.5*( sw(i) - dudx(i) - dudz(i) )
      enddo                                                
c
      do i=1,ntot
       u(i) = un(i) + alf*dt*su(i)
       w(i) = wn(i) + alf*dt*sw(i)
      enddo
c
c
c  fix boundary condition
c
      call wallbc(u,w,ubc,wbc)
c
      enddo           
c
      call copy(u,su,ntot)
      call copy(w,sw,ntot)
      call copy(un,u,ntot)
      call copy(wn,w,ntot)
c      
      return
      end             
c
c*********************************************************************  
      subroutine helmv(u,su,ubc1,ubc2,
     &                 w,sw,wbc1,wbc2,
     &                 temp,st,xw,
     &                 xsq,facv,fact,diff,amass,
     &                 dt,ndim)
c*********************************************************************  
      include 'dim_2d.h'

      common  /bc  /  ubc(nxpp),wbc(nxpp)
      common  /leng/  zpts(0:nzm),wg(0:nzm),d(0:nzm,0:nzm)                                 
      common /scratch/ ox(ntot),oz(ntot)
      common /width/   ybar
      common /solver/ cgstol

      dimension ubc1(*),wbc1(*)
      dimension ubc2(*),wbc2(*)
      dimension u(*),su(*),w(*),sw(*)
      dimension temp(*),st(*)
      dimension diff(ndim,*), amass(*)
      dimension xw(*),xsq(*)
      dimension apv(nz,nz),
     1          fur(nz),fui(nz),
     3          fwr(nz),fwi(nz),
     4          fr(nz),fi(nz),xr(nz),xi(nz)
c
      dimension furxy(nz),fuixy(nz),
     2          fwrxy(nz),fwixy(nz)  
c
c
      dimension ub(ntot),wb(ntot)
c
c
c sweep through the horizontal plane
c in spectral space
c
      ii = 0
      do 100 i=1,nxpp,2
      ii = ii + 1
c
       ak = xsq(ii) 
c
       do k2=1,nz
        do k1=1,nz
         apv(k1,k2) = diff(k1,k2) 
        enddo                 
       enddo

c
       do k1=1,nz
         apv(k1,k1) = apv(k1,k1) + (facv+ak)*amass(k1)
       enddo


c
c  impose velocity boundary condition
c
        do k5=1,nz
         furxy(k5) =  apv(k5,1) *ubc1(i)  + 
     &                apv(k5,nz)*ubc(i)
         fuixy(k5) =  apv(k5,1) *ubc1(i+1)  + 
     &                apv(k5,nz)*ubc(i+1)

        enddo
c

       do k=1,nz
         apv(1,k) = 0.0
         apv(k,1) = 0.0
         apv(nz,k) = 0.0
         apv(k,nz) = 0.0
       enddo
         apv(1,1) = 1.0
         apv(nz,nz) = 1.0
c
c
       do k=1,nz
         ik = i + (k-1)*nxpp
         ipk = ik + 1
         fur(k) = su(ik)*0.5*ybar*wg(k-1)   - furxy(k)
         fui(k) = su(ipk)*0.5*ybar*wg(k-1)  - fuixy(k)
      enddo
c
        fur(1)  = ubc1(i)
        fur(nz) = ubc(i)
        fui(1)  = ubc1(i+1)
        fui(nz) = ubc(i+1)

c
c  update the u array
c

        do k=1,nz-1
         xr(k) = 0.0
         xi(k) = 0.0 
        enddo
 
        xr(nz) = ubc(i)
        xi(nz) = ubc(i+1)


      call jacgs(apv,xr,fur,500,nz,nz,nz,0,cgstol)
      call jacgs(apv,xi,fui,500,nz,nz,nz,0,cgstol)

        do k=1,nz
         ik = i + (k-1)*nxpp
         ipk = ik + 1
         ub(ik)   = xr(k)
         ub(ipk)  = xi(k)
        enddo


c
c  update the w array
c  
c  no-slip boundary condition along the bottom wall
c  and specified u-velocity along the top wall
c  with the gradient dv/dy derived from continuity condition
c

      do k2=1,nz
        do k1=1,nz
         apv(k1,k2) = diff(k1,k2) 
        enddo                 
       enddo

c
       do k1=1,nz
         apv(k1,k1) = apv(k1,k1) + (facv+ak)*amass(k1)
       enddo


        do k5=1,nz
         fwrxy(k5) =  apv(k5,1) *wbc1(i) 
         fwixy(k5) =  apv(k5,1) *wbc1(i+1)  
        enddo



       do k=1,nz
         apv(1,k) = 0.0
         apv(k,1) = 0.0
       enddo
         apv(1,1) = 1.0
c
c
       do k=1,nz
         ik = i + (k-1)*nxpp
         ipk = ik + 1
         fwr(k) = sw(ik)*0.5*ybar*wg(k-1)   - fwrxy(k)  
         fwi(k) = sw(ipk)*0.5*ybar*wg(k-1)  - fwixy(k)
      enddo
c
        fwr(1)  = wbc1(i)
        fwi(1)  = wbc1(i+1)

        fwr(nz) = fwr(nz) + wbc(i)
        fwi(nz) = fwi(nz) + wbc(i+1)

        do k=1,nz
         xr(k) = 0.0 
         xi(k) = 0.0 
        enddo




      call jacgs(apv,xr,fwr,500,nz,nz,nz,0,cgstol)
      call jacgs(apv,xi,fwi,500,nz,nz,nz,0,cgstol)

        do k=1,nz
         ik = i + (k-1)*nxpp
         ipk = ik + 1
         wb(ik)   = xr(k)
         wb(ipk)  = xi(k)
        enddo


c
100     continue
c
      do ijk=1,ntot
       ub(ijk) = 2.*ub(ijk) - u(ijk)
       wb(ijk) = 2.*wb(ijk) - w(ijk)
      enddo
c        
c  pressure step
c
      do m=1,ntot
       su(m) = 0.
      enddo
c
      call div(ub,wb,xw,su)
c
      do m=1,ntot
       su(m) = -su(m)/dt
       sw(m) = 0.
      enddo
c
c sweep through the horizontal plane
c
      ii = 0
      do 110 i=1,nxpp,2
      ii = ii + 1
c
        ak = xsq(ii) 
c
       do k2=1,nz
        do k1=1,nz
         apv(k1,k2) = diff(k1,k2) 
        enddo                 
       enddo
c
       do k1=1,nz
         apv(k1,k1) = apv(k1,k1) + ak*amass(k1)
       enddo


        do k=1,nz
         ik = i + (k-1)*nxpp
         ipk = ik + 1
         fr(k) = su(ik)*0.5*ybar*wg(k-1)
         fi(k) = su(ipk)*0.5*ybar*wg(k-1)
        enddo
c

      if( ak .eq. 0) then

       do k=1,nz
         apv(1,k) = 0.0
         apv(k,1) = 0.0
       enddo
         apv(1,1) = 1.0

          fr(1) = 0.0
          fi(1) = 0.0

       endif

c
c  update the pressure array
c

        do k=1,nz
         xr(k) = 0.0  
         xi(k) = 0.0 
        enddo


       call jacgs(apv,xr,fr,500,nz,nz,nz,0,cgstol)
       call jacgs(apv,xi,fi,500,nz,nz,nz,0,cgstol)

        do k=1,nz
         ik = i + (k-1)*nxpp
         ipk = ik + 1
         su(ik)  = xr(k)
         su(ipk) = xi(k)
        enddo
c
110     continue
c
c correct the velocity
c                     
      call update(u,w,ub,wb,xw,su,dt)
c
c update velocity boundary conditions
c
      call exbc(su,ubc1,wbc1,ubc2,wbc2)
c
        return
        end
c
c
c*********************************************************************
      subroutine edge(u,w)
c*********************************************************************
      include 'dim_2d.h'
c
      complex u(nxhp,nz),w(nxhp,nz)
      complex zero
      zero = cmplx(0.0,0.0)
c
      do k=1,nz
        u(nxhp,k) = zero
        w(nxhp,k) = zero
       enddo
c
      return
      end
c                                                                       
c*********************************************************************  
      subroutine setup                                                
c*********************************************************************  
      include 'dim_2d.h'
c
      common  /wave/  xw(nxhp),xsq(nxhp),alpha,u00
      common  /leng/  zpts(0:nzm),wg(0:nzm),d(0:nzm,0:nzm)                                 
      common /coeff/  diff(nz,nz), amass(nz)
      common /index/  li(ntot),lj(ntot),lk(ntot)
      common /width/  ybar
c
      common /profile/ ubar(ntot)

      dimension z(nz)
c                                                                       
c  determine the collocation and quadrature points in the z-direction
c  transform coordinates are used, i.e. z = [-1,1]
c
      call jacobl(nzm,0.,0.,zpts,nzm)
      call quad(nzm,zpts,wg,nzm)
      call derv(nzm,zpts,d,nzm)
c
      pi2 = 8.0*atan(1.0)
      y0 = -1.0
      xl = pi2/alpha
c
       dx = xl/nx
c
      do k=1,nz
       z(k) = y0 + 0.5*ybar*(zpts(k-1)+1.0) 
      enddo
c
c  set the discretized differential operators in the z-direction
c
c  diffusion term
c
      do i=0,nzm
       ii = i + 1
       do j=0,nzm
        jj = j + 1
        sum = 0
        do n=0,nzm
         sum = sum + (2./ybar)*wg(n)*d(n,i)*d(n,j)
        enddo                                    
        diff(ii,jj) = sum
       enddo
      enddo
c
c mass matrix term
c
      do i=0,nzm
       ii = i + 1

       amass(ii) = 0.5*ybar*wg(i)
      enddo
c
c  setting the wave number
c
      do 1 i=1,nxhp
         xw(i) = (i - 1)*alpha
         xsq(i) = xw(i)*xw(i) 
    1 continue
      xw(nxhp) = 0.0
      xsq(nxhp) = 0.0


c
c set up the mean velocity profile
c

      z0 = 0.0
c     
      do k=1,nz
         z5 = z0 + 0.5*ybar*(zpts(k-1)+1.0)
         kk = (k-1)*nxpp
         do i=1,nxpp
            ik = kk + i
            ubar(ik) = z5
         enddo
      enddo
      




      return                                                            
      end                                                            

c*********************************************************************  
      subroutine init(u,w,temp)                                       
c*********************************************************************  
      include 'dim_2d.h'
c     
      common  /wave/  xw(nxhp),xsq(nxhp),alpha,u00
      common  /time5/ t,dt,nsteps
      common  /leng/ zpts(0:nzm),wg(0:nzm),d(0:nzm,0:nzm)  
      common  /bc/    ubc(nxpp),wbc(nxpp)
      common  /width/ ybar
      common /scratch/ ox(ntot),oz(ntot)
      common /source/ body(ntot)
      common /exforce/ wavlen
      dimension u(*),w(*),temp(*)
c     
c ASSUMING nz IS ODD
c
      y0 = -1.0
      tol = 1.0e-6
c     
      do i=1,ntot  
         u(i) = 0.  
         w(i) = 0. 
         temp(i) = 0. 
         ox(i) = 0.
      enddo
c     
c     boundary condition along the upper wall
c     
      pi2 = 8.0*atan(1.0)
      xl = pi2/alpha
      dx = xl/nx
      fac5 = pi2/wavlen
      
      
      xcl = 0.5*xl
      hafwav = 0.5*wavlen
      xs = xcl - hafwav
      
      
      print *,'upper boundary condition'
      
      do i=1,nx
         xx = (i-1)*dx
         dis = abs(xx-xcl)
         if( dis .gt. hafwav) then
            ubc(i) = 0.0
         else
            x6 = xx - xs
            x5 = x6*fac5 - pi2/2.
            tan5 = tanh(x5)
            sec5 = 1.0/cosh(x5)
            ubc(i) = u00*(tan5*tan5 - 1.0)
c
c wbc = -du/dx
c
            wbc(i) = -2.0*u00*fac5*tan5*sec5*sec5
c
c     ubc(i)=0.1      
         endif
         print *,xx,ubc(i),wbc(i)
      enddo
      
      do k=2,nz-1
         kk = (k-1)*nxpp
         do i=1,nx
            ik = kk + i
            u(ik) = 0.
            w(ik) = 0.1*random() 
c     
         enddo
      enddo
      
      
      kk = (nz-1)*nxpp
      do i=1,nx
         ik = kk + i
         u(ik) = ubc(i)
      enddo
c     
c--------------------------------------------------
c     transform to spectral space 
c--------------------------------------------------
c     
      call  horfft (u,-1,nz)
      call  horfft (w,-1,nz)
      
c     
      
      call  horfft (ubc,-1,1)
      call  horfft (wbc,-1,1)
      
      
      return
      
      end
      
c*********************************************************************  
      subroutine cwbc(u,dw,xw)
c*********************************************************************
      include 'dim_2d.h'
c
      dimension xw(*)
      complex u(nxhp),dw(nxhp)
c
       do i = 1,nxhp 
         dw(i) = -(0.,1.)*xw(i)*u(i)
       enddo
c
      return
      end      
      
      
c*********************************************************************
      subroutine copy(a,b,n)                                              
c*********************************************************************  
      dimension a(*), b(*)                    
c                                                                       
      do i = 1,n                                                   
       b(i) = a(i)                                               
      enddo
      return                                                            
      end                                                               
c                                                                       
c**********************************************************************
      subroutine restart(nhist,to,uw1,ww1,uw2,ww2)
c**********************************************************************
c                                                                       
      include 'dim_2d.h'
c
      common  /wave/  xw(nxhp),xsq(nxhp),alpha,u00
      common  /leng/  zpts(0:nzm),wg(0:nzm),d(0:nzm,0:nzm)                                 
      common  /time5/  t,dt,nsteps                                                 
      common  /flow/  u(ntot),w(ntot),temp(ntot),                 
     1                un(ntot),wn(ntot),tn(ntot)

      common /param/  re
      common /coeff/  diff(nz,nz), amass(nz)
      common  /bc/    ubc(nxpp),wbc(nxpp)
      dimension uw1(*),ww1(*),uw2(*),ww2(*)
c
      read(nhist) to
      read(nhist) re
      read(nhist) li,lj,lk,lkm,lipp,lijk
c
      if( (li.ne.nx) .or. (lj.ne.ny) .or. (lk.ne.nz) ) then
         print *,'warning !!!!'
         print *,'dimensions are inconsistent'
         print *,'li,lj,lk: ',li,lj,lk
         print *,'nx,ny,nz: ',nx,ny,nz
      endif
c
      li5 = li/2 + 1
c
      read(nhist) alpha                                            
      read(nhist) (xw(i),xsq(i),i=1,li5)                             
      read(nhist) (zpts(i),i=0,lkm)                                     
      read(nhist) (wg(i),i=0,lkm)
      read(nhist) ((d(i,j),i=0,lkm),j=0,lkm)
      read(nhist) ((diff(i,j),i=1,lk),j=1,lk)
      read(nhist) (amass(i),i=1,lk)
      read(nhist) (ubc(i),i=1,lipp)
      read(nhist) (wbc(i),i=1,lipp)
      read(nhist) (uw1(i),i=1,lipp),
     &            (ww1(i),i=1,lipp)
      read(nhist) (uw2(i),i=1,lipp),
     &            (ww2(i),i=1,lipp)
c
c velocity and temperature at (n) time level in spectral-physical space
c
      read(nhist) ( u(i),i=1,lijk)
      read(nhist) ( w(i),i=1,lijk)
c
      return                                                            
      end                                                               
c
c                                                                       
c**********************************************************************
      subroutine output(nhist,t,uw1,ww1,uw2,ww2)                      
c**********************************************************************
c                                                                       
      include 'dim_2d.h'
c
      common  /wave/  xw(nxhp),xsq(nxhp),alpha,u00
      common  /leng/  zpts(0:nzm),wg(0:nzm),d(0:nzm,0:nzm)                                 
      common  /flow/  u(ntot),w(ntot),temp(ntot),                 
     1                un(ntot),wn(ntot),tn(ntot)

      common /coeff/  diff(nz,nz), amass(nz)
      common /param/  re
      common  /bc/    ubc(nxpp),wbc(nxpp)
      dimension uw1(*),ww1(*),uw2(*),ww2(*)
      write(nhist) t
      write(nhist) re
      write(nhist) nx,ny,nz,nzm,nxpp,ntot
      write(nhist) alpha                                            
      write(nhist) (xw(i),xsq(i),i=1,nxhp)                             
      write(nhist) (zpts(i),i=0,nzm)                                     
      write(nhist) (wg(i),i=0,nzm)
      write(nhist) ((d(i,j),i=0,nzm),j=0,nzm)
      write(nhist) ((diff(i,j),i=1,nz),j=1,nz)
      write(nhist) (amass(i),i=1,nz)
      write(nhist) (ubc(i),i=1,nxpp)
      write(nhist) (wbc(i),i=1,nxpp)
      write(nhist) (uw1(i),i=1,nxpp),
     &             (ww1(i),i=1,nxpp)
      write(nhist) (uw2(i),i=1,nxpp),
     &             (ww2(i),i=1,nxpp)
c
c velocity and temperature at (n) time level in spectral-physical space
c
      write(nhist) ( u(i),i=1,ntot)
      write(nhist) ( w(i),i=1,ntot)
c
      write(6,*) 'finished writing to tape ',nhist 
      write(6,*) 'at time ',t
c
      return                                                            
      end  
c
c
c*********************************************************************
      subroutine  horfft (b,is,nk)
c*********************************************************************
      include 'dim_2d.h'
c
      common  /fft/   trigsx(nxf),work(ntot),ifaxx(13)

      dimension b(*),a(nxpp)
c
      if (is.eq.-1) then
c
c  forward transform, i.e. from physcial to spectral space
c
      do 1000 kz=1,nk
       kk = (kz-1)*nxpp
c
       do i=1,nxpp
        a(i) = b(i+kk)
       enddo
c
         call fft991 (a,work,trigsx,ifaxx,1,nxpp,nx,1,is)
c
       do i=1,nxpp
        b(i+kk) = a(i)
       enddo
c
1000   continue
c
      else
c
c  inverse transform, from spectral to physical space
c
      do 2000 kz=1,nk
       kk = (kz-1)*nxpp
c
       do i=1,nxpp
        a(i) = b(i+kk)
       enddo
c
         call fft991 (a,work,trigsx,ifaxx,1,nxpp,nx,1,is)
c
       do i=1,nxpp
        b(i+kk) = a(i)
       enddo
c
2000   continue
c
      endif
c
      return
      end
c
c
c*********************************************************************
      subroutine  dpdz (a,b,d,ndim,ybar)
c*********************************************************************
      include 'dim_2d.h'
c
      dimension a(*),b(*),d(0:ndim,0:ndim)
      dimension sum(nxpp)
c
c finding the derivative of a function in the z-direction 
c
      fac = 2./ybar
c
      do k1=1,nz
       kk1 = (k1-1)*nxpp
c
        do ij=1,nxpp
          sum(ij) = 0.0
        enddo
c
        do k2=1,nz 
         kk2 = (k2-1)*nxpp
c
         do ij=1,nxpp
          sum(ij) = sum(ij) + d(k1-1,k2-1)*a(ij+kk2)
         enddo
c
        enddo
c
         do ij=1,nxpp
          b(ij+kk1) = sum(ij)*fac
         enddo
c
      enddo
c
      return
      end
c
c**********************************************************************
      subroutine legen(al,alp,n,xc,ndim)
c**********************************************************************
      dimension al(0:ndim),alp(0:ndim)
c
      al(0) = 1.
      al(1) = xc
      alp(0) = 0
      alp(1) = 1.
c
      do k=1,n
        kp = k + 1
        km = k - 1
        al(kp) = (2.*k+1.)*xc*al(k)/kp - k*al(km)/kp
      enddo
c
      do k=1,n
        kp = k + 1
        km = k - 1
        alp(kp) = (2.*k+1.)*(al(k)+xc*alp(k))/kp - 
     &            k*alp(km)/kp
      enddo
c
      return
      end
c
c**********************************************************************
      subroutine jacobl(n,alpha,beta,xcol,ndim)
c**********************************************************************
c
c  computes the gauss-lobatto collocation points
c
c   N:              degree of approximation (order of polynomials=n+1)
c   ALPHA:          parameter in Jacobi weight
c   BETA:           parameter in Jacobi weight
c   XJAC:           roots from largest to smallest
c
c
c  for Chebyshev-Gauss-Lobatto points use alpha=-0.5 and beta=-0.5
c  for Legendre-Gauss-Lobatto points use           0             0
c
      dimension xjac(200),xcol(0:ndim)
      common /jacpar/alp,bet,rv
      data kstop/10/
      data eps/1.0e-12/
c
      alp = alpha
      bet = beta
      rv = 1. + alp
      np = n + 1
c
      call jacobf(np,pnp1p,pdnp1p,pnp,pdnp,pnm1p,pdnm1,1.0)
      call jacobf(np,pnp1m,pdnp1m,pnm,pdnm,pnm1m,pdnm1,-1.0)
      det = pnp*pnm1m - pnm*pnm1p
      rp = -pnp1p
      rm = -pnp1m
      a = (rp*pnm1m - rm*pnm1p)/det
      b = (rm*pnp   - rp*pnm)/det
      xjac(1) = 1.0
      nh = (n+1)/2
      dth = 3.14159265/(2*n+1)
      cd = cos(2.*dth)
      sd = sin(2.*dth)
      cs = cos(dth)
      ss = sin(dth)
c
      do 39 j=2,nh
       x = cs
       do 29 k=1,kstop
        call jacobf(np,pnp1,pdnp1,pn,pdn,pnm1,pdnm1,x)
        poly = pnp1 + a*pn + b*pnm1
        pder = pdnp1 + a*pdn + b*pdnm1
        recsum = 0.0
        jm = j - 1
        do 27 i=1,jm
          recsum = recsum + 1.0/(x-xjac(i))
27      continue
28      continue
        delx = -poly/(pder-recsum*poly)
        x = x +delx
        if(abs(delx) .lt. eps) go to 30
29      continue
30      continue
        xjac(j) = x
        cssave = cs*cd - ss*sd
        ss = cs*sd + ss*cd
        cs = cssave
39      continue
        xjac(np) = -1.0
        npp = n + 2
        do 49 i=2,nh
          xjac(npp-i) = -xjac(i)
49      continue
        if(n.ne. 2*(n/2)) then 
        do k=0,n
         kk = n - k + 1
         xcol(k) = xjac(kk)
         enddo
         go to 56
        else
        xjac(nh+1) = 0.0
        do k=0,n
         kk = n - k + 1
         xcol(k) = xjac(kk)
         enddo
         return
        endif
c
56      return
        end
c
c**********************************************************************
        subroutine jacobf(n,poly,pder,polym1,pderm1,polym2,pderm2,x)
c**********************************************************************
        common /jacpar/ alp,bet,rv
        apb = alp + bet
        poly = 1.0
        pder = 0.0
        if(n.eq.0) return
        polylst = poly
        pderlst = pder
        poly = rv*x
        pder = rv
        if(n.eq.1) return
        do 19 k=2,n
          a1 = 2.*k*(k+apb)*(2.*k+apb-2.)
          a2 = (2.*k+apb-1.)*(alp**2-bet**2)
          b3 = (2.*k+apb-2.)
          a3 = b3*(b3+1.)*(b3+2.)
          a4 = 2.*(K+alp-1.)*(k+bet-1.)*(2.*k+apb)
          polyn = ((a2+a3*x)*poly - a4*polylst) / a1
          pdern = ((a2+a3*x)*pder - a4*pderlst + a3*poly) / a1
          psave = polylst
          pdsave = pderlst
          polylst = poly
          poly = polyn
          pderlst = pder
          pder = pdern
19      continue
        polym1 = polylst
        pderm1 = pderlst
        polym2 = psave
        pderm2 = pdsave
        return
        end
c
c**********************************************************************
      subroutine quad(n,x,w,ndim)
c**********************************************************************
      parameter (nn=200)
      dimension x(0:ndim),w(0:ndim),alp1(0:nn),al1(0:nn)
c
c  determine the Gauss Quadrature weighting factors
c
      if( ndim .gt. nn) print *,'dimensioning error in Quad'

      small = 1.0e-30
      do k=0,n
       xc = x(k)
       call legen(al1,alp1,n,xc,nn)  
       w(k) = 2. / 
     &         ( n*(n+1)*al1(n)*al1(n) + small )
      enddo
      return
      end
c
c**********************************************************************
      subroutine derv(nterm,x,d,ndim)
c**********************************************************************
      parameter (nn=200)
      dimension x(0:ndim),d(0:ndim,0:ndim),al1(0:nn),alp1(0:nn),
     &          al2(0:nn),alp2(0:nn)
c
      if( ndim .gt. nn) print *,'dimensioning error in Derv'
c
c  determine the derivative at the collocation points
c
      do i=0,nterm
        xi = x(i)
        call legen(al1,alp1,nterm,xi,nn)  
       do j=0,nterm
        xj = x(j)
        call legen(al2,alp2,nterm,xj,nn)  
        if(i.eq.j) then
         d(i,j) = 0
        else
         d(i,j) = al1(nterm)/(al2(nterm)*(xi-xj))
        endif
       enddo
      enddo
c
      ann = 0.25*nterm*(nterm+1)
      d(0,0) = -ann
      d(nterm,nterm) =  ann
      return
      end
c
c**********************************************************************
      subroutine jacgs(a,x,b,niter,ndx,nx,ny,iprt,tol)
c**********************************************************************
c
c  Conjugate Gradient Method with Jacobi Conditioning
c
c  [a][x] = [b]
c
      parameter (ndim=129)
      dimension a(ndx,*),x(*),b(*)
      dimension r(ndim),p(ndim),q(ndim),apn(ndim),diag(ndim),sum(ndim)
      data small,factor/1.0e-30,1.0e-3/
c
      nxy = nx*ny
c
      if(nx.gt.ndim .or. ny.gt.ndim)
     & print *,'error in CGS, increase array sizes'
c
c
CDIR@ IVDEP
      do i = 1, nx
         diag(i) = 1.0/a(i*(ndx+1)-ndx,1)
         sum(i) = 0.
      enddo
c
      do j=1,ny
       do i=1,nx
        sum(i) = sum(i) + a(i,j)*x(j)
       enddo
      enddo
c
CDIR@ IVDEP
      do i=1,nx
       r(i) = b(i) - sum(i)
       q(i) = r(i)*diag(i)
       p(i) = q(i)
      enddo
c
      res0 = 0.0

      do i=1,nx
        res0 = res0 + r(i)*r(i)
      enddo
c
      res0 = sqrt(res0)
c
      if(res0.le.tol) return
c
      do 1000 iter=1,niter
c
CDIR@ IVDEP
      do i=1,nx
       apn(i) = 0.
      enddo
c
      do j=1,ny
       do i=1,nx
        apn(i) = apn(i) + a(i,j)*p(j)
       enddo
      enddo
c
      qdr = 0
      pdapn = 0
c
CDIR@ IVDEP
      do i=1,nx
       qdr = qdr + q(i)*r(i)
       pdapn = pdapn + p(i)*apn(i)
      enddo
c
      alfa = qdr/(pdapn+small)
c
CDIR@ IVDEP
      do i=1,nx
       x(i) = x(i) + alfa*p(i)
       r(i) = r(i) - alfa*apn(i)
       q(i) = r(i)*diag(i)
      enddo
c
      qdrnew = 0 
c
CDIR@ IVDEP
      do i=1,nx
       qdrnew = qdrnew + q(i)*r(i)
      enddo
c
      beta = qdrnew/(qdr+small)
c 
CDIR@ IVDEP
      do i=1,nx
       p(i) = r(i)*diag(i) + beta*p(i)
      enddo
c    
      res1 = 0
      do i=1,nx
        res1 = res1 + r(i)*r(i)  
      enddo
c
      res1 = sqrt(res1)
      resfac = res1/(res0+small)
c
        if (iprt.eq.1) then
         print 100,iter,resfac 
        endif
c
100   format('*** iter = ',i5,5x,'residual= ',e14.6)
      if( res1.le.tol ) return
1000  continue
c
      print *,'non-convergent in cgs ','residual= ',res1
c
      return
      end
c
c**********************************************************************
      subroutine plotf(iunit,iform,t,ud,wd,td)                                         
c**********************************************************************
c                                                                       
      include 'dim_2d.h'
c
      common  /wave/  xw(nxhp),xsq(nxhp),alpha,u00
      common  /leng/  zpts(0:nzm),wg(0:nzm),d(0:nzm,0:nzm)  
      common  /width/ ybar                               
      dimension  ud(*),wd(*),td(*)                 
      dimension  x(nx),z(nz),u(ntot),
     &           w(ntot),temp(ntot)
c
      pi2 = 8.0*atan(1.0)
      y0 = 0.0
      xl = pi2/alpha
c
      do i=1,nx
       x(i) = (i-1)*xl/nx
      enddo
c
                
c
      do k=1,nz
       z(k) = y0 + 0.5*ybar*(zpts(k-1)+1.0) 
      enddo
c
      call copy(ud,u,ntot)
      call copy(wd,w,ntot)
      call copy(td,temp,ntot)
c
      call horfft(u,1,nz)
      call horfft(w,1,nz)
      call horfft(temp,1,nz)
c
      if(iform.eq.0) then
      write(iunit) t            
      write(iunit) nx,nz,ntot
      write(iunit) (x(i),i=1,nx)
      write(iunit) (z(i),i=1,nz)
      write(iunit) (u(i),i=1,ntot),
     1            (w(i),i=1,ntot),(temp(i),i=1,ntot)
c
      else
c
      write(iunit,100) t            
      write(iunit,110) nx,nz,ntot
      write(iunit,100) (x(i),i=1,nx)
      write(iunit,100) (z(i),i=1,nz)
      write(iunit,100) (u(i),i=1,ntot),
     1            (w(i),i=1,ntot),(temp(i),i=1,ntot)
      endif
c
100   format(5e14.6)
110   format(5i10)
c
      return
      end
c
c-----------------------------------------------------
      function random()
c-----------------------------------------------------
c  Routine returns a pseudo-random number between 0-1. 
c-----------------------------------------------------
      integer m, i, md, seed
c      double precision fmd

c      data m/25173/,i/13849/,md/65536/,fmd/65536.d0/,seed/17/
      data m/25173/,i/13849/,md/65536/,fmd/65536.0/,seed/17/

      save seed

      seed   = mod(m*seed+i,md)
      random = seed/fmd
      return
      end


c*********************************************************************  
      subroutine cdpdx(u,du,xw)
c*********************************************************************
      include 'dim_2d.h'
c
      dimension xw(*)
      complex du(nxhp,nz),u(nxhp,nz)
c
      do k = 1,nz 
        do i = 1,nxhp 
         du(i,k) = (0.,1.)*xw(i)*u(i,k)
        enddo
       enddo
c
      return
      end

c*********************************************************************  
      subroutine horvis(phi,su,xsq,re)
c*********************************************************************  
      include 'dim_2d.h'
c
      dimension xsq(*)
      complex phi(nxhp,nz),su(nxhp,nz)
c
      rei = 1.0/re
c
      do k = 1,nz 
        do i = 1,nxhp                                                   
         su(i,k) = su(i,k) - 
     &             xsq(i)*phi(i,k)*rei
        enddo
      enddo
c
      return
      end
c
c*********************************************************************  
      subroutine stavg(u,w,time,dudz)
c*********************************************************************  
      include 'dim_2d.h'
      dimension u(*),w(*),dudz(*)

c
c compute horizontal averaged shear stress along
c both bottom and top walls as well as channel center
c
c assuming odd number of collocation points
c
        kl = 0
        ku = (nz-1)*nxpp

        suml = 0.
        sumu = 0.

        navg = nx

          do i=1,nx
           ijkl = i 
           ijku = i + ku
           suml = suml + dudz(ijkl)
           sumu = sumu + dudz(ijku)
          enddo
c
      write(10,10) time,suml/navg,sumu/navg
      write(6,10) time,suml/navg,sumu/navg
10    format(3e14.6)
c

      return
      end


c
c*********************************************************************
      subroutine wallbc(u,w,ubc,wbc)
c*********************************************************************
      include 'dim_2d.h'
c
      complex u(nxhp,nz),w(nxhp,nz),ubc(nxhp),wbc(nxhp)
      complex zero
      zero = cmplx(0.0,0.0)
c
       do i=1,nxhp
        u(i,1) = zero
        w(i,1) = zero
        u(i,nz) = ubc(i)
       enddo
c
      return
      end








