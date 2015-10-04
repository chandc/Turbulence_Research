c
c  least squares  discretization
c  using tensor product
c  dual time stepping
c
      parameter (norder=15,nem=10,ndepp=4,
     &           ntdof=norder*norder*ndepp,npm=norder-1,
     &           ndim=norder*norder)
      dimension xp(ntdof,nem),yp(ntdof,nem)
      dimension f(ntdof,nem),fn(ntdof,nem),fnn(ntdof,nem)
      dimension rms(8),temp(ntdof,nem)
      dimension wg(norder),d(norder,norder),zpts(norder)
      dimension wht(nem),wid(nem)
      dimension inorth(nem),isouth(nem),iwest(nem),ieast(nem)
      dimension ibcw(nem),ibce(nem),ibcs(nem),ibcn(nem)
      dimension mask(ntdof,nem),res(ntdof,nem)
      dimension p(ntdof,nem),q(ntdof,nem),apn(ntdof,nem),rin(ntdof,nem)
      dimension diag(ntdof,nem)
      dimension u_rel(ntdof,nem),u_img(ntdof,nem),
     &          v_rel(ntdof,nem),v_img(ntdof,nem)
      dimension u(ndim,nem),v(ndim,nem),pp(ndim,nem),
     &          om(ndim,nem),
     &          un(ndim,nem),vn(ndim,nem),
     &          unn(ndim,nem),vnn(ndim,nem),
     &          fu(ndim,nem),fv(ndim,nem),
     &          u_res(ndim,nem),v_res(ndim,nem),
     &          p_res(ndim,nem),om_res(ndim,nem)
      common /index/ iu,iv,ip,iom
      character*80 fin,fout,frun
      data iu,iv,ip,iom/1,2,3,4/
      data small/1.0e-30/
c
      data re,dt,ntime,nsub,iprt,tol,nitcgs,istart,iform,cgsfac/
     &     7500.,0.01,10,2,0,1.0e-14,1000,0,0,0.01/
c
      namelist /input/fin,fout,re,dt,ntime,nsub,iprt,
     &                tol,nitcgs,istart,frun,iform,cgsfac
c
      fac1 = 1.0
      fac2 =-1.0
      fac3 = 0.0
c
      read(5,input)
      write(6,input)
c
      open(9,file=fout,status='unknown')
c
      open(2,file=fin,form='unformatted')
c
      read(2) nelem,nterm
      read(2) (wht(ne),ne=1,nelem),(wid(ne),ne=1,nelem)
      do ne=1,nelem
       read(2) (xp(k,ne),k=1,nterm),(yp(k,ne),k=1,nterm)
       read(2) iwest(ne),ieast(ne),isouth(ne),inorth(ne)
       read(2) ibcw(ne),ibce(ne),ibcs(ne),ibcn(ne)
      enddo
c
      npoly = nterm - 1
c
      print *,'order of polynomials: ',npoly
      print *,'number of elements: ',nelem
c   
      if(nterm.gt.norder  .or. nelem.gt.nem) then
       print *,'Increase Dimension Size !!'
       print *,'Stop the program'
      endif
c
      call jacobl(npoly,0.,0.,zpts(1),npm)
      call quad(npoly,zpts(1),wg(1),npm)
      call derv(npoly,zpts(1),d(1,1),npm)
c 
      ndep = 4
      time = 0.0
c
      pr = 1./re
c
       neig = nterm*nterm
       nee = neig*ndep
c
c
      if(istart.eq.1) then
       fac1 =  1.5
       fac2 = -2.0
       fac3 =  0.5
       if(iform.eq.1) then
         open(1,file='rstart.dat',status='unknown')
         read(1,144)time
         read(1,143)nelem,neig,nterm,ndep,nee
         do ne=1,nelem
          read(1,144) (xp(k,ne),k=1,nterm),(yp(k,ne),k=1,nterm)
          read(1,144) (fn(i,ne),i=1,nee),(fnn(i,ne),i=1,nee)
          read(1,143) iwest(ne),ieast(ne),isouth(ne),inorth(ne)
          read(1,143) ibcw(ne),ibce(ne),ibcs(ne),ibcn(ne)
         enddo
       else
        open(1,file='rstart.dat',form='unformatted')
c
         read(1)time
         read(1)nelem,neig,nterm,ndep,nee
         do ne=1,nelem
          read(1) (xp(k,ne),k=1,nterm),(yp(k,ne),k=1,nterm)
          read(1) (fn(i,ne),i=1,nee),(fnn(i,ne),i=1,nee)
          read(1) iwest(ne),ieast(ne),isouth(ne),inorth(ne)
          read(1) ibcw(ne),ibce(ne),ibcs(ne),ibcn(ne)
         enddo
c
       endif
       close(1)
       go to 932
       endif
c
      iter = 0
c
      do ne=1,nelem
       do i=1,nee
        f(i,ne) = 0.
        fn(i,ne) = 0.
        fnn(i,ne) = 0.
       enddo 
      enddo
c
c  initial condition
c
      do ne=1,nelem
       do i=1,nterm
        ii = (i-1)*nterm
        do j=1,nterm
         ij = ii + j
         ij = (ij-1)*ndep
         ue = 1.0 - yp(j,ne)**2
         fn(ij+1,ne) = 1.0 
         fn(ij+2,ne) = 0.0
        enddo
       enddo
      enddo
c
c
143      format(10i5)
144      format(5e14.6)
100   format(5e14.6)
c
932   continue
c
      do ne=1,nelem
       do i=1,nee
        f(i,ne) = fn(i,ne)
       enddo
      enddo
c
      do ne=1,nelem
       do j=1,nee
        mask(j,ne) = 1
       enddo
      enddo    
c
c setting up connectivity and b.c. on different elements
c
c  fix pressure at one point
c
      ij = (nterm/2-1)*nterm + 1
      ipc = (ij-1)*ndep + ip
      mask(ipc,1) = 0
c
      do ne=1,nelem
c
      if( iwest(ne).eq.0 ) then  
       do j=1,nterm
        ij = (j-1)*ndep
        mask(ij+iu,ne) = 0
        mask(ij+iv,ne) = 0
       enddo
      endif
c
      if( ieast(ne).eq.0 ) then  
       do j=1,nterm
        ii = (nterm-1)*nterm
        ij = (ii+j-1)*ndep
        mask(ij+iu,ne) = 0
        mask(ij+iv,ne) = 0
       enddo
      endif
c
      if( isouth(ne).eq.0 ) then  
c
       ibg = 1
       iend = nterm
       if(iwest(ne).eq.0) ibg = 2
       if(ieast(ne).eq.0) iend = nterm - 1
c
        do i=ibg,iend
         ii = (i-1)*nterm
         ij = ii*ndep
         mask(ij+iu,ne) = 0
         mask(ij+iv,ne) = 0
        enddo
c
      endif 
c
      if( inorth(ne).eq.0 ) then  
c
       ibg = 1
       iend = nterm
       if(iwest(ne).eq.0) ibg = 2
       if(ieast(ne).eq.0) iend = nterm - 1
c
        do i=ibg,iend
         ii = (i-1)*nterm
         ij = (ii+nterm-1)*ndep
         mask(ij+iu,ne) = 0
         mask(ij+iv,ne) = 0
        enddo
c
      endif   
c
      enddo
c
c
      do 5000 it=1,ntime
c
      time = time + dt
c
c sub-iteration to resolve linearization error
c
      do 1200 im=1,nsub
c
      f(ipc,1) = 0.0
c
      do 10 ne=1,nelem
c
      if( iwest(ne).eq.0 ) then  
       do j=1,nterm
        ij = (j-1)*ndep
        f(ij+iu,ne)  =  0.0
        f(ij+iv,ne)  =  0.0
       enddo
      endif
c
      if( ieast(ne).eq.0 ) then  
       do j=1,nterm
        ii = (nterm-1)*nterm
        ij = (ii+j-1)*ndep
        f(ij+iu,ne)  = 0.0
        f(ij+iv,ne)  = 0.0
       enddo
      endif
c
      if( isouth(ne).eq.0 ) then  
c
       ibg = 1
       iend = nterm
       if(iwest(ne).eq.0) ibg = 2
       if(ieast(ne).eq.0) iend = nterm - 1
c
        do i=ibg,iend
         ii = (i-1)*nterm
         ij = ii*ndep
        f(ij+iu,ne)  = 0.0
        f(ij+iv,ne)  = 0.0
        enddo
c
      endif 
c
      if( inorth(ne).eq.0 ) then  
c
       ibg = 1
       iend = nterm
       if(iwest(ne).eq.0) ibg = 2
       if(ieast(ne).eq.0) iend = nterm - 1
c
        do i=ibg,iend
         ii = (i-1)*nterm
         ij = (ii+nterm-1)*ndep
        f(ij+iu,ne)  = 0.0
        f(ij+iv,ne)  = 0.0
        enddo
c
      endif   
c
10    continue
c
      call dge(nelem,nterm,ndep,ntdof,norder,
     &         fac1,dt,pr,
     &         diag,wid,wht,wg,
     &         f,d)
c
      small = 1.0e-30
c
c calculate incomplete residuals
c
      do ne=1,nelem
      do ij=1,neig
       kk = (ij-1)*ndep
       u(ij,ne) = f(kk+1,ne)
       v(ij,ne) = f(kk+2,ne)
       pp(ij,ne) = f(kk+3,ne)
       om(ij,ne) = f(kk+4,ne)
       un(ij,ne) = fn(kk+1,ne)
       vn(ij,ne) = fn(kk+2,ne)
       unn(ij,ne) = fnn(kk+1,ne)
       vnn(ij,ne) = fnn(kk+2,ne)
       fu(ij,ne) = f(kk+1,ne)
       fv(ij,ne) = f(kk+2,ne)
      enddo
      enddo
c
      call rhs(u,v,pp,om,un,vn,unn,vnn,
     &         fu,fv,
     &         dt,pr,nelem,nterm,neig,norder,
     &         fac1,fac2,fac3,
     &         wid,wht,wg,d,
     &         u_res,v_res,p_res,om_res,ndim)
c
      do ne=1,nelem
       do i=1,neig
        ij=(i-1)*ndep
        res(ij+1,ne) = u_res(i,ne)*mask(ij+1,ne)
        res(ij+2,ne) = v_res(i,ne)*mask(ij+2,ne)
        res(ij+3,ne) = p_res(i,ne)*mask(ij+3,ne)
        res(ij+4,ne) = om_res(i,ne)*mask(ij+4,ne)
       enddo
c
       do i=1,nee
        rin(i,ne) = res(i,ne)
       enddo
c
       enddo  
c
c  forming the complete residuals
c
      call collect(nelem,nterm,ndep,ntdof,
     &             isouth,inorth,iwest,ieast,
     &             res)
c
c  forming the complete diagonals
c
      call collect(nelem,nterm,ndep,ntdof,
     &             isouth,inorth,iwest,ieast,
     &             diag)
c
      qdrold = 1.0e30
c
c L2-Norm of residual
c
      res0 = 0.0
      do ne=1,nelem
       do i=1,nee
        res0 = res0 + res(i,ne)*res(i,ne)
       enddo
      enddo
c
      res0 = sqrt(res0)
c
      print *,im,res0
c
      if(res0 .le. tol) go to 1250
c
c  Preconditioned Conjugate Gradient Method
c
      do 1009 itercg=1,nitcgs
c
c pre-conditioning step
c
      qdr = 0.0
      do ne=1,nelem
c
      do i=1,nee
       q(i,ne) = res(i,ne)*mask(i,ne)/(diag(i,ne)+small)
       qdr = qdr + q(i,ne)*rin(i,ne)*mask(i,ne)
      enddo
      enddo
c
      beta = qdr/qdrold
c
      pdapn = 0
      do ne=1,nelem
      do i=1,nee
       p(i,ne) = (q(i,ne) + beta*p(i,ne))*mask(i,ne)
      enddo
      enddo
c
c forming incomplete [A]{p}
c
      do ne=1,nelem
      do ij=1,neig
       kk = (ij-1)*ndep
       u(ij,ne) = p(kk+1,ne)
       v(ij,ne) = p(kk+2,ne)
       pp(ij,ne) = p(kk+3,ne)
       om(ij,ne) = p(kk+4,ne)
      enddo
      enddo
c
      call lhs(u,v,pp,om,fu,fv,
     &         fac1,
     &         dt,pr,nelem,nterm,neig,norder,
     &         wid,wht,wg,d,
     &         u_res,v_res,p_res,om_res,ndim)
c
      pdapn = 0
c
      do ne=1,nelem
       do i=1,neig
        ij=(i-1)*ndep
        apn(ij+1,ne) = u_res(i,ne)
        apn(ij+2,ne) = v_res(i,ne)
        apn(ij+3,ne) = p_res(i,ne)
        apn(ij+4,ne) = om_res(i,ne)
       enddo
c
      do i=1,nee
       pdapn = pdapn + p(i,ne)*apn(i,ne)*mask(i,ne)
      enddo
      enddo
c
      alfa = qdr/(pdapn+small)
c
      do ne=1,nelem
       do i=1,nee
        f(i,ne)   = f(i,ne) + alfa*p(i,ne)*mask(i,ne)
       enddo                   
      enddo
c
c calculate incomplete residuals
c
      do ne=1,nelem
      do ij=1,neig
       kk = (ij-1)*ndep
       u(ij,ne)   = f(kk+1,ne)
       v(ij,ne)   = f(kk+2,ne)
       pp(ij,ne)  = f(kk+3,ne)
       om(ij,ne)  = f(kk+4,ne)
      enddo
      enddo
c
      call rhs(u,v,pp,om,un,vn,unn,vnn,
     &         fu,fv,
     &         dt,pr,nelem,nterm,neig,norder,
     &         fac1,fac2,fac3,
     &         wid,wht,wg,d,
     &         u_res,v_res,p_res,om_res,ndim)
c
      do ne=1,nelem
       do i=1,neig
        ij=(i-1)*ndep
        rin(ij+1,ne) = u_res(i,ne)*mask(ij+1,ne)
        rin(ij+2,ne) = v_res(i,ne)*mask(ij+2,ne)
        rin(ij+3,ne) = p_res(i,ne)*mask(ij+3,ne)
        rin(ij+4,ne) = om_res(i,ne)*mask(ij+4,ne)
       enddo
c
       do i=1,nee
        res(i,ne) = rin(i,ne)
       enddo
c
       enddo  
c
c  forming the complete residuals
c
      call collect(nelem,nterm,ndep,ntdof,
     &             isouth,inorth,iwest,ieast,
     &             res)
c
      qdrold = qdr
c 
      res1 = 0.0
      do ne=1,nelem
       do i=1,nee
         res1 = res1 + res(i,ne)*res(i,ne)
       enddo
      enddo
c
      res1 = sqrt(res1)
      resfac = res1/res0
c
c
        if (iprt.eq.1) print 105,itercg,res1 
c
105   format('*** iter = ',i5,5x,'residual= ',e14.6)
c
      if(resfac .le. cgsfac  .or. res1 .le. tol) go to 1200
c
1009  continue
c
      print *,'non-convergent in cgs ','residual= ',reso
c
1200   continue
c
1250  continue
c
      do ne=1,nelem
       do ij=1,nee
         fnn(ij,ne) = fn(ij,ne)
         fn(ij,ne) = f(ij,ne)
        enddo
      enddo
c
      fac1 =  1.5
      fac2 = -2.0
      fac3 =  0.5
c
5000  continue
c
       print *,'finished at time= ',time
c
       if(iform.eq.1) then
         open(1,file=frun,status='unknown')
         write(1,144)time
         write(1,143)nelem,neig,nterm,ndep,nee
         do ne=1,nelem
          write(1,144) (xp(k,ne),k=1,nterm),(yp(k,ne),k=1,nterm)
          write(1,144) (fn(i,ne),i=1,nee),(fnn(i,ne),i=1,nee)
          write(1,143) iwest(ne),ieast(ne),isouth(ne),inorth(ne)
          write(1,143) ibcw(ne),ibce(ne),ibcs(ne),ibcn(ne)
         enddo
       else
        open(1,file=frun  ,form='unformatted')
         write(1)time
         write(1)nelem,neig,nterm,ndep,nee
         do ne=1,nelem
          write(1) (xp(k,ne),k=1,nterm),(yp(k,ne),k=1,nterm)
          write(1) (fn(i,ne),i=1,nee),(fnn(i,ne),i=1,nee)
          write(1) iwest(ne),ieast(ne),isouth(ne),inorth(ne)
          write(1) ibcw(ne),ibce(ne),ibcs(ne),ibcn(ne)
         enddo
       endif
c
      stop
      end                      
c
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
      subroutine quad(n,x,w,ndim)
c**********************************************************************
      parameter (nn=100)
      dimension x(0:ndim),w(0:ndim),alp1(0:nn),al1(0:nn)
c
c  determine the Gauss Quadrature weighting factors
c
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
      parameter (nn=100)
      dimension x(0:ndim),d(0:ndim,0:ndim),al1(0:nn),alp1(0:nn),
     &          al2(0:nn),alp2(0:nn)
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
      subroutine rhs(u,v,p,om,un,vn,unn,vnn,
     &               fu,fv,
     &               dt,pr,nelem,nterm,neig,norder,
     &               fac1,fac2,fac3,
     &               wid,wht,wg,d,
     &               c1,c2,c3,c4,ndim)
c**********************************************************************
c
      parameter (ndeg=20, ntdof=ndeg*ndeg)
      dimension dudx(ntdof),dudy(ntdof),dvdx(ntdof),dvdy(ntdof),
     &          dpdx(ntdof),dpdy(ntdof),domdx(ntdof),domdy(ntdof),
     &          dfudx(ntdof),dfudy(ntdof),dfvdx(ntdof),dfvdy(ntdof)
      dimension su(ntdof,4)
      dimension wid(*),wht(*),wg(*),d(norder,*)
      dimension c1(ndim,*),c2(ndim,*),c3(ndim,*),c4(ndim,*)
      dimension u(ndim,*),v(ndim,*),p(ndim,*),om(ndim,*),
     &          un(ndim,*),vn(ndim,*),unn(ndim,*),vnn(ndim,*),
     &          fu(ndim,*),fv(ndim,*)
      dimension li(ndeg)
c
      do n=1,nterm
       li(n) = (n-1)*nterm
      enddo
c
      do 100 ne=1,nelem
c
       ajac = 0.25*wid(ne)*wht(ne)
       facx = 2./wid(ne)
       facy = 2./wht(ne)
c
c derivatives in the x-direction
c
      do i=1,neig
       dudx(i) = 0.0
       dvdx(i) = 0.0
       dpdx(i) = 0.0
       domdx(i) = 0.0
       dfudx(i) = 0.0
       dfvdx(i) = 0.0
      enddo
c
      do n=1,nterm
       do i=1,nterm
        do j=1,nterm
          ij = li(i) + j
          kk = li(n) + j
          dx = d(i,n)*facx
          dudx(ij) = dudx(ij) + dx*u(kk,ne)                       
          dvdx(ij) = dvdx(ij) + dx*v(kk,ne)                       
          dpdx(ij) = dpdx(ij) + dx*p(kk,ne)                       
          domdx(ij) = domdx(ij) + dx*om(kk,ne)
          dfudx(ij) = dfudx(ij) + dx*fu(kk,ne)
          dfvdx(ij) = dfvdx(ij) + dx*fv(kk,ne)
        enddo
       enddo
      enddo                       
c
c derivatives in the y-direction
c
      do i=1,neig
       dudy(i) = 0.0
       dvdy(i) = 0.0
       dpdy(i) = 0.0
       domdy(i) = 0.0
       dfudy(i) = 0.0
       dfvdy(i) = 0.0
      enddo
c
      do n=1,nterm
       do i=1,nterm
        do j=1,nterm
         ij = li(i) + j
         kk = li(i) + n
         dy = d(j,n)*facy
         dudy(ij) = dudy(ij) + dy*u(kk,ne)                       
         dvdy(ij) = dvdy(ij) + dy*v(kk,ne)                       
         dpdy(ij) = dpdy(ij) + dy*p(kk,ne)                       
         domdy(ij) = domdy(ij) + dy*om(kk,ne)
         dfudy(ij) = dfudy(ij) + dy*fu(kk,ne)
         dfvdy(ij) = dfvdy(ij) + dy*fv(kk,ne)
        enddo
       enddo
      enddo
c
c
       do i=1,nterm
        do j=1,nterm
         ij = li(i) + j
         facem = ajac*wg(i)*wg(j)
c
         su(ij,1) =   fac1*u(ij,ne)*facem + dt*( u(ij,ne) +
     &                
c >>> convective term
     &                    fu(ij,ne)*dudx(ij)  + fv(ij,ne)*dudy(ij) +
     &                     u(ij,ne)*dfudx(ij) + v(ij,ne)*dfudy(ij) +
c <<<
     &                     dpdx(ij) + pr*domdy(ij)  )*facem
c
         su(ij,2) =  fac1*v(ij,ne)*facem + dt*( v(ij,ne) +
     &                 
c >>> convective term
     &                    fu(ij,ne)*dvdx(ij)  + fv(ij,ne)*dvdy(ij) +
     &                     u(ij,ne)*dfvdx(ij) + v(ij,ne)*dfvdy(ij) +
c <<<
     &                     dpdy(ij) - pr*domdx(ij)  )*facem
c
         su(ij,3) = ( dudx(ij) + dvdy(ij) )*facem
         su(ij,4) = ( om(ij,ne) + dudy(ij) - dvdx(ij) )*facem
        enddo
       enddo
c
       do i=1,nterm
        do j=1,nterm
         ij = li(i) + j
         facem = ajac*wg(i)*wg(j) 
         su(ij,1) = (   dt*fu(ij,ne) 
     &                 - fac2*un(ij,ne) - fac3*unn(ij,ne)
c >>> convective term
     &                + ( fu(ij,ne)*dfudx(ij) + fv(ij,ne)*dfudy(ij) )*dt
c <<<
c
c mean pressure term added
c
     &               + 2.0*pr*dt
c
     &               )*facem -
     &               su(ij,1)
         su(ij,2) = ( dt*fv(ij,ne)
     &                - fac2*vn(ij,ne) - fac3*vnn(ij,ne)
c >>> convective term
c23456789012345678901234567890123456789012345678901234567890123456789012
     &                + (fu(ij,ne)*dfvdx(ij) + fv(ij,ne)*dfvdy(ij))*dt
c <<<
     &                )*facem -
     &               su(ij,2)
         su(ij,3) = -su(ij,3)
         su(ij,4) = -su(ij,4)
        enddo
       enddo     
c
c Multiply by the transpose
c
c at collocation point
c
       do i=1,nterm
        do j=1,nterm
         ij = li(i) + j
         c1(ij,ne) =   (dt+fac1)                  *su(ij,1) 
c >>> convective term
     &               + dt*dfudx(ij)           *su(ij,1) +
     &                 dt*dfvdx(ij)              *su(ij,2) 
c <<<
c
         c2(ij,ne) =   (dt+fac1)                  *su(ij,2) 
c >>> convective term
     &               + dfudy(ij)*dt           *su(ij,1) +
     &                 dt*dfvdy(ij)           *su(ij,2)
c <<<
c
         c3(ij,ne) = 0.0
c
         c4(ij,ne) = su(ij,4)
        enddo
       enddo 
c
c along the x-direction
c
       do n=1,nterm
        do i=1,nterm
         do j=1,nterm
          ij = li(i) + j
          kk = li(n) + j
          dx = d(n,i)*facx

          c1(ij,ne) =  dx                  *su(kk,3)  +
c >>> convective term
     &                 dt*fu(kk,ne)*dx      *su(kk,1)  +    
c <<<
     &                 c1(ij,ne)

          c2(ij,ne) = -dx                 *su(kk,4)  +
c >>> convective term
     &                dt*fu(kk,ne)*dx      *su(kk,2)  +   
c <<<
     &                c2(ij,ne)

          c3(ij,ne) = dx*dt               *su(kk,1)  +
     &                c3(ij,ne)

          c4(ij,ne) = -pr*dt*dx           *su(kk,2)  +
     &                c4(ij,ne)
         enddo
        enddo
       enddo
c
       do n=1,nterm
        do i=1,nterm
         do j=1,nterm
          ij = li(i) + j
          kk = li(i) + n
          dy = d(n,j)*facy
          c1(ij,ne) =             dy     *su(kk,4)  +
c >>> convective term
     &                dt*fv(kk,ne)*dy     *su(kk,1)  +   
c <<<
     &                c1(ij,ne)
          c2(ij,ne) =             dy     *su(kk,3)  +
c >>> convective term
     &                dt*fv(kk,ne)*dy     *su(kk,2)  +
c <<<     
     &                c2(ij,ne)
          c3(ij,ne) = dt*dy              *su(kk,2)  +
     &                c3(ij,ne)
          c4(ij,ne) = dt*pr*dy           *su(kk,1)  +
     &                c4(ij,ne)
         enddo
        enddo
       enddo
c
100    continue
       return
       end      
c
c
c 
c**********************************************************************
      subroutine lhs(u,v,p,om,
     &               fu,fv,
     &               fac1,
     &               dt,pr,nelem,nterm,neig,norder,
     &               wid,wht,wg,d,
     &               c1,c2,c3,c4,ndim)
c**********************************************************************
c
      parameter (ndeg=20, ntdof=ndeg*ndeg)
      dimension dudx(ntdof),dudy(ntdof),dvdx(ntdof),dvdy(ntdof),
     &          dpdx(ntdof),dpdy(ntdof),domdx(ntdof),domdy(ntdof),
     &          dfudx(ntdof),dfudy(ntdof),dfvdx(ntdof),dfvdy(ntdof)
      dimension su(ntdof,4)
      dimension wid(*),wht(*),wg(*),d(norder,*)
      dimension c1(ndim,*),c2(ndim,*),c3(ndim,*),c4(ndim,*)
      dimension u(ndim,*),v(ndim,*),p(ndim,*),om(ndim,*),
     &          fu(ndim,*),fv(ndim,*)
      dimension li(ndeg)
c
      do n=1,nterm
       li(n) = (n-1)*nterm
      enddo
c
      do 100 ne=1,nelem
c
       ajac = 0.25*wid(ne)*wht(ne)
       facx = 2./wid(ne)
       facy = 2./wht(ne)
c
c derivatives in the x-direction
c
      do i=1,neig
       dudx(i) = 0.0
       dvdx(i) = 0.0
       dpdx(i) = 0.0
       domdx(i) = 0.0
       dfudx(i) = 0.0
       dfvdx(i) = 0.0
      enddo
c
      do n=1,nterm
       do i=1,nterm
        do j=1,nterm
          ij = li(i) + j
          kk = li(n) + j
          dx = d(i,n)*facx
          dudx(ij) = dudx(ij) + dx*u(kk,ne)                       
          dvdx(ij) = dvdx(ij) + dx*v(kk,ne)                       
          dpdx(ij) = dpdx(ij) + dx*p(kk,ne)                       
          domdx(ij) = domdx(ij) + dx*om(kk,ne)
          dfudx(ij) = dfudx(ij) + dx*fu(kk,ne)
          dfvdx(ij) = dfvdx(ij) + dx*fv(kk,ne)
        enddo
       enddo
      enddo                       
c
c derivatives in the y-direction
c
      do i=1,neig
       dudy(i) = 0.0
       dvdy(i) = 0.0
       dpdy(i) = 0.0
       domdy(i) = 0.0
       dfudy(i) = 0.0
       dfvdy(i) = 0.0
      enddo
c
      do n=1,nterm
       do i=1,nterm
        do j=1,nterm
         ij = li(i) + j
         kk = li(i) + n
         dy = d(j,n)*facy
         dudy(ij) = dudy(ij) + dy*u(kk,ne)                       
         dvdy(ij) = dvdy(ij) + dy*v(kk,ne)                       
         dpdy(ij) = dpdy(ij) + dy*p(kk,ne)                       
         domdy(ij) = domdy(ij) + dy*om(kk,ne)
         dfudy(ij) = dfudy(ij) + dy*fu(kk,ne)
         dfvdy(ij) = dfvdy(ij) + dy*fv(kk,ne)
        enddo
       enddo
      enddo
c
c
       do i=1,nterm
        do j=1,nterm
         ij = li(i) + j
         facem = ajac*wg(i)*wg(j)
c
         su(ij,1) =   fac1*u(ij,ne)*facem + dt*( u(ij,ne) +
     &                
c >>> convective term
     &                    fu(ij,ne)*dudx(ij)  + fv(ij,ne)*dudy(ij) +
     &                     u(ij,ne)*dfudx(ij) + v(ij,ne)*dfudy(ij) +
c <<<
     &                     dpdx(ij) + pr*domdy(ij)  )*facem
c
         su(ij,2) =  fac1*v(ij,ne)*facem + dt*( v(ij,ne) +
     &                 
c >>> convective term
     &                    fu(ij,ne)*dvdx(ij)  + fv(ij,ne)*dvdy(ij) +
     &                     u(ij,ne)*dfvdx(ij) + v(ij,ne)*dfvdy(ij) +
c <<<
     &                     dpdy(ij) - pr*domdx(ij)  )*facem
c
         su(ij,3) = ( dudx(ij) + dvdy(ij) )*facem
         su(ij,4) = ( om(ij,ne) + dudy(ij) - dvdx(ij) )*facem
        enddo
       enddo
c
c
c Multiply by the transpose
c
c at collocation point
c
       do i=1,nterm
        do j=1,nterm
         ij = li(i) + j
         c1(ij,ne) =   (dt+fac1)                  *su(ij,1) 
c >>> convective term
     &               + dt*dfudx(ij)           *su(ij,1) +
     &                 dt*dfvdx(ij)              *su(ij,2) 
c <<<
c
         c2(ij,ne) =   (dt+fac1)                  *su(ij,2) 
c >>> convective term
     &               + dfudy(ij)*dt           *su(ij,1) +
     &                 dt*dfvdy(ij)           *su(ij,2)
c <<<
c
         c3(ij,ne) = 0.0
c
         c4(ij,ne) = su(ij,4)
        enddo
       enddo 
c
c along the x-direction
c
       do n=1,nterm
        do i=1,nterm
         do j=1,nterm
          ij = li(i) + j
          kk = li(n) + j
          dx = d(n,i)*facx

          c1(ij,ne) =  dx                  *su(kk,3)  +
c >>> convective term
     &                 dt*fu(kk,ne)*dx      *su(kk,1)  +    
c <<<
     &                 c1(ij,ne)

          c2(ij,ne) = -dx                 *su(kk,4)  +
c >>> convective term
     &                dt*fu(kk,ne)*dx      *su(kk,2)  +   
c <<<
     &                c2(ij,ne)

          c3(ij,ne) = dx*dt               *su(kk,1)  +
     &                c3(ij,ne)

          c4(ij,ne) = -pr*dt*dx           *su(kk,2)  +
     &                c4(ij,ne)
         enddo
        enddo
       enddo
c
       do n=1,nterm
        do i=1,nterm
         do j=1,nterm
          ij = li(i) + j
          kk = li(i) + n
          dy = d(n,j)*facy
          c1(ij,ne) =             dy     *su(kk,4)  +
c >>> convective term
     &                dt*fv(kk,ne)*dy     *su(kk,1)  +   
c <<<
     &                c1(ij,ne)
          c2(ij,ne) =             dy     *su(kk,3)  +
c >>> convective term
     &                dt*fv(kk,ne)*dy     *su(kk,2)  +
c <<<     
     &                c2(ij,ne)
          c3(ij,ne) = dt*dy              *su(kk,2)  +
     &                c3(ij,ne)
          c4(ij,ne) = dt*pr*dy           *su(kk,1)  +
     &                c4(ij,ne)
         enddo
        enddo
       enddo
c
100    continue
       return
       end  

  

c
c**********************************************************
      subroutine dge(nelem,nterm,ndep,ntdof,norder,
     &               fac1,dt,pr,
     &               diag,wid,wht,wg,
     &               f,d)
c**********************************************************
      dimension diag(ntdof,*),f(ntdof,*),d(norder,*),
     &          wid(*),wht(*),wg(*)
      dimension aa(4,4)
      common/index/iu,iv,ip,iom
c
       neig = nterm*nterm
       nee = neig*ndep
c
      do i=1,ndep
       do j=1,ndep
        aa(i,j) = 0.0
       enddo
      enddo
c
      do 1000 ne=1,nelem
c
      do i=1,nee
       diag(i,ne) = 0.0
      enddo
c
      facx = 2./wid(ne)
      facy = 2./wht(ne)
      ajac = 0.25*wid(ne)*wht(ne)
c
      do i=1,nterm
       ii = (i-1)*nterm
       do j=1,nterm    
       ij = ii + j
c
      do k1=1,nterm  
      do k2=1,nterm
c
      facem = ajac*wg(k1)*wg(k2)
c
       kk1 = (k1-1)*nterm + k2
       kk1 = (kk1-1)*ndep
       lu = kk1 + iu
       lv = kk1 + iv
       uo = f(lu,ne)
       vo = f(lv,ne)
c
      dudx = 0.0
      dudy = 0.0
      dvdx = 0.0
      dvdy = 0.0
c
      do m=1,nterm
       mm = (m-1)*nterm
       do n=1,nterm    
       mn = (mm+n-1)*ndep
       lu = mn + iu
       lv = mn + iv
       shapx = 0.0
       shapy = 0.0
       if(k1.eq.m) shapx = 1.0
       if(k2.eq.n) shapy = 1.0
       dx = d(k1,m)*facx*shapy         
       dy = d(k2,n)*facy*shapx 
       dudx = dudx + dx*f(lu,ne)
       dudy = dudy + dy*f(lu,ne)
       dvdx = dvdx + dx*f(lv,ne)
       dvdy = dvdy + dy*f(lv,ne)
       enddo
      enddo
c
          shapx = 0.0
          shapy = 0.0
          if(k1.eq.i) shapx = 1.0
          if(k2.eq.j) shapy = 1.0
          dx = d(k1,i)*facx*shapy         
          dy = d(k2,j)*facy*shapx 
          shape = shapx*shapy
c
          aa(1,1) =  fac1*shape + (shape  
     &               + uo*dx + vo*dy + dudx*shape)*dt
          aa(1,2) =  shape*dt*dudy
          aa(1,3) =  dx*dt
          aa(1,4) =  pr*dy*dt

          aa(2,1) =  shape*dvdx*dt
          aa(2,2) =  fac1*shape + (shape 
     &               + uo*dx + vo*dy + dvdy*shape)*dt
          aa(2,3) =  dy*dt  
          aa(2,4) = -pr*dx*dt

          aa(3,1) = dx
          aa(3,2) = dy

          aa(4,1) =   dy
          aa(4,2) =  -dx
          aa(4,4) =  shape  

c
        ip1 = (ij-1)*ndep
        diag(ip1+1,ne) = diag(ip1+1,ne) + (
     &                  aa(1,1)*aa(1,1) + aa(2,1)*aa(2,1) +
     &                  aa(3,1)*aa(3,1) + aa(4,1)*aa(4,1) )*facem
        diag(ip1+2,ne) = diag(ip1+2,ne) + (
     &                  aa(1,2)*aa(1,2) + aa(2,2)*aa(2,2) +
     &                  aa(3,2)*aa(3,2) + aa(4,2)*aa(4,2) )*facem
        diag(ip1+3,ne) = diag(ip1+3,ne) + (
     &                  aa(1,3)*aa(1,3) + aa(2,3)*aa(2,3) +
     &                  aa(3,3)*aa(3,3) + aa(4,3)*aa(4,3) )*facem
        diag(ip1+4,ne) = diag(ip1+4,ne) + (
     &                  aa(1,4)*aa(1,4) + aa(2,4)*aa(2,4) +
     &                  aa(3,4)*aa(3,4) + aa(4,4)*aa(4,4) )*facem

c
      enddo
      enddo
      enddo
      enddo
c
1000  continue
c
      return
      end
c
c************************************************************************
      subroutine collect(nelem,nterm,ndep,ntdof,
     &                   isouth,inorth,iwest,ieast,
     &                   res)
c************************************************************************
      common /index/ iu,iv,ip,iom
      dimension res(ntdof,*),isouth(*),inorth(*),iwest(*),ieast(*)
c
c bidirectional exchange in the y-direction
c
      do ne=1,nelem
       if(isouth(ne).ne.0) then
       do i=1,nterm
        ii = (i-1)*nterm
        ijs = ii*ndep                          
        ijn = (ii+nterm-1)*ndep
        resu = res(ijs+iu,ne) + res(ijn+iu,isouth(ne))
        resv = res(ijs+iv,ne) + res(ijn+iv,isouth(ne))
        resp = res(ijs+ip,ne) + res(ijn+ip,isouth(ne))
        resm = res(ijs+iom,ne) + res(ijn+iom,isouth(ne))
        res(ijs+iu,ne) = resu
        res(ijs+iv,ne) = resv
        res(ijs+ip,ne) = resp
        res(ijs+iom,ne) = resm
        res(ijn+iu,isouth(ne)) = resu
        res(ijn+iv,isouth(ne)) = resv
        res(ijn+ip,isouth(ne)) = resp
        res(ijn+iom,isouth(ne)) = resm
       enddo
       endif
      enddo
c
c bidirectional exchange in the x-direction
c
      do ne=1,nelem
       if(iwest(ne).ne.0) then
       do j=1,nterm
        ijw = (j-1)*ndep                          
        ii = (nterm-1)*nterm
        ije = (ii+j-1)*ndep
        resu = res(ijw+iu,ne) + res(ije+iu,iwest(ne))
        resv = res(ijw+iv,ne) + res(ije+iv,iwest(ne))
        resp = res(ijw+ip,ne) + res(ije+ip,iwest(ne))
        resm = res(ijw+iom,ne) + res(ije+iom,iwest(ne))


        res(ijw+iu,ne) = resu
        res(ijw+iv,ne) = resv
        res(ijw+ip,ne) = resp
        res(ijw+iom,ne) = resm
        res(ije+iu, iwest(ne)) = resu
        res(ije+iv, iwest(ne)) = resv
        res(ije+ip, iwest(ne)) = resp
        res(ije+iom,iwest(ne)) = resm


       enddo
       endif
      enddo
c
      return
      end
