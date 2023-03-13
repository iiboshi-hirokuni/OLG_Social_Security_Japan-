module mpivars

  implicit none
  !integer                        :: nproc=1,id=0,ierr,n  ! 2020/03/15 by iiboshi
  integer   nproc,id,ierr,n  ! 2020/03/15 by iiboshi

end module mpivars

program trantf

!
! OLG MODEL WITH MULTIPLE TYPES, TRANSITIONAL DYNAMICS, BEQUESTS
!
!
!    Dynamic program of individuals age j, productivity levels eps(j):
!
!       V(j,a) =   max   u(c,1-l) + beta ps(j,j')* V(j',a')
!                 c,l,a'
!
!            s.t. (1+gz) a' = (1 + returns on capital) a + w[t] eps(j) l + pension income(j,type,time)
!                          + health transfers (j,type,time)   - c   -wtax(w[t]*eps(j,k)*l)+ rtra(i[t]*a)
!       see the papers for more details on the specific shape of the households' problem (as health is adjusted
!       for minimum levels of health consumption
!       Also the interest rate is weighted by the fraction held in assets that pay i (endogenous) and a fractions that pays
!       debtRt , which is exogenous
!       where growing variables are detrended relative to technology
!       trend (1+gz)^t. Types differ by age and productivity. If
!       m=number of ages and n=number of productivities, then
!
!           type = 1            age=1, eps=1
!                = 2            age=2, eps=1
!                :
!                = m            age=m, eps=1
!                = m+1          age=1, eps=2
!                = m+2          age=2, eps=2
!                :
!                = 2*m          age=m, eps=2
!                :
!          cd       :
!                = (n-1)*m+1    age=1, eps=n
!                = (n-1)*m+2    age=2, eps=n
!                :
!                = n*m          age=m, eps=n
!
!    Functional forms:
!       u(c,1-l)   = log c + gam log(1-l)
!       F(KT,KI,L) = KT^thetT *KI^thetI * (ZL)^{1-thetT-thetI}, 2 sectors
!
!    Input file requires values for
!
!       alb      = lower bound on financial asset holdings
!       aub      = upper bound on financial asset holdings
!       beta     = discount factor
!       delt1    = depreciation rate on tangible capital, sector 1
!       deli1    = depreciation rate on intangible capital, sector 1
!       delt2    = depreciation rate on tangible capital, sector 2
!       deli2    = depreciation rate on intangible capital, sector 2
!       gam      = parameter in utility
!       gz       = growth rate of technology
!       lam      = degree of annuitization
!       thett1   = tangible capital share, sector 1
!       theti1   = intangible capital share, sector 1
!       thett2   = tangible capital share, sector 2
!       theti2   = intangible capital share, sector 2
!       debt0    = initial debt
!
!       Number of types x 3 (read off steady-state code):
!         hpii   = initial indices of state corresponding to asset holdings
!         popi   = initial population, by type
!         eps    = labor productivity by type
!
!       T x Number of elements in workers' net tax function (wtax)
!         wtax   = (y,a,b) triplets, tax=a+b*y paid on income in [y,y(+1)]
!
!       T x Number of elements in retirees' transfer function (rtra)
!         rtra   = (y,a,b) triplets, tra=a+b*y is transferred for [y,y(+1)]
!
!       T x 9:
!         debt   = government borrowing
!         gspend = government spending
!         ctran  = per capita transfers that are not type-dependent
!         taud   = tax rate on distributions
!         taul   = tax rate on labor (set to 0 to use tax rate schedule)
!         taupc  = tax rate on corporate profits
!         taupu  = tax rate on unincorporated profits (distributed in full)
!         tfp    = market TFP parameter
!         gn     = growth rate of new borns
!
!       T x number of types:
!         ps     = probability of survival from age j to j+1
!         xi     = type-dependent transfers
!         trns   =  age-type transfers
!       T x 4:
!         x(:,1) = initial guess for interest rate
!         x(:,2) = initial guess for wage rate
!         x(:,3) = initial guess for GDP
!         x(:,4) = initial guess for tax rate on consumption
!         x(:,5) = initial guess for debt holdings fraction of total assets
!
!       T x number of bequests:
!         beq    = per capita amount bequested
!         jbeq   = indices for bequesters
!         jinh   = indices for inheritors
!
!    Output file is loaded into Matlab -- run plottranm to view results.
!
!    NOTE: edit this file to change maximum age (nj), the number of
!          productivity levels (ne), the number of households (nh=
!          nj*ne), the number of rates in the tax/transfer schedules
!          (mtx), the maximum years of transition (tmax), the maximum
!          number of years times number of households (hmax=nh*tmax),
!          the number of grid points in the asset grid (na), and the
!          needed storage for optimal decisions (omax=na*nh*tmax/nproc).
!          Note that each processor is assigned a certain number of
!          cohorts (tmax/nproc) and computes dynamic programs for all
!          households in these cohorts (that is, ne*tmax/nproc in all).
!


! by Adrian Peralta-Alva December 2018, building on the program by
! Ellen McGrattan, 3-12-14


  use mpivars
  use mpi
  implicit none
!  include 'mpif.h'
! 以下na = 500 → na = 50 に変更
  integer, parameter             :: ne=4,nj=101,nh=ne*nj,mtx=14,npar=16,     &
                                    na=1000,tmax=240,hmax=nh*tmax,           &
                                    omax=na*nh*10,tbeq=1,nbeq=4,             &
                                    maxit1=15, maxit2=15
!  integer, parameter             :: ne=4,nj=101,nh=ne*nj,mtx=12,npar=16,     &
!                                    na=500,tmax=240,hmax=nh*tmax,           &
!                                    omax=na*nh*8,tbeq=1,nbeq=4

  integer                :: h,i,j,l,t,idone,it,icomp,ntx,ntr,iteracion,&
                                    indic, maxit
  integer, dimension(nh)         :: hpii
  integer, dimension(omax)       :: myo
  integer, dimension(2)          :: ipar
  integer, dimension(tbeq,nbeq)  :: jbeq,jinh
  real(8), dimension(npar)          :: par
  real(8), dimension(nh)            :: eps,popi
  real(8), dimension(tmax,3*mtx)    :: wtaxt,rtrat
  real(8), dimension(tmax,nh)       :: pst,xit,pens,health,income1,income2, &
                                       medrate, ltcare, ltcrate    ! Add 2020/09/27                                       
  real(8), dimension(tmax,13)        :: exog
  real(8), dimension(tmax,5)        :: x,x1,res
  real(8), dimension(hmax,4)        :: allv
  real(8), dimension(tmax,19)       :: ts
  real(8), dimension(tmax,ne)       :: ltaxj,wpsij,rpsij
  real(8), dimension(tmax)         :: agedep,debtH,debtR
  real(8), dimension(tbeq,nbeq)     :: beq
  real(8)                           :: crit,sum1,nrm, tfp_growth_rate
  real,    dimension(ne,3)        ::  prop_agent     ! Add 2020/09/27
  real(8), dimension(tmax,ne)       :: pen_rate      ! Add 2020/10/25
  external trans



   call MPI_INIT(ierr)
   call mpi_COMM_RANK(MPI_COMM_WORLD, id, ierr)
   call mpi_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

  if (id==0) then

   !open(unit=5,  file='trantfWPHGhighertfp.inp')

!   open(unit=5,  file='trantf_JP.inp')
   open(unit=5,  file='./file_in/trantf_XX.inp')
   
   open(unit=7,  file='./file_out/trantf_JP.nxt')
   open(unit=8,  file='./file_out/trantf_JP.dat')
   open(unit=10, file='./file_out/trantf_JP.ts')
   open(unit=11, file='./file_out/trantf_JP.allv')
   
   open(unit=12,  file='./tfp_growth_rate.txt')

!   open(unit=5,  file='trantfWPHGhighertfp.inp')
!   open(unit=7,  file='trantfWPHGhighertfp.nxt')
!   open(unit=8,  file='trantfWPHGhighertfp.dat')
!   open(unit=10, file='trantfWPHGhighertfp.ts')
!   open(unit=11, file='trantfWPHGhighertfp.allv')


  endif


  if (na*hmax/nproc>omax) write(*,*) 'WARNING: Code requires larger omax'
  !
  ! On processor 0, read inputs
  !
  if (id==0) then

!  The following interest rates are the real rates on government debt and can vary over time (for a few years)
!  blelow we set them constant


     do t=1, tmax
     debtR(t)=1.01_8
     end do



    !
    ! Read in parameters
    !
    do i=1,npar
      read(5,*) par(i)
    enddo
    read(5,*) ipar(1)
    read(5,*) ipar(2)
    
    
    !
    !  read  tpf_growth_rate    add 2020/12/05
    !     
    read(12,*) tfp_growth_rate
    

   !
   ! Read in zeta and replacement rate tables
   !
   do i=1,ne
      read(5,*)  (prop_agent(i,j), j=1,3)
   enddo
   
    !
    ! Read in TFP and fiscal exogenous variables
    !
    do t=1,tmax
      read(5,*) exog(t,:)
    enddo
    
    !
    !  set tfp growth rate  
    !
       do t=1, tmax-1
            if (t < 121) then
               exog(t+1,8) = exog(t,8)*(1+tfp_growth_rate)           
           else
               exog(t+1,8) = exog(t,8)    
           endif
           ! pupulation growth rate
             if (t > 16) then
                exog(t,9) = exog(t,9)+0.05     
             elseif (t > 121) then
                exog(t,9) = exog(t,9)                     
             endif
           !  gov spending by logit  
           exog(t,2) = 0.2*exp(0.3*t)/(1+0.2*exp(0.3*t))*0.19                    
       enddo   
    !
    
    
    do i=1,ne
      do t = 1,tmax
          pen_rate(t,i) =  exog(t,i+9)   
      enddo
    enddo

    !
    ! Read in initial asset distributions, populations, productivities
    !
    do h=1,nh
      read(5,*) hpii(h),popi(h),eps(h)
    enddo

    !
    ! Read in the net tax tables
    !
    read(5,*) ntx
    if (ntx>mtx) then
      write(*,*) 'Update trantf.f90 with parameter mtx>=',ntx
      stop
    endif
    do t=1,tmax
      read(5,*) (wtaxt(t,i),i=1,3*ntx)
    enddo

    !
    ! Read in the transfer tables
    !
    read(5,*) ntr
    if (ntr>mtx) then
      write(*,*) 'Update trantf.f90 with parameter mtx>=',ntr
      stop
    endif
    do t=1,tmax
      read(5,*) (rtrat(t,i),i=1,3*ntr)
    enddo

  

    !
    ! Read in survival probabilities and age and time dependent lump sum transfers (including pensions and health)
    !
    write(*,*) 'read transfers:'
    do l=1,ne
      do j=1,nj
        h  = (l-1)*nj+j
        do t=1,tmax
          read(5,*) pst(t,h),xit(t,h),pens(t,h),health(t,h), medrate(t,h), ltcare(t,h), ltcrate(t,h)

        enddo
      enddo
      pst(tmax,l*nj) = 0.
    enddo
    
    
    !  
    !  pension 所得代替率の変更による年金支給額の計算
    !
    do l=1,ne
      do j=1,nj
        h  = (l-1)*nj+j
        do t=1,tmax
                pens(t,h) = pens(t,h) * pen_rate(t,l)/ pen_rate(2,l)
        enddo
      enddo 
    enddo  

    !
    ! Read in bequest information
    !
    if (nbeq>0) then
      do t=1,tbeq
        do j=1,nbeq
          read(5,*) beq(t,j),jbeq(t,j),jinh(t,j)
        enddo
      enddo
    endif

    !
    ! Read in initial guess for the interest rate, wage rate, and transfers
    !
    do t=1,tmax
      read(5,*) x(t,:)
    enddo

    indic=0
    iteracion=1
    debtH(:)=x(:,5)

  endif

  if (nproc>1) then
     call mpi_BCAST(par,npar,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(ipar,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(wtaxt,3*mtx*tmax,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(rtrat,3*mtx*tmax,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(ntx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(ntr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(hpii,nh,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(popi,nh,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(exog,9*tmax,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(eps,nh,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(pst,tmax*nh,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(xit,tmax*nh,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(pens,tmax*nh,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(health,tmax*nh,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(beq,tbeq*nbeq,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(jbeq,tbeq*nbeq,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(jinh,tbeq*nbeq,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(x,5*tmax,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(debtH,tmax,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(debtR,tmax,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(iteracion,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(indic,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     !  Add 2020/09/27 medrate(t,h), ltcare(t,h), ltcrate(t,h)
     call mpi_BCAST(medrate,tmax*nh,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(ltcare,tmax*nh,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(ltcrate,tmax*nh,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     call mpi_BCAST(prop_agent,ne*3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  endif

   call mpi_BARRIER(MPI_COMM_WORLD,ierr)

  !
  ! Compute transitions using nproc processors
  !
  if ((id==0).and.(mod(tmax,nproc).ne.0))                                    &
     write(*,*) 'WARNING: Code assumes tmax is divisible by nproc'

  ! The following do while recomputes interest rates as debtHt affect it
  ! in the case below it runs only once because we have a close enough solution
  ! but it may have to be run several times



  do while(iteracion<2)


   if (id==0) write(*,*) 'to restart, iteracion', iteracion
  !
  ! 0. Restart and update debt levels (after one iteration)
  !
     maxit    = maxit2   ! = 20
!    crit     = 1.e-5
!    idone    = 0
!    it       = 0
!    myo      = -1
    if (id==0) write(*,*) 'Residuals in transition:'


      if(iteracion>1) then

      indic=1

      if(id==0) then
       write(*,*) 'to update debt holdings'

      end if

      call trans(x,nh,ne,nj,tmax,omax,hmax,na,eps,popi,pst,xit,pens,health,tbeq,nbeq,beq,  &
      jbeq,jinh,exog,par,npar,wtaxt,rtrat,mtx,ntx,ntr,ipar,hpii,myo,      &
      res,allv,ts,ltaxj,wpsij,rpsij,agedep,debtH,debtR,income1,income2,indic, &
      medrate, ltcare, ltcrate,prop_agent)
      call mpi_BARRIER(MPI_COMM_WORLD,ierr)

        x1(:,1:3) = .9_8*x(:,1:3)+.1_8*(x(:,1:3)-res(:,1:3))
        nrm = maxval(abs(x1(:,1:3)-x(:,1:3))/(1+abs(x(:,1:3))))
       write(*,'(5x,e10.4)') nrm


      if (id==0) then
            write(*,*) 'after trans'
       end if

      end if


      if(id==0) then
      debtH(:)=x(:,5)
      end if


      call mpi_BCAST(x,5*tmax,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call mpi_BCAST(debtH,tmax,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call mpi_BARRIER(MPI_COMM_WORLD,ierr)


  !
  ! 1. Hold tauc fixed and compute irate
  !

    indic=0
    maxit    = maxit1  ! maxit    = 15
    crit     = 1.e-6 ! = 1.e-6
    idone    = 0
    it       = 0
    myo      = -1

    if (id==0) write(*,*) 'to enter res3 loop'

    do while ((idone==0).and.(it<maxit))

      if (id==0) then
            write(*,*) 'to enter trans'
            write(*,*) id,  'it_1 =',  it
      end if




     call trans(x,nh,ne,nj,tmax,omax,hmax,na,eps,popi,pst,xit,pens,health,tbeq,nbeq,beq,  &
         jbeq,jinh,exog,par,npar,wtaxt,rtrat,mtx,ntx,ntr,ipar,hpii,myo,      &
         res,allv,ts,ltaxj,wpsij,rpsij,agedep,debtH,debtR,income1,income2,indic, &
         medrate, ltcare, ltcrate,prop_agent)

     if (id==0) then
       write(*,*) id, 'resx=', res(2,:)   ! 2020/03/15 by iiboshi
     end if

     x1(:,1:3) = .9_8*x(:,1:3)+.1_8*(x(:,1:3)-res(:,1:3))
     if (id==0) then
        nrm = maxval(abs(x1(:,1:3)-x(:,1:3))/(1+abs(x(:,1:3))))
       write(*,'(5x,e10.4)') nrm

       if (nrm<crit) idone=1
       x(:,1:3)  = x1(:,1:3)
     endif
     it = it+1
    call mpi_BCAST(x,5*tmax,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_BCAST(idone,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_BARRIER(MPI_COMM_WORLD,ierr)

    enddo

    if(id==0) then
    write(*,*) 'debtH,X(5) after irate'
    end if


  !
  ! 2. Take solution for irate, w and compute tauc and debt holdings
  !
     maxit    = maxit2  ! 15
     idone    = 0
     it       = 0
     crit     = 1.e-5
     myo      = -1

    do while ((idone==0).and.(it<maxit))

     if (id==0) then
            write(*,*) id,  'it_2 =',  it
      end if

     call trans(x,nh,ne,nj,tmax,omax,hmax,na,eps,popi,pst,xit,pens,health,tbeq,nbeq,beq,  &
         jbeq,jinh,exog,par,npar,wtaxt,rtrat,mtx,ntx,ntr,ipar,hpii,myo,      &
         res,allv,ts,ltaxj,wpsij,rpsij,agedep,debtH,debtR,income1,income2,indic, &
         medrate, ltcare, ltcrate,prop_agent)

     if (id==0) then
         write(*,*) id, 'resx =', res(1,:)   ! 2020/03/15 by iiboshi
     end if
!      if (it.lt.20) then
!      x1(:,4) = x(:,4)-res(:,4)
!      else
       x1(:,4) = .9_8*x(:,4)+.1_8*(x(:,4)-res(:,4))
       ! x1(:,5) = .99998_8*x(:,5)+.00002_8*(x(:,5)-res(:,5))
!      end if

     if (id==0) then
       ! nrm = maxval(abs(x1(:,4:5)-x(:,4:5))/(1+abs(x(:,4:5))))
        nrm = maxval(abs(x1(:,4)-x(:,4))/(1+abs(x(:,4))))

       write(*,'(5x,e10.4)') nrm
       if (nrm<crit) idone=1
       ! x(:,4:5)  = x1(:,4:5)
       x(:,4)  = x1(:,4)
     endif
     it = it+1
    call mpi_BCAST(x,5*tmax,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_BCAST(idone,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_BARRIER(MPI_COMM_WORLD,ierr)
    enddo

    if( id==0) then
       write(*,*) id, 'debtH,X(5) after tauc'

       do i=1, tmax
        ! write(*,*) debtH(i),x(i,5)
       end do

       iteracion=iteracion+1
       write(*,*) 'iteracion', iteracion
    end if



    call mpi_BCAST(iteracion,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_BARRIER(MPI_COMM_WORLD,ierr)
  enddo


    call mpi_BARRIER(MPI_COMM_WORLD,ierr)
  !
  ! Write out new input file
  !



   if(id==0) then
      write(*,*) 'about to write output'


    write(*,*) 'to write new input file'
    write(7,'(1x,f14.8,6x,a3)') par(1),'alb'
    write(7,'(1x,f14.8,6x,a3)') par(2),'aub'
    write(7,'(1x,f14.8,6x,a4)') par(3),'alpha'
    write(7,'(1x,f14.8,6x,a4)') par(4),'beta'
    write(7,'(1x,f14.8,6x,a5)') par(5),'delt1'
    write(7,'(1x,f14.8,6x,a5)') par(6),'deli1'
    write(7,'(1x,f14.8,6x,a5)') par(7),'delt2'
    write(7,'(1x,f14.8,6x,a5)') par(8),'deli2'
    write(7,'(1x,f14.8,6x,a3)') par(9),'gam'
    write(7,'(1x,f14.8,6x,a2)') par(10),'gz'
    write(7,'(1x,f14.8,6x,a3)') par(11),'lam'
    write(7,'(1x,f14.8,6x,a6)') par(12),'thett1'
    write(7,'(1x,f14.8,6x,a6)') par(13),'theti1'
    write(7,'(1x,f14.8,6x,a6)') par(14),'thett2'
    write(7,'(1x,f14.8,6x,a6)') par(15),'theti2'
    write(7,'(1x,f14.8,6x,a5)') par(16),'debt0'
    write(7,'(1x,i14,6x,a5)')   ipar(1),'idebt'
    write(7,'(1x,i14,6x,a6)')   ipar(2),'ispend'
    write(7,*)
    do h=1,nh
      write(7,'(1x,i14,1x,e14.6,1x,e14.6)') hpii(h),popi(h),eps(h)
    enddo
    write(7,*)
    write(7,'(1x,i14,6x,a19)')  ntx,'Number of tax rates'
    do t=1,tmax
      write(7,'(48(1x,e14.6))') (wtaxt(t,i),i=1,3*ntx)
    enddo
    write(7,'(1x,i14,6x,a19)')  ntr,'Number of transfers'
    do t=1,tmax
      write(7,'(48(1x,e14.6))') (rtrat(t,i),i=1,3*ntr)
    enddo
    write(7,*)
    do t=1,tmax
      write(7,'(10(1x,e14.6))') exog(t,:)
    enddo
    write(7,*)
    do h=1,nh
      do t=1,tmax
        write(7,'(4(1x,e14.6))') pst(t,h),xit(t,h),pens(t,h),health(t,h)
      enddo
    enddo
    write(7,*)
    if (nbeq>0) then
      do t=1,tbeq
        do j=1,nbeq
          write(7,'(1x,e16.6,1x,i3,1x,i3)') beq(t,j),jbeq(t,j),jinh(t,j)
        enddo
      enddo
    endif
    write(7,*)
    do t=1,tmax
      write(7,'(5(1x,e14.6))') x(t,:)
    enddo

    !
    ! Write out results for plottranm.m
    !
    do i=1,npar
      write(8,'(f13.6)') par(i)
    enddo
    write(8,'(i13)') ipar(1)
    write(8,'(i13)') ipar(2)
    write(8,'(i13)') nj
    write(8,'(i13)') ne
    write(8,'(i13)') tmax
    write(8,'(i13)') na
    write(8,'(i13)') ntx
    write(8,'(i13)') ntr
    write(8,'(i13)') nbeq
    write(8,'(i13)') tbeq
    do h=1,nh
      write(8,'(e13.6)') eps(h)
    enddo
    do t=1,tmax
      write(8,'(e13.6)') debtH(t)
    enddo
    do h=1,nh
      do t=1,tmax
        write(8,'(f13.6)') pst(t,h)
      enddo
    enddo
    do h=1,nh
      do t=1,tmax
        write(8,'(f13.6)') xit(t,h)
      enddo
    enddo
    do h=1,nh
      do t=1,tmax
        write(8,'(f13.6)') pens(t,h)
      enddo
    enddo
    do h=1,nh
      do t=1,tmax
        ! write(8,'(f13.6)') health(t,h)
        write(8,'(f13.6)') health(t,h)*medrate(t,h)+ltcare(t,h)*ltcrate(t,h)  ! add 2021/02/20
      enddo
    enddo


    do h=1,nh
      write(8,'(i13)') hpii(h)
    enddo
    do h=1,nh
      write(8,'(f13.6)') popi(h)
    enddo
    do i=1,2
      do t=1,tmax
        write(8,'(f13.6)') exog(t,i)
      enddo
    enddo
    do t=1,tmax
      write(8,'(f13.6)') x(t,4)
    enddo
    do i=4,9
      do t=1,tmax
        write(8,'(f13.6)') exog(t,i)
      enddo
    enddo
    do i=1,3*ntx
      do t=1,tmax
        write(8,'(f13.6)') wtaxt(t,i)
      enddo
    enddo
    do i=1,3*ntr
      do t=1,tmax
        write(8,'(f13.6)') rtrat(t,i)
      enddo
    enddo
    do t=1,tmax
      write(10,'(19(1x,f13.6))') (ts(t,i),i=1,19)
    enddo
    do t=1,hmax
      write(11,'(4(1x,f13.6))') (allv(t,i),i=1,4)
    enddo
    do i=1,19
      do t=1,tmax
        write(8,'(f13.6)') ts(t,i)
      enddo
    enddo
    do i=1,2
      do t=1,tmax
        write(8,'(f13.6)') x(t,i)
      enddo
    enddo
    do t=1,tmax
      write(8,'(f13.6)') exog(t,3)
    enddo
    do t=1,tmax
      write(8,'(f13.6)') x(t,3)
    enddo
    do i=1,4
      do t=1,hmax
        write(8,'(f13.6)') allv(t,i)
      enddo
    enddo
    if (nbeq>0) then
      do t=1,tbeq
        do j=1,nbeq
          write(8,'(f13.6)') beq(t,j)
          write(8,'(i13)')   jbeq(t,j)
          write(8,'(i13)')   jinh(t,j)
        enddo
      enddo
    endif

    do l=1,ne
      do t=1,tmax
        write(8,'(f13.6)') ltaxj(t,l)
      enddo
    enddo
    do l=1,ne
      do t=1,tmax
        write(8,'(f13.6)') wpsij(t,l)
      enddo
    enddo

    do l=1,ne
      do t=1,tmax
        write(8,'(f13.6)') rpsij(t,l)
      enddo
    enddo


  do t=1, tmax
  write(8,'(f13.6)') agedep(t)
    end do


    do h=1,nh
      do t=1,tmax
        write(8,'(f13.6)') income1(t,h)
      enddo
    enddo
   do h=1,nh
      do t=1,tmax
        write(8,'(f13.6)') income2(t,h)
      enddo
    enddo
    end if

  call mpi_FINALIZE(ierr)

end program trantf

subroutine trans(x,nh,ne,nj,tmax,omax,hmax,na,eps,popi,pst,xit,pens,health,tbeq,nbeq,    &
         beq,jbeq,jinh,exog,par,npar,wtaxt,rtrat,mtx,ntx,ntr,ipar,hpii,myo,  &
         resx,allv,ts,ltaxj,wpsij,rpsij,agedep,debtH,debtR,income1,income2,indic, &
         medrate, ltcare, ltcrate,prop_agent)
  !
  ! Starting point is the household overlapping generations matrix (HHOGM):
  !
  !                       Households -->
  !      ------------------------------------------------------------>
  !     |     ages 1 to nj,        |     ages 1 to nj,         |
  !     | productivity level 1     | productivity level 2      | ...
  !      ------------------------------------------------------------>
  !     | 1  2  3  ...        nj   | 1  2  3  ...        nj    |  <== nh in
  !  T  | 2  1  2  ...       nj-1  | 2  1  2  ...       nj-1   |      total
  !  i  | 3  2  1  ...             | 3  2  1  ...              |
  !  m  | 4  3  2  ...             | 4  3  2  ...              |
  !  e  | 5  4  3  ...        :    | 5  4  3  ...        :     |
  !     | 6  5  4             3    | 6  5  4             3     |
  !  |  | :  :  :             2    | :  :  :             2     |
  !  |  |                     1    |                     1     |
  !  v  |                     2    |                     2     |
  !     |                     3    |                     3     |
  !     |                     4    |                     4     |
  !     |                     5    |                     5     |
  !     |                          |                           |
  !     |                          |                           |
  !     |                          |                           |
  !     |tmax-nj   ...             |                           |
  !     | :                        |                           |
  !     |tmax      ...             |                           |
  !      ------------------------------------------------------------>
  !
  !      Note: particular cohorts run along diagonals within each subbox
  !
  !
  ! Store results for households as follows:
  !
  !   Rows:
  !     1:nh           t=1 results
  !        1:              for age=1, eps=1--------|
  !        2:              for age=2, eps=1 .......|.........
  !         :                                      | same   :
  !        nh:             for age=nj, eps=ne      | people :
  !     nh+1:2*nh      t=2 results                 |        :
  !        nh+1:           for age=1, eps=1        |        :
  !        nh+2:           for age=2, eps=1    <---|        :
  !        nh+3:           for age=3, eps=1  <.....|........:
  !        :                                       |        :
  !        2*nh:           for age=nj,eps=ne       |        :
  !     2*nh+1:3*nh    t=3 results                 |        :
  !        nh+1:           for age=1, eps=1        |        :
  !        nh+2:           for age=2, eps=1        |        :
  !        nh+3:           for age=3, eps=1    <---|        :
  !        nh+4:           for age=4, eps=1  <..............:
  !        :
  !        3*nh:           for age=nj,eps=ne
  !        :
  !
  !     hmax-nh*nj
  !        :hmax       Final steady state
  !

  use mpivars
  use mpi
  implicit none
!  include 'mpif.h'
  integer, intent(in)                        :: nh,ne,nj,tmax,omax,hmax,na,  &
                                                mtx,npar,tbeq,nbeq,ntx,ntr,&
                                                indic
  integer, dimension(nh),intent(in)          :: hpii
  integer, dimension(2),intent(in)           :: ipar
  integer, dimension(tbeq,nbeq),intent(in)   :: jbeq,jinh
  integer, dimension(nh)                     :: hpi
  integer, dimension(omax),intent(inout)     :: myo
  integer, dimension(2*ne*tmax/nproc)        :: myi
  integer, dimension(2*ne*tmax/nproc,nproc)  :: itotal
  integer, dimension(1)                      :: kk
  integer                                    :: h,i,j,k,l,s,t,i1,i2,j1,j2,   &
                                                k1,k2,imax,vmax,ic,ic1,      &
                                                ncohort,mopt,kzero,icomw
  real(8), dimension(tmax,5)         :: x
  real(8), dimension(npar),intent(in)           :: par
  real(8), dimension(nh),intent(in)             :: eps,popi
  real(8), dimension(tmax,3*mtx),intent(in)     :: wtaxt,rtrat
  real(8), dimension(tmax,nh),intent(in)        :: pst,xit,pens,health,  &
                                                   medrate, ltcare, ltcrate
  real(8), dimension(tbeq,nbeq),intent(in)      :: beq
  real(8), dimension(tmax,10),intent(in)         :: exog
  real(8), dimension(tmax,5),intent(out)        :: resx
  real(8), dimension(hmax,4),intent(out)        :: allv
  real(8), dimension(tmax,19),intent(out)       :: ts
  real(8), dimension(tmax,ne),intent(out)       :: ltaxj,wpsij,rpsij
  real(8), dimension(tmax)                      :: iratet,wt,ctrant,tkap,rt1t,  &
                                                rt2t,ri1t,ri2t,debtt,        &
                                                gspendt,tauct,taudt,         &
                                                taupct,tauput,tault,tfpt,    &
                                                gnt,ldivl1,ldivl2,kdivlt1,   &
                                                kdivlt2,kdivli1,kdivli2,     &
                                                ydivl,krat1,krat2,kaprat,    &
                                                tranj,taudlt,taupclt,  &
                                                taupult,debttn,popr,gdpt,    &
                                                coef1,tottrns, debtH, debtR, &
                                                hlthtrs
  real(8), dimension(tmax)          :: agedep,old, young
  real(8), dimension(tmax,ne)                   :: fwrk,fret
  real(8), dimension(tmax,nh)                   :: popt,bxit,income1,income2
  real(8), dimension(nj)                        :: irate,w,ctran,pc,taul
  real(8), dimension(nj,3*mtx)                  :: wtax,rtra
  real(8), dimension(3*mtx)                     :: wsched,rsched
  real(8), dimension(nh)                        :: pop,ps,xi
  real(8), dimension(4*hmax/nproc)              :: myv
  real(8), dimension(4*hmax/nproc,nproc)        :: vtotal
  real(8), dimension(nh,4)                      :: tem
  real(8), dimension(12)                         :: inputs
  real(8), dimension(3)                         :: outputs
  real(8), dimension(na)                        :: v,vnext,rhs
  real(8)                                       :: alb,aub,alpha,beta,delt1,    &
                                                deli1,delt2,deli2,gam,gn,    &
                                                gz,thett1,theti1,thett2,     &
                                                theti2,thetl1,thetl2,sum1,   &
                                                a,at,walltime,lab,alph1,     &
                                                lam,debt0,taxinc,ltax,psi,   &
                                                mrate
  real, dimension(ne,3)        ::  prop_agent
  external foc

  walltime = -MPI_WTIME()

  !
  ! Initialize parameters
  !
  alb                 = par(1)
  aub                 = par(2)
  alpha               = par(3)
  beta                = par(4)
  delt1               = par(5)
  deli1               = par(6)
  delt2               = par(7)
  deli2               = par(8)
  gam                 = par(9)
  gz                  = par(10) ! Omega growth rate of labor tech level
  lam                 = par(11)
  thett1              = par(12)
  theti1              = par(13)
  thett2              = par(14)
  theti2              = par(15)
  debt0               = par(16)
  alph1               = 1.-alpha
  thetl1              = 1.-thett1-theti1
  thetl2              = 1.-thett2-theti2
  inputs(1)           = gam
  imax                = 2*ne*tmax/nproc
  vmax                = 4*hmax/nproc
  icomw               = MPI_COMM_WORLD
  kzero               = floor((aub-alb*na)/(aub-alb))

  !
  ! Exogenous variables, prices, and transfers
  !
  debttn              = 0.33*exog(:,1)
  gspendt             = exog(:,2)
  ctrant              = exog(:,3)
  taudt               = exog(:,4)
  tault               = exog(:,5)
  taupct              = exog(:,6)
  tauput              = exog(:,7)
  tfpt                = exog(:,8)
  gnt                 = exog(:,9) ! gnowth rate of population
  iratet              = x(:,1)
  wt                  = x(:,2)
  gdpt                = x(:,3)
  tauct               = x(:,4)
  debtH               = x(:,5)

  taudlt(1)           = taudt(1)
  taupclt(1)          = taupct(1)
  taupult(1)          = tauput(1)
  do t=2,tmax
    taudlt(t)         = taudt(t-1)
    taupclt(t)        = taupct(t-1)
    taupult(t)        = tauput(t-1)
  enddo

  !
  ! Compute fraction of population, by type, relative to population trend
  !
  popt                = 0.
  popt(1,:)           = popi
  do l=1,ne
    do t=2,tmax
      h               = (l-1)*nj+1
      popt(t,h)       = (1.+gnt(t-1))*popt(t-1,h)
    enddo
    do t=2,tmax
      do j=2,nj
        h             = (l-1)*nj+j
        popt(t,h)     = popt(t-1,h-1)*pst(t-1,h-1)
      enddo
    enddo
  enddo

  do t=1,tmax
  agedep(t)=0.0_8
  old(t)=0.0_8
  young(t)=0.0_8
  end do



  do l=1,ne
    do t=1,tmax
      fwrk(t,l) = 0.
      fret(t,l) = 0.

      old(t)   =old(t)+sum(popt(t,(l-1)*nj+47:(l-1)*nj+nj))
      young(t) =young(t)+sum(popt(t,(l-1)*nj+1:(l-1)*nj+46))

      do j=1,nj
        h             = (l-1)*nj+j
        if (eps(h)>0.) then
          fwrk(t,l)   = fwrk(t,l)+popt(t,h)/sum(popt(t,:))
        else
          fret(t,l)   = fret(t,l)+popt(t,h)/sum(popt(t,:))
        endif
      enddo
    enddo
  enddo
        do t=1, tmax
        agedep(t)=old(t)/young(t)
        end do

  !
  ! Compute intergenerational transfers
  !
  bxit         = xit
  do t=1,tmax
    do l=1,nbeq
      if (tbeq==1) then
        bxit(t,jbeq(1,l))  = bxit(t,jbeq(1,l))-beq(1,l)
        bxit(t,jinh(1,l))  = bxit(t,jinh(1,l))+popt(t,jbeq(1,l))*beq(1,l)/popt(t,jinh(1,l))
      else
        bxit(t,jbeq(t,l))  = bxit(t,jbeq(t,l))-beq(t,l)
        bxit(t,jinh(t,l))  = bxit(t,jinh(t,l))+popt(t,jbeq(t,l))*beq(t,l)/popt(t,jinh(t,l))
      endif
    enddo
  enddo

  !
  ! Split up the task of computing the problems of each cohort
  !
  !   Each faces prices t to t+nj-s
  !     * If born before t=1, s>1, t=1
  !     * If born at or after, s=1, t>=1
  !
  !   Processors/Cohorts (ncohort=tmax/nproc):
  !     * Proc 0:  computes dynamic program for ne*(ncohort-1) households
  !     * Proc >0: computes dynamic program for ne*ncohort households
  !
  !   Example: nj = 101, ncohort = 5, then the tasks are:
  !     * Proc 0:  ages 98-101, t=1
  !     * Proc 1:  ages 93-97, t=1
  !        :
  !     * Proc 19: ages 3-7, t=1
  !     * Proc 20: ages 1-2, t=1, ages 2-4, t=2:4
  !     * Proc 21: ages 1-5, t=5:9
  !       :
  !     * Proc 47: ages 1-5, t=135:139
  !   with each processor doing this for all productivity levels (ne)
  !

  ! ncohort             =  6  ! = tmax/nproc ! 2020/03/15 by iiboshi
  ncohort             = tmax/nproc ! 2020/03/15 by iiboshi
  myv                 = 0.
  myi                 = 0
  ic1                 = 1
  if (id==0) ic1      = 2
  do ic=ic1,ncohort
    s                 = max(1,nj-id*ncohort-ic+2)
    t                 = max(1,id*ncohort+ic-nj)
    irate(s:nj)       = iratet(t:t+nj-s)
    ctran(s:nj)       = ctrant(t:t+nj-s)
    w(s:nj)           = wt(t:t+nj-s)
    taul(s:nj)        = tault(t:t+nj-s)
    pc(s:nj)          = 1.+tauct(t:t+nj-s)
    wtax(s:nj,:)      = wtaxt(t:t+nj-s,:)
    rtra(s:nj,:)      = rtrat(t:t+nj-s,:)
    do l=1,ne
      ! write(*,*) 'l =', l
      j1              = 2*(l-1)*ncohort +2*(ic-1)
      myi(j1+1)       = s
      myi(j1+2)       = t
      do j=0,nj-s
        h             = (l-1)*nj+s+j
        ps(h)         = pst(t+j,h)
        if (eps(h)>0.) then
          xi(h)       = bxit(t+j,h)
        else
          xi(h)       = bxit(t+j,h)
        endif
      enddo

      !
      ! Iterate backwards for ages j=nj to s, solving
      ! dynamic programs recursively
      !
      h             = (l-1)*nj+s !+nj ! 2020/03/15 by iiboshi
      ! write(*,*) 'h =', h
      ! write(*,*) 's =', s
      inputs(2)       = pc(nj)
      inputs(3)       = 1.-taul(nj)
      inputs(4)       = w(nj)*eps(l*nj)
      !inputs(5)       = ctran(nj)+pens(nj,h)+health(nj,h)+xi(l*nj)
      inputs(5)        = ctran(nj) + pens(nj,h) + medrate(nj,h)*health(nj,h)  &       !add 2020/07/19
                     + ltcrate(nj,h)*ltcare(nj,h) + xi(l*nj)              !add 2020/08/18
      !inputs(8)       = health(nj,h)
      inputs(8)     = health(nj,h) &            !add 2020/05/20
                     + ltcare(nj,h)            !add 2020/08/18
      inputs(10)  = prop_agent(l,1) ! zeta add 2020/08/22
      inputs(11)  = prop_agent(l,3) ! gam add 2020/08/22
      wsched          = wtax(nj,1:3*mtx)
      rsched          = rtra(nj,1:3*mtx)

      do i=1,na
        !write(*,*) 'i=', i  ! 2020/03/15 by iiboshi
        a             = aub*float(i-1)/float(na-1)
        inputs(6)     = ((1.-debtH(nj))*irate(nj)+debtH(nj)*debtR(nj)-1.)*a
        inputs(7)     = a
        call foc(inputs,wsched,rsched,mtx,ntx,ntr,  outputs)
        v(i)          = outputs(1)
        j1            = (l-1)*ncohort*nj*na+(ic-1)*nj*na+(nj-1)*na+i
        !write(*,*) 'jl =', j1  ! 2020/03/15 by iiboshi
        myo(j1)       = kzero
      enddo
      do j=nj-1,s,-1
        !
        ! For age j, find value functions
        !
        h             = (l-1)*nj+j
        vnext         = v
        inputs(2)     = pc(j)
        inputs(3)     = 1.-taul(nj)
        inputs(4)     = w(j)*eps(h)
        !inputs(5)     = ctran(j)+xi(h)+pens(j,h)+health(j,h)
        inputs(5)        = ctran(j) +xi(h)+ pens(j,h) + medrate(j,h)*health(j,h)  &       !add 2020/07/19
                     + ltcrate(j,h)*ltcare(j,h)               !add 2020/08/18
        !inputs(8)     = health(j,h)
        inputs(8)     = health(j,h) &            !add 2020/05/20
                     + ltcare(j,h)            !add 2020/08/18
        wsched        = wtax(j,1:3*mtx)
        rsched        = rtra(j,1:3*mtx)

        do i=1,na
          a           = (aub*float(i-1)+alb*float(na-i))/float(na-1)
          j1          = (l-1)*ncohort*nj*na+(ic-1)*nj*na+(j-1)*na+i
          mopt        = myo(j1)
          if (mopt>0) then
            k1        = max(1,floor(float(mopt)-max(2.,.1*na)))
            k2        = min(floor(float(mopt)+max(2.,.1*na)),na)
          else
            k1        = 1
            k2        = na
          endif
          rhs         = -1.0e+9
          do k=k1,k2
            at        = (aub*float(k-1)+alb*float(na-k))/float(na-1)
            inputs(6) = ((1.-debtH(j))*irate(j)+debtH(j)*debtR(j)-1.)*a
            inputs(7) = a-(1.+gz)*at*(1.-lam*(1.-ps(h)))
            call foc(inputs,wsched,rsched,mtx,ntx,ntr, outputs)
            rhs(k)    = outputs(1)+beta*ps(j)*vnext(k)
          enddo
          kk          = maxloc(rhs)
          k           = kk(1)
          v(i)        = rhs(k)
          myo(j1)     = k
        enddo
      enddo
    enddo

    do l=1,ne
      h               = (l-1)*nj+s
      hpi(h)          = hpii(h)

      !
      ! Update distributions of asset holdings
      !

      do j=s+1,nj
        h             = (l-1)*nj+j
        j1            = (l-1)*ncohort*nj*na +(ic-1)*nj*na+(j-2)*na+hpi(h-1)
        hpi(h)        = myo(j1)
      enddo
      !
      ! Store optimal decisions in myv(.):
      !   |                                |
      !   | data for eps = 1               |
      !   |   data for first cohort        |
      !   |     outputs for age s          |
      !   |     outputs for age s+1        |
      !   |            :                   |
      !   |     outputs for age nj         |
      !   |   data for second cohort       |
      !   |     outputs for age s          |
      !   |     outputs for age s+1        |
      !   |            :                   |
      !   |     outputs for age nj         |
      !   |   data for third cohort        |
      !   |      :                         |
      !   | data for eps = 2               |
      !   |      :                         |
      !
      !
     ! write(*,*) 's 1 =', s
      do j=s,nj
        h             = (l-1)*nj+j
        i             = hpi(h)
        a             = (aub*float(i-1)+alb*float(na-i))/float(na-1)
        if (j<nj) then
          k           = hpi(h+1)
          at          = (aub*float(k-1)+alb*float(na-k))/float(na-1)
        else
          at          = 0.
        endif

        wsched        = wtax(j,1:3*mtx)
        rsched        = rtra(j,1:3*mtx)
        inputs(2)     = pc(j)
        inputs(3)     = 1.-taul(j)
        inputs(4)     = w(j)*eps(h)
        !inputs(5)     = ctran(j)+xi(h)+pens(j,h)+health(j,h)
         inputs(5)        = ctran(j) + pens(j,h) + medrate(j,h)*health(j,h)  &       !add 2020/07/19
                     + ltcrate(j,h)*ltcare(j,h) +xi(h)              !add 2020/08/18
        !inputs(8)     = health(j,h)
        inputs(8)     = health(j,h) &            !add 2020/05/20
                     + ltcare(j,h)            !add 2020/08/18
        inputs(6)     = ((1.-debtH(j))*irate(j)+debtH(j)*debtR(j)-1.)*a
        inputs(7)     = a-(1.+gz)*at*(1.-lam*(1.-ps(h)))

        inputs(10)  = prop_agent(l,1) ! zeta add 2020/08/22
        inputs(11)  = prop_agent(l,3) ! gam add 2020/08/22

        call foc(inputs,wsched,rsched,mtx,ntx,ntr, outputs)

        j1            = 4*(l-1)*ncohort*nj+ 4*((ic-1)*nj+j-1)
        myv(j1+1)     = a
        myv(j1+2)     = outputs(2)
        myv(j1+3)     = outputs(3)
        myv(j1+4)     = at
      enddo
    enddo
  enddo


  walltime = walltime+MPI_WTIME()
  write(*,*) 'my id/walltime  ',id,walltime
  !call flush

  call mpi_GATHER(myv,vmax,MPI_REAL8,  vtotal,vmax,MPI_REAL8,  0,icomw,ierr)
  call mpi_GATHER(myi,imax,MPI_INTEGER,itotal,imax,MPI_INTEGER,0,icomw,ierr)
!        vtotal(:,1) = myv ! 2020/03/15 by iiboshi
!        itotal(:,1) = myi ! 2020/03/15 by iiboshi

  !
  ! Collect results into matrix allv(.,i), which for variable i, where
  !
  !       allv(:,i) = [HHOGM row 1; HHOGM row 2; ...HHOGM row tmax], i=1:4
  !
  !write(*,*) 's = 2', s
  if (id==0) then
    do n=0,nproc-1
      ic1             = 1
      if (n==0) ic1   = 2
      do ic=ic1,ncohort
        do l=1,ne
          j1          = 2*(l-1)*ncohort +2*(ic-1)
          s           = itotal(j1+1,n+1)
          t           = itotal(j1+2,n+1)
          do j=1,nj
            j1        = 4*(l-1)*ncohort*nj+ 4*((ic-1)*nj+j-1)+1
            j2        = 4*(l-1)*ncohort*nj+ 4*((ic-1)*nj+j)
            tem(j,:)  = vtotal(j1:j2,n+1)
          enddo
          j1            = 1
          do j=(t-1)*nh+(l-1)*nj+s,(t+nj-s-1)*nh+l*nj,nh+1
           ! write(*,*) 's =', s
            allv(j,:)   = tem(j1+s-1,:)
            j1          = j1+1
          enddo
        enddo
      enddo
    enddo
  endif
  call mpi_BARRIER(MPI_COMM_WORLD,ierr)

  !
  ! Fill in steady state results for generations near the end
  !
  resx                  = 0.
  if (id==0) then
    do l=1,ne
      k = 0
      do t=tmax-nj+1,tmax
        k               = k+1
        i1              = (t-1)*nh+(l-1)*nj+1
        i2              = (t-1)*nh+(l-1)*nj+k
        j1              = (tmax-nj-1)*nh+(l-1)*nj+1
        j2              = (tmax-nj-1)*nh+(l-1)*nj+k
        allv(i1:i2,:)   = allv(j1:j2,:)
      enddo
    enddo
    !
    ! Adding everything up for each date and put results in ts(.,.)
    !   ts(.,1)       = KT, sector 1 (detrended by gz,pop)
    !   ts(.,2)       = KT, sector 2 (")
    !   ts(.,3)       = KI, sector 1 (")
    !   ts(.,4)       = KI, sector 2 (")
    !   ts(.,5)       = L, sector 1  (detrended by pop)
    !   ts(.,6)       = L, sector 2  (")
    !   ts(.,7)       = C (detrended by gz,pop)
    !   ts(.,8)       = Y (")
    !   ts(.,9)       = XT, sector 1 (")
    !   ts(.,10)      = XT, sector 2 (")
    !   ts(.,11)      = XI, sector 1 (")
    !   ts(.,12)      = XI, sector 2 (")
    !   ts(.,13)      = Corporate profits (")
    !   ts(.,14)      = Noncorporate profits (")
    !   ts(.,15)      = Corporate dividends (")
    !   ts(.,16)      = GDP (")
    !   ts(.,17)      = B'-(1+i)*B (")
    !   ts(.,18)      = Labor tax
    !   ts(.,19)      = Transfers associated with tax/transfer functions
    !
    rt1t              = ((1.-taudlt)*iratet/(1.-taudt)-1.)/(1.-taupct)+delt1
    rt2t              = (iratet-1.)/(1.-tauput)+delt2
    ri1t              = (1.-taudlt)*(1.-taupclt)*iratet/((1.-taudt)*(1.-     &
                        taupct))-1.+deli1
    ri2t              = (1.-taupult)*iratet/(1.-tauput)-1.+deli2
    ldivl1            = thetl1*alpha/(thetl1*alpha+thetl2*alph1)
    ldivl2            = thetl2*alph1/(thetl1*alpha+thetl2*alph1)
    kdivlt1           = thett1*wt/(thetl1*rt1t)
    kdivlt2           = thett2*wt/(thetl2*rt2t)
    kdivli1           = theti1*wt/(thetl1*ri1t)
    kdivli2           = theti2*wt/(thetl2*ri2t)
    krat1             = kdivli1/kdivlt1
    krat2             = kdivli2/kdivlt2
    kaprat            = (kdivlt2/kdivlt1)*(ldivl2/ldivl1)
    coef1             = (1.-taudlt)*(1.+(1.-taupclt)*krat1)+(1.+(1.-         &
                        taupult)*krat2)*kaprat
    ydivl             = ((tfpt*kdivlt1**thett1*kdivli1**theti1*ldivl1)**     &
                        alpha*(tfpt*kdivlt2**thett2*kdivli2**theti2*         &
                        ldivl2)**alph1)*2.
    debtt(2:tmax)     = debttn(1:tmax-1)
    popr(1:tmax-1)    = sum(popt(2:tmax,:),dim=2)/sum(popt(1:tmax-1,:),dim=2)
    popr(tmax)        = popr(tmax-1)


    do t=1,tmax
      i1              = (t-1)*nh+1
      i2              = t*nh
      pop             = popt(t,:)
      ps              = pst(t,:)
      lab             = sum(pop*allv(i1:i2,3)*eps)/sum(pop)
      ts(t,18)        = 0.
      ts(t,19)        = 0.
      ltaxj(t,:)      = 0.
      wpsij(t,:)      = 0.
      rpsij(t,:)      = 0.
      tottrns     = 0.
      hlthtrs     = 0.
      wsched          = wtaxt(t,1:3*mtx)
      rsched          = rtrat(t,1:3*mtx)

      if (tault(t)>0.) then
        ts(t,18)      = tault(t)*wt(t)*lab

      else
        h             = 1
        l             = 1

        do i=i1,i2

          if (eps(h)>0.) then
            taxinc      = wt(t)*allv(i,3)*eps(h)
            call comptax(taxinc,wsched,mtx,ntx, ltax,mrate)
            ts(t,18)    = ts(t,18)+pop(h)*mrate*taxinc/sum(pop)
            ts(t,19)    = ts(t,19)+pop(h)*(mrate*taxinc-ltax)/sum(pop)
            ltaxj(t,l)  = ltaxj(t,l)+pop(h)*mrate*taxinc/sum(pop)
            wpsij(t,l)  = wpsij(t,l)+pop(h)*(mrate*taxinc-ltax)/sum(pop)
            psi=0.0_8
          else
            ! taxinc      = ((1.-debtH(t))*irate(t)+debtH(t)*debtR(t)-1.)*allv(i,1) ! 2020/03/15 irate-> iratet
            taxinc      = ((1.-debtH(t))*iratet(t)+debtH(t)*debtR(t)-1.)*allv(i,1)  ! 2020/03/15 irate-> iratet
            call comptax(taxinc,rsched,mtx,ntr, psi,mrate)
            ts(t,19)    = ts(t,19)+pop(h)*psi/sum(pop)
            rpsij(t,l)  = rpsij(t,l)+pop(h)*psi/sum(pop)
          endif

          income1(t,h)=wt(t)*allv(i,3)*eps(h)
!+((1.-debtH(t))*irate(t)+debtH(t)*debtR(t)-1.)*&
!          allv(i,1)+allv(i,1)
         ! write(*,*) 't = ', t  ! 2020/03/15 ctran->ctrant
!          income2(t,h)=income1(t,h)+ctrant(t)+pens(t,h)+health(t,h)+psi-ltax   ! 2020/03/15 ctran->ctrant
         income2(t,h)=income1(t,h)+ctran(t)+pens(t,h)+health(t,h)+psi-ltax

          if (h==nj*l) then
            l         = l+1
          endif
          h           = h+1
        enddo
      endif
      ts(t,8)         = ydivl(t)*lab


      if (t==1) then
!        ts(t,1)       = (sum(pop*allv(i1:i2,1))/sum(pop)-debt0)/coef1(1)
       ts(t,1)       = (sum(pop*allv(i1:i2,1))/sum(pop)-debt0)/coef1(1)
      else
        if (ipar(1)==0) then
          ts(t,1)     = (sum(pop*allv(i1:i2,1))/sum(pop)-debtt(t))/coef1(t)
        else
!          ts(t,1)     = (sum(pop*allv(i1:i2,1))/sum(pop)-debtt(t)*gdpt(t))/  &
!                        coef1(t)
ts(t,1)     = (sum(pop*allv(i1:i2,1))/sum(pop)-debtt(t)*gdpt(t))/  &
                        coef1(t)

       endif
      endif
      ts(t,5)         = ldivl1(t)*lab
      ts(t,6)         = ldivl2(t)*lab
      ts(t,7)         = sum(pop*allv(i1:i2,2))/sum(pop)
      tkap(t)         = sum(pop*allv(i1:i2,4)*(1.-ps)*(1.-lam))/sum(pop)

    enddo
    ts(:,2)           = kaprat*ts(:,1)
    ts(:,3)           = krat1*ts(:,1)
    ts(:,4)           = krat2*ts(:,2)
    ts(1:tmax-1,9)    = (1.+gz)*popr(1:tmax-1)*ts(2:tmax,1)-(1.-delt1)*      &
                        ts(1:tmax-1,1)
    ts(tmax,9)        = ts(tmax-1,9)
    ts(1:tmax-1,10)   = (1.+gz)*popr(1:tmax-1)*ts(2:tmax,2)-(1.-delt2)*      &
                        ts(1:tmax-1,2)
    ts(tmax,10)       = ts(tmax-1,10)
    ts(1:tmax-1,11)   = (1.+gz)*popr(1:tmax-1)*ts(2:tmax,3)-(1.-deli1)*      &
                        ts(1:tmax-1,3)
    ts(tmax,11)       = ts(tmax-1,11)
    ts(1:tmax-1,12)   = (1.+gz)*popr(1:tmax-1)*ts(2:tmax,4)-(1.-deli2)*      &
                        ts(1:tmax-1,4)
    ts(tmax,12)       = ts(tmax-1,12)
    ts(:,13)          = alpha*ts(:,8)-wt*ts(:,5)-delt1*ts(:,1)-ts(:,11)
    ts(:,14)          = alph1*ts(:,8)-wt*ts(:,6)-delt2*ts(:,2)-ts(:,12)
    ts(:,15)          = alpha*ts(:,8)-wt*ts(:,5)-ts(:,9)-ts(:,11)-taupct*    &
                        ts(:,13)
    ts(:,16)          = ts(:,8)-ts(:,11)-ts(:,12)
    if (ipar(1)==0) then
      debtt(1)        = debt0
      ts(:,17)        = (1.+gz)*popr*debttn-debtR*debtt
    else
      debtt(1)        = debt0/ts(1,16)
      ts(1:tmax-1,17) = (1.+gz)*popr(1:tmax-1)*debttn(1:tmax-1)*             &
                        ts(2:tmax,16)-debtR(1:tmax-1)*debtt(1:tmax-1)*      &
                        ts(1:tmax-1,16)
      ts(tmax,17)     = ts(tmax-1,17)
      gspendt         = gspendt*ts(:,16)
    endif
    tranj             = sum(popt*bxit,dim=2)/sum(popt,dim=2)
    ! tottrns          = sum(popt*pens,dim=2)/sum(popt,dim=2)+sum(popt*health,dim=2)/sum(popt,dim=2)

     tottrns          = sum(popt*pens,dim=2)/sum(popt,dim=2)+sum(popt*medrate*health,dim=2)/sum(popt,dim=2) &
                      + sum(popt*ltcrate*ltcare,dim=2)/sum(popt,dim=2)  ! Rev 2020/09/27

    ! Below computes all heath transfers which are part of consumption but
    ! do not generate revenues
    hlthtrs          = sum(popt*medrate*health,dim=2)/sum(popt,dim=2)  &
                     + sum(popt*ltcrate*ltcare,dim=2)/sum(popt,dim=2)  ! Rev 2020/09/27
    !
    ! Residual equations
    !
    resx(:,1)         = rt1t-thett1*alpha*ts(:,8)/ts(:,1)
    resx(:,2)         = wt-thetl1*alpha*ydivl*(ts(:,5)+ts(:,6))/ts(:,5)
    resx(:,3)         = gdpt-ts(:,8)+ts(:,11)+ts(:,12)

    ! relative to original program we take health transfers from consumption in terms of the tax base
    resx(:,4)         = tauct-(ctrant+tranj+tottrns(:)+gspendt-(1.+gz)*tkap-ts(:,18)+   &
                        ts(:,19)-taudt*ts(:,15)-taupct*ts(:,13)-tauput*ts(:,14)-ts(:,17)) &
                        /(ts(:,7)-hlthtrs(:))
    resx(2:tmax,5)    = debtt(2:tmax)/(debtt(2:tmax)+ts(2:tmax,1)+ts(2:tmax,2)+ts(2:tmax,3)+ts(2:tmax,4)) &
                       -debtH(2:tmax)
    ! from resx(:,5) ! 2020/03/15 by iiboshi

  endif

    if(indic.eq.1) then
    x(:,5)=debtt(:)/(debtt(:)+ts(:,1)+ts(:,2)+ts(:,3)+ts(:,4))
    end if

end subroutine trans

subroutine foc(input,wsched,rsched,mtx,ntx,ntr, output)
  !
  ! Two equations in two unknowns (c,l)
  !
  !   -- if taxes are proportional:
  !
  !      gam*p*c/(1-l) = lwedge * w
  !      p*c = nlinc+ lwedge * w*l
  !
  !   -- if taxes are proportional:
  !
  !      gam*p*c/(1-l) = w-dtax(w*l)/dl
  !      p*c = nlinc+w*l-tax(w*l)
  !
  ! where
  !    u(c,1-l) = log(c)+gam*log(1-l)
  !    nlinc    = nonlabor income
  !    l in [0,1] is checked
  !

  implicit none
  integer                             :: i,maxit
  integer, intent(in)                 :: mtx,ntx,ntr
  real(8), dimension(12),intent(in)       :: input
  real(8), dimension(3*mtx),intent(in)   :: wsched,rsched
  real(8), dimension(3),intent(out)      :: output
  real(8)                                :: gam,p,lwedge,w,otran,capinc,newa,   &
                                         nlinc,atw,del,res,dres,u,c,l,l0,lp, &
                                         taxinc,tax,mrate,psi, hlth, zeta

  gam        = input(11)
  p          = input(2)
  lwedge     = input(3)   ! Turn this on for proportional taxes
  w          = input(4)
  otran      = input(5)
  capinc     = input(6)
  newa       = input(7)
  hlth       = input(8)   !health expenditure(added in Japan code)
  zeta       = input(10)
  nlinc      = otran+capinc+newa

  if (w<=0.) then
    !
    ! Retiree
    !
    call comptax(capinc,rsched,mtx,ntr, psi,mrate)
    l        = 0.
    ! c          = hlth+(nlinc+psi-hlth)/p !消費にhealth支出を追加（以下も同様に）
      c          = (nlinc+psi-hlth)/p   ! revise 2020/09/22

  else
    !
    ! Worker
    !
    if (lwedge<.99) then
      !
      ! Facing proportional taxes (tax tableで計算する際には使っていない)
      !
      if ( zeta == 1 ) then
         atw    = lwedge*w
         l      = (atw-gam*(nlinc-hlth))/(atw+atw*gam)

     else
         atw    = lwedge*w
         l      = (.6**zeta*atw+.4*zeta*(.6)**(zeta-1.)-gam*(nlinc-hlth))/             &
               (zeta*(.6)**(zeta-1.)*atw+gam*atw)
        do i=1,10
           res  = atw*(1.-l)**zeta-gam*(atw*l+nlinc-hlth)
           dres = -zeta*atw*(1.-l)**(zeta-1.)-gam*atw
           l    = l-res/dres
        enddo
     endif

         ! c      = hlth+(nlinc+atw*l-hlth)/p
          c      = (nlinc+atw*l-hlth)/p
      if (l>1.) then
        l    = 1.
        ! c    = hlth+(nlinc+lwedge*w-hlth)/p
        c    = (nlinc+lwedge*w-hlth)/p
      endif
      if (l<0.) then
        l    = 0.
        !c    = hlth+(nlinc-hlth)/p
        c    = (nlinc-hlth)/p
      endif
    else
    !
        ! Facing net tax function wtax(.)：ここ以下で計算
        !
      maxit  = 5
      l0     = .3
!      maxit  = 8
!      l0     = .1


      del    = .001
      do i=1,maxit
        taxinc = w*l0
        call comptax(taxinc,wsched,mtx,ntx, tax,mrate)
        res    = gam*(nlinc+taxinc-tax)-w*(1.-mrate)*(1.-l0)**zeta
        lp     = l0+del
        taxinc = w*lp
        call comptax(taxinc,wsched,mtx,ntx, tax,mrate)
        dres   = (gam*(nlinc+taxinc-tax)-w*(1-mrate)*(1.-lp)**zeta-res)/.001
        l0     = l0-res/dres
      enddo
      l        = l0
      taxinc   = w*l
      call comptax(taxinc,wsched,mtx,ntx, tax,mrate)
      !  c        = hlth+(nlinc+w*l-tax-hlth)/p
          c          = (nlinc+w*l-tax-hlth)/p   ! revise 2020/09/22
      if (l>1.) then
        l      = 1.
        call comptax(taxinc,wsched,mtx,ntx, tax,mrate)
        ! c      = hlth+(nlinc+w*l-tax-hlth)/p
          c      = (nlinc+w*l-tax-hlth)/p ! rev 2020/09/22
      endif
      if (l<0.) then
        l      = 0.
         !  c      = hlth+(nlinc-hlth)/p
          c      = (nlinc-hlth)/p ! rev 2020/09/22
      endif
    endif
  endif
  if (c<=0.) then
    u        = -1.0e+9
    c        = 0.
  else
     if ( zeta == 1 ) then
       u        = log(c) + gam * log(1.-l)
!      u        = log(c - hlth) + gam * log(1.-l)  ! adj 2020/06/07 効用関数を消費から健康支出を引いた値を最大化する形に修正　→ 2020/08/30 上で消費からhlthを引いているのでこれは使わない。
     else
        if (l>=1.) then
           u     = log(c)
        else
           u     = log(c)+gam*(1.-l)**(1.-zeta)/(1.-zeta)
        endif
     endif
  endif

  output     = (/ u,c,l /)

end subroutine foc


subroutine comptax(taxinc,taxsched,mtx,ltx, tax,mrate)
  !
  ! Compute taxes paid and the marginal rate given
  !   (1) taxable income
  !   (2) the tax schedule: (x,y) pairs, where x is the rate for income over y
  !
  implicit none
  integer                            :: i,j,k,l
  integer, intent(in)                :: mtx,ltx
  real(8), dimension(3*mtx),intent(in)  :: taxsched
  real(8), intent(in)                   :: taxinc
  real(8), intent(out)                  :: tax,mrate

  tax    = 0.
  mrate  = 0.
  do i=1,ltx
    j       = 3*(i-1)+1
    k       = 3*(i-1)+2
    l       = 3*(i-1)+3
    if (taxinc>=taxsched(j)) then
      tax    = taxsched(k)+taxsched(l)*taxinc
      mrate  = taxsched(l)
    endif
  enddo


end subroutine comptax
