program trantf

implicit none

  integer, parameter             :: ne=4,nj=101,nh=ne*nj,mtx=14,npar=16,     &
                                    na=800,tmax=240,hmax=nh*tmax,           &
                                    omax=na*nh*40,tbeq=1,nbeq=4


integer                :: h,i,j,l,t,idone,it,maxit,icomp,ntx,ntr,iteracion,&
                                    indic
  integer, dimension(nh)         :: hpii
  integer, dimension(omax)       :: myo
  integer, dimension(2)          :: ipar
  integer, dimension(tbeq,nbeq)  :: jbeq,jinh
  real(8), dimension(npar)          :: par
  real(8), dimension(nh)            :: eps,popi
  real(8), dimension(tmax,3*mtx)    :: wtaxt,rtrat
  real(8), dimension(tmax,nh)       :: pst,xit,pens,health,income1,income2,  &
                                       medrate, ltcare, ltcrate
  real(8), dimension(tmax,13)        :: exog
  real(8), dimension(tmax,5)        :: x,x1,res
  real(8), dimension(hmax,4)        :: allv
  real(8), dimension(tmax,19)       :: ts
  real(8), dimension(tmax,ne)       :: ltaxj,wpsij,rpsij
  real(8), dimension(tmax)         :: agedep,debtH,debtR
  real(8), dimension(tbeq,nbeq)     :: beq
  real(8)                           :: crit,sum1,nrm
   real, dimension(ne,3)        ::  prop_agent !add 2020/08/18
  external trans

   ! Input File
!   open(unit=5,  file='./file/trantf_JP.inp')
!   open(unit=5,  file='./file/trantf_JP_debt15.inp')
   open(unit=5,  file='./file_in/trantf_JP_2sec_g1500_gs199_baseline.inp')
   
   !open(unit=5,  file='trantf_JP_ex3.inp')

   ! open(unit=5,  file='trantfWPHGbaseline.inp')
   ! open(unit=5,  file='trantfWPDEBT.inp')
   ! open(unit=5,  file='trantfWPfavorabledemog.inp')
   ! open(unit=5,  file='trantfWPPITfinhealth.inp')
   ! open(unit=5,  file='trantfWPHGhighertfp.inp')


   ! output File
   open(unit=71,  file='./file_in/trantf_XX1_para.txt')
   open(unit=72,  file='./file_in/trantf_XX2_pop.txt')
   open(unit=73,  file='./file_in/trantf_XX3_wtax.txt')
   open(unit=74,  file='./file_in/trantf_XX4_transf.txt')
   open(unit=75,  file='./file_in/trantf_XX5_exog.txt')
   open(unit=76,  file='./file_in/trantf_XX6_pst.txt')
   open(unit=77,  file='./file_in/trantf_XX7_beq.txt')
   open(unit=78,  file='./file_in/trantf_XX8_ini_guess.txt')
  ! open(unit=79,  file='trantf_XX.csv')

    !
    ! Read in parameters
    !
    do i=1,npar
      read(5,*) par(i)
    enddo
    read(5,*) ipar(1)
    read(5,*) ipar(2)

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
          !  read(5,*) pst(t,h),xit(t,h),pens(t,h),health(t,h)
           read(5,*) pst(t,h),xit(t,h),pens(t,h),health(t,h),medrate(t,h), ltcare(t,h), ltcrate(t,h)
        enddo
      enddo
      pst(tmax,l*nj) = 0.
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

    do t=1,tmax
      read(5,*) x(t,:)
    enddo




      write(*,*) 'about to write output'


    write(*,*) 'to write new input file'

    ! File No 1
    write(71,'(1x,f14.8,6x,a3)') par(1),'alb'
    write(71,'(1x,f14.8,6x,a3)') par(2),'aub'
    write(71,'(1x,f14.8,6x,a4)') par(3),'alpha'
    write(71,'(1x,f14.8,6x,a4)') par(4),'beta'
    write(71,'(1x,f14.8,6x,a5)') par(5),'delt1'
    write(71,'(1x,f14.8,6x,a5)') par(6),'deli1'
    write(71,'(1x,f14.8,6x,a5)') par(7),'delt2'
    write(71,'(1x,f14.8,6x,a5)') par(8),'deli2'
    write(71,'(1x,f14.8,6x,a3)') par(9),'gam'
    write(71,'(1x,f14.8,6x,a2)') par(10),'gz'
    write(71,'(1x,f14.8,6x,a3)') par(11),'lam'
    write(71,'(1x,f14.8,6x,a6)') par(12),'thett1'
    write(71,'(1x,f14.8,6x,a6)') par(13),'theti1'
    write(71,'(1x,f14.8,6x,a6)') par(14),'thett2'
    write(71,'(1x,f14.8,6x,a6)') par(15),'theti2'
    write(71,'(1x,f14.8,6x,a5)') par(16),'debt0'
    write(71,'(1x,i14,6x,a5)')   ipar(1),'idebt'
    write(71,'(1x,i14,6x,a6)')   ipar(2),'ispend'
    write(71,*)


  do i=1,ne
     write(71,'(3(2x,e10.4))') prop_agent(i,1),  prop_agent(i,2), prop_agent(i,3)
  enddo



    !  File No 2
    do h=1,nh
      write(72,'(1x,i14,1x,e14.6,1x,e14.6)') hpii(h),popi(h),eps(h)
    enddo
    write(72,*)

    !  File No 3
    write(73,'(1x,i14,6x,a19)')  ntx,'Number of tax rates'
    do t=1,tmax
      write(73,'(48(1x,e14.6))') (wtaxt(t,i),i=1,3*ntx)
    enddo
    write(73,*)

    !  File No 4
    write(74,'(1x,i14,6x,a19)')  ntr,'Number of transfers'
    do t=1,tmax
      write(74,'(48(1x,e14.6))') (rtrat(t,i),i=1,3*ntr)
    enddo
    write(74,*)

    !  File No 5
    do t=1,tmax
      write(75,'(10(1x,e14.6))') exog(t,:)
    enddo
    write(75,*)

    !  File No 6
    do h=1,nh
      do t=1,tmax
        write(76,'(7(1x,e14.6))') pst(t,h), xit(t,h), pens(t,h), health(t,h), medrate(t,h), ltcare(t,h), ltcrate(t,h)
       ! write(76,'(4(1x,e14.6))') pst(t,h),xit(t,h),pens(t,h),health(t,h)
      enddo
    enddo
    write(76,*)

    !  File No 7
    if (nbeq>0) then
      do t=1,tbeq
        do j=1,nbeq
          write(77,'(1x,e16.6,1x,i3,1x,i3)') beq(t,j),jbeq(t,j),jinh(t,j)
        enddo
      enddo
    endif
    write(77,*)

    do t=1,tmax
      write(78,'(5(1x,e14.6))') x(t,:)
    enddo



end program trantf

