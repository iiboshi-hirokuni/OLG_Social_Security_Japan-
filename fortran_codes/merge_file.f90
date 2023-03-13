program trantf

implicit none

  integer, parameter             :: ne=4,nj=101,nh=ne*nj,mtx=14,npar=16,     &
                                    na=500,tmax=240,hmax=nh*tmax,           &
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
  real(8), dimension(tmax,nh)       :: pst,xit,pens,health,income1,income2, &
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
   

   ! Input File
   !open(unit=5,  file='trantfWPHGbaseline.inp')
   ! open(unit=5,  file='trantfWPHGhighertfp.inp')
   open(unit=51,  file='./file_in/trantf_XX1_para.txt')
   open(unit=52,  file='./file_in/trantf_XX2_pop.txt')
   open(unit=53,  file='./file_in/trantf_XX3_wtax.txt')
   open(unit=54,  file='./file_in/trantf_XX4_transf.txt')  
   open(unit=55,  file='./file_in/trantf_XX5_exog.txt')
   ! open(unit=56,  file='./file_in/trantf_XX6_pst.txt')
   open(unit=56,  file='./file_in/trantf_XX6_pst_new.txt')
   open(unit=57,  file='./file_in/trantf_XX7_beq.txt')
   open(unit=58,  file='./file_in/trantf_XX8_ini_guess.txt')

    ! Output File
   open(unit=7,  file='./file_in/trantf_XX.inp')
 

  !============================================================================  
  !============================================================================  

    !
    ! Read in parameters
    !
    do i=1,npar
      read(51,*) par(i)
    enddo
    read(51,*) ipar(1)
    read(51,*) ipar(2)
    
   !
   ! Read in zeta and replacement rate tables
   !
   do i=1,ne
      read(51,*)  (prop_agent(i,j), j=1,3)
   enddo
   
   !
    ! Read in TFP and fiscal exogenous variables
    !
    do t=1,tmax
      read(55,*) exog(t,:)
    enddo



    !
    ! Read in initial asset distributions, populations, productivities
    !
    do h=1,nh
      read(52,*) hpii(h),popi(h),eps(h)
    enddo

    !
    ! Read in the net tax tables
    !
    read(53,*) ntx
    if (ntx>mtx) then
      write(*,*) 'Update trantf.f90 with parameter mtx>=',ntx
      stop
    endif
    do t=1,tmax
      read(53,*) (wtaxt(t,i),i=1,3*ntx)
    enddo

    !
    ! Read in the transfer tables
    !
    read(54,*) ntr
    if (ntr>mtx) then
      write(*,*) 'Update trantf.f90 with parameter mtx>=',ntr
      stop
    endif
    do t=1,tmax
      read(54,*) (rtrat(t,i),i=1,3*ntr)
    enddo

    !
    ! Read in survival probabilities and age and time dependent lump sum transfers (including pensions and health)
    !
    write(*,*) 'read transfers:'
    do l=1,ne
      do j=1,nj
        h  = (l-1)*nj+j
        do t=1,tmax
         !  read(56,*) pst(t,h),xit(t,h),pens(t,h),health(t,h)
            read(56,*) pst(t,h),xit(t,h),pens(t,h),health(t,h),medrate(t,h), ltcare(t,h), ltcrate(t,h)
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
          read(57,*) beq(t,j),jbeq(t,j),jinh(t,j)
        enddo
      enddo
    endif
    
    
    !  File No 8 
    do t=1,tmax
      read(58,*) x(t,:)
    enddo   
          
  !============================================================================  
  !============================================================================  
 
    
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
    
    do i=1,ne
     write(7,'(3(2x,e10.4))') prop_agent(i,1),  prop_agent(i,2), prop_agent(i,3)
   enddo
  
   do t=1,tmax
      write(7,'(10(1x,e14.6))') exog(t,:)
    enddo
    
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
    
    
    write(7,*)
    
    do h=1,nh
      do t=1,tmax
        ! write(7,'(4(1x,e14.6))') pst(t,h),xit(t,h),pens(t,h),health(t,h)
        write(7,'(7(1x,e14.6))') pst(t,h), xit(t,h), pens(t,h), health(t,h), medrate(t,h), ltcare(t,h), ltcrate(t,h)
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

    
  
end program trantf

