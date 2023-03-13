program bg_jp_ex
!
! 2020.04.12 balanced groth： 定常状態を計算して、そのまま拡大・縮小して伸ばす過程
!
! OLG MODEL WITH MULTIPLE PRODUCTIVITY (MP) LEVELS
!
!    Dynamic program of individuals (j=type,a=asset)
!
!       V(j,a) = max  u(c,1-l) + beta ps(j,j')* V(j',a')
!               c,l,a'
!
!                s.t. (1+gz) a' = (1+i) a + w eps(j,k) l -c
!                                   -wtax(w*eps(j,k)*l)+ rtra(i*a)
!                (*labor productivity for household type h)
!
!       where growing variables are detrended relative to technology
!       trend (1+gz)^t.  Types differ by age and productivity.
!       If m=number of ages and n=number of productivities, then
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
!                :
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
!       alb        = lower bound on asset holdings
!       aub        = upper bound on asset holdings
!       alpha      = sector 1 share
!       beta       = discount factor
!       delt1      = depreciation rate on tangible capital, sector 1
!       deli1      = depreciation rate on intangible capital, sector 1
!       delt2      = depreciation rate on tangible capital, sector 2
!       deli2      = depreciation rate on intangible capital, sector 2
!       gam        = parameter in utility
!       gn         = growth rate of population
!       gz         = growth rate of technology
!       lam        = degree of annuitization
!       ctran      = common transfers to all types
!       taud       = common tax rate on distributions
!       taul       = common tax rate on labor
!       taupc      = tax rate on corporate profits
!       taupu      = tax rate on unincorporated profits (distributed in full)
!       tfp        = TFP parameter = Z^{1-theta}
!       thett1     = tangible capital share, sector 1
!       theti1     = intangible capital share, sector 1
!       thett2     = tangible capital share, sector 2
!       theti2     = intangible capital share, sector 2
!       debt       = government borrowing (/GDP)
!       gspend     = government spending(/GDP))
!       idebt      = 0 if debt in levels, 1 if share of GDP
!       ispend     = 0 if gspend in levels, 1 if share of GDP
!       irate      = initial guess for interest rate
!       tauc       = initial guess for the tax rate on consumption
!
!       Net tax and transfer functions for workers and retirees
!         wtax     = (y,a,b) triplets, tax=a+b*y paid on income in [y,y(+1)]
!         rtra     = (y,a,b) triplets, tra=a+b*y is transfered for [y,y(+1)]
!
!       Number of ages * Number of productivities x 4:
!         eps      = labor productivity for household type h
!         ps       = probability of survival from age j to j'
!         xi       = type-dependent transfers
!         mu       = distribution of household types
!
!       nbeq x 3
!         beq      = nonaccidental bequest
!         jbeq     = index for bequester
!         jinh     = index for inheritor
!
!       NOTE:        edit file to change maximum age (nj), the number of
!                    productivity levels (ne), the dimension of the tax/
!                    transfer function (mtx), or the number of grid points
!                    in the asset grid (n).
!
!    Output file is loaded into Matlab -- run plotbgmp to view results.
!
!    The algorithm used for solving the fixed point is Gauss-Seidel.
!

! Ellen McGrattan, 2-26-14
! Revised, 1-15-15

  use mod_agg_jp_ex
  use mod_func_jp

  implicit none
  integer, parameter         :: ne=4,nj=101,nh=ne*nj,n=1500,mtx=14,npar=26,  &
                                 tmax=240,nx=3,mbeq=4, &
                                 gam_flag= 2 ! gam_flag = 1  -> constant,
                                              ! gam_flag = 2  -> depend on agents
                                              ! gam_flag = 3  -> depend on age
                                 ! nxは3に変更（元は2）
                                 ! nはグリッド数（手軽に400程度は必要、もとは1000）
  integer                     :: h,i,j,l,t,idone,it,maxit,ieven,iprod,nbeq, &
                                 ntx,ntr
  integer, dimension(nh)     :: hpi
  integer, dimension(nh,n)   :: opt
  integer, dimension(mbeq)   :: jbeq,jinh
  integer, dimension(2)      :: ipar
  real, dimension(npar)      :: par
  real, dimension(nh)        :: cj,lj,aj,apj,ltaxj,wpsij,rpsij,mtrj
  real, dimension(nh)        :: eps,ps,xi,mu,                               &
                                 pens, health, dislab,medrate, &    ! add 2020/05/29, dislab追加 2020/06/07
                                 ltcare, ltcrate, health_in,  ltcare_in    ! medrate 追加 2020/07/19 lng_care, crate 2020/08/18
  real, dimension(mtx,3)     :: wtax,rtra
  real, dimension(nx)        :: x,x1,res
  real, dimension(mbeq)      :: beq
  real, dimension(ne,4)      :: pscale
  real                        :: kapt1,kapt2,kapi1,kapi2,lab1,lab2,c,y,ltax,  &
                                 wpsi,rpsi,crit,sum1,xt1,xt2,xi1,xi2,          &
                                 debtR,                                   &
                                 nrm, asset, re_npar   ! add 2020/07/19,  add 2020/08/18
   real, dimension(ne,3)        ::  prop_agent

!  external agg   !【疑問　これなんですか？？】

!  open(unit=5, file='bgt_JP.txt')
!  open(unit=5, file='bgt_JP_ex_gam.txt')
!  open(unit=5, file='bgt_JP_ex_gam_eff.txt') !労働生産性を賃構2019で入れ替えた

!  open(unit=5, file='bgt_JP_ex_gam_eff_mu.txt')   !労働生産性を賃構2019＋資産クラスの構成比変更＋人口分布2015
!  open(unit=5, file='bgt_JP_ex_gam_eff_mu2115.txt') !労働生産性を賃構2019＋資産クラスの構成比変更＋人口分布2115
!  open(unit=5, file='bgt_JP_ex_gam_eff_munon.txt') !労働生産性を賃構2019＋資産クラスの構成比変更＋人口分布non-aging

!  open(unit=5, file='bgt_JP_ex_type_ltc.txt') ! 男女正規非正規＋人口分布2015
!  open(unit=5, file='bgt_JP_ex_type_2115.txt') ! 男女正規非正規＋人口分布2015
!  open(unit=5, file='bgt_JP_ex_type_nonaging.txt') ! 男女正規非正規＋人口分布non-aging

!!!!!2020/08/29 data adj Ozaki
!  open(unit=5, file='bgt_JP_ex_type_ltc_new2.txt') !医療・介護分離、フリッシュ弾力性など設定
!  open(unit=5, file='bgt_JP_ex_type_ltc_new2_tfp1.txt') !TFPパラメターを1に設定
!  open(unit=5, file='bgt_JP_ex_type_ltc_new_muprop.txt') !タイプ人口比例
  open(unit=5, file='bgt_JP_ex_type_ltc_new_muprop_tfp1.txt') !タイプ人口比例+TFPパラメターを1に設定


  open(unit=7, file='bgt_JP.nxt')         !inputファイルとoutputファイルにpensとhealthを足す。
  open(unit=8, file='bgt_JP.dat')
  open(unit=9, file='trantf_JP.inp')
  open(unit=10, file='tran_init_8.txt')


   debtR = 1.01_8  ! add 2020/05/29 【これどういう指定？】

  !
  ! Read in parameters
  !
  do i=1,npar
    read(5,*) par(i)
  enddo

  read(5,*) ipar(1)
  read(5,*) ipar(2)

  !
  ! Read in guess for the interest rate and the tax rate on consumption
  !
  read(5,*) x(1)
  read(5,*) x(2)
  read(5,*) x(3)   ! add 2020/05/30
  x1       = x

  !
  ! Check for consistency on number of productivity levels
  !
  read(5,*) i
  if (i/=ne) then
    write(*,*) 'Update bgtcmp.f90 with parameter ne=',i
    stop
  endif

  !
  ! Read in zeta and replacement rate tables
  !
  do i=1,ne
    read(*,*)  (prop_agent(i,j), j=1,3)
  enddo

  !
  ! Read in the tax and transfer tables
  !
  read(5,*) ntx
  if (ntx>mtx) then
    write(*,*) 'Update taxfunc.f90 with parameter mtx >',ntx
    stop
  endif
  do i=1,ntx
    read(5,*) (wtax(i,j),j=1,3)
  enddo
  read(5,*) ntr
  if (ntr>mtx) then
    write(*,*) 'Update taxfunc.f90 with parameter mtx >',ntr
    stop
  endif
  do i=1,ntr
    read(5,*) (rtra(i,j),j=1,3)
  enddo

  !
  ! Read in labor productivities, survival probabilities, and transfers
  !   if iprod = 0, read in parameters for all productivity levels
  !   if iprod = 1, read in for one productivity level and scale for others
  !

  read(5,*) iprod
  if (iprod==0) then
    do h=1,nh
!      read(5,*) eps(h),ps(h),xi(h),mu(h)
      read(5,*) eps(h),ps(h),xi(h),mu(h), pens(h), health_in(h), dislab(h) ,medrate(h), ltcare_in(h), ltcrate(h)
                  !2020/05/20  add pens, health, dislab追加 2020/06/07
                  !2020/07/19  add medrate  !2020/08/18  ltcare(h), ltcrate(h)
    enddo
  else
    do l=1,ne
      read(5,*) pscale(l,:)
    enddo
    do j=1,nj
      read(5,*) eps(j),ps(j),xi(j),mu(j)
    enddo
    do l=ne,1,-1
      do j=1,nj
        i      = (l-1)*nj+j
        eps(i) = pscale(l,1)*eps(j)
        ps(i)  = pscale(l,2)*ps(j)
        xi(i)  = pscale(l,3)*xi(j)
        mu(i)  = pscale(l,4)*mu(j)
      enddo
    enddo
  endif

  !
  ! Read in parameters for bequests
  !
  read(5,*) nbeq
  if (nbeq>0) then
    do l=1,nbeq
      read(5,*) beq(l),jbeq(l),jinh(l)
    enddo
  endif
  !
  ! Compute steady state
  !
  maxit   = 10 !収束回数：実際の計算では100くらいは必要。!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  crit    = 1.e-8
  idone   = 0
  it      = 0
  opt     = -1

  !
  ! 1. Hold tauc fixed and compute irate
  !

  do while ((idone==0).and.(it<maxit))
    call agg(x,nx,nh,ne,nj,n,eps,ps,xi,mu,par,npar,ipar,wtax,rtra,mtx,ntx,   &
             ntr,beq,jbeq,jinh,nbeq,mbeq,opt,  res,hpi,cj,lj,aj,apj,ltaxj,   &
             wpsij,rpsij,mtrj,kapt1,kapt2,kapi1,kapi2,lab1,lab2,c,y,ltax,    &
             wpsi,rpsi,pens, health, debtR, dislab, medrate, asset, &
             ltcare, ltcrate , prop_agent, gam_flag, health_in,  ltcare_in ) ! 2020/05/29 add pens, health, debtR, add dislab 2020/06/07
    x1(1:2) = .975*x(1:2)+.025*(x(1:2)-res(1:2))      !残差、technical Appsの(4.1)式に当たる？
                                            !!新しい候補の値に更新

    write(*,*) 'it_1 = ', it  ! 2020/03/15 by iiboshi
    write(*,'(10(1x,e15.9))') x1,maxval(abs(x1-x)/(1+abs(x)))
!    if (abs(x1(1)-x(1))/(1+abs(x(1)))<crit) idone=1
    nrm = (abs(x1(1)-x(1))/(1+abs(x(1))))
    if (nrm<crit) idone=1
    it    = it+1
    x(1:2)  = x1(1:2)
!    if (isnan(x1(1)) .or. isnan(x1(2))) then
!       x(2)  = x1(2)
!    else
!     x(1)  = x1(1)
 !   endif

     write(*,*)  'asset=',asset,'c=',c, 'y=',y,  'k1=', kapt1, 'k2=',kapt2, 'ki1=',kapi1,'ki2=',kapi2, &
                 'l1=',lab1, 'l2=',lab2, 'x1=',x(1), 'x2=',x(2), 'x3=',x(3)   ! add 2020/07/21
  enddo


  !
  ! 2. Take solution for irate and compute tauc
  !
  maxit   = 10 !収束回数：実際の計算では100くらいは必要。!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  idone   = 0
  it      = 0
  opt     = -1
  do while ((idone==0).and.(it<maxit))
    call agg(x,nx,nh,ne,nj,n,eps,ps,xi,mu,par,npar,ipar,wtax,rtra,mtx,ntx,   &  !【疑問　ここでcallしてるのが「agg」でOKなのはなぜ？
             ntr,beq,jbeq,jinh,nbeq,mbeq,opt,  res,hpi,cj,lj,aj,apj,ltaxj,   &
             wpsij,rpsij,mtrj,kapt1,kapt2,kapi1,kapi2,lab1,lab2,c,y,ltax,    &
             wpsi,rpsi, pens, health, debtR, dislab, medrate, asset, &
             ltcare, ltcrate , prop_agent, gam_flag, health_in,  ltcare_in  )  ! 2020/05/29 add pens, health, debtR, add dislab 2020/06/07
    x1(3) = x(3)-res(3) !!新しい候補の値に更新

    write(*,*) 'it_2 = ', it  ! 2020/03/15 by iiboshi
    write(*,'(10(1x,e15.9))') x1,maxval(abs(x1-x)/(1+abs(x)))
    if (maxval(abs(x1-x)/(1+abs(x)))<crit) idone=1
    it    = it+1
    x(3)  = x1(3)

     write(*,*)  'asset=',asset,'c=',c, 'y=',y,  'k1=', kapt1, 'k2=',kapt2, 'ki1=',kapi1,'ki2=',kapi2, &
                 'l1=',lab1, 'l2=',lab2, 'x1=',x(1), 'x2=',x(2), 'x3=',x(3)   ! add 2020/07/21
  enddo

  if (maxit==0) then
    write(*,*) 'Results based on user-supplied inputs...'
    call agg(x1,nx,nh,ne,nj,n,eps,ps,xi,mu,par,npar,ipar,wtax,rtra,mtx,ntx,  &
             ntr,beq,jbeq,jinh,nbeq,mbeq,opt,  res,hpi,cj,lj,aj,apj,ltaxj,   &
             wpsij,rpsij,mtrj,kapt1,kapt2,kapi1,kapi2,lab1,lab2,c,y,ltax,    &
             wpsi,rpsi, pens, health, debtR, dislab, medrate, asset , &
             ltcare, ltcrate  , prop_agent, gam_flag , health_in,  ltcare_in )  ! 2020/05/29 add pens, health, debtR, add dislab 2020/06/07
  endif
  xt1   = ((1.+par(10))*(1.+par(11))-1.+par(5))*kapt1
  xt2   = ((1.+par(10))*(1.+par(11))-1.+par(7))*kapt2
  xi1   = ((1.+par(10))*(1.+par(11))-1.+par(6))*kapi1
  xi2   = ((1.+par(10))*(1.+par(11))-1.+par(8))*kapi2


  !
  ! Write out new input file for steady state calculations
  !

   write(*,*) 'about to write output'
   write(*,*) 'to write new input file'


  write(7,'(1x,f14.8,5x,a3)') par(1),'alb'
  write(7,'(1x,f14.8,5x,a3)') par(2),'aub'
  write(7,'(1x,f14.8,5x,a4)') par(3),'alpha'
  write(7,'(1x,f14.8,5x,a4)') par(4),'beta'
  write(7,'(1x,f14.8,5x,a5)') par(5),'delt1'
  write(7,'(1x,f14.8,5x,a5)') par(6),'deli1'
  write(7,'(1x,f14.8,5x,a5)') par(7),'delt2'
  write(7,'(1x,f14.8,5x,a5)') par(8),'deli2'
  write(7,'(1x,f14.8,5x,a3)') par(9),'gam'
  write(7,'(1x,f14.8,5x,a2)') par(10),'gn'
  write(7,'(1x,f14.8,5x,a2)') par(11),'gz'
  write(7,'(1x,f14.8,5x,a3)') par(12),'lam'
  write(7,'(1x,f14.8,5x,a5)') par(13),'ctran'
  write(7,'(1x,f14.8,5x,a4)') par(14),'taud'
  write(7,'(1x,f14.8,5x,a4)') par(15),'taul'
  write(7,'(1x,f14.8,5x,a5)') par(16),'taupc'
  write(7,'(1x,f14.8,5x,a5)') par(17),'taupu'
  write(7,'(1x,f14.8,5x,a3)') par(18),'tfp'
  write(7,'(1x,f14.8,5x,a6)') par(19),'thett1'
  write(7,'(1x,f14.8,5x,a6)') par(20),'theti1'
  write(7,'(1x,f14.8,5x,a6)') par(21),'thett2'
  write(7,'(1x,f14.8,5x,a6)') par(22),'theti2'
  write(7,'(1x,f14.8,5x,a4)') par(23),'debt'
  write(7,'(1x,f14.8,5x,a6)') par(24),'gspend'
  write(7,'(1x,i14,5x,a5)')   ipar(1),'idebt'
  write(7,'(1x,i14,5x,a6)')   ipar(2),'ispend'
  write(7,'(1x,e14.8,5x,a5)') x1(1),'irate'
  write(7,'(1x,e14.8,5x,a4)') x1(2),'tauc'
  write(7,*)
  write(7,'(1x,i14,5x,a29)') ne,'Number of productivity levels'
  write(7,*)
  write(7,'(1x,i14,5x,a24)') ntx,'Tax function for workers'
  do i=1,ntx
    write(7,'(3(1x,e14.8))') (wtax(i,j),j=1,3)
  enddo
  write(7,'(1x,i14,5x,a30)') ntr,'Transfer function for retirees'
  do i=1,ntr
    write(7,'(3(1x,e14.8))') (rtra(i,j),j=1,3)
  enddo
  write(7,*)
  write(7,'(1x,i14,5x,a25)')  iprod,'Scale eps,ps,xi,mu below?'
  if (iprod==0) then
    write(7,'(4(1x,e14.6),5x,a2)') eps(1),ps(1),xi(1),mu(1),'No'
    do h=2,nh
      write(7,'(4(1x,e14.6))') eps(h),ps(h),xi(h),mu(h)
    enddo
  else
    write(7,'(4(1x,e14.6),5x,a3)') pscale(1,:),'Yes'
    do l=2,ne
      write(7,'(4(1x,e14.6))') pscale(l,:)
    enddo
    write(7,*)
    do j=1,nj
      write(7,'(4(1x,e14.6))') eps(j)/pscale(1,1),ps(j)/pscale(1,2),         &
                               xi(j)/pscale(1,3),mu(j)/pscale(1,4)
    enddo
  endif
  write(7,'(1x,i14,5x,a18)') nbeq,'Number of bequests'
  do l=1,nbeq
    write(7,'(1x,e14.8,2(1x,i3))') beq(l),jbeq(l),jinh(l)
  enddo

  !
  ! Write out new input file for transition calculations
  ! (注)ここで移行過程用のデータを計算（移行過程は、定常状態からのズレをプロット。ここのデータがないと移行過程を再現できない）
  !
  write(9,'(1x,f14.8,5x,a3)') par(1),'alb'
  write(9,'(1x,f14.8,5x,a3)') par(2),'aub'
  write(9,'(1x,f14.8,5x,a5)') par(3),'alpha'
  write(9,'(1x,f14.8,5x,a4)') par(4),'beta'
  write(9,'(1x,f14.8,5x,a5)') par(5),'delt1'
  write(9,'(1x,f14.8,5x,a5)') par(6),'deli1'
  write(9,'(1x,f14.8,5x,a5)') par(7),'delt2'
  write(9,'(1x,f14.8,5x,a5)') par(8),'deli2'
  write(9,'(1x,f14.8,5x,a3)') par(9),'gam'
  write(9,'(1x,f14.8,5x,a2)') par(11),'gz'
  write(9,'(1x,f14.8,5x,a3)') par(12),'lam'
  write(9,'(1x,f14.8,5x,a6)') par(19),'thett1'
  write(9,'(1x,f14.8,5x,a6)') par(20),'theti1'
  write(9,'(1x,f14.8,5x,a6)') par(21),'thett2'
  write(9,'(1x,f14.8,5x,a6)') par(22),'theti2'
  if (ipar(1)==0) then
    write(9,'(1x,f14.8,5x,a5)') par(23),'debt0'
  else
    write(9,'(1x,f14.8,5x,a5)') par(23)*(y-xi1-xi2),'debt0'
  endif
  write(9,'(1x,i14,5x,a5)')   ipar(1),'idebt'
  write(9,'(1x,i14,5x,a6)')   ipar(2),'ispend'
  write(9,*)

  do i=1,ne
     write(9,'(3(2x,e10.4))') prop_agent(i,1),  prop_agent(i,2), prop_agent(i,3)
  enddo

  write(9,*)
  do t=1,tmax
    write(9,'(13(1x,e14.6))') par(23),par(24),par(13:18),par(10), prop_agent(1,2), prop_agent(2,2), &
                              prop_agent(3,2), prop_agent(4,2)
  enddo


  write(9,*)
  do h=1,nh
    write(9,'(1x,i14,1x,e14.8,1x,e14.8)') hpi(h),mu(h),eps(h)
  enddo
  write(9,*)
  write(9,'(1x,i14,5x,a23)')  ntx,'Tax functions over time'
  do t=1,tmax
    do i=1,ntx
      write(9,'(3(2x,e10.4))',advance='no') wtax(i,:)
    enddo
    write(9,*)
  enddo
  write(9,'(1x,i14,5x,a28)')  ntr,'Transfer functions over time'
  do t=1,tmax
    do i=1,ntr
      write(9,'(3(2x,e10.4))',advance='no') rtra(i,:)
    enddo
    write(9,*)
  enddo

  write(9,*)
  do h=1,nh
    do t=1,tmax
      ! write(9,'(2(1x,e14.6))') ps(h),0.
      write(9,'(7(1x,e14.6))') ps(h), xi(h), pens(h), health(h), medrate(h), ltcare(h), ltcrate(h)
      !add 2020/05/29 adj 2020/06/06
    enddo
  enddo
  write(9,*)
  do l=1,nbeq
    write(9,'(1x,e14.6,2(1x,i4))') beq(l),jbeq(l),jinh(l)
  enddo

  write(9,*)
  do t=1,tmax
!    write(9,'(4(1x,e14.6))') x1(1),(1.-par(19)-par(20))*par(3)*y/lab1,       & !2020/06/06
!                             y-xi1-xi2,x1(2)
    write(9,'(5(1x,e14.6))') x1(1),(1.-par(19)-par(20))*par(3)*y/lab1,       &
                             y-xi1-xi2,x1(2),                               &
                             par(23)*(y-xi1-xi2)/(kapt1+kapt2+kapi1+kapi2+par(23)*(y-xi1-xi2)) ! 2020/06/06
  enddo

  write(10,*)
  do t=1,tmax
    write(10,'(5(1x,e14.6))') x1(1),(1.-par(19)-par(20))*par(3)*y/lab1,       &
                             y-xi1-xi2,x1(2), par(23)/(kapt1+kapt2+kapi1+kapi2+par(23))
  enddo


  !
  ! Write out results for plotbgmp.m
  ! 描画用のMATLABコードに渡す部分　→ここにはpens, health, debtRは不要？→Matlab側で取り出す
  !
  re_npar = npar - 2  ! 2020/08/18

  do i=1,12
    write(8,'(f13.8)') par(i)
  enddo
  write(8,'(f13.8)') x1(2)
  do i=14,re_npar-2
    write(8,'(f13.8)') par(i)
  enddo
  if (ipar(1)==0) then
    write(8,'(f13.8)') par(re_npar-1)
  else
    write(8,'(f13.8)') par(re_npar-1)*(y-xi1-xi2)
  endif
  if (ipar(2)==0) then
    write(8,'(f13.8)') par(re_npar)
  else
    write(8,'(f13.8)') par(re_npar)*(y-xi1-xi2)
  endif
  write(8,'(f13.8)') kapt1
  write(8,'(f13.8)') kapt2
  write(8,'(f13.8)') kapi1
  write(8,'(f13.8)') kapi2
  write(8,'(f13.8)') lab1
  write(8,'(f13.8)') lab2
  write(8,'(f13.8)') c
  write(8,'(f13.8)') xt1
  write(8,'(f13.8)') xt2
  write(8,'(f13.8)') xi1
  write(8,'(f13.8)') xi2
  write(8,'(f13.8)') y
  write(8,'(f13.8)') ltax
  write(8,'(f13.8)') wpsi
  write(8,'(f13.8)') rpsi
  write(8,'(f13.8)') x1(1)
  write(8,'(f13.8)') par(13)
  write(8,'(i13)')   nj
  write(8,'(i13)')   ne
  write(8,'(i13)')   n
  do h=1,nh
    write(8,'(f13.8)') eps(h)
  enddo
  do h=1,nh
    write(8,'(f13.8)') ps(h)
  enddo
  do h=1,nh
    write(8,'(f13.8)') xi(h)
  enddo
  do h=1,nh
    write(8,'(f13.8)') mu(h)
  enddo
  do h=1,nh
    write(8,'(i13)') hpi(h)
  enddo
  do h=1,nh
    write(8,'(f13.8)') cj(h)
  enddo
  do h=1,nh
    write(8,'(f13.8)') lj(h)
  enddo
  do h=1,nh
    write(8,'(f13.8)') aj(h)
  enddo
  do h=1,nh
    write(8,'(f13.8)') apj(h)
  enddo
  do h=1,nh
    write(8,'(f13.8)') ltaxj(h)
  enddo
  do h=1,nh
    write(8,'(f13.8)') wpsij(h)+ medrate(h)*health(h) + ltcrate(h)*ltcare(h) ! revise 2020/08/18
  enddo
  do h=1,nh
    write(8,'(f13.8)') rpsij(h)+pens(h)  ! revise 2020/07/19
  enddo
  do h=1,nh
    write(8,'(f13.8)') mtrj(h)
  enddo
  do i=1,n
    do h=1,nh
      write(8,'(i13)') opt(h,i)
    enddo
  enddo
  write(8,'(i13)') ntx
  do j=1,3
    do i=1,ntx
      write(8,'(f13.8)') wtax(i,j)
    enddo
  enddo
  write(8,'(i13)') ntr
  do j=1,3
    do i=1,ntr
      write(8,'(f13.8)') rtra(i,j)
    enddo
  enddo
  write(8,'(i13)') nbeq
  do l=1,nbeq
    write(8,'(f13.8)') beq(l)
  enddo
  do l=1,nbeq
    write(8,'(i13)') jbeq(l)
  enddo
  do l=1,nbeq
    write(8,'(i13)') jinh(l)
  enddo

  do j=1,3
    do i=1,ne
      write(8,'(f13.8)') prop_agent(i,j)
    enddo
  enddo

  do h=1,nh                        ! add 2020/08/29
    write(8,'(f13.8)') ltcare(h)     ! add 2020/08/29
  enddo                            ! add 2020/08/29
  do h=1,nh                        ! add 2020/05/30
    write(8,'(f13.8)') pens(h)     ! add 2020/05/30
  enddo                            ! add 2020/05/30
  do h=1,nh                        ! add 2020/05/30
    write(8,'(f13.8)') health(h)   ! add 2020/05/30
  enddo


end program bg_jp_ex

