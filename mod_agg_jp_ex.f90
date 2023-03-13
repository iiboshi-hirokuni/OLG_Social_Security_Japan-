module mod_agg_jp_ex

   implicit none

   contains

subroutine agg(x,nx,nh,ne,nj,n,eps,ps,xi,mu,par,npar,ipar,wtax,rtra,mtx,ntx, &
               ntr,beq,jbeq,jinh,nbeq,mbeq,opt, resx,hpi,cj,lj,aj,apj,ltaxj, &
               wpsij,rpsij,mtrj,kapt1,kapt2,kapi1,kapi2,lab1,lab2,c,y,ltax,  &
               wpsi,rpsi,pens, health, debtR, dislab,medrate, asset, &
               ltcare, ltcrate, prop_agent, gam_flag, health_in,  ltcare_in  )
               !年金pens、医療費health, debtR 追加 2020/05/17, add dislab 2020/06/07
               !介護費ltcare, 介護負担率ltrate, 医療負担率mdrate

  !
  ! Dynamic programs for the OLG model
  ! inputはdimension(8)／outputはdimension(7)
  !

  use mod_func_jp   ! add 2020/05/21

!  implicit none
  integer, intent(in)                      :: nh,ne,nj,n,npar,nx,nbeq,mbeq,  &
                                              mtx,ntx,ntr, gam_flag
  integer                                  :: h,i,j,k,k1,k2,kzero,l
  integer, dimension(1)                    :: kk
  integer, dimension(2),intent(in)         :: ipar
  integer, dimension(nh),intent(out)       :: hpi
  integer, dimension(nh,n),intent(inout)   :: opt
  integer, dimension(mbeq),intent(in)      :: jbeq,jinh
  real, dimension(nx),intent(in)           :: x
  real, dimension(nx),intent(out)          :: resx
  real, dimension(npar),intent(in)         :: par
  real, dimension(mtx,3),intent(in)        :: wtax,rtra
  real, dimension(mbeq),intent(in)         :: beq
  real, dimension(nh),intent(inout)      :: eps,ps,xi,mu, &
                                               pens,health, dislab, & ! add 2020/05/20, add dislab 2020/06/07
                                               medrate,  &   ! add 2020/07/19
                                               ltcare, ltcrate, health_in,  ltcare_in ! add 2020/08/18
  real, dimension(nh),intent(out)        :: cj,lj,aj,apj,ltaxj,wpsij,rpsij,mtrj
  real, intent(out)                       :: kapt1,kapt2,kapi1,kapi2,lab1,  &
                                              lab2,c,y,ltax,wpsi,rpsi,asset
  real, dimension(nh)                     :: bxi
  real, dimension(n)                      :: vnext,v,rhs
  real, dimension(12)                     :: input
  real, dimension(7)                      :: output
  real                                     :: alb,aub,alpha,alph1,beta,      &
                                              delt1,deli1,delt2,deli2,gn,    &
                                              gz,gam,tfp,tauc,taud,taul,     &
                                              taupc,taupu,thett1,theti1,     &
                                              thetl1,thett2,theti2,thetl2,   &
                                              w,sum1,ctran,tkap,a,at,        &
                                              irate,cdiv,cprof,uprof,        &
                                              dshare,gshare,gdp,rt1,rt2,ri1, &
                                              ri2,tranj,invt1,invt2,invi1,   &
                                              invi2,p1,p2,y1,y2,krat1,krat2, &
                                              kaprat,kdivlt1,kdivlt2,        &
                                              kdivli1,kdivli2,ydivl,ldivl1,  &
                                              ldivl2,lab,wtmp,coef1,lam,     &
                                              debtH,                         &
                                              tottrns, hlthtrs,              & ! add 2020/05/23
                                              wj, wn, ltctrs    ! add 2020/08/18
  real, dimension(ne,2),intent(in)        :: prop_agent            ! add 2020/08/18
  real, intent(in)                         :: debtR                  ! add 2020/05/20


  input(12) = gam_flag

  !
  ! Initialize parameters
  !
  alb            = par(1)
  aub            = par(2)
  alpha          = par(3)
  beta           = par(4)
  delt1          = par(5)
  deli1          = par(6)
  delt2          = par(7)
  deli2          = par(8)
  !gam            = par(9)
  gn             = par(10)
  gz             = par(11)
  lam            = par(12)
  ctran          = par(13)
  taud           = par(14)
  taul           = par(15)
  taupc          = par(16)
  taupu          = par(17)
  tfp            = par(18)
  thett1         = par(19)
  theti1         = par(20)
  thett2         = par(21)
  theti2         = par(22)

  alph1          = 1.-alpha
  thetl1         = 1.-thett1-theti1
  thetl2         = 1.-thett2-theti2
  irate          = x(1)
  tauc           = x(2)
  rt1            = (irate-1.)/(1.-taupc)+delt1
  rt2            = (irate-1.)/(1.-taupu)+delt2
  ri1            = irate-1.+deli1
  ri2            = irate-1.+deli2
  ldivl1         = thetl1*alpha/(thetl1*alpha+thetl2*alph1)
  ldivl2         = thetl2*alph1/(thetl1*alpha+thetl2*alph1)
  kdivlt1        = thett1/(thetl1*rt1)
  kdivlt2        = thett2/(thetl2*rt2)
  kdivli1        = theti1/(thetl1*ri1)
  kdivli2        = theti2/(thetl2*ri2)
  ydivl          = 2.*(tfp*kdivlt1**thett1*kdivli1**theti1*ldivl1)**alpha*   &
                      (tfp*kdivlt2**thett2*kdivli2**theti2*ldivl2)**alph1
  w              = ((thetl1*alpha+thetl2*alph1)*ydivl)**                     &
                   (1./(1.-((1.-thetl1)*alpha+(1.-thetl2)*alph1)))
  kdivlt1        = kdivlt1*w
  kdivlt2        = kdivlt2*w
  kdivli1        = kdivli1*w
  kdivli2        = kdivli2*w
  ydivl          = ydivl*w**((1.-thetl1)*alpha+(1.-thetl2)*alph1)

  bxi            = xi
  do l=1,nbeq
    bxi(jbeq(l)) = bxi(jbeq(l))-beq(l)
    bxi(jinh(l)) = bxi(jinh(l))+mu(jbeq(l))*beq(l)/mu(jinh(l))
  enddo

  kzero          = (aub-alb*n)/(aub-alb)
  input(1)       = par(9) ! gam
  input(2)       = 1.+tauc
  input(3)       = 1.-taul

  ! 年金計算  ! 2020/08/18
  write(*,*) ' '

  do l = 1,ne
    wj = 0
    wn = 0   ! working year
    do j = 1,nj
       h  = (l-1)*nj+j
       if ( eps(h) > 0 ) then
           wj = wj + w * eps(h) * 1/3
           wn = wn + 1
           pens(h) = 0
           health(h) = health_in(h)*w/3
           ltcare(h) = ltcare_in(h)*w/3
       else
          pens(h) =wj/wn * prop_agent(l,2)
           health(h) = health_in(h)*w/3
           ltcare(h) = ltcare_in(h)*w/3
       endif
    enddo
    write(*,*) ' type = ',l, ',  average income =',  wj/wn,  'pens=', pens(h)
    write(*,*) 'wk yrs =', wn, 'zeta =', prop_agent(l,1), 'rep rate of pens =', prop_agent(l,2)
   enddo
  ! For each productivity level, iterate backwards for ages nj to 1,
  ! solving dynamic programs recursively (Nishiyama-Smetter 2014 Step 2)
  !

  do l=1,ne
    h            = l*nj
    input(4)     = w*eps(h)
    do i=1,n
      a          = (aub*float(i-1)+alb*float(n-i))/float(n-1)
      ! input(5)   = ctran+bxi(h) !rem 2020/05/20
      ! input(5)   = ctran + pens(h) + health(h) + xi(h)   !add 2020/05/20

      input(5)   = ctran + pens(h) + medrate(h)*health(h) &  !add 2020/07/19
                 + ltcrate(h)*ltcare(h) + xi(h)            !add 2020/08/18
      !input(6)   = (irate-1.)*a
      input(6)   = ((1.-debtH)*irate+debtH*debtR-1.)*a !add 2020/05/20
      input(7)   = a
      input(8)   = (health(h)+ ltcare(h)) !add 2020/05/20  2020/08/18 for ltcare(h)
      input(9)   = dislab(h) !add 2020/06/07 労働不効用を年齢で可変に
      input(10)  = prop_agent(l,1) ! zeta add 2020/08/22
      input(11)  = prop_agent(l,3) ! gam add 2020/08/22
      call static(input,wtax,rtra,mtx,ntx,ntr, output)
      vnext(i)   = output(1)
    enddo
    opt(h,:)     = kzero


    do j=nj-1,1,-1
      h          = (l-1)*nj+j
      input(4)   = w*eps(h)
      do i=1,n
        a          = (aub*float(i-1)+alb*float(n-i))/float(n-1)
        if (opt(h,i)>0) then
          k1       = max(1,floor(float(opt(h,i))-max(2.,.1*n)))
          k2       = min(floor(float(opt(h,i))+max(2.,.1*n)),n)
        else
          k1       = 1
          k2       = n
        endif
        rhs        = -1.e+9
        do k=k1,k2
          at       = (aub*float(k-1)+alb*float(n-k))/float(n-1)
          !input(5) = ctran+bxi(h) ! rem 2020/05/20
          ! input(5) = ctran + pens(h) + health(h) + xi(h) ! add 2020/05/20
          input(5)   = ctran + pens(h) + medrate(h)*health(h)  &   !add 2020/07/19
                     + ltcrate(h)*ltcare(h) + xi(h)              !add 2020/08/18
          !input(6) = (irate-1.)*a   !rem 2020/05/20  !grossの金利だから1引いて利率に。日本版ではphiを反映する必要
          input(6)   = ((1.-debtH)*irate+debtH*debtR-1.)*a !add 2020/05/20
          input(7) = a-(1.+gz)*at *(1.-lam*(1.-ps(h)))
          input(8) = health(h) + ltcare(h)  !add 2020/05/20  2020/08/18 for ltcare(h)
          input(9) = dislab(h) !add 2020/06/07 労働不効用を年齢で可変に
          input(10)  = prop_agent(l,1) ! zeta add 2020/08/22
          input(11)  = prop_agent(l,3) ! gam add 2020/08/22
          call static(input,wtax,rtra,mtx,ntx,ntr, output)
          rhs(k)   = output(1)+beta*ps(h)*vnext(k)
        enddo
        kk         = maxloc(rhs)
        k          = kk(1)
        opt(h,i)   = k
        v(i)       = rhs(k)
      enddo
      vnext        = v
    enddo
  enddo

  !
  ! Update distributions of asset holdings
  ! (Nishiyama-Smetter 2014 Step 3)
  !
  do l=1,ne
    h              = (l-1)*nj+1
    hpi(h)         = kzero
    do j=1,nj-1
      h            = (l-1)*nj+j
      hpi(h+1)     = opt(h,hpi(h))
    enddo
  enddo
  !
  ! Add everything up
  !
  asset          = 0.
  tkap           = 0.
  lab            = 0.
  c              = 0.
  y              = 0.
  ltax           = 0.
  wpsi           = 0.
  rpsi           = 0.
  cj             = 0.
  lj             = 0.
  aj             = 0.
  apj            = 0.
  mtrj           = 0.
  ltaxj          = 0.
  wpsij          = 0.
  rpsij          = 0.

  do l=1,ne
    do j=1,nj
      h            = (l-1)*nj+j
      input(4)     = w*eps(h)
      i            = hpi(h)
      a            = (aub*float(i-1)+alb*float(n-i))/float(n-1)
      if (j<nj) then
        k          = hpi(h+1)
        at         = (aub*float(k-1)+alb*float(n-k))/float(n-1)
      else
        at         = 0.
      endif
      !input(5)     = ctran+bxi(h) !rem 2020/05/20
      ! input(5)     = ctran + pens(h) + health(h) + xi(h) !add 2020/05/20
      input(5)   = ctran + pens(h) + medrate(h)*health(h)  &       !add 2020/07/19
                     + ltcrate(h)*ltcare(h) + xi(h)              !add 2020/08/18

      !input(6)     = (irate-1.)*a  !rem 2020/05/20
      input(6)     = ((1.-debtH)*irate+debtH*debtR-1.)*a !add 2020/05/20
      input(7)     = a-(1.+gz)*at *(1.-lam*(1.-ps(h)))
      input(8)     = health(h) &            !add 2020/05/20
                     + ltcare(h)            !add 2020/08/18
      input(9)     = dislab(h) !add 2020/06/07 労働不効用を年齢で可変に
       input(10)  = prop_agent(l,1) ! zeta add 2020/08/22
       input(11)  = prop_agent(l,3) ! gam add 2020/08/22
      call static(input,wtax,rtra,mtx,ntx,ntr, output)
      asset        = asset+ mu(h)*a
      tkap         = tkap + mu(h)*(1.-ps(h))*(1.-lam)*at
      c            = c    + mu(h)*output(2)
      lab          = lab  + mu(h)*eps(h)*output(3)
      ltax         = ltax + mu(h)*output(5)
      wpsi         = wpsi + mu(h)*output(6)
      rpsi         = rpsi + mu(h)*output(7)
      cj(h)        = output(2)
      lj(h)        = output(3)
      aj(h)        = a
      apj(h)       = at
      mtrj(h)      = output(4)
      ltaxj(h)     = output(5)
      wpsij(h)     = output(6)
      rpsij(h)     = output(7)

       ! write(*,*)  output
    enddo
  enddo

  tranj          = sum(mu*bxi)
  hlthtrs        = sum(mu*medrate*health)         !追加 2020/07/19
  ltctrs         = sum(mu*ltcrate*ltcare)         !追加 2020/08/18  long-term-care transfer
  tottrns        = sum(mu*(pens))+  hlthtrs  +  ltctrs    !修正 2020/08/19
  y              = ydivl*lab
  krat1          = kdivli1/kdivlt1
  krat2          = kdivli2/kdivlt2
  kaprat         = (kdivlt2/kdivlt1)*(ldivl2/ldivl1)
  coef1          = (1.-taud)*(1.+(1.-taupc)*krat1)+(1.+(1.-taupu)*krat2)*kaprat
  if (ipar(1)==0) then
    kapt1        = (asset-par(23))/coef1
    kapt2        = kaprat*kapt1
    kapi1        = krat1*kapt1
    kapi2        = krat2*kapt2
    invi1        = ((1.+gz)*(1.+gn)-1.+deli1)*kapi1
    invi2        = ((1.+gz)*(1.+gn)-1.+deli2)*kapi2
    gdp          = y-invi1-invi2
    dshare       = par(23)/gdp
  else
    dshare       = par(23)
    kapt1        = (asset-dshare*y)/(coef1-dshare*((1.+gz)*(1.+gn)-1.+       &
                   deli1)*krat1-dshare*((1.+gz)*(1.+gn)-1.+deli2)*krat2*     &
                   kaprat)
    kapt2        = kaprat*kapt1
    kapi1        = krat1*kapt1
    kapi2        = krat2*kapt2
    invi1        = ((1.+gz)*(1.+gn)-1.+deli1)*kapi1
    invi2        = ((1.+gz)*(1.+gn)-1.+deli2)*kapi2
    gdp          = y-invi1-invi2
  endif
  invt1          = ((1.+gz)*(1.+gn)-1.+delt1)*kapt1
  invt2          = ((1.+gz)*(1.+gn)-1.+delt2)*kapt2
  lab1           = ldivl1*lab
  lab2           = ldivl2*lab
  y1             = tfp*(kapt1)**thett1*(kapi1)**theti1*(lab1)**thetl1
  y2             = tfp*(kapt2)**thett2*(kapi2)**theti2*(lab2)**thetl2
  p1             = alpha*y/y1
  p2             = alph1*y/y2
  cprof          = p1*y1-w*lab1-delt1*kapt1-invi1
  uprof          = p2*y2-w*lab2-delt2*kapt2-invi2
  cdiv           = p1*y1-w*lab1-invt1-invi1-taupc*cprof

  !
  ! Update guess for capital/labor ratio and transfers
  !
  if (ipar(2)==0) then
    gshare       = par(24)/gdp
  else
    gshare       = par(24)
  endif
  resx(1)        = rt1-thett1*alpha*y/kapt1
!  resx(2)        = tauc-(ctran+tranj+gshare*gdp-(1.+gz)*tkap-ltax+wpsi+      &
!                   rpsi-taud*cdiv-taupc*cprof-taupu*uprof-((1.+gz)*          &
!                   (1.+gn)-irate)*dshare*gdp)/c                              ! 分母からhealthを引く。分子にpenとhealthを足す cf sub.trans
 resx(2)        = tauc-(ctran+tottrns+tranj+gshare*gdp-(1.+gz)*tkap-ltax+wpsi+      &  ! tranjはゼロ
                   rpsi-taud*cdiv-taupc*cprof-taupu*uprof-((1.+gz)*          &
                   (1.+gn)-debtR)*dshare*gdp)/(c-hlthtrs)                              ! 分母からhealthを引く。分子にpenとhealthを足す cf sub.trans


 resx(3)        = dshare*gdp/(dshare*gdp+kapt1+kapt2+kapi1+kapi2)-debtH  !add 2020/05/20

end subroutine agg

end module mod_agg_jp_ex
