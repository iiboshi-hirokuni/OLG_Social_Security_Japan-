module mod_func_jp

   implicit none
   contains


subroutine static(input, wtax, rtra, mtx, ntx, ntr, output)
    !
    ! Two equations in two unknowns (c, l)
    !
    ! -- if taxes are proportional:
    !
    !    gam*p*c/ (1-l) = lwedge * w
    !    p * c = nlinc + lwedge * w * l
    !
    ! -- if taxes aren't proportional
    !
    !    gam*p*c/ (1-l) = w- dtax (w*l) / dl
    !    p*c = nlinc+w*l-tax(w*l)
    !
    ! where
    !    u(c, 1-l) = log(c) + gam*log(1-l)
    !    nlinc = nonlabor income
    !    l in [0, 1] is checked

    implicit none
    integer                             :: i, maxit, gam_flag
    integer, intent(in)                 :: mtx, ntx, ntr
    real, dimension(12), intent(in)      :: input
    real, dimension(mtx,3), intent(in)  ::wtax,rtra
    real, dimension(7), intent(out)     :: output
    real                                :: gam, p, lwedge, w, otran, capinc, newa,  &
                                        hlth, nlinc, atw, del, res, dres, u, c, l,  &
                                        l0, lp, taxinc, tax, mrate, ltax,           &
                                        wpsi, rpsi, zeta
  gam_flag = input(12)

   if     ( gam_flag == 1 ) then
                                 gam  = input(1)
   elseif ( gam_flag == 2 ) then
                                 gam  = input(11)
   else
                                 gam  = input(9) !gamを年齢で可変に
   endif

   ! write(*,*) 'gam = ' , gam

    p         = input(2)
    lwedge    = input(3) ! Turn this on for proportional taxex
    w         = input(4)
    otran     = input(5)
    capinc    = input(6)
    newa      = input(7)
    hlth      = input(8) !health expenditure(added in Japan code)
    zeta      = input(10)
   ! zeta      = 1
    nlinc     = otran + capinc + newa !non-labor income
    wpsi      = 0.
    rpsi      = 0.

    if (w<=0.) then
      !
      ! Retiree
      !
      call comptax(capinc, rtra, mtx, ntr, rpsi, mrate)
      l          = 0.
      ! c          = hlth+(nlinc+rpsi-hlth)/p !消費にhealth支出を追加（以下も同様に）
      c          = (nlinc+rpsi-hlth)/p   ! revise 2020/08/19
      ! c          = (nlinc+rpsi)/p
      mrate      = 0.
      ltax       = 0.
      ! rpsi       = otran ! add 2020/07/19
     else
      !
      ! Worker
      !
      if (lwedge<.99) then
        !
        ! Facing proportional taxes (tax tableで計算する際には使っていない)
        !
        if ( zeta == 1 ) then
            atw      = lwedge*w
            l        = (atw-gam*(nlinc-hlth))/(atw+atw*gam) !!!労働供給を調整するか 2020/06/06
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

            c        = (nlinc+atw*l-hlth)/p
            ! c        = (nlinc+atw*l)/p
            mrate    = 1.-lwedge
            ltax     = (1.-lwedge)*w*l
            wpsi     = 0.
        if (l>1.) then
          l      = 1.
          atw    = lwedge*w
          c      = (nlinc+atw-hlth)/p
          ltax   = (1.-lwedge)*w
        endif
        if (l<0.) then
          l      = 0.
         ! c      = hlth +(nlinc-hlth)/p
          c      = (nlinc-hlth)/p
          ltax   = 0.
          mrate  = 0.
        endif
      else
        !
        ! Facing net tax function wtax(.)：ここ以下で計算
        !
        maxit    = 20
        l0       = .3 !2020/09/20　zeta, gamのサイズによって初期値を変える必要

        del      = .001
        do i=1, maxit
          taxinc = w*l0
          call comptax(taxinc,wtax,mtx,ntx,tax,mrate)
          ltax   = mrate*taxinc
          wpsi   = ltax-tax ! add hlth 20/07/19
          res    = gam*(nlinc+taxinc-ltax+wpsi)-w*(1.-mrate)*(1.-l0)**zeta
          lp     = l0+del
          taxinc = w*lp
          call comptax(taxinc,wtax,mtx,ntx,tax,mrate)
          ltax   = mrate*taxinc
          wpsi   = ltax-tax ! add hlth 20/07/19
          dres   = (gam*(nlinc+taxinc-ltax+wpsi)-w*(1-mrate)*(1.-lp)**zeta-res)/.001
          l0     = l0-res/dres
        enddo
        l         = l0

       ! l        = 0.8  !add 2020/06/07 labor supplyを固定（上のdo i=1　～　l = l0までをコメントしてこれを生かす

        taxinc   = w*l
        call comptax(taxinc,wtax,mtx,ntx,tax,mrate)
!        ltax     = mrate*taxinc   ! del 2020/06/07
!        wpsi     = ltax-tax      ! del 2020/06/07
        ltax     = mrate*taxinc
        wpsi     = ltax-tax
        !  c        = hlth+(nlinc+w*l-ltax+wpsi-hlth)/p  ! rev 2020/07/19
          c          = (nlinc+w*l-ltax+wpsi-hlth)/p   ! revise 2020/08/19
        ! c        = (nlinc+w*l-ltax+wpsi)/p
        if (l<0.) then
          l      = 0.
        !  c      = hlth+(nlinc-hlth)/p
          c      = (nlinc-hlth)/p ! rev 2020/09/22
          mrate  = 0.
          ltax   = 0.
          wpsi   = 0.
        endif
        if (l>=1.) then
         l         = 1.
          call comptax(w,wtax,mtx,ntx,tax,mrate)
!          ltax     = mrate*w    ! del 2020/06/07
!          wpsi     = ltax-tax   ! del 2020/06/07
          ltax     = mrate*taxinc
          wpsi     = ltax-tax+hlth
         !  c        = hlth+(nlinc+w*l-ltax+wpsi-hlth)/p  ! rev 2020/07/19
           c          = (nlinc+w*l-ltax+wpsi-hlth)/p   ! revise 2020/08/19
          ! c        = (nlinc+w*l-ltax+wpsi)/p
        endif
      endif
    endif

    if (c<=0.) then
      u        = -1.e+9
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

    output     =(/ u,c,l,mrate,ltax,wpsi,rpsi /)

end subroutine static


subroutine comptax(taxinc,taxtab,mtx,ltx, tax,mrate)
  !
  ! Compute taxes (or transfers) paid and the marginal rate given
  !   (1) taxable income
  !   (2) the tax table
  !
  implicit none
  integer                            :: i,idone
  integer, intent(in)                :: mtx,ltx
  real, dimension(mtx,3),intent(in)  :: taxtab
  real, intent(in)                   :: taxinc
  real, intent(out)                  :: tax,mrate

  tax    = 0.
  mrate  = 0.
  do i=1,ltx
    if (taxinc>=taxtab(i,1)) then
      tax    = taxtab(i,2)+taxtab(i,3)*taxinc
      mrate  = taxtab(i,3)
    endif
  enddo

end subroutine comptax


end module mod_func_jp


