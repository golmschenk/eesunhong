module polyroots_cmplx_roots_gen

    use iso_fortran_env
    use ieee_arithmetic

    implicit none

    private

    integer, parameter, public :: polyroots_module_rk = real64
    integer, parameter :: wp = polyroots_module_rk

    real(wp), parameter :: eps = epsilon(1.0_wp) !! machine epsilon
    real(wp), parameter :: pi = acos(-1.0_wp)
    real(wp), parameter :: deg2rad = pi / 180.0_wp

    public :: cmplx_roots_gen_laguerre

contains

    !*****************************************************************************************
    !>
    !  Modified from `cmplx_roots_gen` to comment in the use of only using Laguerre.
    !  See `cmplx_roots_gen` for other details.

    subroutine cmplx_roots_gen_laguerre(degree, poly, roots, polish_roots_after, use_roots_as_starting_points)

        implicit none

        integer, intent(in) :: degree !! degree of the polynomial and size of 'roots' array
        complex(wp), dimension(degree + 1), intent(in) :: poly !! coeffs of the polynomial, in order of increasing powers.
        complex(wp), dimension(degree), intent(inout) :: roots !! array which will hold all roots that had been found.
        !! If the flag 'use_roots_as_starting_points' is set to
        !! .true., then instead of point (0,0) we use value from
        !! this array as starting point for cmplx_laguerre
        logical, intent(in), optional :: polish_roots_after !! after all roots have been found by dividing
        !! original polynomial by each root found,
        !! you can opt in to polish all roots using full
        !! polynomial. [default is false]
        logical, intent(in), optional :: use_roots_as_starting_points !! usually we start Laguerre's
        !! method from point (0,0), but you can decide to use the
        !! values of 'roots' array as starting point for each new
        !! root that is searched for. This is useful if you have
        !! very rough idea where some of the roots can be.
        !! [default is false]

        complex(wp), dimension(:), allocatable :: poly2 !! `degree+1` array
        integer :: i, n, iter
        logical :: success
        complex(wp) :: coef, prev

        integer, parameter :: MAX_ITERS = 50
        ! constants needed to break cycles in the scheme
        integer, parameter :: FRAC_JUMP_EVERY = 10
        integer, parameter :: FRAC_JUMP_LEN = 10
        real(wp), dimension(FRAC_JUMP_LEN), parameter :: FRAC_JUMPS = &
                [0.64109297_wp, 0.91577881_wp, 0.25921289_wp, 0.50487203_wp, 0.08177045_wp, &
                        0.13653241_wp, 0.306162_wp, 0.37794326_wp, 0.04618805_wp, 0.75132137_wp] !! some random numbers
        real(wp), parameter :: FRAC_ERR = 10.0_wp * eps  !! fractional error
        !! (see. Adams 1967 Eqs 9 and 10)
        !! [2.0d-15 in original code]
        complex(wp), parameter :: zero = cmplx(0.0_wp, 0.0_wp, wp)
        complex(wp), parameter :: c_one = cmplx(1.0_wp, 0.0_wp, wp)

        ! initialize starting points
        if (present(use_roots_as_starting_points)) then
            if (.not.use_roots_as_starting_points) roots = zero
        else
            roots = zero
        end if

        ! skip small degree polynomials from doing Laguerre's method
        if (degree<=1) then
            if (degree==1) roots(1) = -poly(1) / poly(2)
            return
        endif

        allocate(poly2(degree + 1))
        poly2 = poly

        do n = degree, 3, -1

            ! find root with Laguerre's method
            call cmplx_laguerre(poly2, n, roots(n), iter, success)
            ! or
            ! find root with (Laguerre's method -> SG method -> Newton's method)
            ! call cmplx_laguerre2newton(poly2, n, roots(n), iter, success, 2)
            if (.not.success) then
                roots(n) = zero
                call cmplx_laguerre(poly2, n, roots(n), iter, success)
            endif

            ! divide the polynomial by this root
            coef = poly2(n + 1)
            do i = n, 1, -1
                prev = poly2(i)
                poly2(i) = coef
                coef = prev + roots(n) * coef
            enddo
            ! variable coef now holds a remainder - should be close to 0

        enddo

        ! find all but last root with Laguerre's method
        call cmplx_laguerre(poly2, 2, roots(2), iter, success)
        ! or
        ! call cmplx_laguerre2newton(poly2, 2, roots(2), iter, success, 2)
        if (.not.success) then
            call solve_quadratic_eq(roots(2), roots(1), poly2)
        else
            ! calculate last root from Viete's formula
            roots(1) = -(roots(2) + poly2(2) / poly2(3))
        endif

        if (present(polish_roots_after)) then
            if (polish_roots_after) then
                do n = 1, degree ! polish roots one-by-one with a full polynomial
                    call cmplx_laguerre(poly, degree, roots(n), iter, success)
                    !call cmplx_newton_spec(poly, degree, roots(n), iter, success)
                enddo
            endif
        end if

    contains

        recursive subroutine cmplx_laguerre(poly, degree, root, iter, success)

            !  Subroutine finds one root of a complex polynomial using
            !  Laguerre's method. In every loop it calculates simplified
            !  Adams' stopping criterion for the value of the polynomial.
            !
            !  For a summary of the method go to:
            !  http://en.wikipedia.org/wiki/Laguerre's_method

            implicit none

            integer, intent(in) :: degree !! a degree of the polynomial
            complex(wp), dimension(degree + 1), intent(in) :: poly !! an array of polynomial cooefs
            !! length = degree+1, poly(1) is constant
            !!```
            !!        1              2             3
            !!   poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
            !!```
            integer, intent(out) :: iter !! number of iterations performed (the number of polynomial
            !! evaluations and stopping criterion evaluation)
            complex(wp), intent(inout) :: root !! input: guess for the value of a root
            !! output: a root of the polynomial
            !!
            !! Uses 'root' value as a starting point (!!!!!)
            !! Remember to initialize 'root' to some initial guess or to
            !! point (0,0) if you have no prior knowledge.

            logical, intent(out) :: success !! is false if routine reaches maximum number of iterations

            real(wp) :: faq !! jump length
            complex(wp) :: p         !! value of polynomial
            complex(wp) :: dp        !! value of 1st derivative
            complex(wp) :: d2p_half  !! value of 2nd derivative
            integer :: i, k
            logical :: good_to_go
            complex(wp) :: denom, denom_sqrt, dx, newroot, fac_netwon, fac_extra, F_half, c_one_nth
            real(wp) :: ek, absroot, abs2p, one_nth, n_1_nth, two_n_div_n_1, stopping_crit2

            iter = 0
            success = .true.

            ! next if-endif block is an EXTREME failsafe, not usually needed, and thus turned off in this version.
            !if (.false.) then ! change false-->true if you would like to use caution about having first coefficient == 0
            if (degree<0) then
                write(*, *) 'Error: cmplx_laguerre: degree<0'
                return
            endif
            if (poly(degree + 1)==zero) then
                if (degree==0) return
                call cmplx_laguerre(poly, degree - 1, root, iter, success)
                return
            endif
            if (degree<=1) then
                if (degree==0) then  ! we know from previous check than poly(1) not equal zero
                    success = .false.
                    write(*, *) 'Warning: cmplx_laguerre: degree=0 and poly(1)/=0, no roots'
                    return
                else
                    root = -poly(1) / poly(2)
                    return
                endif
            endif
            !endif
            !  end EXTREME failsafe

            good_to_go = .false.
            one_nth = 1.0_wp / degree
            n_1_nth = (degree - 1.0_wp) * one_nth
            two_n_div_n_1 = 2.0_wp / n_1_nth
            c_one_nth = cmplx(one_nth, 0.0_wp, wp)

            do i = 1, MAX_ITERS
                ! prepare stoping criterion
                ek = abs(poly(degree + 1))
                absroot = abs(root)
                ! calculate value of polynomial and its first two derivatives
                p = poly(degree + 1)
                dp = zero
                d2p_half = zero
                do k = degree, 1, -1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                    d2p_half = dp + d2p_half * root
                    dp = p + dp * root
                    p = poly(k) + p * root    ! b_k
                    ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
                    ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
                    ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
                    ! Eq 8.
                    ek = absroot * ek + abs(p)
                enddo
                iter = iter + 1

                abs2p = real(conjg(p) * p)
                if (abs2p==0.0_wp) return
                stopping_crit2 = (FRAC_ERR * ek)**2
                if (abs2p<stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967)
                    ! do additional iteration if we are less than 10x from stopping criterion
                    if (abs2p<0.01d0 * stopping_crit2) then
                        return ! return immediately, because we are at very good place
                    else
                        good_to_go = .true. ! do one iteration more
                    endif
                else
                    good_to_go = .false.  ! reset if we are outside the zone of the root
                endif

                faq = 1.0_wp
                denom = zero
                if (dp/=zero) then
                    fac_netwon = p / dp
                    fac_extra = d2p_half / dp
                    F_half = fac_netwon * fac_extra

                    denom_sqrt = sqrt(c_one - two_n_div_n_1 * F_half)

                    !G=dp/p  ! gradient of ln(p)
                    !G2=G*G
                    !H=G2-2.0_wp*d2p_half/p  ! second derivative of ln(p)
                    !denom_sqrt=sqrt( (degree-1)*(degree*H-G2) )

                    ! NEXT LINE PROBABLY CAN BE COMMENTED OUT
                    if (real(denom_sqrt, wp)>=0.0_wp) then
                        ! real part of a square root is positive for probably all compilers. You can
                        ! test this on your compiler and if so, you can omit this check
                        denom = c_one_nth + n_1_nth * denom_sqrt
                    else
                        denom = c_one_nth - n_1_nth * denom_sqrt
                    endif
                endif
                if (denom==zero) then !test if demoninators are > 0.0 not to divide by zero
                    dx = (absroot + 1.0_wp) * exp(cmplx(0.0_wp, FRAC_JUMPS(mod(i, FRAC_JUMP_LEN) + 1) * 2 * pi, wp)) ! make some random jump
                else
                    dx = fac_netwon / denom
                    !dx=degree/denom
                endif

                newroot = root - dx
                if (newroot==root) return ! nothing changes -> return
                if (good_to_go) then       ! this was jump already after stopping criterion was met
                    root = newroot
                    return
                endif

                if (mod(i, FRAC_JUMP_EVERY)==0) then ! decide whether to do a jump of modified length (to break cycles)
                    faq = FRAC_JUMPS(mod(i / FRAC_JUMP_EVERY - 1, FRAC_JUMP_LEN) + 1)
                    newroot = root - faq * dx ! do jump of some semi-random length (0<faq<1)
                endif
                root = newroot
            enddo
            success = .false.
            ! too many iterations here
        end subroutine cmplx_laguerre

        subroutine solve_quadratic_eq(x0, x1, poly)

            ! Quadratic equation solver for complex polynomial (degree=2)

            implicit none

            complex(wp), intent(out) :: x0, x1
            complex(wp), dimension(*), intent(in) :: poly !! coeffs of the polynomial
            !! an array of polynomial cooefs,
            !! length = degree+1, poly(1) is constant
            !!```
            !!        1              2             3
            !!   poly(1) x^0 + poly(2) x^1 + poly(3) x^2
            !!```
            complex(wp) :: a, b, c, b2, delta

            a = poly(3)
            b = poly(2)
            c = poly(1)
            ! quadratic equation: a z^2 + b z + c = 0

            b2 = b * b
            delta = sqrt(b2 - 4.0_wp * (a * c))
            if (real(conjg(b) * delta, wp)>=0.0_wp) then  ! scalar product to decide the sign yielding bigger magnitude
                x0 = -0.5_wp * (b + delta)
            else
                x0 = -0.5_wp * (b - delta)
            endif
            if (x0==cmplx(0.0_wp, 0.0_wp, wp)) then
                x1 = cmplx(0.0_wp, 0.0_wp, wp)
            else ! Viete's formula
                x1 = c / x0
                x0 = x0 / a
            endif

        end subroutine solve_quadratic_eq

        recursive subroutine cmplx_laguerre2newton(poly, degree, root, iter, success, starting_mode)

            !  Subroutine finds one root of a complex polynomial using
            !  Laguerre's method, Second-order General method and Newton's
            !  method - depending on the value of function F, which is a
            !  combination of second derivative, first derivative and
            !  value of polynomial [F=-(p"*p)/(p'p')].
            !
            !  Subroutine has 3 modes of operation. It starts with mode=2
            !  which is the Laguerre's method, and continues until F
            !  becames F<0.50, at which point, it switches to mode=1,
            !  i.e., SG method (see paper). While in the first two
            !  modes, routine calculates stopping criterion once per every
            !  iteration. Switch to the last mode, Newton's method, (mode=0)
            !  happens when becomes F<0.05. In this mode, routine calculates
            !  stopping criterion only once, at the beginning, under an
            !  assumption that we are already very close to the root.
            !  If there are more than 10 iterations in Newton's mode,
            !  it means that in fact we were far from the root, and
            !  routine goes back to Laguerre's method (mode=2).
            !
            !  For a summary of the method see the paper: Skowron & Gould (2012)

            implicit none

            integer, intent(in) :: degree !! a degree of the polynomial
            complex(wp), dimension(degree + 1), intent(in) :: poly !! is an array of polynomial cooefs
            !! length = degree+1, poly(1) is constant
            !!```
            !!        1              2             3
            !!   poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
            !!```
            complex(wp), intent(inout) :: root !! input: guess for the value of a root
            !! output: a root of the polynomial
            !!
            !! Uses 'root' value as a starting point (!!!!!)
            !! Remember to initialize 'root' to some initial guess or to
            !! point (0,0) if you have no prior knowledge.
            integer, intent(in) :: starting_mode !! this should be by default = 2. However if you
            !! choose to start with SG method put 1 instead.
            !! Zero will cause the routine to
            !! start with Newton for first 10 iterations, and
            !! then go back to mode 2.
            integer, intent(out) :: iter !! number of iterations performed (the number of polynomial
            !! evaluations and stopping criterion evaluation)
            logical, intent(out) :: success !! is false if routine reaches maximum number of iterations

            real(wp) :: faq ! jump length
            complex(wp) :: p         ! value of polynomial
            complex(wp) :: dp        ! value of 1st derivative
            complex(wp) :: d2p_half  ! value of 2nd derivative
            integer :: i, j, k
            logical :: good_to_go
            complex(wp) :: denom, denom_sqrt, dx, newroot
            real(wp) :: ek, absroot, abs2p, abs2_F_half
            complex(wp) :: fac_netwon, fac_extra, F_half, c_one_nth
            real(wp) :: one_nth, n_1_nth, two_n_div_n_1
            integer :: mode
            real(wp) :: stopping_crit2

            iter = 0
            success = .true.
            stopping_crit2 = 0.0_wp  !  value not important, will be initialized anyway on the first loop (because mod(1,10)==1)

            ! next if-endif block is an EXTREME failsafe, not usually needed, and thus turned off in this version.
            !if (.false.)then ! change false-->true if you would like to use caution about having first coefficient == 0
            if (degree<0) then
                write(*, *) 'Error: cmplx_laguerre2newton: degree<0'
                return
            endif
            if (poly(degree + 1)==zero) then
                if (degree==0) return
                call cmplx_laguerre2newton(poly, degree - 1, root, iter, success, starting_mode)
                return
            endif
            if (degree<=1) then
                if (degree==0) then  ! we know from previous check than poly(1) not equal zero
                    success = .false.
                    write(*, *) 'Warning: cmplx_laguerre2newton: degree=0 and poly(1)/=0, no roots'
                    return
                else
                    root = -poly(1) / poly(2)
                    return
                endif
            endif
            !endif
            !  end EXTREME failsafe

            j = 1
            good_to_go = .false.
            mode = starting_mode  ! mode=2 full laguerre, mode=1 SG, mode=0 newton

            do ! infinite loop, just to be able to come back from newton, if more than 10 iteration there

                !------------------------------------------------------------- mode 2
                if (mode>=2) then  ! LAGUERRE'S METHOD
                    one_nth = 1.0_wp / degree
                    n_1_nth = (degree - 1.0_wp) * one_nth
                    two_n_div_n_1 = 2.0_wp / n_1_nth
                    c_one_nth = cmplx(one_nth, 0.0_wp, wp)

                    do i = 1, MAX_ITERS  !
                        faq = 1.0_wp

                        ! prepare stoping criterion
                        ek = abs(poly(degree + 1))
                        absroot = abs(root)
                        ! calculate value of polynomial and its first two derivatives
                        p = poly(degree + 1)
                        dp = zero
                        d2p_half = zero
                        do k = degree, 1, -1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                            d2p_half = dp + d2p_half * root
                            dp = p + dp * root
                            p = poly(k) + p * root    ! b_k
                            ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
                            ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
                            ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
                            ! Eq 8.
                            ek = absroot * ek + abs(p)
                        enddo
                        abs2p = real(conjg(p) * p, wp) !abs(p)
                        iter = iter + 1
                        if (abs2p==0.0_wp) return

                        stopping_crit2 = (FRAC_ERR * ek)**2
                        if (abs2p<stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967)
                            ! do additional iteration if we are less than 10x from stopping criterion
                            if (abs2p<0.01_wp * stopping_crit2) then ! ten times better than stopping criterion
                                return ! return immediately, because we are at very good place
                            else
                                good_to_go = .true. ! do one iteration more
                            endif
                        else
                            good_to_go = .false. ! reset if we are outside the zone of the root
                        endif

                        denom = zero
                        if (dp/=zero) then
                            fac_netwon = p / dp
                            fac_extra = d2p_half / dp
                            F_half = fac_netwon * fac_extra

                            abs2_F_half = real(conjg(F_half) * F_half, wp)
                            if (abs2_F_half<=0.0625_wp) then     ! F<0.50, F/2<0.25
                                ! go to SG method
                                if (abs2_F_half<=0.000625_wp) then ! F<0.05, F/2<0.025
                                    mode = 0 ! go to Newton's
                                else
                                    mode = 1 ! go to SG
                                endif
                            endif

                            denom_sqrt = sqrt(c_one - two_n_div_n_1 * F_half)

                            ! NEXT LINE PROBABLY CAN BE COMMENTED OUT
                            if (real(denom_sqrt, wp)>=0.0_wp) then
                                ! real part of a square root is positive for probably all compilers. You can
                                ! test this on your compiler and if so, you can omit this check
                                denom = c_one_nth + n_1_nth * denom_sqrt
                            else
                                denom = c_one_nth - n_1_nth * denom_sqrt
                            endif
                        endif
                        if (denom==zero) then !test if demoninators are > 0.0 not to divide by zero
                            dx = (abs(root) + 1.0_wp) * exp(cmplx(0.0_wp, FRAC_JUMPS(mod(i, FRAC_JUMP_LEN) + 1) * 2 * pi, wp)) ! make some random jump
                        else
                            dx = fac_netwon / denom
                        endif

                        newroot = root - dx
                        if (newroot==root) return ! nothing changes -> return
                        if (good_to_go) then       ! this was jump already after stopping criterion was met
                            root = newroot
                            return
                        endif

                        if (mode/=2) then
                            root = newroot
                            j = i + 1    ! remember iteration index
                            exit     ! go to Newton's or SG
                        endif

                        if (mod(i, FRAC_JUMP_EVERY)==0) then ! decide whether to do a jump of modified length (to break cycles)
                            faq = FRAC_JUMPS(mod(i / FRAC_JUMP_EVERY - 1, FRAC_JUMP_LEN) + 1)
                            newroot = root - faq * dx ! do jump of some semi-random length (0<faq<1)
                        endif
                        root = newroot
                    enddo ! do mode 2

                    if (i>=MAX_ITERS) then
                        success = .false.
                        return
                    endif

                endif ! if mode 2

                !------------------------------------------------------------- mode 1
                if (mode==1) then  ! SECOND-ORDER GENERAL METHOD (SG)

                    do i = j, MAX_ITERS  !
                        faq = 1.0_wp

                        ! calculate value of polynomial and its first two derivatives
                        p = poly(degree + 1)
                        dp = zero
                        d2p_half = zero
                        if (mod(i - j, 10)==0) then
                            ! prepare stoping criterion
                            ek = abs(poly(degree + 1))
                            absroot = abs(root)
                            do k = degree, 1, -1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                                d2p_half = dp + d2p_half * root
                                dp = p + dp * root
                                p = poly(k) + p * root    ! b_k
                                ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
                                ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
                                ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
                                ! Eq 8.
                                ek = absroot * ek + abs(p)
                            enddo
                            stopping_crit2 = (FRAC_ERR * ek)**2
                        else
                            do k = degree, 1, -1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                                d2p_half = dp + d2p_half * root
                                dp = p + dp * root
                                p = poly(k) + p * root    ! b_k
                            enddo
                        endif

                        abs2p = real(conjg(p) * p, wp) !abs(p)**2
                        iter = iter + 1
                        if (abs2p==0.0_wp) return

                        if (abs2p<stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967)
                            if (dp==zero) return
                            ! do additional iteration if we are less than 10x from stopping criterion
                            if (abs2p<0.01_wp * stopping_crit2) then ! ten times better than stopping criterion
                                return ! return immediately, because we are at very good place
                            else
                                good_to_go = .true. ! do one iteration more
                            endif
                        else
                            good_to_go = .false. ! reset if we are outside the zone of the root
                        endif

                        if (dp==zero) then !test if demoninators are > 0.0 not to divide by zero
                            dx = (abs(root) + 1.0_wp) * exp(cmplx(0.0_wp, FRAC_JUMPS(mod(i, FRAC_JUMP_LEN) + 1) * 2 * pi, wp)) ! make some random jump
                        else
                            fac_netwon = p / dp
                            fac_extra = d2p_half / dp
                            F_half = fac_netwon * fac_extra

                            abs2_F_half = real(conjg(F_half) * F_half, wp)
                            if (abs2_F_half<=0.000625_wp) then ! F<0.05, F/2<0.025
                                mode = 0 ! set Newton's, go there after jump
                            endif

                            dx = fac_netwon * (c_one + F_half)  ! SG
                        endif

                        newroot = root - dx
                        if (newroot==root) return ! nothing changes -> return
                        if (good_to_go) then       ! this was jump already after stopping criterion was met
                            root = newroot
                            return
                        endif

                        if (mode/=1) then
                            root = newroot
                            j = i + 1    ! remember iteration number
                            exit     ! go to Newton's
                        endif

                        if (mod(i, FRAC_JUMP_EVERY)==0) then ! decide whether to do a jump of modified length (to break cycles)
                            faq = FRAC_JUMPS(mod(i / FRAC_JUMP_EVERY - 1, FRAC_JUMP_LEN) + 1)
                            newroot = root - faq * dx ! do jump of some semi-random length (0<faq<1)
                        endif
                        root = newroot

                    enddo ! do mode 1

                    if (i>=MAX_ITERS) then
                        success = .false.
                        return
                    endif

                endif ! if mode 1

                !------------------------------------------------------------- mode 0
                if (mode==0) then  ! NEWTON'S METHOD

                    do i = j, j + 10  ! do only 10 iterations the most, then go back to full Laguerre's
                        faq = 1.0_wp

                        ! calculate value of polynomial and its first two derivatives
                        p = poly(degree + 1)
                        dp = zero
                        if (i==j) then ! calculate stopping crit only once at the begining
                            ! prepare stoping criterion
                            ek = abs(poly(degree + 1))
                            absroot = abs(root)
                            do k = degree, 1, -1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                                dp = p + dp * root
                                p = poly(k) + p * root    ! b_k
                                ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
                                ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
                                ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
                                ! Eq 8.
                                ek = absroot * ek + abs(p)
                            enddo
                            stopping_crit2 = (FRAC_ERR * ek)**2
                        else        !
                            do k = degree, 1, -1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                                dp = p + dp * root
                                p = poly(k) + p * root    ! b_k
                            enddo
                        endif
                        abs2p = real(conjg(p) * p, wp) !abs(p)**2
                        iter = iter + 1
                        if (abs2p==0.0_wp) return

                        if (abs2p<stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967)
                            if (dp==zero) return
                            ! do additional iteration if we are less than 10x from stopping criterion
                            if (abs2p<0.01_wp * stopping_crit2) then ! ten times better than stopping criterion
                                return ! return immediately, because we are at very good place
                            else
                                good_to_go = .true. ! do one iteration more
                            endif
                        else
                            good_to_go = .false. ! reset if we are outside the zone of the root
                        endif

                        if (dp==zero) then ! test if demoninators are > 0.0 not to divide by zero
                            dx = (abs(root) + 1.0_wp) * exp(cmplx(0.0_wp, FRAC_JUMPS(mod(i, FRAC_JUMP_LEN) + 1) * 2 * pi, wp)) ! make some random jump
                        else
                            dx = p / dp
                        endif

                        newroot = root - dx
                        if (newroot==root) return ! nothing changes -> return
                        if (good_to_go) then
                            root = newroot
                            return
                        endif

                        ! this loop is done only 10 times. So skip this check
                        !if (mod(i,FRAC_JUMP_EVERY)==0) then ! decide whether to do a jump of modified length (to break cycles)
                        !  faq=FRAC_JUMPS(mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1)
                        !  newroot=root-faq*dx ! do jump of some semi-random length (0<faq<1)
                        !endif
                        root = newroot

                    enddo ! do mode 0 10 times

                    if (iter>=MAX_ITERS) then
                        ! too many iterations here
                        success = .false.
                        return
                    endif
                    mode = 2 ! go back to Laguerre's. This happens when we were unable to converge in 10 iterations with Newton's

                endif ! if mode 0

            enddo ! end of infinite loop

            success = .false.

        end subroutine cmplx_laguerre2newton

    end subroutine cmplx_roots_gen_laguerre

end module polyroots_cmplx_roots_gen