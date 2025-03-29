program test
    use iso_fortran_env, only : wp=>real64
    use numerical_integration
    implicit none

    real(wp), parameter :: analytical_result = 0.84270079295, sqrt_pi = sqrt(3.14159265358979323846_wp)
    real(wp) :: a = 0.0_wp, b = 1.0_wp
    integer :: n = 20, m = 3
    procedure(f_x), pointer :: f => null()
    real(wp) :: I

    m = min(n, m) ! n >= m
    if (a >= b) then ! b > a
        b = a + b
    end if

    !print *, "(N, M) : ", n, ",", m

    f => gaussian

    print *, "ANALYTICAL RESULT : ", analytical_result
    print *, ""

    ! I = romberg(a, b, n, m, f)
    ! print *, "ROMBERG : ", I
    ! print *, "RELATIVE ERROR : ", (I - analytical_result) * 100.0_wp / analytical_result
    ! print *, ""

    I = trapezoid(a, b, 1000, f)
    print *, "TRAPEZOID : ", I
    print *, "RELATIVE ERROR : ", (I - analytical_result) * 100.0_wp / analytical_result
    print *, ""

    I = simpson(a, b, 1000, f)
    print *, "SIMPSON : ", I
    print *, "RELATIVE ERROR : ", (I - analytical_result) * 100.0_wp / analytical_result
    print *, ""

    contains
        pure function gaussian(x) result(y)
            implicit none
            real(8), intent(in) :: x
            real(8) :: y
            y = exp(-x**2) * 2.0_wp / sqrt_pi
        end function
end program