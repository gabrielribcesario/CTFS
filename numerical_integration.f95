module numerical_integration
    implicit none

    private
    public :: romberg, f_x, simpson, trapezoid

    abstract interface
        pure function f_x(x) result(y)
            use iso_fortran_env, only : wp=>real64
            implicit none
            real(wp), intent(in) :: x
            real(wp) :: y
        end function
    end interface

    contains
        recursive pure function romberg(a, b, n, m, f) result(area)
            use iso_fortran_env, only : wp=>real64
            implicit none
            procedure(f_x), pointer, intent(in) :: f
            real(wp), intent(in), value :: a, b
            integer, intent(in), value :: n, m
            real(wp) :: h_n, area
            integer :: k

            if (m == 0 .and. n == 0) then
                area = (f(b) + f(a)) * (b - a)
            else if (m == 0) then
                h_n = (b - a) * (0.5_wp ** n)
                area = 0.0_wp
                do k = 1, 2 ** (n - 1)
                    area = area + f(a + (2.0_wp * k - 1.0_wp) * h_n)
                end do
                area = area * h_n + 0.5_wp * romberg(a, b, n - 1, 0, f)
            else ! m != 0 .and. n != 0
                area = ( 4.0_wp ** m * romberg(a, b, n, m - 1, f) - romberg(a, b, n - 1, m - 1, f) ) / (4.0_wp ** m - 1.0_wp)
            end if
        end function

        pure function simpson(a, b, n, f) result(area)
            use iso_fortran_env, only : wp=>real64
            implicit none
            procedure(f_x), pointer, intent(in) :: f
            real(wp), intent(in) :: a, b
            integer, intent(in) :: n
            real(wp) :: h_n, area
            integer :: i

            h_n = (b - a) / n
            if (mod(n, 2) /= 0) then ! If n is odd
                ! Apply the 1/3 rule to the first n - 3 intervals
                area = f(a) + 4.0_wp * f(a + h_n)
                do i = 2, n - 4, 2
                    area = area + 2.0_wp * f(a + i * h_n) + 4.0_wp * f(a + (i + 1) * h_n)
                end do
                area = (area + f(a + i * h_n)) * h_n / 3.0_wp
                ! Apply the 3/8 rule to the remaining 3 intervals
                area = area + ( f(a + i * h_n) &
                                + 3.0_wp * ( f(a + (i + 1) * h_n) + f(a + (i + 2) * h_n) ) &
                                + f(b) ) * h_n * 3.0_wp / 8.0_wp
            else ! If n is even
                ! Apply the 1/3 rule to all intervals
                area = f(a) + 4.0_wp * f(a + h_n) + f(b)
                do i = 2, n - 1, 2
                    area = area + 2.0_wp * f(a + i * h_n) + 4.0_wp * f(a + (i + 1) * h_n)
                end do
                area = area * h_n / 3.0_wp
            end if
        end function

        pure function trapezoid(a, b, n, f) result(area)
            use iso_fortran_env, only : wp=>real64
            implicit none
            procedure(f_x), pointer, intent(in) :: f
            real(wp), intent(in) :: a, b
            integer, intent(in) :: n
            real(wp) :: h_n, area
            integer :: i

            h_n = (b - a) / n
            area = (f(a) + f(b)) / 2.0_wp
            do i = 1, n - 1
                area = area + f(a + i * h_n)
            end do
            area = area * h_n
        end function
end module