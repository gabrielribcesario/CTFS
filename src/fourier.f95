module fourier
    use iso_fortran_env, only : wp=>real64
    implicit none
    real(wp), parameter :: tau = 2.0_wp * 3.14159265358979323846_wp
    complex(wp), parameter :: imag_tau = (0.0_wp, tau)

    private
    public :: fourier_transform, f_x

    abstract interface
        pure function f_x(x) result(y)
            use iso_fortran_env, only : wp=>real64
            implicit none
            real(wp), intent(in) :: x
            real(wp) :: y
        end function
    end interface

    contains
        function fourier_transform(f, a, t0, ncoeff) result(coeff)
            use iso_fortran_env, only : wp=>real64
            implicit none
            procedure(f_x), pointer :: f ! Function to be approximated
            integer, intent(in) :: ncoeff ! Number of Fourier coefficients
            real(wp), intent(in) :: t0 ! Period [s]
            real(wp), intent(in) :: a ! Initial t
            complex(wp) :: coeff(ncoeff)
            complex(wp) :: area ! Helper variable for trapezoidal rule integration
            real(wp) :: w0, h_n
            integer :: i, n, nintervals

            nintervals = ceiling(t0 / 5.0E-5_wp)
            h_n = t0 / nintervals
            ! Fundamental frequency [rad/s], w0 = 2pi / T0 
            w0 = tau / t0
            do n = 0, ncoeff - 1
                ! Trapezoidal rule integration
                area = (integrand(a) + integrand(a + t0)) / 2.0_wp
                do i = 1, nintervals - 1
                    area = area + integrand(a + i * h_n)
                end do
                ! n-th coefficient
                coeff(1 + n) = area * h_n / t0
            end do

            contains 
                function integrand(t) result(y)
                    implicit none
                    real(wp), intent(in) :: t
                    complex(wp) :: y
                    y = f(t) * exp( -imag_tau * n * w0 * t)
                end function
        end function
end module