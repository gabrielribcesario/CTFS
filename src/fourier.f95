MODULE fourier
    USE iso_fortran_env, ONLY : wp=>REAL64
    IMPLICIT NONE
    REAL(wp), PARAMETER :: tau = 2.0_wp * 3.14159265358979323846_wp
    COMPLEX(wp), PARAMETER :: imag_tau = (0.0_wp, tau)

    PRIVATE
    PUBLIC :: ctft_coefficients, f_x

    ABSTRACT INTERFACE
        PURE FUNCTION f_x(x) RESULT(y)
            USE iso_fortran_env, ONLY : wp=>REAL64
            IMPLICIT NONE
            REAL(wp), INTENT(IN) :: x
            REAL(wp) :: y
        END FUNCTION
    END INTERFACE

    CONTAINS
        FUNCTION ctft_coefficients(f, a, t0, ncoeff) RESULT(coeff)
            USE iso_fortran_env, ONLY : wp=>REAL64
            IMPLICIT NONE
            PROCEDURE(f_x), POINTER :: f ! Function to be approximated
            INTEGER, INTENT(IN) :: ncoeff ! Number of Fourier coefficients
            REAL(wp), INTENT(IN) :: t0 ! Period [s]
            REAL(wp), INTENT(IN) :: a ! Initial t
            COMPLEX(wp) :: coeff(ncoeff)
            COMPLEX(wp) :: area, complex_freq ! Helper variables for trapezoidal rule integration
            REAL(wp) :: w0, h_n
            INTEGER :: i, n, nintervals

            w0 = tau / t0 ! Fundamental frequency [rad/s], w0 = 2pi/T0
            nintervals = ceiling(t0 / 1.0E-4_wp) 
            h_n = t0 / nintervals 
            DO n = 1, ncoeff
                complex_freq = imag_tau * (n - 1) * w0
                ! Trapezoidal rule integration
                area = (integrand(a, complex_freq) + integrand(a + t0, complex_freq)) / 2.0_wp
                DO i = 1, nintervals - 1
                    area = area + integrand(a + i * h_n, complex_freq)
                END DO
                coeff(n) = area * h_n / t0 ! n-th coefficient
            END DO

            CONTAINS 
                FUNCTION integrand(t, jw) RESULT(y)
                    IMPLICIT NONE
                    COMPLEX(wp), INTENT(IN) :: jw ! Complex frequency [Hz]
                    REAL(wp), INTENT(IN) :: t ! Time [s]
                    COMPLEX(wp) :: y
                    y = f(t) * exp( -jw * t )
                END FUNCTION
        END FUNCTION
END MODULE