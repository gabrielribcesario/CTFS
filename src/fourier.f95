MODULE fourier
    USE iso_fortran_env, ONLY : wp=>REAL64
    USE functions, ONLY : f_t
    USE omp_lib
    IMPLICIT NONE
    REAL(wp), PARAMETER :: M_TAU = 2.0_wp * 3.14159265358979323846_wp

    PRIVATE
    PUBLIC :: ctft_coefficients

    CONTAINS
        COMPLEX(wp) FUNCTION ctft_coefficients(f, p, ncoeff) RESULT(coeff)
            USE iso_fortran_env, ONLY : wp=>REAL64
            IMPLICIT NONE
            PROCEDURE(f_t), POINTER, INTENT(IN) :: f ! Function to be approximated
            INTEGER, INTENT(IN) :: ncoeff ! Number of Fourier coefficients
            REAL(wp), INTENT(IN) :: p(4) ! [t0, DC, A, phi]
            DIMENSION coeff(ncoeff)

            COMPLEX(wp) :: area, jw ! Helper variables for trapezoidal rule integration
            REAL(wp) :: t0 ! Period [s]
            REAL(wp) :: w0, h_n, t
            INTEGER :: i, n
            INTEGER :: nintervals

            t0 = p(1)
            nintervals = ceiling(t0 / 1.0E-5_wp) ! ~1E-5 width intervals should be more than accurate enough
            h_n = t0 / nintervals
            w0 = M_TAU / t0 ! Fundamental frequency [rad/s], w0 = 2pi/T0
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jw, area, t) 
            DO n = 1, ncoeff
                jw = COMPLEX(0.0_wp, (n - 1) * w0) ! Complex frequency [rad/s]
                ! Trapezoidal rule integration
                area = (f(0.0_wp, p) + f(t0, p) * exp( -jw * t0 )) / 2.0_wp
                DO i = 1, nintervals - 1
                    t = i * h_n
                    area = area + f(t, p) * exp( -jw * t )
                END DO
                area = area / nintervals ! area * h_n / t0
                coeff(n) = area ! n-th coefficient
            END DO
            !$OMP END PARALLEL DO
        END FUNCTION
END MODULE