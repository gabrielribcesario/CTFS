PROGRAM ctft
    USE iso_fortran_env, ONLY : wp=>REAL64
    USE fourier
    IMPLICIT NONE

    ! gfortran dft.f95 fourier.f95 -o dft -O2

    ! For a causal signal x(t):
    ! x[n] = ∫ x(t) * δ(t - nT) * dt
    ! x[n] = x(nT), n = 0, 1, 2 ...
    ! I = ∫ Σ x(t) * δ(t - nT) * dt
    ! I = Σ x[n]

    PROCEDURE(f_x), POINTER :: func => null()
    INTEGER, PARAMETER :: ncoeff = 100 ! Number of Fourier coefficients
    REAL(wp), PARAMETER :: pi = 3.14159265358979323846_wp
    REAL(WP) :: trig_coeff(2, ncoeff)
    COMPLEX(wp) :: exp_coeff(ncoeff)

    func => wrapper

    ! Exponential series
    exp_coeff = ctft_coefficients(func, 0.0_wp, 2.0_wp*pi, ncoeff)
    ! Write exp series to file
    OPEN(1, file="exp_coeff.dat", status="replace", action="write")
    WRITE(1, '(A)') "Real,Imag"
    WRITE(1, '(SP,G0.12,",",G0.12)') exp_coeff
    CLOSE(1)

    ! Trigonometric series
    trig_coeff(:, 1) = [real(exp_coeff(1)), 0.0_wp]
    trig_coeff(1, 2:) = 2.0_wp * abs(exp_coeff(2:))
    trig_coeff(2, 2:) = atan2(imag(exp_coeff(2:)), &
                              real(exp_coeff(2:)))
    ! Write trig series to file
    OPEN(1, file="trig_coeff.dat", status="replace", action="write")
    WRITE(1, '(A)') "Amplitude,Phase"
    WRITE(1, '(SP,G0.12,",",G0.12)') trig_coeff
    CLOSE(1)

    CONTAINS
        PURE FUNCTION wrapper(x) RESULT(y)
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: x
            REAL(8) :: y
            y = cos(x)
        END FUNCTION

        PURE FUNCTION mse(y_true, y_false) RESULT(out)
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: y_true(:), y_false(size(y_true))
            REAL(8) :: out
            INTEGER :: i, len
            len = size(y_true)
            out = 0.0_wp
            DO i = 1, len
                out = out + (y_true(i) - y_false(i))**2
            END DO
            out = out / len
        END FUNCTION

        PURE FUNCTION rmse(y_true, y_false) RESULT(out)
            IMPLICIT NONE
            REAL(8), INTENT(IN) :: y_true(:), y_false(size(y_true))
            REAL(8) :: out
            out = sqrt(mse(y_true, y_false))
        END FUNCTION
END PROGRAM