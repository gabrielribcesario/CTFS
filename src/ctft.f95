PROGRAM ctft
    USE iso_fortran_env, ONLY : wp=>REAL64
    USE fourier
    USE functions
    IMPLICIT NONE

    ! gfortran ctft.f95 fourier.f95 functions.f95 -o ctft -O2

    PROCEDURE(f_t), POINTER :: func => null()
    INTEGER, PARAMETER :: ncoeff = 1000 ! Number of Fourier coefficients
    REAL(wp), PARAMETER :: M_PI = 3.14159265358979323846_wp
    REAL(wp) :: tol = 0.1_wp ! Stop calculation when RMSE < tol
    REAL(WP) :: t0, DC, A, phi
    REAL(WP) :: p(4)
    REAL(WP) :: trig_coeff(2, ncoeff) ! Change to pointer
    COMPLEX(wp) :: exp_coeff(ncoeff) ! Change to pointer

    t0 = 1.0_wp ! Fundamental period [s]
    DC = 0.0_wp ! DC offset
    A = 1.0_wp ! Amplitude
    phi = 0.0_wp!M_PI * 0.5_wp ! Phase shift
    p = [t0, DC, A, phi] ! Parameter vector
    func => sawtooth ! Function to be approximated

    ! Exponential series
    exp_coeff = ctft_coefficients(func, p, ncoeff)

    ! Write exp series to file
    OPEN(1, file="exp_coeff.dat", status="replace", action="write")
    WRITE(1, '(A)') "Real,Imag"
    WRITE(1, '(SP,G0.12,",",G0.12)') exp_coeff
    CLOSE(1)

    ! Compact trigonometric series
    trig_coeff(:, 1) = [real(exp_coeff(1)), 0.0_wp]
    trig_coeff(1, 2:) = 2.0_wp * abs(exp_coeff(2:))
    trig_coeff(2, 2:) = atan2(imag(exp_coeff(2:)), &
                              real(exp_coeff(2:)))

    ! Write trig series to file
    OPEN(1, file="trig_coeff.dat", status="replace", action="write")
    WRITE(1, '(A)') "Amplitude,Phase"
    WRITE(1, '(SP,G0.12,",",G0.12)') trig_coeff
    CLOSE(1)
END PROGRAM