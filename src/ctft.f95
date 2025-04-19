PROGRAM ctft
    USE iso_fortran_env, ONLY : dp=>REAL64
    USE fourier
    USE functions
    IMPLICIT NONE

    ! gfortran ctft.f95 fourier.f95 functions.f95 -o ctft -O2

    PROCEDURE(f_t), POINTER :: func => null()
    INTEGER, PARAMETER :: max_iter = 1000 ! Number of Fourier coefficients
    REAL(dp), PARAMETER :: M_PI = 3.14159265358979323846_dp
    REAL(dp) :: tol = 0.1_dp ! Stop calculation when RMSE < tol
    REAL(dp) :: eps = 1E-5_dp ! Trapezoid integration maximum interval width
    REAL(dp) :: t0, DC, A, phi
    REAL(dp) :: p(4)
    REAL(dp) :: trig_coeff(2, max_iter) ! Change to pointer
    COMPLEX(dp) :: exp_coeff(max_iter) ! Change to pointer
    INTEGER :: nrun

    t0 = 1.0_dp ! Fundamental period [s]
    DC = 0.0_dp ! DC offset
    A = 1.0_dp ! Amplitude
    phi = M_PI * 0.5_dp ! Phase shift
    p = [t0, DC, A, phi] ! Parameter vector
    func => sawtooth ! Function to be approximated

    ! Exponential series
    CALL ctft_coefficients(func, p, eps, tol, max_iter, exp_coeff, nrun)

    ! Write exp series to file
    OPEN(1, file="exp_coeff.dat", status="replace", action="write")
    WRITE(1, '(A)') "Real,Imag"
    WRITE(1, '(SP,G0.12,",",G0.12)') exp_coeff(1:nrun)
    CLOSE(1)

    ! Compact trigonometric series
    trig_coeff(:, 1) = [exp_coeff(1)%re, 0.0_dp]
    trig_coeff(1, 2:nrun) = 2.0_dp * abs(exp_coeff(2:nrun))
    trig_coeff(2, 2:nrun) = atan2(exp_coeff(2:nrun)%im, exp_coeff(2:nrun)%re)

    ! Write trig series to file
    OPEN(1, file="trig_coeff.dat", status="replace", action="write")
    WRITE(1, '(A)') "Amplitude,Phase"
    WRITE(1, '(SP,G0.12,",",G0.12)') trig_coeff(:, 1:nrun)
    CLOSE(1)
END PROGRAM