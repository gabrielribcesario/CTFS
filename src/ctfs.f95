PROGRAM ctft
    USE iso_fortran_env, ONLY : dp=>REAL64
    USE fourier
    USE functions
    IMPLICIT NONE

    ! gfortran -o ctfs functions.f95 fourier.f95 ctfs.f95 -O2

    PROCEDURE(f_t), POINTER :: func => null()
    INTEGER, PARAMETER :: max_iter = 1001 ! Maximum number of Fourier coefficients
    REAL(dp), PARAMETER :: M_PI = 3.14159265358979323846_dp
    REAL(dp) :: tol ! Stopping criteria (RMSE < tol)
    REAL(dp) :: eps = 1.0E-5_dp ! Integration accuracy
    REAL(dp) :: t0, DC, A, phi
    REAL(dp) :: p(4)
    REAL(dp) :: trig_coeff(2, max_iter) 
    COMPLEX(dp) :: exp_coeff(max_iter) 
    REAL(dp) :: history(2, max_iter) ! [MSE, RMSE] history
    CHARACTER(len=32) :: arg ! Command line arguments
    INTEGER :: nrun

    !! Parse arguments
    ! Fundamental period [s]
    CALL get_command_argument(1, arg)
    READ(arg, fmt=*) t0
    ! DC offset
    CALL get_command_argument(2, arg)
    READ(arg, fmt=*) DC
    ! Amplitude
    CALL get_command_argument(3, arg)
    READ(arg, fmt=*) A
    ! Phase shift
    CALL get_command_argument(4, arg)
    READ(arg, fmt=*) phi
    ! Function to be approximated
    CALL get_command_argument(5, arg)
    IF (ARG == "SINE") THEN
        func => sine
    ELSE IF (ARG == "SQUARE") THEN
        func => square 
    ELSE IF (ARG == "SAWTOOTH") THEN
        func => sawtooth 
    ELSE IF (ARG == "TRIANGLE") THEN
        func => triangle
    ELSE
        CALL exit(1)
    END IF
    ! Tolerance
    CALL get_command_argument(6, arg)
    READ(arg, fmt=*) tol

    p = [t0, DC, A, phi] ! Parameter vector

    ! Exponential series
    CALL ctfs_coefficients(func, p, eps, tol, max_iter, exp_coeff, history, nrun)

    IF (nrun == max_iter) THEN
        PRINT *, "SERIES DID NOT CONVERGE (RMSE >= TOL)"
        PRINT *, "RMSE = ", history(2, max_iter)
    END IF

    ! Write MSE/RMSE history
    OPEN(1, file="hist.dat", status="replace", action="write")
    WRITE(1, '(A)') "MSE,RMSE"
    WRITE(1, '(G0.12,",",G0.12)') history(:, 1:nrun)
    CLOSE(1)

    ! Write exp series to file
    OPEN(1, file="exp.dat", status="replace", action="write")
    WRITE(1, '(A)') "Real,Imag"
    WRITE(1, '(SP,G0.12,",",G0.12)') exp_coeff(1:nrun)
    CLOSE(1)

    ! Compact trigonometric series
    trig_coeff(:, 1) = [exp_coeff(1)%re, 0.0_dp] ! C_0 = D_0; Im(D_0) = 0
    trig_coeff(1, 2:nrun) = 2.0_dp * abs(exp_coeff(2:nrun)) ! C_n = 2*|D_n|; n > 0
    trig_coeff(2, 2:nrun) = atan2(exp_coeff(2:nrun)%im, &
                                  exp_coeff(2:nrun)%re) ! θ_n = ∠D_n; n > 0

    ! Write trig series to file
    OPEN(1, file="trig.dat", status="replace", action="write")
    WRITE(1, '(A)') "Amplitude,Phase"
    WRITE(1, '(SP,G0.12,",",G0.12)') trig_coeff(:, 1:nrun)
    CLOSE(1)

    CALL exit(0) ! EXIT_SUCESS
END PROGRAM