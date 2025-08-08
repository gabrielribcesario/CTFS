program ctfs
    use iso_fortran_env, ONLY : dp=>real64, error_unit
    use fourier
    use functions
    implicit none

    procedure(f_t), pointer :: func => null()
    integer, parameter :: max_iter = 1001 ! Maximum number of Fourier coefficients
    real(dp), parameter :: M_PI = 3.14159265358979323846_dp
    real(dp) :: tol ! Stopping criteria (RMSE < tol)
    real(dp) :: eps = 1.0E-5_dp ! Integration accuracy
    real(dp) :: t0, DC, A, phi
    real(dp) :: p(4)
    real(dp) :: trig_coeff(2, max_iter) 
    complex(dp) :: exp_coeff(max_iter) 
    real(dp) :: history(2, max_iter) ! [MSE, RMSE] history
    character(len=32) :: arg ! Command line arguments
    integer :: nrun

    !! Parse arguments
    ! Fundamental period [s]
    call get_command_argument(1, arg)
    read(arg, fmt=*) t0
    ! DC offset
    call get_command_argument(2, arg)
    read(arg, fmt=*) DC
    ! Amplitude
    call get_command_argument(3, arg)
    read(arg, fmt=*) A
    ! Phase shift
    call get_command_argument(4, arg)
    read(arg, fmt=*) phi
    ! Function to be approximated
    call get_command_argument(5, arg)
    if (ARG == "SINE") then
        func => sine
    else if (ARG == "SQUARE") then
        func => square 
    else if (ARG == "SAWTOOTH") then
        func => sawtooth 
    else if (ARG == "TRIANGLE") then
        func => triangle
    else
        call exit(1)
    end if
    ! Tolerance
    call get_command_argument(6, arg)
    read(arg, fmt=*) tol

    p = [t0, DC, A, phi] ! parameter vector

    ! Exponential series
    call ctfs_coefficients(func, p, eps, tol, max_iter, exp_coeff, history, nrun)

    if (nrun == max_iter) then
        write(error_unit, '(A)') "SERIES DID NOT CONVERGE (RMSE >= TOL)"
        write(error_unit, '(A,G0.12)') "RMSE = ", history(2, max_iter)
    end if

    ! Write MSE/RMSE history
    open(1, file="hist.dat", status="replace", action="write")
    write(1, '(A)') "MSE,RMSE"
    write(1, '(G0.12,",",G0.12)') history(:, 1:nrun)
    close(1)

    ! Write exp series to file
    open(1, file="exp.dat", status="replace", action="write")
    write(1, '(A)') "Real,Imag"
    write(1, '(SP,G0.12,",",G0.12)') exp_coeff(1:nrun)
    close(1)

    ! Compact trigonometric series
    trig_coeff(:, 1) = [exp_coeff(1)%re, 0.0_dp] ! C_0 = D_0; Im(D_0) = 0
    trig_coeff(1, 2:nrun) = 2.0_dp * abs(exp_coeff(2:nrun)) ! C_n = 2*|D_n|; n > 0
    trig_coeff(2, 2:nrun) = atan2(exp_coeff(2:nrun)%im, &
                                  exp_coeff(2:nrun)%re) ! θ_n = ∠D_n; n > 0

    ! Write trig series to file
    open(1, file="trig.dat", status="replace", action="write")
    write(1, '(A)') "Amplitude,Phase"
    write(1, '(SP,G0.12,",",G0.12)') trig_coeff(:, 1:nrun)
    close(1)

    call exit(0) ! EXIT_SUCESS
end program