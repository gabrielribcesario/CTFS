module fourier
    use iso_fortran_env, only : dp=>real64
    use functions, only : f_t
    implicit none
    real(dp), parameter :: M_TAU = 2.0_dp * 3.14159265358979323846_dp

    private
    public :: ctfs_coefficients

    contains
    subroutine ctfs_coefficients(f, p, eps, tol, max_iter, coeff, history, nrun)
        use iso_fortran_env, only : dp=>real64
        implicit none
        procedure(f_t), POINTER, intent(in) :: f ! Function to be approximated
        integer, intent(in) :: max_iter ! Maximum number of Fourier coefficients
        real(dp), intent(in) :: p(4) ! [t0, DC, A, phi]
        real(dp), intent(in) :: eps ! Trapezoid integration maximum interval width
        real(dp), intent(in) :: tol ! Stopping criteria tolerance
        complex(dp), intent(out) :: coeff(max_iter) ! Fourier coefficients
        real(dp), intent(out) :: history(2, max_iter) ! [MSE, mse] history
        integer, intent(out) :: nrun ! Less than or equal to max_iter

        complex(dp) :: area, jw ! Helper variables for the trapezoidal rule integration step
        real(dp) :: w0, h_n, mse, kk ! More helper variables
        real(dp) :: t0 ! Period [s]
        real(dp), ALLOCATABLE :: signal(:) ! Function to be approximated, evaluated at the points that define the integration subintervals
        real(dp), ALLOCATABLE :: err(:) ! Residual, i.e. err(t) = signal(t) - Σ D_n*exp(jnw*t)
        integer :: i, n
        integer :: nsteps

        t0 = p(1) ! Period [s]. Integration bounds are always [0, t0]
        nsteps = ceiling(t0 / eps) ! ~1E-5 width intervals should be accurate enough
        h_n = t0 / nsteps ! Integration interval width
        w0 = M_TAU / t0 ! Fundamental frequency [rad/s], w0 = 2pi/T0

        ! Evaluate function at the points that define the integration subintervals
        allocate(signal(nsteps + 1))
        allocate(err(nsteps + 1))
        do concurrent (i = 1 : nsteps + 1)
            signal(i) = f((i - 1) * h_n, p)
            err(i) = signal(i)
        end do

        ! Calculate coefficients (complex exponential form)
        do n = 1, max_iter
            ! complex frequency [rad/s]
            jw = complex(0.0_dp, (n - 1) * w0) 

            ! Trapezoidal rule integration
            area = signal(1) ! f(0) = f(t0), exp(-jwt0) = 1
            do i = 2, nsteps
                area = area + signal(i) * exp( -jw * (i - 1) * h_n )
            end do
            area = area / nsteps ! area = area * h_n / t0; h_n = t0 / nsteps
            coeff(n) = area ! n-th coefficient

            kk = merge(1.0_dp, 2.0_dp, n == 1) ! If true, k=1; If false, k=2

            ! Update the residual
            err(1) = err(1) - kk * real(coeff(n))
            err(nsteps + 1) = err(nsteps + 1) - kk * real(coeff(n) * exp( jw * t0 ))

            ! Check stopping criterion: MSE = ∫[(f(t) - g(t))^2]dt / ∫dt
            mse = err(1)**2 ! f(0) = f(t0), exp(-jwt0) = 1
            do i = 2, nsteps
                err(i) = err(i) - kk * real(coeff(n) * exp( jw * (i - 1) * h_n )) ! Update the residual
                mse = mse + err(i)**2
            end do
            mse = mse / (nsteps * t0) ! MSE = MSE * h_n / t0^2; h_n = t0 / nsteps
            history(1, n) = mse ! MSE
            history(2, n) = sqrt( mse ) ! RMSE

            ! Stop if RMSE < tol
            if (history(2, n) < tol) exit
        end do
        nrun = min(n, max_iter)
    end subroutine
end module