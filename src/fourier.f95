MODULE fourier
    USE iso_fortran_env, ONLY : dp=>REAL64
    USE functions, ONLY : f_t
    IMPLICIT NONE
    REAL(dp), PARAMETER :: M_TAU = 2.0_dp * 3.14159265358979323846_dp

    PRIVATE
    PUBLIC :: ctfs_coefficients

    CONTAINS
        SUBROUTINE ctfs_coefficients(f, p, eps, tol, max_iter, coeff, history, nrun)
            USE iso_fortran_env, ONLY : dp=>REAL64
            IMPLICIT NONE
            PROCEDURE(f_t), POINTER, INTENT(IN) :: f ! Function to be approximated
            INTEGER, INTENT(IN) :: max_iter ! Maximum number of Fourier coefficients
            REAL(dp), INTENT(IN) :: p(4) ! [t0, DC, A, phi]
            REAL(dp), INTENT(IN) :: eps ! Trapezoid integration maximum interval width
            REAL(dp), INTENT(IN) :: tol ! Stopping criteria tolerance
            COMPLEX(dp), INTENT(OUT) :: coeff(max_iter) ! Fourier coefficients
            REAL(dp), INTENT(OUT) :: history(2, max_iter) ! [MSE, mse] history
            INTEGER, INTENT(OUT) :: nrun ! Less than or equal to max_iter

            COMPLEX(dp) :: area, jw ! Helper variables for the trapezoidal rule integration step
            REAL(dp) :: w0, h_n, mse, kk ! More helper variables
            REAL(dp) :: t0 ! Period [s]
            REAL(dp), ALLOCATABLE :: signal(:) ! Function to be approximated, evaluated at the points that define the integration subintervals
            REAL(dp), ALLOCATABLE :: err(:) ! Residual, i.e. err(t) = signal(t) - Σ D_n*exp(jnw*t)
            INTEGER :: i, n
            INTEGER :: nsteps

            t0 = p(1) ! Period [s]. Integration bounds are always [0, t0]
            nsteps = ceiling(t0 / eps) ! ~1E-5-1E-6 width intervals should be more than accurate enough
            h_n = t0 / nsteps ! Integration interval width
            w0 = M_TAU / t0 ! Fundamental frequency [rad/s], w0 = 2pi/T0

            ! Evaluate function at the points that define the integration subintervals
            ALLOCATE(signal(nsteps + 1))
            ALLOCATE(err(nsteps + 1))
            DO CONCURRENT (i = 1 : nsteps + 1)
                signal(i) = f((i - 1) * h_n, p)
                err(i) = signal(i)
            END DO

            ! Calculate coefficients (complex exponential form)
            DO n = 1, max_iter
                jw = COMPLEX(0.0_dp, (n - 1) * w0) ! Complex frequency [rad/s]
                ! Trapezoidal rule integration
                area = (signal(1) + signal(nsteps + 1) * exp( -jw * t0 )) / 2.0_dp
                DO i = 2, nsteps
                    area = area + signal(i) * exp( -jw * (i - 1) * h_n )
                END DO
                area = area / nsteps ! area = area * h_n / t0; h_n = t0 / nsteps
                coeff(n) = area ! n-th coefficient

                kk = merge(1.0_dp, 2.0_dp, n == 1) ! If true, k=1; If false, k=2
                ! Update the residual
                err(1) = err(1) - kk * REAL(coeff(n))
                err(nsteps + 1) = err(nsteps + 1) - kk * REAL(coeff(n) * exp( jw * t0 ))
                ! Check stopping criteria: MSE = ∫[(f(t) - g(t))^2]dt / ∫dt. 
                ! Trapezoidal rule integration
                mse = (err(1)**2 + err(nsteps + 1)**2) / 2.0_dp
                DO i = 2, nsteps
                    err(i) = err(i) - kk * REAL(coeff(n) * exp( jw * (i - 1) * h_n )) ! Update the residual
                    mse = mse + err(i)**2
                END DO
                mse = mse / (nsteps * t0) ! MSE = MSE * h_n / t0^2; h_n = t0 / nsteps
                history(1, n) = mse ! MSE
                history(2, n) = sqrt( mse ) ! RMSE

                ! Stop if RMSE < tol
                IF (history(2, n) < tol) EXIT
            END DO
            nrun = min(n, max_iter)

            DEALLOCATE(signal)
            DEALLOCATE(err)
        END SUBROUTINE
END MODULE