program dft
    use iso_fortran_env, only : wp=>real64
    use fourier
    implicit none

    ! gfortran dft.f95 fourier.f95 -o dft -O2

    ! For a causal signal x(t):
    ! x[n] = ∫ x(t) * δ(t - nT) * dt
    ! x[n] = x(nT); n = 0, 1, 2 ...
    ! I = ∫ Σ x(t) * δ(t - nT) * dt
    ! I = Σ x[n]

    integer, parameter :: ncoeff = 100 ! Number of Fourier coefficients
    real(wp), parameter :: pi = 3.14159265358979323846_wp
    complex(wp) :: coeff(ncoeff)
    procedure(f_x), pointer :: func => null()

    func => wrapper

    coeff = fourier_transform(func, 0.0_wp, pi, ncoeff)
    open(1, file="c.dat", status="replace", action="write")
    write(1, '(G0.12," ",SP,G0.12,"i")') coeff
    close(1)

    contains
        pure function wrapper(x) result(y)
            implicit none
            real(8), intent(in) :: x
            real(8) :: y
            y = exp(-x / 2.0_8)
        end function

        pure function mse(y_true, y_false) result(out)
            implicit none
            real(8), intent(in) :: y_true(:), y_false(size(y_true))
            real(8) :: out, aux
            integer :: i, len

            len = size(y_true)
            do i = 1, len
                aux = y_true(i) - y_false(i)
                out = out + aux * aux
            end do
            out = sqrt(out) / len
        end function

        pure function rmse(y_true, y_false) result(out)
            implicit none
            real(8), intent(in) :: y_true(:), y_false(size(y_true))
            real(8) :: out
            integer :: len

            len = size(y_true)
            out = mse(y_true, y_false)
            out = sqrt(out * len) / len
        end function
end program