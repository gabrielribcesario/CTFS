module functions
    use iso_fortran_env, only : dp=>real64

    real(dp), parameter :: M_PI = 3.14159265358979323846_dp 
    real(dp), parameter :: M_1_PI = 1.0_dp / M_PI
    real(dp), parameter :: M_TAU = 2.0_dp * M_PI
    real(dp), parameter :: M_1_TAU = 1.0_dp / M_TAU

    private
    public :: f_t, sine, square, sawtooth, triangle

    abstract interface
        real(dp) pure function f_t(t, p) result(y)
            use iso_fortran_env, only : dp=>real64
            implicit none
            real(dp), intent(in) :: t
            real(dp), intent(in) :: p(4) ! [t0, DC, A, phi]
        end function
    end interface

    contains
    real(dp) pure function sine(t, p) result(y)
        ! Periodic sine wave
        ! y = DC + A*sin[ t/(2pi) + phi ]
        implicit none
        real(dp), intent(in) :: t
        real(dp), intent(in) :: p(4) ! [t0, DC, A, phi]
        y = p(2) + p(3) * sin(M_TAU * t / p(1) + p(4))
    end function

    real(dp) pure function square(t, p) result(y)
        ! Periodic square wave
        ! y = DC + A*(-1)^floor( 2*t/t0 + phi/pi )
        implicit none
        real(dp), intent(in) :: t
        real(dp), intent(in) :: p(4) ! [t0, DC, A, phi]
        integer :: rad
        rad = floor(2.0_dp * t / p(1) + p(4) * M_1_PI)
        if (mod(rad, 2) == 0) then
            y = p(2) + p(3)
        else
            y = p(2) - p(3)
        end if
    end function

    real(dp) pure function sawtooth(t, p) result(y)
        ! Periodic sawtooth wave
        ! y = DC + 2*A * [ t/t0 + phi/(2pi) - floor( t/t0 + phi/(2pi) + 1/2 ) ] 
        implicit none
        real(dp), intent(in) :: t
        real(dp), intent(in) :: p(4) ! [t0, DC, A, phi]
        real(dp) :: rad
        rad = t / p(1) + p(4) * M_1_TAU
        y = p(2) + p(3) * 2.0_dp * (rad - floor(rad + 0.5_dp))
    end function

    real(dp) pure function triangle(t, p) result(y)
        ! Periodic triangle wave
        ! y = DC + A*( 4 * | t0/2 - [t - t0/4 + phi*t0/(2pi)] mod(t0) | / t0 - 1 )
        implicit none
        real(dp), intent(in) :: t
        real(dp), intent(in) :: p(4) ! [t0, DC, A, phi]
        real(dp) :: rad
        rad = t - 0.25_dp * p(1) + p(4) * p(1) * M_1_TAU
        y = p(2) + p(3) * (4.0_dp * abs((0.5_dp * p(1) - modulo(rad, p(1))) / p(1)) - 1.0_dp)
    end function
end module