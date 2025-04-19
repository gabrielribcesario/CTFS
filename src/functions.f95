MODULE functions
    use iso_fortran_env, only : dp=>REAL64

    REAL(dp), PARAMETER :: M_PI = 3.14159265358979323846_dp 
    REAL(dp), PARAMETER :: M_1_PI = 1.0_dp / M_PI
    REAL(dp), PARAMETER :: M_TAU = 2.0_dp * M_PI
    REAL(dp), PARAMETER :: M_1_TAU = 1.0_dp / M_TAU

    PRIVATE
    PUBLIC :: f_t, sinusoidal, square, sawtooth, triangle

    ABSTRACT INTERFACE
        REAL(dp) PURE FUNCTION f_t(t, p) RESULT(y)
            USE iso_fortran_env, ONLY : dp=>REAL64
            IMPLICIT NONE
            REAL(dp), INTENT(IN) :: t
            REAL(dp), INTENT(IN) :: p(4) ! [t0, DC, A, phi]
        END FUNCTION
    END INTERFACE

    CONTAINS
        REAL(dp) PURE FUNCTION sinusoidal(t, p) RESULT(y)
            ! Periodic sine wave
            ! y = DC + A*sin[ t/(2pi) + phi ]
            IMPLICIT NONE
            REAL(dp), INTENT(IN) :: t
            REAL(dp), INTENT(IN) :: p(4) ! [t0, DC, A, phi]
            y = p(2) + p(3) * sin(M_TAU * t / p(1) + p(4))
        END FUNCTION

        REAL(dp) PURE FUNCTION square(t, p) RESULT(y)
            ! Periodic square wave
            ! y = DC + A*(-1)^floor( 2*t/t0 + phi/pi )
            IMPLICIT NONE
            REAL(dp), INTENT(IN) :: t
            REAL(dp), INTENT(IN) :: p(4) ! [t0, DC, A, phi]
            INTEGER :: rad
            rad = floor(2.0_dp * t / p(1) + p(4) * M_1_PI)
            IF (mod(rad, 2) == 0) THEN
                y = p(2) + p(3)
            ELSE
                y = p(2) - p(3)
            END IF
        END FUNCTION

        REAL(dp) PURE FUNCTION sawtooth(t, p) RESULT(y)
            ! Periodic sawtooth wave
            ! y = DC + 2*A * [ t/t0 + phi/(2pi) - floor( t/t0 + phi/(2pi) + 1/2 ) ] 
            IMPLICIT NONE
            REAL(dp), INTENT(IN) :: t
            REAL(dp), INTENT(IN) :: p(4) ! [t0, DC, A, phi]
            REAL(dp) :: rad
            rad = t / p(1) + p(4) * M_1_TAU
            y = p(2) + p(3) * 2.0_dp * (rad - floor(rad + 0.5_dp))
        END FUNCTION

        REAL(dp) PURE FUNCTION triangle(t, p) RESULT(y)
            ! Periodic triangle wave
            ! y = DC + A*( 4 * | t0/2 - [t - t0/4 + phi*t0/(2pi)] mod(t0) | / t0 - 1 )
            IMPLICIT NONE
            REAL(dp), INTENT(IN) :: t
            REAL(dp), INTENT(IN) :: p(4) ! [t0, DC, A, phi]
            y = p(2) + p(3) * (4.0_dp * abs(0.5_dp * p(1) - mod(t - 0.25_dp * p(1) + p(4) * p(1) * M_1_TAU, p(1))) / p(1) - 1.0_dp)
        END FUNCTION
END MODULE