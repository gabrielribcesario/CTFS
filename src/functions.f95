MODULE functions
    use iso_fortran_env, only : wp=>REAL64

    REAL(wp), PARAMETER :: M_PI = 3.14159265358979323846_wp 
    REAL(wp), PARAMETER :: M_1_PI = 1.0_wp / M_PI
    REAL(wp), PARAMETER :: M_TAU = 2.0_wp * M_PI
    REAL(wp), PARAMETER :: M_1_TAU = 1.0_wp / M_TAU

    PRIVATE
    PUBLIC :: f_t, sinusoidal, square, sawtooth, triangle

    ABSTRACT INTERFACE
        REAL(wp) PURE FUNCTION f_t(t, p) RESULT(y)
            USE iso_fortran_env, ONLY : wp=>REAL64
            IMPLICIT NONE
            REAL(wp), INTENT(IN) :: t
            REAL(wp), INTENT(IN) :: p(4) ! [t0, DC, A, phi]
        END FUNCTION
    END INTERFACE

    CONTAINS
        REAL(wp) PURE FUNCTION sinusoidal(t, p) RESULT(y)
            ! Periodic sine wave
            IMPLICIT NONE
            REAL(wp), INTENT(IN) :: t
            REAL(wp), INTENT(IN) :: p(4) ! [t0, DC, A, phi]
            y = p(2) + p(3) * sin(M_TAU * t / p(1) + p(4))
        END FUNCTION

        REAL(wp) PURE FUNCTION square(t, p) RESULT(y)
            ! Periodic square wave
            IMPLICIT NONE
            REAL(wp), INTENT(IN) :: t
            REAL(wp), INTENT(IN) :: p(4) ! [t0, DC, A, phi]
            y = p(2) + p(3) * (-1.0_wp)**(floor(2.0_wp * t / p(1) + p(4) * M_1_PI))
        END FUNCTION

        REAL(wp) PURE FUNCTION sawtooth(t, p) RESULT(y)
            ! Periodic sawtooth wave
            IMPLICIT NONE
            REAL(wp), INTENT(IN) :: t
            REAL(wp), INTENT(IN) :: p(4) ! [t0, DC, A, phi]
            REAL(wp) :: rad
            rad = t / p(1) + p(4) * 0.5_wp * M_1_PI
            y = p(2) + p(3) * 2.0_wp * (rad - floor(rad + 0.5_wp))
        END FUNCTION

        REAL(wp) PURE FUNCTION triangle(t, p) RESULT(y)
            ! Periodic triangle wave
            IMPLICIT NONE
            REAL(wp), INTENT(IN) :: t
            REAL(wp), INTENT(IN) :: p(4) ! [t0, DC, A, phi]
            y = p(2) + p(3) * (4.0_wp * abs(0.5_wp * p(1) - mod(t - 0.25_wp * p(1) + p(4) * p(1) * M_1_TAU, p(1))) / p(1) - 1.0_wp)
        END FUNCTION

        REAL(wp) PURE FUNCTION mse(y_true, y_false) RESULT(out)
            IMPLICIT NONE
            REAL(wp), INTENT(IN) :: y_true(:), y_false(size(y_true))
            INTEGER :: i, len
            len = size(y_true)
            out = 0.0_wp
            DO i = 1, len
                out = out + (y_true(i) - y_false(i))**2
            END DO
            out = out / len
        END FUNCTION

        REAL(wp) PURE FUNCTION rmse(y_true, y_false) RESULT(out)
            IMPLICIT NONE
            REAL(wp), INTENT(IN) :: y_true(:), y_false(size(y_true))
            out = sqrt(mse(y_true, y_false))
        END FUNCTION
END MODULE