MODULE functions
    use iso_fortran_env, only : wp=>REAL64

    ABSTRACT INTERFACE
        PURE FUNCTION f_x(x) RESULT(y)
            USE iso_fortran_env, ONLY : wp=>REAL64
            IMPLICIT NONE
            REAL(wp), INTENT(IN) :: x
            REAL(wp) :: y
        END FUNCTION
    END INTERFACE

    CONTAINS
        PURE FUNCTION sinusoidal(x) RESULT(y)
        ! Periodic sine wave
        IMPLICIT NONE
        REAL(wp), INTENT(IN) :: x
        REAL(wp) :: y
        y = cos(x)
        END FUNCTION

        PURE FUNCTION square(x) RESULT(y)
        ! Periodic square wave
        IMPLICIT NONE
        REAL(wp), INTENT(IN) :: x
        REAL(wp) :: y
        y = cos(x)
        END FUNCTION

        PURE FUNCTION saw(x) RESULT(y)
        ! Periodic sawtooth wave
        IMPLICIT NONE
        REAL(wp), INTENT(IN) :: x
        REAL(wp) :: y
        y = cos(x)
        END FUNCTION

        PURE FUNCTION triangle(x) RESULT(y)
        ! Periodic triangle wave
        IMPLICIT NONE
        REAL(wp), INTENT(IN) :: x
        REAL(wp) :: y
        y = cos(x)
        END FUNCTION

        PURE FUNCTION mse(y_true, y_false) RESULT(out)
            IMPLICIT NONE
            REAL(wp), INTENT(IN) :: y_true(:), y_false(size(y_true))
            REAL(wp) :: out
            INTEGER :: i, len
            len = size(y_true)
            out = 0.0_wp
            DO i = 1, len
                out = out + (y_true(i) - y_false(i))**2
            END DO
            out = out / len
        END FUNCTION

        PURE FUNCTION rmse(y_true, y_false) RESULT(out)
            IMPLICIT NONE
            REAL(wp), INTENT(IN) :: y_true(:), y_false(size(y_true))
            REAL(wp) :: out
            out = sqrt(mse(y_true, y_false))
        END FUNCTION
END MODULE