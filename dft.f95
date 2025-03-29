program dft
    use iso_fortran_env, only : wp=>real64
    implicit none

    ! For a causal signal x(t):
    ! x[n] = ∫ x(t) * δ(t - nT) * dt
    ! x[n] = x(nT); n = 0, 1, 2 ...
    ! I = ∫ Σ x(t) * δ(t - nT) * dt
    ! I = Σ x[n]

end program