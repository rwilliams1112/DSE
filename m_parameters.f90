module m_parameters

    implicit none

    integer, parameter  :: r8 = kind(1.d0)
    real(r8), parameter :: PI =                                         &
            3.1415926535897932384626433832795028841971693993751_r8
    real(r8), parameter :: N_C = 3._r8, C_F = 4._r8/3._r8, C_A = 3._r8

    complex(r8), parameter :: i_  = cmplx(0._r8, 1._r8, kind=r8)
    complex(r8), parameter :: C0_ = cmplx(0._r8, 0._r8, kind=r8)
    complex(r8), parameter :: C1_ = cmplx(1._r8, 0._r8, kind=r8)

end module
