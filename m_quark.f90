module m_quark

    ! Global parameters
    use m_parameters

    ! Global libraries
    use m_library
    use m_dirac_matrices

    implicit none
    private

    type quark

        logical  :: is_allocated = .false.    ! allocataion guard

        ! Renormalization constants
        real(r8) :: mass
        real(r8) :: mu_sq
        real(r8) :: Z_2 = 1._r8
        real(r8) :: Z_M = 1._r8

        ! Grid
        integer  :: n
        real(r8) :: eps_sq
        real(r8) :: lam_sq

        ! Spline interpolator
        real(r8), allocatable :: log_p_sq(:)
        real(r8), allocatable ::     p_sq(:)

        real(r8), allocatable ::     A   (:)
        real(r8), allocatable ::     A2  (:)
        real(r8), allocatable ::     B   (:)
        real(r8), allocatable ::     B2  (:)

    contains

        procedure :: create    => quark_create
        procedure :: delete    => quark_delete
        procedure :: update    => quark_update

        ! Access
        procedure :: get_S     => quark_get_S
        procedure :: get_A_B   => quark_get_A_B
        procedure :: get_sv_ss => quark_get_sv_ss

    end type

    public quark

contains
!_______________________________________________________________________
!
subroutine quark_create(self)

    ! passed variables
    class(quark)         :: self

    ! local variables

    if(self%is_allocated.eqv..true.)then
        write(*,'(A)')'# ... (warning) quark already allocated. Overwriting object.'
        call self%delete
    end if

    ! Allocate grid
    allocate(self%log_p_sq(1:self%n))
    allocate(self%p_sq    (1:self%n))
    allocate(self%A       (1:self%n))
    allocate(self%B       (1:self%n))
    allocate(self%A2      (1:self%n))
    allocate(self%B2      (1:self%n))

    ! Log linear grid.
    call gauss_legendre(log(self%eps_sq),log(self%lam_sq),self%log_p_sq)
    self%p_sq(:) = exp(self%log_p_sq)

    ! Initialise A and B functions
    self%A(:) = 1._r8
    self%B(:) = 1._r8

    call self%update() ! updates spline

    ! Mark the object as being allocated
    self%is_allocated = .true.

    write(*,'(A)')'# ... created quark object.'

end subroutine
!_______________________________________________________________________
!
subroutine quark_delete(self)

    ! passed variables
    class(quark)         :: self

    if(self%is_allocated.eqv..false.)then
        write(*,'(A)')'# ... (warning) quark already deallocated. '
        return
    end if

    ! Deallocate grid
    deallocate(self%log_p_sq)
    deallocate(self%p_sq    )
    deallocate(self%A       )
    deallocate(self%B       )
    deallocate(self%A2      )
    deallocate(self%B2      )

    ! Mark the object as being deallocated
    self%is_allocated = .false.

    write(*,'(A)')'# ... deleted quark object.'

end subroutine
!_______________________________________________________________________
!
subroutine quark_update(self)

    ! passed variables
    class(quark)         :: self

    call spline(self%log_p_sq, self%A, self%A2)
    call spline(self%log_p_sq, self%B, self%B2)

end subroutine
!_______________________________________________________________________
!
subroutine quark_get_A_B(self, p_sq, A, B)

    ! passed variables
    class(quark)          :: self
    real(r8), intent(in)  :: p_sq
    real(r8), intent(out) :: A, B

    A = splint(self%log_p_sq, self%A, self%A2, log(p_sq))
    B = splint(self%log_p_sq, self%B, self%B2, log(p_sq))

end subroutine
!_______________________________________________________________________
!
subroutine quark_get_sv_ss(self, p_sq, sv, ss)

    ! passed variables
    class(quark)          :: self
    real(r8), intent(in)  :: p_sq
    real(r8), intent(out) :: sv, ss

    real(r8) :: tmp

    call quark_get_A_B(self, p_sq, sv, ss)
    tmp = sv**2 * p_sq + ss**2
    sv = sv / tmp
    ss = ss / tmp

end subroutine
!_______________________________________________________________________
!
function quark_get_S(self, p) result(S)

    ! passed variables
    class(quark)         :: self
    real(r8), intent(in) :: p(4)

    ! returned variables
    complex(r8)          :: S(16)

    ! local variables
    real(r8) :: p_sq, sA, sB

    p_sq = sum(p*p)
    call quark_get_sv_ss(self, p_sq, sA, sB)

    S(:) = gV_(:,1)*(i_*p(1)*sA) + gV_(:,2)*(i_*p(2)*sA) &
         + gV_(:,3)*(i_*p(3)*sA) + gV_(:,4)*(i_*p(4)*sA) &
         + idV_(:) * sB

end function
!_______________________________________________________________________
!
end module
