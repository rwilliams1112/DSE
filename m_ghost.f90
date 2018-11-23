module m_ghost

    ! Global parameters
    use m_parameters

    use m_library

    implicit none
    private

    type ghost

        ! (de)allocated flag
        logical  :: is_allocated = .false.

        ! Renormalization constants
        real(r8) :: xm   ! subtraction point
        real(r8) :: xs   ! renormalisation point
        real(r8) :: Z_3t = 1._r8

        ! Grid
        integer  :: n
        real(r8) :: eps_sq
        real(r8) :: lam_sq

        ! Boundary condition for the ghost
        real(r8) :: one_over_Gxs

        ! Spline interpolator
        real(r8), allocatable :: log_p_sq(:)
        real(r8), allocatable ::     p_sq(:)
        real(r8), allocatable ::     G   (:)
        real(r8), allocatable ::     G2  (:)

        ! Good idea to have IR and UV fits
        real(r8) :: kappa, UV_A, UV_B, IR_A, IR_B

    contains

        ! Initialisation
        procedure :: create => ghost_create
        procedure :: delete => ghost_delete

        ! Access
        procedure :: update => ghost_update
        procedure :: get_G  => ghost_get_G
        procedure :: get_DG => ghost_get_DG
        procedure :: get_D  => ghost_get_D   ! return ghost propagator
        procedure :: set_G  => ghost_set_G

        ! Input/Output
        procedure :: output     => ghost_output
        procedure :: output_raw => ghost_output_raw
        procedure :: input_raw  => ghost_input_raw

    end type

    public ghost

contains
!_______________________________________________________________________
!
subroutine ghost_input_raw( self, filename )

    class(ghost)             :: self
    character(*), intent(in) :: filename

    integer :: unit_number

    open(newunit = unit_number,file = trim(filename),                   &
         status  = 'old',      form = 'unformatted')

    read(unit_number)self%n
    read(unit_number)self%eps_sq
    read(unit_number)self%lam_sq

    call self%create()

    read(unit_number)self%xm
    read(unit_number)self%xs
    read(unit_number)self%one_over_Gxs
    read(unit_number)self%Z_3t

    read(unit_number)self%log_p_sq(:)
    read(unit_number)self%p_sq    (:)
    read(unit_number)self%G       (:)
    read(unit_number)self%G2      (:)

    close(unit_number)

end subroutine
!_______________________________________________________________________
!
subroutine ghost_output_raw( self, filename )

    class(ghost)             :: self
    character(*), intent(in) :: filename

    integer :: unit_number

    open(newunit = unit_number,file = trim(filename),                   &
         status  = 'replace',  form = 'unformatted')

    write(unit_number)self%n
    write(unit_number)self%eps_sq
    write(unit_number)self%lam_sq

    write(unit_number)self%xm
    write(unit_number)self%xs
    write(unit_number)self%one_over_Gxs
    write(unit_number)self%Z_3t

    write(unit_number)self%log_p_sq(:)
    write(unit_number)self%p_sq    (:)
    write(unit_number)self%G       (:)
    write(unit_number)self%G2      (:)

    close(unit_number)

end subroutine
!_______________________________________________________________________
!
subroutine ghost_create(self)

    ! passed variables
    class(ghost)         :: self

    ! local variables
    integer :: i

    if(self%is_allocated.eqv..true.)then
        write(*,'(A)')'# ... (warning) ghost already allocated. Overwriting object.'
        call self%delete
    end if

    ! Allocate grid
    allocate(self%log_p_sq(1:self%n))
    allocate(self%p_sq    (1:self%n))
    allocate(self%G       (1:self%n))
    allocate(self%G2      (1:self%n))

    ! Log linear grid.
    do i = 1, self%n
        self%log_p_sq(i) = log(self%eps_sq) + log(self%lam_sq/self%eps_sq)*(i-1._r8)/(self%n-1._r8)
    end do

    self%p_sq(:) = exp(self%log_p_sq(:))

    ! Mark the object as being allocated
    self%is_allocated = .true.

    write(*,'(A)')'# ... created ghost object.'

end subroutine
!_______________________________________________________________________
!
subroutine ghost_delete(self)

    ! passed variables
    class(ghost)         :: self

    if(self%is_allocated.eqv..false.)then
        write(*,'(A)')'# ... (warning) ghost already deallocated.'
        return
    end if

    ! Deallocate grid
    deallocate(self%log_p_sq)
    deallocate(self%p_sq    )
    deallocate(self%G       )
    deallocate(self%G2      )

    ! Mark the object as being deallocated
    self%is_allocated = .false.

    write(*,'(A)')'# ... deleted ghost object.'

end subroutine
!_______________________________________________________________________
!
function ghost_get_G(self, p_sq) result (erg)

    ! passed variables
    class(ghost)         :: self
    real(r8), intent(in) :: p_sq

    ! returned variables
    real(r8) :: erg

    erg = splint(self%log_p_sq, self%G, self%G2, log(p_sq))

end function
!_______________________________________________________________________
!
function ghost_get_DG(self, p_sq) result(erg)

    ! passed variables
    class(ghost)         :: self
    real(r8), intent(in) :: p_sq

    ! returned variables
    real(r8)             :: erg

    erg = self%get_G(p_sq)/p_sq

end function
!_______________________________________________________________________
!
function ghost_get_D(self, p) result(erg)

    ! passed variables
    class(ghost)         :: self
    real(r8), intent(in) :: p(4)

    ! returned variables
    real(r8)             :: erg

    ! local variables
    real(r8) :: p_sq
    p_sq = sum(p*p)

    erg = -self%get_G(p_sq)/p_sq

end function
!_______________________________________________________________________
!
subroutine ghost_set_G(self, G, chisq_)

    ! passed variables
    class(ghost)         :: self
    real(r8), intent(in) :: G(:)
    real(r8), intent(out), optional :: chisq_

    ! local variables
    integer  :: i
    real(r8) :: chisq

    chisq = 0._r8
    do i = 1, self%n
        chisq = max(abs(self%G(i)/G(i) - 1._r8), chisq)
        self%G(i) = G(i)
    end do

    call self%update()

    if(present(chisq_))chisq_ = chisq

end subroutine
!_______________________________________________________________________
!
subroutine ghost_update(self)

    ! passed variables
    class(ghost)         :: self

    ! local variables
    real(r8), parameter  :: delta = - 9._r8/44._r8
    integer              :: n
    real(r8)             :: dy, x1, x2, y1, y2

    integer, parameter   :: ir_1 = 2, ir_2 = 4

    call spline(self%log_p_sq, self%G, self%G2)

    ! update the IR fits
    call fit(self%log_p_sq(ir_1:ir_2),log(self%G(ir_1:ir_2)),self%IR_A,self%IR_B)

    ! update UV extrapolation
    n  = self%n
    x1 = self%p_sq(n  )
    x2 = self%p_sq(n-1)
    y1 = self%G(n  )
    y2 = self%G(n-1)

    dy = (y1 - y2)/(x1 - x2)

    self%UV_A = y1*   ( delta*y1/(x1*dy))**(-delta)
    self%UV_B = x1*exp(-delta*y1/(x1*dy))

contains

    subroutine fit(x,y,a,b)

        real(r8), dimension(:), intent(in) :: x, y
        real(r8), intent(out)              :: a, b

        ! local
        integer  :: ndata
        real(r8) :: ss, sx, sxoss, sy, st2
        real(r8) :: t(1:size(x))

        ndata=size(x)
        ss=real(size(x))
        sx=sum(x)
        sy=sum(y)

        sxoss=sx/ss
        t(:)=x(:)-sxoss
        b=dot_product(t,y)

        st2=dot_product(t,t)
        b=b/st2
        a=(sy-sx*b)/ss

    end subroutine

end subroutine
!_______________________________________________________________________
!
subroutine ghost_output(self, filename)

    ! passed variables
    class(ghost)     :: self
    character(len=*) :: filename

    ! local variables
    integer  :: i, unit_number

    write(*,'(A)')'# ... writing ghost output to file: '//trim(adjustl(filename))

    open(newunit = unit_number, file = trim(adjustl(filename)) )

    write(unit_number,'(A)')'#               p_sq                   G'
    do i = 1, self%n
        write(unit_number,'(20es20.10e3)')self%p_sq(i), self%G(i)
    end do

    close(unit_number)

end subroutine
!_______________________________________________________________________
!
end module
