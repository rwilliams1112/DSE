module m_gluon

    ! Global parameters
    use m_parameters
    use m_library
    use m_dirac_matrices

    implicit none
    private

    type gluon

        ! (de)allocated flag
        logical  :: is_allocated = .false.

        ! Renormalization constants
        real(r8) :: xm   ! subtraction point
        real(r8) :: xs   ! renormalisation point
        real(r8) :: Z_3 = 1._r8

        ! Grid
        integer  :: n
        real(r8) :: eps_sq  ! infrared cutoff
        real(r8) :: lam_sq  ! ultraviolet cutoff

        ! Boundary condition for the gluon
        real(r8) :: Z_at_m
        real(r8) :: Z_at_s

        ! Good idea to have IR and UV fits
        real(r8)     :: kappa, UV_A, UV_B, IR_A, IR_B

        ! Spline interpolator
        real(r8), allocatable :: log_p_sq(:)
        real(r8), allocatable ::     p_sq(:)
        real(r8), allocatable ::     Z   (:)
        real(r8), allocatable ::     Z2  (:)

    contains

        ! Initialisation
        procedure :: create => gluon_create
        procedure :: delete => gluon_delete

        ! Access
        procedure :: update => gluon_update
        procedure :: get_Z  => gluon_get_Z   ! return scalar dressing
        procedure :: get_D  => gluon_get_D   ! return gluon propagator
        procedure :: set_Z  => gluon_set_Z

        ! return a fit function for the gluon
        procedure, nopass :: fit_Z  => gluon_fit_Z

        ! Input/Output
        procedure :: output     => gluon_output
        procedure :: output_raw => gluon_output_raw
        procedure :: input_raw  => gluon_input_raw

    end type

    public gluon

contains
!_______________________________________________________________________
!
subroutine gluon_input_raw( self, filename )

    class(gluon)             :: self
    character(*), intent(in) :: filename

    integer :: unit_number

    open(newunit = unit_number,file = trim(filename),                   &
         status  = 'old',      form = 'unformatted')

    read(unit_number)self%n
    read(unit_number)self%eps_sq
    read(unit_number)self%lam_sq

    call self%create()

    read(unit_number)self%xs
    read(unit_number)self%xm
    read(unit_number)self%Z_at_s
    read(unit_number)self%Z_at_m
    read(unit_number)self%Z_3

    read(unit_number)self%log_p_sq(:)
    read(unit_number)self%p_sq    (:)
    read(unit_number)self%Z       (:)
    read(unit_number)self%Z2      (:)

    close(unit_number)

end subroutine
!_______________________________________________________________________
!
subroutine gluon_output_raw( self, filename )

    class(gluon)             :: self
    character(*), intent(in) :: filename

    integer :: unit_number

    open(newunit = unit_number,file = trim(filename),                   &
         status  = 'replace',  form = 'unformatted')

    write(unit_number)self%n
    write(unit_number)self%eps_sq
    write(unit_number)self%lam_sq

    write(unit_number)self%xs
    write(unit_number)self%xm
    write(unit_number)self%Z_at_s
    write(unit_number)self%Z_at_m
    write(unit_number)self%Z_3

    write(unit_number)self%log_p_sq(:)
    write(unit_number)self%p_sq    (:)
    write(unit_number)self%Z       (:)
    write(unit_number)self%Z2      (:)

    close(unit_number)

end subroutine
!_______________________________________________________________________
!
subroutine gluon_create(self)

    ! passed variables
    class(gluon)         :: self

    ! local variables
    integer :: i

    if(self%is_allocated.eqv..true.)then
        write(*,'(A)')'# ... (warning) gluon already allocated. Overwriting object.'
        call self%delete
    end if

    ! Allocate grid
    allocate(self%log_p_sq(1:self%n))
    allocate(self%p_sq    (1:self%n))
    allocate(self%Z       (1:self%n))
    allocate(self%Z2      (1:self%n))

    ! Log linear grid.
    do i = 1, self%n
        self%log_p_sq(i) = log(self%eps_sq) + log(self%lam_sq/self%eps_sq)*(i-1._r8)/(self%n-1._r8)
    end do

    self%p_sq(:) = exp(self%log_p_sq(:))

    ! Mark the object as being allocated
    self%is_allocated = .true.

    write(*,'(A)')'# ... created gluon object.'

end subroutine
!_______________________________________________________________________
!
subroutine gluon_delete(self)

    ! passed variables
    class(gluon)         :: self

    if(self%is_allocated.eqv..false.)then
        write(*,'(A)')'# ... (warning) gluon already deallocated. '
        return
    end if

    ! Deallocate grid
    deallocate(self%log_p_sq)
    deallocate(self%p_sq    )
    deallocate(self%Z       )
    deallocate(self%Z2      )

    ! Mark the object as being deallocated
    self%is_allocated = .false.

    write(*,'(A)')'# ... deleted gluon object.'

end subroutine
!_______________________________________________________________________
!
function gluon_fit_Z(p_sq) result (erg)

    ! passed variables
    !class(gluon)         :: self
    real(r8), intent(in) :: p_sq

    ! returned variables
    real(r8) :: erg

    ! local vairables
    real(r8), parameter :: a      =  0.595_r8
    real(r8), parameter :: b      =  1.355_r8
    real(r8), parameter :: c      = 11.5_r8
    real(r8), parameter :: a_mu   = 0.3_r8
    real(r8), parameter :: gamma  = -13._r8/22._r8
    real(r8), parameter :: beta_0 = 11._r8 !- 2._r8 * nF / 3._r8

    real(r8), parameter :: PI =                                &
            3.1415926535897932384626433832795028841971693993751_r8
    real(r8) :: x

    x = p_sq / 1.96_r8

    erg = x/(x + 1._r8)**2 * ( (c/(x+a))**b + x * (a_mu*beta_0 /4._r8/PI * log(1._r8+x))**gamma  )

end function
!_______________________________________________________________________
!
function gluon_get_Z(self, p_sq) result (erg)

    ! passed variables
    class(gluon)         :: self
    real(r8), intent(in) :: p_sq

    ! returned variables
    real(r8) :: erg

    erg = splint(self%log_p_sq, self%Z, self%Z2, log(p_sq))

    ! Pauli-Villars
    erg = erg! / (1._r8 + p_sq / 1e3_r8)

end function
!_______________________________________________________________________
!
function gluon_get_D(self, p) result (D)

    ! passed variables
    class(gluon)         :: self
    real(r8), intent(in) :: p(4)

    ! returned variables
    real(r8) :: D(4,4)

    integer  :: i_nu
    real(r8) :: DZ, p_sq

    p_sq = sum(p*p)
    DZ   = self%get_Z(p_sq) / p_sq

    do i_nu = 1, 4
        D(:,i_nu) = (d_(:,i_nu) - p(:)*p(i_nu) / p_sq) * DZ
    end do

end function
!_______________________________________________________________________
!
subroutine gluon_set_Z(self, Z)

    ! passed variables
    class(gluon)         :: self
    real(r8), intent(in) :: Z(:)

    self%Z(:) = Z(:)

    call self%update()

end subroutine
!_______________________________________________________________________
!
subroutine gluon_update(self)

    ! passed variables
    class(gluon)         :: self

    ! local variables
    real(r8), parameter  :: gamma = -13._r8/22._r8
    integer              :: n
    real(r8)             :: dy, x1, x2, y1, y2

    integer, parameter   :: ir_1 = 2, ir_2 = 4

    call spline(self%log_p_sq, self%Z, self%Z2)

    ! update the IR fits    ! must be log(x) and log(Z), as stored in splines
    call fit(self%log_p_sq(ir_1:ir_2),log(self%Z(ir_1:ir_2)),self%IR_A,self%IR_B)

    ! update UV extrapolation
    n  = self%n
    x1 = self%p_sq(n  )
    x2 = self%p_sq(n-1)
    y1 = self%Z(n  )
    y2 = self%Z(n-1)

    dy = (y1 - y2)/(x1 - x2)

    self%UV_A = y1*   ( gamma*y1/(x1*dy))**(-gamma)
    self%UV_B = x1*exp(-gamma*y1/(x1*dy))

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
subroutine gluon_output(self, filename)

    ! passed variables
    class(gluon)     :: self
    character(len=*) :: filename

    ! local variables
    integer  :: i, unit_number

    write(*,'(A)')'# ... writing gluon output to file: '//trim(adjustl(filename))

    open(newunit = unit_number, file = trim(adjustl(filename)) )

    write(unit_number,'(A)')'#               p_sq                   Z                   D'
    do i = 1, self%n
        write(unit_number,'(20es20.10e3)')self%p_sq(i), self%Z(i), self%Z(i)/self%p_sq(i)
    end do

    close(unit_number)

end subroutine
!_______________________________________________________________________
!
end module
