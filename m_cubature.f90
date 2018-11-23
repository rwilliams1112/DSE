module m_hypercube

    ! Global parameters
    use m_parameters

    implicit none
    private

    type hypercube

        integer               :: n_dim
        real(r8), allocatable :: hwidth(:)
        real(r8), allocatable :: center(:)
        real(r8)              :: volume

    contains

        procedure :: cut_left  => hypercube_cut_left
        procedure :: cut_right => hypercube_cut_right

        procedure :: delete => hypercube_delete
        ! copy/assign
        !procedure :: hypercube_assign
        !generic   :: assignment(=) => hypercube_assign

    end type

    ! Constructors for the hypercube
    interface hypercube
        module procedure :: constructor_n
        module procedure :: constructor_nab
        module procedure :: constructor_hsplit
    end interface

    public :: hypercube

contains

subroutine hypercube_delete(this)
    class(hypercube)     :: this

    deallocate(this%hwidth)
    deallocate(this%center)

end subroutine
!_______________________________________________________________________
!
function constructor_n(n_dim) result(this)

    type(hypercube)     :: this
    integer, intent(in) :: n_dim

    ! default [0,1]^d
    this%n_dim = n_dim
    allocate(this%hwidth(1:n_dim), source = 0.5_r8)
    allocate(this%center(1:n_dim), source = 0.0_r8)
    call hypercube_calculate_volume(this)

end function
!_______________________________________________________________________
!
function constructor_nab(n_dim, a, b) result(this)

    type(hypercube)      :: this
    integer, intent(in)  :: n_dim
    real(r8), intent(in) :: a(1:n_dim),b(1:n_dim)

    ! scaled [a,b]^d
    this%n_dim = n_dim
    allocate(this%hwidth(1:n_dim), source = (b(1:n_dim)-a(1:n_dim))/2._r8)
    allocate(this%center(1:n_dim), source = (b(1:n_dim)+a(1:n_dim))/2._r8)
    call hypercube_calculate_volume(this)

end function
!_______________________________________________________________________
!
function constructor_hsplit(h, split) result(this)

    type(hypercube)                 :: this
    type(hypercube), intent(in out) :: h
    integer, intent(in)             :: split

    this = h                   ! copy the passed hypercube
    call h%cut_left(split)     ! cut passed hypercube in left direction of split
    call this%cut_right(split) ! cut new hypercube in right direction of split
end function
!_______________________________________________________________________
!
subroutine hypercube_assign(to, from)

    class(hypercube), intent(out) :: to
    type (hypercube), intent(in)  :: from

    to%n_dim  = from%n_dim
    to%volume = from%volume
    if(allocated(from%hwidth))then
        allocate(to%hwidth(1:ubound(from%hwidth, 1)), source=from%hwidth)
    else
        allocate(to%hwidth(1:from%n_dim))
    end if
    if(allocated(from%center))then
        allocate(to%center(1:ubound(from%center, 1)), source=from%center)
    else
        allocate(to%center(1:from%n_dim))
    end if
end subroutine
!_______________________________________________________________________
!
subroutine hypercube_calculate_volume(this)

    type(hypercube) :: this

    this%volume = product(2._r8*this%hwidth(:))
    return

end subroutine
!_______________________________________________________________________
!
subroutine hypercube_cut_left(this, split)

    class(hypercube) :: this
    integer, intent(in) :: split

    this%hwidth(split) = this%hwidth(split) / 2._r8
    this%center(split) = this%center(split) - this%hwidth(split)
    this%volume = this%volume / 2._r8
    return

end subroutine
!_______________________________________________________________________
!
subroutine hypercube_cut_right(this, split)

    class(hypercube) :: this
    integer, intent(in) :: split

    this%hwidth(split) = this%hwidth(split) / 2._r8
    this%center(split) = this%center(split) + this%hwidth(split)
    this%volume = this%volume / 2._r8
    return

end subroutine
!_______________________________________________________________________
!
subroutine hypercube_set(this, i, a, b)

    type(hypercube)      :: this
    integer, intent(in)  :: i
    real(r8), intent(in) :: a, b

    this%center(i) = (b+a)/2._r8
    this%hwidth(i) = (b-a)/2._r8
    return

end subroutine
!_______________________________________________________________________
!
subroutine hypercube_move(this, i, distance)

    type(hypercube)      :: this
    integer, intent(in)  :: i
    real(r8), intent(in) :: distance

    this%center(i) = this%center(i) + distance
    return

end subroutine
!_______________________________________________________________________
!
end module

module m_cubature_buffer

    implicit none
    private

    integer, parameter :: r8 = kind(1.d0)

    type cubature_buffer

        integer :: n_dim
        integer :: n_fun
        integer :: n_max
        integer :: n_pts

        real(r8), allocatable :: x(:,:)
        real(r8), allocatable :: f(:,:)

    end type

    public :: cubature_buffer

end module

module m_heap

    use m_hypercube

    implicit none
    private

    integer, parameter :: r8 = kind(1.d0)

    type element_t

        type(hypercube), allocatable :: h
        integer                      :: n_fun
        integer                      :: splitDim
        real(r8), allocatable        :: result(:)
        real(r8), allocatable        :: error (:)
        real(r8)                     :: max_error
    contains
        ! copy/assign
        procedure :: element_t_assign
        generic   :: assignment(=) => element_t_assign
    end type

    type heap

        integer :: cur_size
        integer :: max_size
        integer :: n_fun

        real(r8), allocatable :: result(:)
        real(r8), allocatable :: error (:)

        type(element_t), allocatable :: leaf(:)

    contains

        ! copy/assign
        !procedure :: heap_assign
        !generic   :: assignment(=) => heap_assign

        procedure :: insert => heap_insert
        procedure :: remove => heap_remove

        procedure :: grow   => heap_grow

    end type
    ! Constructor for type tuple
    interface element_t
        module procedure :: element_t_constructor_n
        module procedure :: element_t_constructor_split
    end interface

    ! Constructor for type tuple
    interface heap
        module procedure :: heap_constructor_nn
        module procedure :: heap_constructor_n0
    end interface

    public :: element_t
    public :: heap

contains
!_______________________________________________________________________
!
subroutine element_t_assign(to, from)

    class(element_t), intent(out) :: to
    type (element_t), intent(in)  :: from
    type (hypercube)              :: aux


    to%splitDim  = from%splitDim
    to%max_error = from%max_error
    to%n_fun     = from%n_fun
    !allocate(to%result(ubound(from%result, 1)), source=from%result)
    !allocate(to%error (ubound(from%error,  1)), source=from%error )
    allocate(to%result, source=from%result)
    allocate(to%error , source=from%error )

    ! Helios fix
    aux = from%h
    allocate(to%h, source= aux)  ! broken allocate(to%h, source= from%h)

end subroutine
!_______________________________________________________________________
!
! The element_t constructor
!_______________________________________________________________________
!
function element_t_constructor_n(n_fun, h) result(this)

    type(element_t)     :: this
    integer, intent(in) :: n_fun
    type(hypercube)     :: h

    allocate(this%result(1:n_fun))
    allocate(this%error (1:n_fun))
    allocate(this%h, source= h)

    this%n_fun = n_fun
    this%max_error = 0._r8
    this%splitDim  = 0

end function
!_______________________________________________________________________
!
function element_t_constructor_split(elem) result(this)

    type(element_t), intent(in out) :: elem
    type(element_t)                 :: this

    allocate(this%h); this%h = hypercube(elem%h,elem%splitDim)
    this%splitDim  = 0
    this%max_error = 0._r8
    elem%splitDim  = 0
    elem%max_error = 0._r8

    this%n_fun     = elem%n_fun

    elem%result = 0._r8
    elem%error  = 0._r8
    allocate(this%result(ubound(elem%result, 1)), source=elem%result)
    allocate(this%error (ubound(elem%error,  1)), source=elem%error )

end function
!_______________________________________________________________________
!
! The heap constructor
!_______________________________________________________________________
!
function heap_constructor_n0(n_fun) result(this)

    type(heap)          :: this
    integer, intent(in) :: n_fun

    this%n_fun    = n_fun
    this%max_size = 1
    this%cur_size = 0           ! heap is empty

    allocate(this%leaf(1:this%max_size ))
    allocate(this%result(1:n_fun ), source = 0._r8)
    allocate(this%error (1:n_fun ), source = 0._r8)

end function
!_______________________________________________________________________
!
! The heap constructor
!_______________________________________________________________________
!
function heap_constructor_nn(n_fun, max_size) result(this)

    type(heap)          :: this
    integer, intent(in) :: n_fun
    integer, intent(in) :: max_size

    this%n_fun    = n_fun
    this%cur_size = 0           ! heap is empty

    if(max_size.le.0)then
        this%max_size = 1
        allocate(this%leaf(1:max_size))
    else
        this%max_size = max_size
        allocate(this%leaf(1:max_size))
    end if
    allocate(this%result(1:n_fun ), source = 0._r8)
    allocate(this%error (1:n_fun ), source = 0._r8)

end function
!_______________________________________________________________________
!
! The heap constructor
!_______________________________________________________________________
!
subroutine heap_assign(to, from)

    class(heap), intent(out) :: to
    type (heap), intent(in)  :: from

    to%n_fun    = from%n_fun
    to%max_size = from%max_size
    to%cur_size = from%cur_size
    allocate(to%leaf(ubound(from%leaf, 1)), source=from%leaf)
    allocate(to%result(ubound(from%result, 1)), source=from%result)
    allocate(to%error (ubound(from%error , 1)), source=from%error )

end subroutine
!_______________________________________________________________________
!
subroutine heap_insert(this, leaf)

    class(heap), intent(in out)   :: this
    class(element_t), intent(in) :: leaf

    ! local variables
    integer :: i_current, i_parent
    class(element_t), allocatable :: temp

    allocate(temp)

    ! Store the result and error in the tree
    this%result(:) = this%result(:) + leaf%result(:)
    this%error (:) = this%error (:) + leaf%error (:)

    ! Increment the number of leaves and track the current location in i_current
    this%cur_size = this%cur_size + 1
    i_current     = this%cur_size

    ! If we have reached the maximum, grow the heap
    if(i_current.ge.this%max_size)call this%grow()

    ! Bubble the new value up the tree so that it is partially sorted.
    do while (i_current.gt.1)
        i_parent = floor(i_current/2._r8)

        ! Move bigger errors to the top
        if(leaf%max_error .lt. this%leaf(i_parent)%max_error )exit
        this%leaf(i_current) = this%leaf(i_parent)
        i_current = i_parent

    end do
    this%leaf(i_current) = leaf

end subroutine
!_______________________________________________________________________
!
subroutine heap_remove(this, leaf)

    ! we assiume that an empty heap will NOT be called
    class(heap), intent(in out)   :: this

    type(element_t) :: leaf
    type(element_t) :: temp

    ! local variables
    integer :: i_current, i_child, i_largest

    !allocate(temp)

    ! We want to return the top node
    leaf = this%leaf(1)

    ! move last leaf to place of first leaf
    if(this%cur_size.gt.1)this%leaf(1)  = this%leaf(this%cur_size)

    this%result(:) = this%result(:) - leaf%result(:)
    this%error (:) = this%error (:) - leaf%error (:)

    ! reduce count by one
    this%cur_size = this%cur_size - 1

    ! We bubble the leaf DOWN the heap
    i_current = 1
    do while (i_current.lt.this%cur_size)

        i_child = 2*i_current

        if(i_child.gt.this%cur_size)exit

        ! check if error of child is greater
        if(this%leaf(i_child)%max_error .lt. this%leaf(i_current)%max_error)then
            ! child has a great error. We want to swap
            i_largest = i_current
        else
            i_largest = i_child
        end if


        ! let's check the next child
        i_child = i_child + 1
        if(i_child.le.this%cur_size)then ! this child exists
            if(this%leaf(i_largest)%max_error .lt. this%leaf(i_child)%max_error)then
                ! child has a greater error. We want to swap
                i_largest = i_child
                ! we swap new largest with current
            end if
        end if

        if(i_largest.ne.i_current)then
            temp = this%leaf(i_largest)
            this%leaf(i_largest) = this%leaf(i_current)
            this%leaf(i_current) = temp
            ! Pick child for next comparison and repeat.
            i_current = i_largest
            cycle
        end if

        exit

    end do

end subroutine
!_______________________________________________________________________
!
subroutine heap_grow(this)

    class(heap), intent(in out)  :: this
    type(element_t), allocatable :: temp(:)

    allocate(temp(1:2*this%max_size+1))
    temp(1:this%max_size) = this%leaf(1:this%max_size)

    this%max_size = 2*this%max_size+1
    call move_alloc(from=temp,to=this%leaf)

end subroutine
!_______________________________________________________________________
!
end module
!_______________________________________________________________________
!
module m_rule

    use m_hypercube

    implicit none
    private

    integer, parameter :: r8 = kind(1.d0)

    type rule

        integer :: n_dim
        integer :: n_pts

        real(r8), allocatable :: x  (:,:)
        real(r8), allocatable :: wgk(:)
        real(r8), allocatable :: wg (:)

        ! dimension independent constants
        real(r8) :: weight2, weight4, weightE2, weightE4, ratio

        ! dimension dependent constants
        real(r8) :: weight1, weight3, weight5, weightE1, weightE3

    contains

        procedure :: apply => rule_apply
        procedure :: delete => rule_delete

        ! type-bound copy/assign
        !procedure :: rule_assign
        !generic   :: assignment(=) => rule_assign

    end type

    ! Constructor for type rule
    interface rule
        module procedure :: constructor_n
    end interface

    public :: rule

contains

subroutine rule_delete(this)

    class(rule) :: this

    deallocate(this%x)
    deallocate(this%wgk)
    deallocate(this%wg)

end subroutine
!_______________________________________________________________________
!
! The rule construtor
!_______________________________________________________________________
!
function constructor_n(n_dim) result (this)

    type(rule) :: this
    integer, intent(in)      :: n_dim

    if(n_dim.lt.1)then
        print '(A)', '# (error) invalid number of dimensions passed'
        stop
    end if

    select case(n_dim)
    case(1)
        this%n_dim = n_dim
        this%n_pts = 15

        allocate(this%x  (1:this%n_dim,this%n_pts))
        allocate(this%wg (             this%n_pts))
        allocate(this%wgk(             this%n_pts))
        call rule15GaussKronrod_create(this)

    case default
        this%n_dim = n_dim
        this%n_pts = 1 + 2*(2*n_dim) + 2*n_dim* (n_dim-1) + 2**n_dim

        ! set dimension independent constants
        this%weight2  = 980._r8 /  6561._r8
        this%weight4  = 200._r8 / 19683._r8
        this%weightE2 = 245._r8 /   486._r8
        this%weightE4 =  25._r8 /   729._r8
        this%ratio    = 1._r8/7._r8

        ! set dimension dependent constants
        this%weight1  = (12824._r8-9120._r8*n_dim+400._r8*n_dim**2)/19683._r8
        this%weight3  = (1820._r8-400._r8*n_dim)/19683._r8
        this%weight5  = (6859._r8/19683._r8)/(2._r8**n_dim)
        this%weightE1 = (729._r8-950._r8*n_dim+50._r8*n_dim**2)/729._r8
        this%weightE3 = (265._r8-100._r8*n_dim)/1458._r8

        allocate(this%x  (1:this%n_dim,this%n_pts))
        allocate(this%wg (             this%n_pts))
        allocate(this%wgk(             this%n_pts))
        call rule75GenzMalik_create(this)

    end select

end function
!_______________________________________________________________________
!
! The rule copy/assign operator
!_______________________________________________________________________
!
subroutine rule_assign(to, from)

    class(rule), intent(out) :: to
    type (rule), intent(in)  :: from

    to%n_dim = from%n_dim
    to%n_pts = from%n_pts

    ! assign and copy points array
    if(allocated(from%x))allocate(to%x(lbound(from%x,1):ubound(from%x,1), &
                                       lbound(from%x,2):ubound(from%x,2)),&
                       source = from%x(lbound(from%x,1):ubound(from%x,1), &
                                       lbound(from%x,2):ubound(from%x,2)))
    if(allocated(from%wg ))allocate(to%wg (lbound(from%wg, 1):ubound(from%wg, 1)), &
                         source = from%wg (lbound(from%wg, 1):ubound(from%wg, 1)))
    if(allocated(from%wgk))allocate(to%wgk(lbound(from%wgk,1):ubound(from%wgk,1)), &
                         source = from%wgk(lbound(from%wgk,1):ubound(from%wgk,1)))

    ! copy dimension dependent weights
    to%weight1  = from%weight1
    to%weight2  = from%weight2
    to%weight3  = from%weight3
    to%weight4  = from%weight4
    to%weight5  = from%weight5
    to%weightE1 = from%weightE1
    to%weightE3 = from%weightE3
    to%weightE2 = from%weightE2
    to%weightE4 = from%weightE4
    to%ratio    = from%ratio

end subroutine
!_______________________________________________________________________
!
! The rule15GaussKronrod creator for unit width
!_______________________________________________________________________
!
subroutine rule15GaussKronrod_create(this)

    type(rule) :: this

    real(r8), parameter :: xgk(1:8) =          &
     [ 0.991455371120812639206854697526329_r8, &
       0.949107912342758524526189684047851_r8, &
       0.864864423359769072789712788640926_r8, &
       0.741531185599394439863864773280788_r8, &
       0.586087235467691130294144838258730_r8, &
       0.405845151377397166906606412076961_r8, &
       0.207784955007898467600689403773245_r8, &
       0.000000000000000000000000000000000_r8 ]

    real(r8), parameter :: wg(1:4) =           &
     [ 0.129484966168869693270611432679082_r8, &
       0.279705391489276667901467771423780_r8, &
       0.381830050505118944950369775488975_r8, &
       0.417959183673469387755102040816327_r8 ]

    real(r8), parameter :: wgk(1:8) =           &
      [ 0.022935322010529224963732008058970_r8, &
        0.063092092629978553290700663189204_r8, &
        0.104790010322250183839876322541518_r8, &
        0.140653259715525918745189590510238_r8, &
        0.169004726639267902826583426598550_r8, &
        0.190350578064785409913256402421014_r8, &
        0.204432940075298892414161999234649_r8, &
        0.209482141084727828012999174891714_r8 ]

    integer  :: i

    this%wg (:) = 0._r8
    ! unit interval of (-0.5,0.5)
    this%x(:,1) = xgk(8) / 2._r8
    this%wgk(1) = wgk(8) / 2._r8
    this%wg (1) = wg (4) / 2._r8

    ! 15 kronrod weights
    do i = 1, 7
        this%x  (:,2*i)   = - xgk(i) / 2._r8
        this%x  (:,2*i+1) = + xgk(i) / 2._r8
        this%wgk(  2*i)   =   wgk(i) / 2._r8
        this%wgk(  2*i+1) =   wgk(i) / 2._r8
    end do

    ! 7 gauss weights
    do i = 1, 3
        this%wg(4*i)   = wg(i) / 2._r8
        this%wg(4*i+1) = wg(i) / 2._r8
    end do

end subroutine
!_______________________________________________________________________
!
! Create a rule for a unit hypercube, center = 0, hwidth = 1/2
!_______________________________________________________________________
!
subroutine rule75GenzMalik_create(this)

    type(rule) :: this


    ! dimension independent constants
    real(r8), parameter :: lambda2  = sqrt(9._r8/70._r8)
    real(r8), parameter :: lambda4  = sqrt(9._r8/10._r8)
    real(r8), parameter :: lambda5  = sqrt(9._r8/19._r8)

    integer  :: i_count, j
    integer  :: i, signs, d, mask ! grey ordering / hamilton cycle
    real(r8) :: x(1:this%n_dim)

    i_count = 0
    x(:)    = 0._r8
    i_count = i_count + 1; this%x(:,i_count) = x(:)

    do i = 1, this%n_dim
        x(i) = - lambda2 / 2._r8
        i_count = i_count + 1; this%x(:,i_count) = x(:)
        x(i) = + lambda2 / 2._r8
        i_count = i_count + 1; this%x(:,i_count) = x(:)
        x(i) = - lambda4 / 2._r8
        i_count = i_count + 1; this%x(:,i_count) = x(:)
        x(i) = + lambda4 / 2._r8
        i_count = i_count + 1; this%x(:,i_count) = x(:)
        x(i) = 0._r8
    end do

    do i = 1, this%n_dim-1
        x(i) = 0._r8 - lambda4 / 2._r8
        do j = i+1, this%n_dim
            x(j) = - lambda4 / 2._r8
            i_count = i_count + 1; this%x(:,i_count) = x(:)
            x(i) = + lambda4 / 2._r8
            i_count = i_count + 1; this%x(:,i_count) = x(:)
            x(j) = + lambda4 / 2._r8
            i_count = i_count + 1; this%x(:,i_count) = x(:)
            x(i) = - lambda4 / 2._r8
            i_count = i_count + 1; this%x(:,i_count) = x(:)
            x(j) = 0._r8
        end do
        x(i) = 0._r8
    end do

    x(:) = lambda5 / 2._r8
    i_count = i_count + 1; this%x(:,i_count) = x(:)

    i = 0 ! Hamilton cycle
    signs = 0
    do while(.true.)
        d = trailz(not(i))+1
        if(d.gt.this%n_dim)exit
        mask  = ishft(1,d-1)
        signs = ieor(signs,mask)
        if(iand(signs,mask).ne.0)then
            x(d) = - lambda5 / 2._r8
        else
            x(d) = + lambda5 / 2._r8
        end if
        i_count = i_count + 1; this%x(:,i_count) = x(:)
        i = i + 1
    end do

end subroutine
!_______________________________________________________________________
!
function rule_apply(this, h, n_fun, ff, result, error) result(dimDiffMax)

    class(rule)     :: this
    type(hypercube) :: h
    integer, intent(in)  :: n_fun
    real(r8), intent(in) :: ff(1:n_fun, 1:this%n_pts)
    integer              :: dimDiffMax

    real(r8), intent(in out):: result(1:n_fun), error(1:n_fun)

    select case(this%n_dim)
    case(1)
        call rule_apply_15GaussKronrod(this, n_fun, ff, h%volume, dimDiffMax, result, error)
    case default
        call rule_apply_75GenzMalik(this, n_fun, ff, h%volume, dimDiffMax, result, error)
    end select

end function
!_______________________________________________________________________
!
subroutine rule_apply_15GaussKronrod(this, n_fun, ff, volume, dimDiffMax, result, err)

    class(rule)              :: this
    integer, intent(in)      :: n_fun
    real(r8), intent(in)     :: ff(1:n_fun,1:this%n_pts), volume
    real(r8), intent(in out) :: result(1:n_fun), err(1:n_fun)
    integer, intent(in out)  :: dimDiffMax

    real(r8), dimension(1:n_fun) :: result_gauss, result_kronrod, result_abs, result_asc, mean
    integer :: i

    dimDiffMax = 1 ! only one dimension to split

    do i = 1, n_fun
        result_kronrod(i) = sum(    ff(i,:) *this%wgk(:))
        result_abs    (i) = sum(abs(ff(i,:))*this%wgk(:))
        result_gauss  (i) = sum(    ff(i,:) *this%wg (:))
    end do

    result(:) = volume * result_kronrod(:)

    mean(:) = result_kronrod(:) / 2._r8
    do i = 1, n_fun
        result_asc(i)     = sum(abs(ff(i,:)-mean(i))*this%wgk(:))
    end do
    err(1:n_fun) = abs(result_kronrod(1:n_fun) - result_gauss(1:n_fun)) * volume
    result_abs(1:n_fun) = result_abs(1:n_fun) * volume
    result_asc(1:n_fun) = result_asc(1:n_fun) * volume

    do i = 1, n_fun
        if((result_asc(i).ne.0._r8).and.(err(i).ne.0._r8))then
            err(i) = err(i) * min((200._r8*err(i) / result_asc(i))**1.5_r8,1._r8)
        end if
        if(result_abs(i).gt. tiny(1._r8) / (50 * epsilon(1._r8)))then
            err(i) = max(50._r8 * epsilon(1._r8) * result_abs(i), err(i))
        end if
    end do

end subroutine
!_______________________________________________________________________
!
subroutine rule_apply_75GenzMalik(this, n_fun, ff, volume, dimDiffMax, result, res5th)

    class(rule)              :: this
    integer, intent(in)      :: n_fun
    real(r8), intent(in)     :: ff(1:n_fun,1:this%n_pts), volume
    real(r8), intent(in out) :: result(1:n_fun), res5th(1:n_fun)
    integer, intent(in out)  :: dimDiffMax

    integer :: i_count, i
    real(r8), dimension(1:n_fun) :: sum0, sum1, sum2, sum3, sum4, sum5, v0, v1, v2, v3
    real(r8) :: maxdiff, frthdf(1:n_fun), dfs(1:this%n_dim)

    ! accumulate the integral
    i_count = 1
    ! (0, 0, 0, ...) generators
    sum0(1:n_fun) = ff(1:n_fun,i_count) ! center point
    sum1(1:n_fun) = 0._r8
    sum2(1:n_fun) = 0._r8
    sum3(1:n_fun) = 0._r8
    sum4(1:n_fun) = 0._r8
    sum5(1:n_fun) = 0._r8

    ! (lambda2, 0, 0, ...), (lambda3=lambda4, 0, 0, ...) generators
    do i = 1, this%n_dim
        i_count = i_count + 1; v0(1:n_fun) = ff(1:n_fun,i_count)
        i_count = i_count + 1; v1(1:n_fun) = ff(1:n_fun,i_count)
        i_count = i_count + 1; v2(1:n_fun) = ff(1:n_fun,i_count)
        i_count = i_count + 1; v3(1:n_fun) = ff(1:n_fun,i_count)
        sum2(1:n_fun) = sum2(1:n_fun) + v0(1:n_fun) + v1(1:n_fun)
        sum3(1:n_fun) = sum3(1:n_fun) + v2(1:n_fun) + v3(1:n_fun)

        frthdf(1:n_fun) = abs( 2*sum0(1:n_fun) - (v2(1:n_fun) + v3(1:n_fun)) &
                        + this%ratio*(v0(1:n_fun) + v1(1:n_fun) - 2._r8*sum0(1:n_fun)))/4._r8

        dfs(i) = SUM( frthdf(:), mask = ABS(sum0(:)) + frthdf(:) .gt. abs(sum0(:)) )

    end do

    ! (lambda4, lambda4, 0, ...)
    do i = 1, 2*(this%n_dim)*(this%n_dim-1)
        i_count = i_count + 1; sum4(1:n_fun) = sum4(1:n_fun) + ff(1:n_fun,i_count)
    end do

    ! (lamda5, lamda5, ..., lamda5)
    do i = 1, 2**this%n_dim
        i_count = i_count + 1; sum5(1:n_fun) = sum5(1:n_fun) + ff(1:n_fun,i_count)
    end do

    ! Calculate the result and the error
    result(1:n_fun) = volume * (this%weight1 *sum0(1:n_fun) + this%weight2 *sum2(1:n_fun) &
                              + this%weight3 *sum3(1:n_fun) + this%weight4 *sum4(1:n_fun) &
                              + this%weight5 *sum5(1:n_fun))
    res5th(1:n_fun) = volume * (this%weightE1*sum0(1:n_fun) + this%weightE2*sum2(1:n_fun) &
                              + this%weightE3*sum3(1:n_fun) + this%weightE4*sum4(1:n_fun))
    res5th(1:n_fun) = abs(res5th(1:n_fun)-result(1:n_fun))

    ! Now decide which dimension to split
    dimDiffMax = 1
    maxdiff    = 0._r8
    do i = 1, this%n_dim
        if(dfs(i).gt.maxdiff)then
            maxdiff = dfs(i)
            dimDiffMax = i
        end if
    end do

end subroutine
!_______________________________________________________________________
!
end module
!_______________________________________________________________________
!
module m_cubature

    use m_heap
    use m_hypercube
    use m_rule

    implicit none
    private

    integer, parameter :: r8 = kind(1.d0)

    type cubature

        integer :: n_dim
        integer :: n_fun
        integer :: n_eval
        integer :: n_regions

        integer :: min_eval
        integer :: max_eval

        integer :: min_regions
        integer :: max_regions

        real(r8) :: AbsError
        real(r8) :: RelError

        type(heap), allocatable :: tree
        type(rule), allocatable :: cubature_rule

    contains

        procedure :: adapt => adaptive_integrate

    end type

    public :: cubature

contains
!_______________________________________________________________________
!
subroutine adaptive_integrate(this, n_dim, x1, x2, n_fun, func, res, err)

    class(cubature) :: this

    integer, intent(in)   :: n_dim, n_fun
    real(r8), intent(in)  :: x1(1:n_dim), x2(1:n_dim)
    real(r8), intent(out) :: res(1:n_fun), err(1:n_fun)

    ! The function to be integrated
    interface
        function func(n_fun, x)
            import
            integer, intent(in)        :: n_fun
            real(r8), intent(in)       :: x(:)
            real(r8), dimension(n_fun) :: func
        end function
    end interface

    type(hypercube), allocatable :: h
    type(element_t), allocatable :: elem, elem1, elem2

    this%n_eval = 0

    ! Create the cubature rule
    allocate(this%cubature_rule); this%cubature_rule = rule(n_dim)
    !allocate(this%cubature_rule, source = rule(n_dim))

    ! Create the heap of default size 1. It will grow logarithmically.
    !allocate(this%tree); this%tree = heap(n_fun, 1)
    allocate(this%tree, source =  heap(n_fun, 1))
    ! Create the initial region
    !allocate(h);    h    = hypercube(n_dim, x1, x2)
    allocate(h,    source = hypercube(n_dim, x1, x2))
    allocate(elem, source = element_t(n_fun, h))
    call eval_region(elem, func, this%cubature_rule)
    this%n_eval = this%n_eval + this%cubature_rule%n_pts

    ! Add this to the heap
    call this%tree%insert(elem)

    allocate(elem1)
    do while((this%n_eval .lt. this%max_eval))

        if(converged(this%tree, this%AbsError, this%RelError)&
           .and.(this%n_eval .ge.this%min_eval))exit

        ! pop the element from the tree
        call this%tree%remove(elem1)

        ! split elem1 into elem1 | elem2
        !allocate(elem2); elem2 = element_t(elem1)
        allocate(elem2, source = element_t(elem1))

        call eval_region(elem1, func, this%cubature_rule)
        call eval_region(elem2, func, this%cubature_rule)
        this%n_eval = this%n_eval + 2*this%cubature_rule%n_pts

        call this%tree%insert(elem1)
        call this%tree%insert(elem2)

        ! traverse the tree
        deallocate(elem2)

    end do
    !print *, '# neval', this%n_eval
    !print *, '# ', this%tree%result(:), this%tree%error(:)

    res(1:n_fun) = this%tree%result(1:n_fun)
    err(1:n_fun)  = this%tree%error(1:n_fun)

        !print *, '# neval', this%n_eval
        !print *, '# ', this%tree%result(:), this%tree%error(:)


        deallocate(elem1)
        deallocate(elem)
        deallocate(h)
        deallocate(this%tree)
        deallocate(this%cubature_rule)

end subroutine
!_______________________________________________________________________
!
logical function converged(tree, AbsError, RelError)

    type(heap) :: tree
    real(r8), intent(in) :: AbsError, RelError
    integer :: i

    converged = .true.
    do i = 1, tree%n_fun
        if((tree%error(i).gt.AbsError).and.(tree%error(i).gt.abs(tree%result(i))*RelError))converged = .false.
    end do

end function
!_______________________________________________________________________
!
subroutine eval_region(elem, func, cubature_rule)

    type(element_t)     :: elem

    interface
        function func(n_fun, x)
            import
            integer, intent(in)        :: n_fun
            real(r8), intent(in)       :: x(:)
            real(r8), dimension(n_fun) :: func
        end function
    end interface
    type(rule)     :: cubature_rule

    integer  :: i_pts
    real(r8) :: x(1:elem%h%n_dim)
    real(r8) :: f(1:elem%n_fun,1:cubature_rule%n_pts)

    !$omp parallel do default(shared) private(x)
    do i_pts = 1, cubature_rule%n_pts
        x(:) = cubature_rule%x(:,i_pts) * (2._r8) * elem%h%hwidth(:) + elem%h%center(:)
        f(1:elem%n_fun,i_pts) = func(elem%n_fun, x)
    end do
    elem%splitDim = cubature_rule%apply(elem%h, elem%n_fun, f, elem%result, elem%error)
    elem%max_error = maxval(elem%error)

end subroutine
!_______________________________________________________________________
!
end module
