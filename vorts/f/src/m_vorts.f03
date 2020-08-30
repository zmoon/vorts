!
! module for advecting N vortons
!   including vorton class and integration schemes which take arrays of vortons as input
!

module m_vorts
  implicit none

  private
  public construct_vorton, count_num_in, FT_step, RK4_step, &
    dp, pi

  integer, parameter :: dp = kind(1.0d0)  ! chosen precision
  real(kind=dp), parameter :: pi = 4*atan(1.0_dp)


  !> Vorton class
  type, public :: Vorton

    !> attributes
    real(dp) :: G   ! vorton strength Gamma
    real(dp), dimension(:), allocatable :: xhist, yhist  ! position history. or 2xN array vhist = [xhist; yhist] ?

  contains  ! bound methods (fortran functions and procedures)
    ! currently no methods. for reference, they are written like this:
    ! procedure :: calc_l => vorton_l  ! "intervortical distance"

  end type Vorton

  !> constructor for Vorton class
  interface Vorton
    procedure construct_vorton
  end interface Vorton


contains

  !====== Vorton class methods =====================================================================

  !> Vorton class constructor
  function construct_vorton(G, xi, yi, nt) result(new_vorton)
    type(Vorton) :: new_vorton
    real(dp), intent(in) :: G, xi, yi
    integer, intent(in) :: nt

    new_vorton%G = G

    allocate(new_vorton%xhist(nt+1))  ! these don't really have to be allocatable either, could be static
    allocate(new_vorton%yhist(nt+1))
    new_vorton%xhist(1) = xi
    new_vorton%yhist(1) = yi

  end function construct_vorton


  !====== Other functions ==========================================================================

  !> distance between one vorton and another, squared
  function calc_lsqd(x1, y1, x2, y2) result(lsqd)
    real(dp), intent(in) :: x1, y1, x2, y2  ! two sets of coords
    real(dp) :: lsqd

    lsqd = (x1 - x2)**2 + (y1 - y2)**2

  end function calc_lsqd


  !> dxdt and dydt calculations
  function calc_tends(G, x0, y0, n) result(tends)
    real(dp), intent(in), dimension(:) :: G, x0, y0  ! arrays of coords and Gamma vals
    integer, intent(in) :: n  ! num vortons + tracers
    real(dp), dimension(2,n) :: tends

    integer :: i, j
    real(dp) :: x0i, y0i, x0j, y0j, Gj, dxdti, dydti, lij_sqd

    do i = 1, n
      x0i = x0(i)
      y0i = y0(i)

      dxdti = 0
      dydti = 0
      do j = 1, n
        Gj  = G(j)

        if ( i /= j .and. Gj /= 0.) then
          x0j = x0(j)
          y0j = y0(j)

          lij_sqd = calc_lsqd(x0i, y0i, x0j, y0j)
          dxdti = dxdti + -1/(2*pi) * Gj * (y0i-y0j) / lij_sqd
          dydti = dydti +  1/(2*pi) * Gj * (x0i-x0j) / lij_sqd

        end if

      end do

      tends(1, i) = dxdti
      tends(2, i) = dydti

    end do

  end function calc_tends


  !> Euler forward integration (1st-order)
  subroutine FT_step(vortons, dt, l, n)
    type(Vorton), intent(inout), dimension(:) :: vortons  ! save values to hists here
    real(dp), intent(in) :: dt
    integer, intent(in)  :: l  ! time step to be calculated
    integer, intent(in)  :: n  ! num vortons + tracers

    real(dp), dimension(n) :: G, x0, y0  ! arrays of coords and Gamma vals
    integer :: i
    real(dp), dimension(2,n) :: tends
    real(dp), dimension(n) :: xtend, ytend, xnew, ynew

    !> prepare arrays for input into calc_tends
    !  could be a separate subroutine as well
    do i = 1, n
      G(i)  = vortons(i)%G
      x0(i) = vortons(i)%xhist(l-1)
      y0(i) = vortons(i)%yhist(l-1)
    end do

    tends = calc_tends(G, x0, y0, n)
    xtend = tends(1,:)
    ytend = tends(2,:)

    xnew = x0 + xtend*dt
    ynew = y0 + ytend*dt

    !> update hist values (could be separate subroutine)
    do i = 1, n
      vortons(i)%xhist(l) = xnew(i)
      vortons(i)%yhist(l) = ynew(i)
    end do

  end subroutine FT_step


  !> add centered time?


  !> RK4 integration
  subroutine RK4_step(vortons, dt, l, n)
    type(Vorton), intent(inout), dimension(:) :: vortons  ! save values to hists here
    real(dp), intent(in) :: dt
    integer, intent(in)  :: l   ! time step to be calculated
    integer, intent(in)  :: n  ! num vortons + tracers

    real(dp), dimension(n) :: G, x0, y0  ! arrays of coords and Gamma vals
    integer :: i
    integer :: num_vortons
    real(dp), dimension(2,n) :: tends
    real(dp), dimension(n) :: xnew, ynew

    real(dp), dimension(n) :: x1, x2, x3, x4, &
                              y1, y2, y3, y4, &
                              k1x, k2x, k3x, k4x, &
                              k1y, k2y, k3y, k4y

    !> prepare arrays for input into calc_tends
    !  this could be a separate subroutine as since done in FT_step also
    do i = 1, n
      G(i)  = vortons(i)%G
      x0(i) = vortons(i)%xhist(l-1)
      y0(i) = vortons(i)%yhist(l-1)
    end do

    x1 = x0
    y1 = y0
    tends = calc_tends(G, x1, y1, n)
    k1x = tends(1,:)
    k1y = tends(2,:)

    x2 = x0 + dt/2*k1x
    y2 = y0 + dt/2*k1y
    tends = calc_tends(G, x2, y2, n)
    k2x = tends(1,:)
    k2y = tends(2,:)

    x3 = x0 + dt/2*k2x
    y3 = y0 + dt/2*k2y
    tends = calc_tends(G, x3, y3, n)
    k3x = tends(1,:)
    k3y = tends(2,:)

    x4 = x0 + dt/1*k3x
    y4 = y0 + dt/1*k3y
    tends = calc_tends(G, x4, y4, n)
    k4x = tends(1,:)
    k4y = tends(2,:)

    xnew = x0 + dt/6*(k1x + 2*k2x + 2*k3x + k4x)
    ynew = y0 + dt/6*(k1y + 2*k2y + 2*k3y + k4y)

    !> update hist values (could be separate subroutine)
    do i = 1, n
      vortons(i)%xhist(l) = xnew(i)
      vortons(i)%yhist(l) = ynew(i)
    end do

  end subroutine RK4_step


  !> determine number of input vortons
  !  or could make subroutine to do this and also load the data for ifname
  function count_num_in(ifname) result(num_in)
    character(len=*), intent(in) :: ifname
    integer :: num_in, iline, ios
    character(len=1) :: firstchar
    integer, parameter :: maxrecs=201

    open(unit=10, file=trim(ifname))
    num_in = 0
    do iline = 1, maxrecs
      read(10, *, iostat=ios) firstchar
      if ( ios /= 0 ) exit
      if ( iline == maxrecs ) then
        print *, 'Error: maximum number of records exceeded...'
        print *, 'exiting program now...'
        stop
      end if
      ! if ( iline > skiprows ) then
      !   num_vortons = num_vortons + 1
      ! else
      !   ! print *, 'reading header'
      ! end if
      if ( firstchar == '#' ) then
        ! print *, 'reading header'
      else
        num_in = num_in + 1
        ! print *, 'reading header'
      end if
    end do
    close(unit=10)

  end function count_num_in

end module m_vorts