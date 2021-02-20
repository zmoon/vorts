!> Time steppers + Vorton type
module m_vorts
  ! use omp_lib
  use m_params, only: rk=>realkind, pi
  implicit none

  private
  public vorton_type, step_FT, step_RK4

  !> Vorton class
  type :: vorton_type
    !> Attributes
    real(rk) :: G   ! vorton strength Gamma
    real(rk), dimension(:), allocatable :: xhist, yhist  ! position history
  end type vorton_type

  !> Specify constructor for Vorton class
  interface vorton_type
    procedure construct_vorton
  end interface vorton_type

contains

  !> Vorton class constructor
  function construct_vorton(G, xi, yi, nt) result(new_vorton)
    type(vorton_type) :: new_vorton
    real(rk), intent(in) :: G, xi, yi
    integer, intent(in) :: nt

    new_vorton%G = G

    allocate(new_vorton%xhist(nt+1))
    allocate(new_vorton%yhist(nt+1))
    new_vorton%xhist(1) = xi
    new_vorton%yhist(1) = yi

  end function construct_vorton


  !> Distance between one vorton and another, squared
  pure real(rk) function lsqd(x1, y1, x2, y2)
    real(rk), intent(in) :: x1, y1, x2, y2  ! two sets of coords

    lsqd = (x1 - x2)**2 + (y1 - y2)**2

  end function lsqd


  !> Tendency calculations: dxdt and dydt for all point vortices (vortons and tracers)
  subroutine calc_tends(G, x0, y0, n_total, n_vortons, dxdt, dydt)
    real(rk), intent(in), dimension(:) :: G, x0, y0  ! arrays of coords and Gamma vals
    integer, intent(in) :: n_total  ! num vortons + tracers
    integer, intent(in) :: n_vortons
    real(rk), intent(out), dimension(:) :: dxdt, dydt

    integer :: i, j
    real(rk) :: x0i, y0i, x0j, y0j, Gj, dxdti, dydti, lij_sqd

    ! !$omp parallel private (i)
    ! !$omp do
    ! do i = 1, n_total
    do concurrent (i = 1:n_total)
      x0i = x0(i)
      y0i = y0(i)

      !> Sum contributions to the tendency of vorton `i`, which may be a normal vorton or tracer
      dxdti = 0
      dydti = 0
      do j = 1, n_vortons  ! only vortons have G values, no tracers
        if ( i /= j ) then
          Gj  = G(j)
          x0j = x0(j)
          y0j = y0(j)

          lij_sqd = lsqd(x0i, y0i, x0j, y0j)

          dxdti = dxdti + (-1)/(2*pi) * Gj * (y0i-y0j) / lij_sqd
          dydti = dydti +  1/(2*pi) * Gj * (x0i-x0j) / lij_sqd
        end if
      end do

      dxdt(i) = dxdti
      dydt(i) = dydti

    end do
    ! !$omp end do
    ! !$omp end parallel

  end subroutine calc_tends


  !> Take an FT time step -- Euler forward integration (1st-order)
  subroutine step_FT(vortons, dt, l, n_total, n_vortons)
    type(vorton_type), intent(inout), dimension(:) :: vortons  ! positions at time `l` will be updated
    real(rk), intent(in) :: dt
    integer, intent(in)  :: l  ! time step to be calculated
    integer, intent(in)  :: n_total  ! num vortons + tracers
    integer, intent(in)  :: n_vortons

    real(rk), dimension(n_total) :: G, x0, y0  ! arrays of coords and Gamma vals
    integer :: i
    real(rk), dimension(n_total) :: dxdt, dydt, xnew, ynew

    !> Prepare arrays for input into `calc_tends` (could be a separate subroutine)
    do i = 1, n_total
      G(i)  = vortons(i)%G
      x0(i) = vortons(i)%xhist(l-1)  ! TODO: just passing `vortons%xhist(l-1)`
      y0(i) = vortons(i)%yhist(l-1)
    end do

    call calc_tends(G, x0, y0, n_total, n_vortons, dxdt, dydt)

    xnew = x0 + dxdt*dt
    ynew = y0 + dydt*dt

    !> Update hist values (could be separate subroutine)
    do i = 1, n_total
      vortons(i)%xhist(l) = xnew(i)
      vortons(i)%yhist(l) = ynew(i)
    end do

  end subroutine step_FT


  !> Take an RK4 time step
  subroutine step_RK4(vortons, dt, l, n_total, n_vortons)
    type(vorton_type), intent(inout), dimension(:) :: vortons  ! positions at time `l` will be updated
    real(rk), intent(in) :: dt
    integer, intent(in)  :: l   ! time step to be calculated
    integer, intent(in)  :: n_total  ! num vortons + tracers
    integer, intent(in)  :: n_vortons

    real(rk), dimension(n_total) :: G, x0, y0  ! arrays of coords and Gamma vals
    integer :: i
    real(rk), dimension(n_total) :: xnew, ynew

    real(rk), dimension(n_total) :: x1, x2, x3, x4, &
                                    y1, y2, y3, y4, &
                                    k1x, k2x, k3x, k4x, &
                                    k1y, k2y, k3y, k4y

    !> Prepare arrays for input into `calc_tends`
    do i = 1, n_total
      G(i)  = vortons(i)%G
      x0(i) = vortons(i)%xhist(l-1)
      y0(i) = vortons(i)%yhist(l-1)
    end do

    !> RK4
    x1 = x0
    y1 = y0
    call calc_tends(G, x1, y1, n_total, n_vortons, k1x, k1y)

    x2 = x0 + dt/2*k1x
    y2 = y0 + dt/2*k1y
    call calc_tends(G, x2, y2, n_total, n_vortons, k2x, k2y)

    x3 = x0 + dt/2*k2x
    y3 = y0 + dt/2*k2y
    call calc_tends(G, x3, y3, n_total, n_vortons, k3x, k3y)

    x4 = x0 + dt/1*k3x
    y4 = y0 + dt/1*k3y
    call calc_tends(G, x4, y4, n_total, n_vortons, k4x, k4y)

    xnew = x0 + dt/6*(k1x + 2*k2x + 2*k3x + k4x)
    ynew = y0 + dt/6*(k1y + 2*k2y + 2*k3y + k4y)

    !> Update hist values
    do i = 1, n_total
      vortons(i)%xhist(l) = xnew(i)
      vortons(i)%yhist(l) = ynew(i)
    end do

  end subroutine step_RK4

end module m_vorts
