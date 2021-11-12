!> C interface
module capi
  use iso_c_binding
  use m_io, only: simsettings_type, time_stepper_from_index
  use m_vorts, only: vorton_type
  implicit none

  private
  public :: run

contains

  subroutine run( &
    nv, x, y, G, &
    dt, nt, &
    imethod, &
    xout, yout &
  ) bind(c, name='asdf')
    ! In
    integer(c_int), intent(in), value :: nv
    real(c_double), intent(in) :: x(nv), y(nv), G(nv)
    real(c_double), intent(in), value :: dt
    integer(c_int), intent(in), value :: nt
    integer(c_int), intent(in), value :: imethod
    ! Out
    real(c_double), intent(inout) :: xout(nv, nt+1), yout(nv, nt+1)
    ! Local
    type(simsettings_type) :: settings
    type(vorton_type), dimension(:), allocatable :: vortons
    integer :: i, l

    ! Initial settings
    settings%dt = dt
    settings%n_timesteps = nt
    settings%write_vortons = .false.
    settings%write_tracers = .false.
    settings%write_ps = .false.
    settings%take_time_step => time_stepper_from_index(imethod)
    settings%n_total = nv

    ! Form vortons array
    settings%n_vortons = 0
    settings%n_tracers = 0
    allocate(vortons(settings%n_total))
    do i = 1, nv
      vortons(i) = vorton_type(G(i), x(i), y(i), settings%n_timesteps)

      if (G(i) == 0) then
        settings%n_tracers = settings%n_tracers + 1
      else
        settings%n_vortons = settings%n_vortons + 1
      end if

      if (settings%n_tracers > 0 .and. G(i) /= 0) stop 'tracers must come after vortons'
    end do

    ! Integrate
    time_loop: do l = 2, settings%n_timesteps+1  ! start at 2 - preserving initial position at l=1
      call settings%take_time_step(vortons, settings%dt, l, settings%n_total, settings%n_vortons)
    end do time_loop

    ! Set outputs
    do i = 1, settings%n_vortons
      xout(i,:) = vortons(i)%xhist
      yout(i,:) = vortons(i)%yhist
    end do

  end subroutine run

end module capi
