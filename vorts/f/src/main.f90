!> This program reads simulation input from text files, runs the simulation,
!> and writes output to text file.
program lets_do_it
  use m_io
  use m_vorts
  implicit none

  type(Vorton), dimension(:), allocatable :: vortons
  type(SimSettings) :: settings

  character(len=50) :: f1
  integer :: i  ! vorton index
  integer :: l  ! current time step

  !> Read in initial vortons, construct their history arrays, read in settings
  call parse_input_txt(vortons, settings)

  !> Check that things were loaded correctly
  print *, ''
  print *, 'The following vortons have been loaded:'
  print *, '  Vorton  Gamma     xi      yi'
  f1 = '(2x, i3, 3x, 3(f8.2))'
  do i = 1, settings%n_vortons
    print f1, i, vortons(i)%G, vortons(i)%xhist(1), vortons(i)%yhist(1)
  end do
  print *, ''
  print *, 'And the first 10 tracers:'
  print *, '  Tracer    xi      yi'
  do i = settings%n_vortons + 1, settings%n_total
    print f1, i - settings%n_vortons, vortons(i)%xhist(1), vortons(i)%yhist(1)
    if ( (i - settings%n_vortons) == 10 ) exit
  end do

  !> Integrate
  print *, ''
  print *, 'Now integrating...'
  time_loop: do l = 2, settings%n_timesteps+1  ! start at 2 to not overwrite values in history
    if ( settings%integration_routine_name == 'FT' ) then  ! is there a way to assign subroutine to variable so as to avoid this if ?

      call FT_step(vortons, settings%dt, l, settings%n_total, settings%n_vortons)

    else if ( settings%integration_routine_name == 'RK4' ) then

      call RK4_step(vortons, settings%dt, l, settings%n_total, settings%n_vortons)

    else

      stop 'Invalid `integration_routine_name`'

    end if
  end do time_loop

  !> Write output
  print *, ''
  print *, 'Now writing output...'
  call write_output_txt(vortons, settings)

  !> Done!
  print *, ''
  print *, 'Done!'

end program lets_do_it
