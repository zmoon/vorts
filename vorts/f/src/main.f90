!> This program reads simulation input from text files, runs the simulation,
!> and writes output to text file.
program lets_do_it
  use m_io
  use m_vorts
  implicit none

  type(vorton_type), dimension(:), allocatable :: vortons
  type(simsettings_type) :: settings

  character(len=50) :: fmt1
  integer :: i  ! vorton index
  integer :: l  ! current time step

  !> Read in initial vortons & construct their history arrays;
  !> Read in settings and construct settings type
  call parse_input_txt(vortons, settings)

  !> Check that things were loaded correctly
  print *, ''
  print *, 'The following vortons have been loaded:'
  print *, '  Vorton  Gamma     xi      yi'
  fmt1 = '(2x, i3, 3x, 3(f8.2))'
  do i = 1, settings%n_vortons
    print fmt1, i, vortons(i)%G, vortons(i)%xhist(1), vortons(i)%yhist(1)
  end do
  print *, ''
  print *, 'And the first 10 tracers:'
  print *, '  Tracer    xi      yi'
  do i = settings%n_vortons + 1, settings%n_total
    print fmt1, i - settings%n_vortons, vortons(i)%xhist(1), vortons(i)%yhist(1)
    if ( (i - settings%n_vortons) == 10 ) exit
  end do

  !> Integrate
  print *, ''
  print *, 'Now integrating...'
  time_loop: do l = 2, settings%n_timesteps+1  ! start at 2 - preserving initial position at l=1
    call settings%take_time_step(vortons, settings%dt, l, settings%n_total, settings%n_vortons)
  end do time_loop

  !> Write output
  print *, ''
  print *, 'Now writing output...'
  call write_output_txt(vortons, settings)

  !> Done!
  print *, ''
  print *, 'Done!'

end program lets_do_it
