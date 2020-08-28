!
! run simulation of interactions between point vortices
! using module m_vorts
!
! loads simulation inputs and writes out the output paths
!

program lets_do_it
  use m_vorts
  implicit none

  !> i/o variables
  real(dp) :: Gamma, xi, yi  ! input strengths and initial locations
  real(dp) :: dxdt_test, dydt_test
  integer :: iline, ios, i, j, l, num_vortons, num_tracers, num_total
  character(len=30) :: f1, f2, f3, f4, f5
  character(len=1) :: firstchar
  character(len=50) :: ofname, vortnum, tracernum
  integer :: ofunit
  integer, parameter :: maxrecs = 100, skiprows = 1

  !> simulation variables
  real(dp) :: dt  ! simulation time step
  integer  :: nt ! number of time steps
  character(len=10) :: integration_routine_name
  type(Vorton), dimension(:), allocatable :: vortons, tracers, vs_and_ts  ! arrays of vortons!
  ! pointer :: integration_routine

  !> read in simulation info: dt and nt
  open(unit=9, file='./in/vorts_sim_in.txt')
  read(9, *) dt
  read(9, *) nt
  read(9, *) integration_routine_name
  close(9)

  !> determine number of inputs so can allocate
  num_vortons = count_num_in('./in/vorts_in.txt')  ! gotta run at the higher level for this to work
  num_tracers = count_num_in('./in/tracers_in.txt')
  allocate(vortons(num_vortons))
  allocate(tracers(num_tracers))
  num_total = num_vortons + num_tracers
  allocate(vs_and_ts(num_total))

  !> load data this time
  open(unit=10, file='./in/vorts_in.txt')
  do iline = 1, num_vortons+skiprows
    if ( iline <= skiprows ) then    ! ultimately don't want to hardcode skiprows. this loading could be separate fn or subroutine?
      read(10, *) firstchar          ! could also check for comment char here too
      ! print *, 'reading header again'
    else
      i = iline - skiprows
      read(10, *) Gamma, xi, yi
      vortons(i) = construct_vorton(Gamma, xi, yi, nt)
    end if
  end do
  close(10)

  open(unit=10, file='./in/tracers_in.txt')
  do iline = 1, num_tracers+skiprows
    if ( iline <= skiprows ) then
      read(10, *) firstchar
    else
      i = iline - skiprows
      read(10, *) xi, yi  ! tracers have no Gamma
      Gamma = 0
      tracers(i) = construct_vorton(Gamma, xi, yi, nt)
    end if
  end do
  close(10)

  !> construct combined array
  print *, 'constructing combined array'
  do i = 1, num_total
    if ( i <= num_vortons ) then
      vs_and_ts(i) = vortons(i)
    else
      vs_and_ts(i) = tracers(i-num_vortons)
    end if
  end do

  !> check that things were loaded correctly
  print *, ''
  print *, 'The following vortons have been loaded:'
  print *, '  Vorton  Gamma     xi      yi'
  f1 = '(2x, i3, 3x, 3(f8.2))'
  do i = 1, num_vortons
    ! print f1, i, vortons(i)%G, vortons(i)%xhist(1), vortons(i)%yhist(1)
    print f1, i, vs_and_ts(i)%G, vs_and_ts(i)%xhist(1), vs_and_ts(i)%yhist(1)
  end do

  print *, ''
  print *, 'And the first 10 tracers'
  print *, '  Tracer   xi      yi'
  do i = 1, num_tracers
    print f1, i, tracers(i)%xhist(1), tracers(i)%yhist(1)
    if ( i == 10 ) exit
  end do

  !> make sure things match up in vs_and_ts
  print *, ''
  print *, 'Everything in vs_and_ts'
  print *, '  Vorton  Gamma     xi      yi'
  f1 = '(2x, i3, 3x, 3(f8.2))'
  do i = 1, num_total
    print f1, i, vs_and_ts(i)%G, vs_and_ts(i)%xhist(1), vs_and_ts(i)%yhist(1)
  end do

  !> do some integrating
  print *, ''
  print *, 'now doing some integrating'
  do l = 2, nt+1  ! start at 2 to not overwrite values in history
    if ( trim(integration_routine_name) == 'FT' ) then  ! is there a way to assign subroutine to variable so as to avoid this if ?

      call FT_step(vs_and_ts, dt, l, num_total)

    else if ( trim(integration_routine_name) == 'RK4' ) then

      call RK4_step(vs_and_ts, dt, l, num_total)

    end if
  end do

  !> save the results
  !  individual files versions
  ! f3 = '(2(1x, f0.8))'  ! could use sci notation instead? also supposedly free-form is much faster
  !
  ! print *, ''
  ! print *, 'Now writing:'
  !
  ! do i = 1, num_vortons+num_tracers
  !   if ( i <= num_vortons ) then
  !
  !     write(vortnum, '(i0)') i
  !     ofunit = i + 10
  !     ofname = './out/vorton' // trim(vortnum) // '.txt'
  !     print *, ofname
  !     open(unit=ofunit, file=ofname)
  !
  !     write(ofunit, fmt='(a)') '# x y'
  !     do j = 1, nt+1
  !       ! write(ofunit, fmt=f3) vortons(i)%xhist(j), vortons(i)%yhist(j)
  !       write(ofunit, fmt=f3) vs_and_ts(i)%xhist(j), vs_and_ts(i)%yhist(j)
  !     end do
  !     close(ofunit)
  !
  !   else
  !
  !     write(tracernum, '(i0)') i-num_vortons
  !     ofunit = i + 10
  !     ofname = './out/tracer' // trim(tracernum) // '.txt'
  !     print *, ofname
  !     open(ofunit, file=ofname)
  !
  !     write(ofunit, fmt='(a)') '# x y'
  !     do j = 1, nt+1
  !       ! write(ofunit, fmt=f3) tracers(i)%xhist(j), tracers(i)%yhist(j)
  !       write(ofunit, fmt=f3) vs_and_ts(i)%xhist(j), vs_and_ts(i)%yhist(j)
  !     end do
  !     close(ofunit)
  !
  !   end if
  ! end do

  !> write all vortons to one csv file
!  write(vortnum, '(i0)') nt+1 - 1
!  f4 = '(g0.5, ' // trim(vortnum) // '(",", g0.5))'
!
!  print *, f4
!
!  open(unit=101, file='./out/vortons.csv')
!  write(101, fmt=*) '# x1(1:nt); y1; x2; y2; ...'
!
!  open(unit=102, file='./out/tracers.csv')
!  write(102, fmt=*) '# x1(1:nt); y1; x2; y2; ...'
!
!  do i = 1, num_total
!
!    if ( i <= num_vortons ) then
!
!      write(101, fmt=f4) vs_and_ts(i)%xhist
!      write(101, fmt=f4) vs_and_ts(i)%yhist
!
!    else
!
!      write(102, fmt=f4) vs_and_ts(i)%xhist
!      write(102, fmt=f4) vs_and_ts(i)%yhist
!
!    end if
!  end do


  !> Poincare section, using 2nd vorton (should be the one at the top, (0, 1)) as ref
  !  use free form output for speed ??
  open(unit=103, file='./out/ps.txt')
  write(103, fmt=*) '# x y'

  do j = 2, nt+1

    do i = num_vortons+1, num_total  ! looping through tracers only

      ! x[j] < xi[i]) & (x[j-1] > xi[i]) & (y[j] > 0
      if ( vs_and_ts(2)%xhist(j) < 0. .and. vs_and_ts(2)%xhist(j-1) > 0. .and. vs_and_ts(2)%yhist(j) > 0 ) then

        write(103, fmt=*) vs_and_ts(i)%xhist(j), vs_and_ts(i)%yhist(j)

      end if

    end do
  end do




end program lets_do_it
