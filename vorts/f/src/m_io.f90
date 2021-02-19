!> I/O routines
module m_io
  use m_vorts, only: dp, Vorton
  implicit none

  private
  public :: SimSettings, parse_input_txt, write_output_txt

  type :: SimSettings
    real(dp) :: dt  ! simulation time step
    integer :: nt  ! number of time steps
    integer :: n_vortons  ! number of vortons (first in the input)
    integer :: n_tracers  ! number of tracers (after vortons in the input)
    integer :: n_total  ! number of vortons + tracers
    character(len=10) :: integration_routine_name
    logical :: write_vortons, write_tracers, write_ps  ! whether to write these output files

  !> Bound methods
  ! contains

  end type SimSettings

  ! interface SimSettings
  !   module procedure :: settings_from_nml
  ! end interface SimSettings

contains

  ! !> Load sim settings data from the nml
  ! type(SimSettings) function settings_from_nml(nml='./in/settings.nml') result(settings)
  !   character(len=*), intent(in) :: nml

  !   real(dp) :: dt
  !   integer :: nt, n_vortons, n_tracers
  !   character(len=*) :: integration_routine_name
  !   logical :: write_vortons, write_tracers, write_ps

  !   namelist /settings/ dt, nt, n_vortons, n_tracers, integration_routine_name, &
  !     write_vortons, write_tracers, write_ps
  !   open(10, file=nml)
  !   read(10, nml=settings)
  !   close(10)

  !   settings%dt = dt
  !   settings%nt = nt
  !   ! settings%n_vortons = n_vortons
  !   ! settings%n_tracers = n_tracers
  !   ! settings%n_total = n_vortons + n_tracers
  !   settings%integration_routine_name = integration_routine_name
  !   settings%write_vortons = write_vortons
  !   settings%write_tracers = write_tracers
  !   settings%write_ps = write_ps

  ! end function settings_from_nml


  !> Load sim settings from normal text file
  ! type(SimSettings) function settings_from_txt(fp) result(settings)
  ! function settings_from_txt(fp='./in/settings.txt') result(settings)
  function settings_from_txt(fp) result(settings)
    character(len=*), intent(in) :: fp
    type(SimSettings) :: settings

    real(dp) :: dt
    integer :: nt, n_vortons, n_tracers
    character(len=50) :: integration_routine_name
    logical :: write_vortons, write_tracers, write_ps

    open(unit=9, file=fp)
    read(9, *) dt
    read(9, *) nt
    ! read(9, *) n_vortons
    ! read(9, *) n_tracers
    read(9, *) integration_routine_name
    read(9, *) write_vortons
    read(9, *) write_tracers
    read(9, *) write_ps
    close(9)

    settings%dt = dt
    settings%nt = nt
    ! settings%n_vortons = n_vortons
    ! settings%n_tracers = n_tracers
    ! settings%n_total = n_vortons + n_tracers
    settings%integration_routine_name = trim(integration_routine_name)
    settings%write_vortons = write_vortons
    settings%write_tracers = write_tracers
    settings%write_ps = write_ps

  end function settings_from_txt


  ! !> Return array of Vorton, reading initial position and G value from text file
  ! function vortons_from_txt(ifname, n, skiprows=1) result(vortons)
  !   character(len=*), intent(in) :: ifname
  !   integer, intent(in) :: n  ! number of vortons (to allocate)
  !   integer, intent(in) :: skiprows
  !   type(Vorton), dimension(:), allocatable :: vortons
  !   integer :: i, iline

  !   allocate(vortons(n))

  !   open(unit=10, file='./in/vorts_in.txt')
  !   do iline = 1, n + skiprows
  !     if ( iline <= skiprows ) then
  !       read(10, *)
  !     else
  !       i = iline - skiprows
  !       read(10, *) Gamma, xi, yi
  !       vortons(i) = Vorton(Gamma, xi, yi, nt)
  !     end if
  !   end do
  !   close(10)

  ! end function vortons_from_txt


  !> Parse input, creating the `Vortons` and `SimSettings` instances
  subroutine parse_input_txt(vortons, settings)
    type(Vorton), dimension(:), allocatable, intent(out) :: vortons
    type(SimSettings), intent(out) :: settings

    integer :: iline, ios, skiprows
    integer :: n_vortons, n_tracers, n_total
    real(dp) :: Gamma, xi, yi
    integer :: i

    !> First read the settings that are in the settings file
    !> This doesn't include the numbers, which we calculate below when reading the vorton input
    settings = settings_from_txt('./in/settings.txt')

    !> Allocate vorton array
    call count_lines_in_txt('./in/vortons.txt', n_total, skiprows)
    allocate(vortons(n_total))

    !> Read the vortons input file
    n_vortons = 0
    n_tracers = 0
    open(unit=10, file='./in/vortons.txt')
    do iline = 1, n_total + skiprows
      if ( iline <= skiprows ) then
        read(10, *)
      else
        !> Read and add to the vortons array
        i = iline - skiprows  ! vorton index
        read(10, *) Gamma, xi, yi
        vortons(i) = Vorton(Gamma, xi, yi, settings%nt)

        !> Increment counts
        if ( Gamma == 0 ) then
          n_tracers = n_tracers + 1
        else
          n_vortons = n_vortons + 1
        end if

        !> Check that vortons are first and then tracers
        if ( n_tracers > 0 .and. Gamma /= 0 ) stop 'Tracers must come after true vortons in the input.'

      end if
    end do
    close(10)

    !> Check that the total number of vorton lines makes sense
    if ( n_total /= n_vortons + n_tracers ) stop 'Vorton input line count inconsistent'

    !> Add vorton counts to the settings
    settings%n_vortons = n_vortons
    settings%n_tracers = n_tracers
    settings%n_total = n_total

  end subroutine parse_input_txt


  !> Count the number of lines, skipping header lines (that start with `#`)
  subroutine count_lines_in_txt(ifname, num_vorton_lines, num_header_lines)
    character(len=*), intent(in) :: ifname
    integer, intent(out) :: num_vorton_lines, num_header_lines
    integer :: iline, ios
    character(len=1) :: firstchar
    integer, parameter :: maxrecs=1001

    open(unit=10, file=trim(ifname))
    num_vorton_lines = 0
    num_header_lines = 0
    do iline = 1, maxrecs
      read(10, *, iostat=ios) firstchar
      if ( ios /= 0 ) exit  ! eof
      if ( iline == maxrecs ) then
        print *, 'Error: maximum number of records exceeded...'
        print *, 'exiting program now...'
        stop
      end if
      if ( firstchar == '#' ) then
        num_header_lines = num_header_lines + 1
      else
        num_vorton_lines = num_vorton_lines + 1
      end if
    end do
    close(unit=10)

  end subroutine count_lines_in_txt


  !> Write text output
  subroutine write_output_txt(vortons, settings)
    type(Vorton), dimension(:), intent(in) :: vortons
    type(SimSettings) :: settings

    character(len=50) :: vortnum, f4
    integer :: i, j

    associate( &
      nt => settings%nt, &
      n_vortons => settings%n_vortons, &
      n_total => settings%n_total &
    )

    !> save the results
    write(vortnum, '(i0)') nt+1 - 1
    f4 = '(g0.5, ' // trim(vortnum) // '(",", g0.5))'
    ! print *, f4

    open(unit=101, file='./out/vortons.csv')
    write(101, fmt=*) '# x1(1:nt); y1; x2; y2; ...'

    open(unit=102, file='./out/tracers.csv')
    write(102, fmt=*) '# x1(1:nt); y1; x2; y2; ...'

    do i = 1, n_total

      if ( i <= n_vortons ) then

        if ( settings%write_vortons ) then
          write(101, fmt=f4) vortons(i)%xhist
          write(101, fmt=f4) vortons(i)%yhist
        end if

      else

        if ( settings%write_tracers ) then
          write(102, fmt=f4) vortons(i)%xhist
          write(102, fmt=f4) vortons(i)%yhist
        end if

      end if
    end do

    close(101)
    close(102)


    !> Poincare section, using 2nd vorton (should be the one at the top, (0, 1)) as ref
    if ( settings%write_ps ) then

      open(unit=103, file='./out/ps.txt')
      write(103, fmt=*) '# x y'

      do j = 2, nt+1

        do i = n_vortons+1, n_total  ! looping through tracers only

          ! x[j] < xi[i]) & (x[j-1] > xi[i]) & (y[j] > 0
          if ( vortons(2)%xhist(j) < 0. .and. vortons(2)%xhist(j-1) > 0. .and. vortons(2)%yhist(j) > 0 ) then

            write(103, fmt=*) vortons(i)%xhist(j), vortons(i)%yhist(j)

          end if

        end do
      end do

    end if

    close(103)

    end associate

  end subroutine write_output_txt

end module m_io
