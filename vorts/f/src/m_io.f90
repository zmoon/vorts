!> I/O routines and simulation settings
module m_io
  use m_params, only: rk=>realkind
  use m_vorts, only: vorton_type
  implicit none

  private
  public :: simsettings_type, parse_input_txt, write_output_txt, time_stepper_from_index

  type :: simsettings_type
    real(rk) :: dt  ! simulation time step
    integer :: n_timesteps  ! number of time steps
    integer :: n_vortons  ! number of vortons (first in the input)
    integer :: n_tracers  ! number of tracers (after vortons in the input)
    integer :: n_total  ! number of vortons + tracers
    character(len=10) :: time_stepper_name
    logical :: write_vortons, write_tracers, write_ps  ! whether to write these output files
    procedure(time_stepper_int), pointer, nopass :: take_time_step
  end type simsettings_type

  ! interface simsettings_type
  !   module procedure :: settings_from_nml
  ! end interface simsettings_type

  !> Interface for stepper selection, in order to allow making a procedure pointer
  !> Can't put in m_vorts source since it uses the vorton type
  abstract interface
    subroutine time_stepper_int(vortons, dt, l, n_total, n_vortons)
      import :: vorton_type, rk  ! in order to access these in the interface

      type(vorton_type), intent(inout), dimension(:) :: vortons  ! save values to hists here
      real(rk), intent(in) :: dt
      integer, intent(in)  :: l   ! time step to be calculated
      integer, intent(in)  :: n_total  ! num vortons + tracers
      integer, intent(in)  :: n_vortons

    end subroutine time_stepper_int
  end Interface

contains

  !> Stepper selector, called when creating the settings instance
  function time_stepper_from_name(name) result(f_ptr)
    use m_vorts, only: step_FT, step_RK4

    character(len=*), intent(in) :: name
    procedure(time_stepper_int), pointer :: f_ptr

    select case ( trim(name) )
      case ('FT')
        f_ptr => step_FT
      case ('RK4')
        f_ptr => step_RK4
      case default
        stop 'invalid integration routine name'
    end select
  end

  function time_stepper_from_index(i) result(f_ptr)
    use m_vorts, only: step_FT, step_RK4

    integer, intent(in) :: i
    procedure(time_stepper_int), pointer :: f_ptr

    select case (i)
      case (1)
        f_ptr => step_FT
      case (2)
        f_ptr => step_RK4
      case default
        stop 'invalid integration routine index'
    end select
  end


  ! !> Load sim settings data from namelist
  ! type(simsettings_type) function settings_from_nml(nml='./in/settings.nml') result(settings)
  !   character(len=*), intent(in) :: nml

  !   real(rk) :: dt
  !   integer :: nt, n_vortons, n_tracers
  !   character(len=*) :: time_stepper_name
  !   logical :: write_vortons, write_tracers, write_ps

  !   namelist /settings/ dt, nt, n_vortons, n_tracers, time_stepper_name, &
  !     write_vortons, write_tracers, write_ps
  !   open(10, file=nml)
  !   read(10, nml=settings)
  !   close(10)

  !   settings%dt = dt
  !   settings%n_timesteps = nt
  !   ! settings%n_vortons = n_vortons
  !   ! settings%n_tracers = n_tracers
  !   ! settings%n_total = n_vortons + n_tracers
  !   settings%time_stepper_name = time_stepper_name
  !   settings%write_vortons = write_vortons
  !   settings%write_tracers = write_tracers
  !   settings%write_ps = write_ps

  ! end function settings_from_nml


  !> Load sim settings from normal text file
  function settings_from_txt(fp) result(settings)
    character(len=*), intent(in) :: fp
    type(simsettings_type) :: settings

    real(rk) :: dt
    integer :: nt
    character(len=50) :: time_stepper_name
    logical :: write_vortons, write_tracers, write_ps

    open(unit=9, file=fp)
    read(9, *) dt
    read(9, *) nt
    read(9, *) time_stepper_name
    read(9, *) write_vortons
    read(9, *) write_tracers
    read(9, *) write_ps
    close(9)

    settings%dt = dt
    settings%n_timesteps = nt
    settings%time_stepper_name = trim(time_stepper_name)
    settings%write_vortons = write_vortons
    settings%write_tracers = write_tracers
    settings%write_ps = write_ps
    settings%take_time_step => time_stepper_from_name(settings%time_stepper_name)

  end function settings_from_txt


  !> Parse input, creating the `Vortons` and `SimSettings` instances
  subroutine parse_input_txt(vortons, settings)
    type(vorton_type), dimension(:), allocatable, intent(out) :: vortons
    type(simsettings_type), intent(out) :: settings

    integer :: iline, ios, skiprows
    integer :: n_vortons, n_tracers, n_total
    real(rk) :: Gamma, xi, yi
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
    open(unit=10, file='./in/vortons.txt', iostat=ios)
    do iline = 1, n_total + skiprows
      if ( iline <= skiprows ) then
        read(10, *)
      else
        !> Read and add to the vortons array
        i = iline - skiprows  ! vorton index
        read(10, *) Gamma, xi, yi
        vortons(i) = vorton_type(Gamma, xi, yi, settings%n_timesteps)

        !> Increment counts
        if ( Gamma == 0 ) then
          n_tracers = n_tracers + 1
        else
          n_vortons = n_vortons + 1
        end if

        !> Check that vortons are first and then tracers only
        if ( n_tracers > 0 .and. Gamma /= 0 ) stop 'Tracers must come after true vortons in the input.'
      end if
      !> If we have gone too far, there is a problem with `n_total` or `skiprows`
      if ( ios /= 0 ) stop 'Problem with vorton input file format'
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
    integer, parameter :: maxrecs = 1001

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
    type(vorton_type), dimension(:), intent(in) :: vortons
    type(simsettings_type) :: settings

    character(len=50) :: s_n_vortons, fmt_csv
    integer :: i, j

    associate( &
      nt => settings%n_timesteps, &
      n_vortons => settings%n_vortons, &
      n_total => settings%n_total &
    )

    !> Construct the format string for vorton/tracer CSV output
    write(s_n_vortons, '(i0)') nt+1 - 1
    fmt_csv = '(g0.5, ' // trim(s_n_vortons) // '(",", g0.5))'

    !> Write vortons?
    if ( settings%write_vortons ) then
      open(unit=101, file='./out/vortons.csv')
      write(101, fmt=*) '# x1(1:nt); y1; x2; y2; ...'
      do i = 1, n_vortons
          write(101, fmt=fmt_csv) vortons(i)%xhist
          write(101, fmt=fmt_csv) vortons(i)%yhist
      end do
      close(101)
    end if

    !> Write tracers?
    if ( settings%write_tracers ) then
      open(unit=102, file='./out/tracers.csv')
      write(102, fmt=*) '# x1(1:nt); y1; x2; y2; ...'
      do i = n_vortons + 1, n_total
        write(102, fmt=fmt_csv) vortons(i)%xhist
        write(102, fmt=fmt_csv) vortons(i)%yhist
      end do
      close(102)
    end if

    !> Write Poincare section?
    !> ! Only works for equilateral triangle with 2nd vorton as top of triangle with x=0 !
    if ( settings%write_ps ) then
      open(unit=103, file='./out/ps.txt')
      write(103, fmt=*) '# x y'
      do j = 2, nt+1
        do i = n_vortons+1, n_total  ! tracers only!
          !> Check that x[j] < xi[i]) & (x[j-1] > xi[i]) & (y[j] > 0
          if ( vortons(2)%xhist(j) < 0. .and. vortons(2)%xhist(j-1) > 0. .and. vortons(2)%yhist(j) > 0 ) then
            write(103, fmt=*) vortons(i)%xhist(j), vortons(i)%yhist(j)
          end if
        end do
      end do
      close(103)
    end if

    end associate

  end subroutine write_output_txt

end module m_io
