!> C interface
module capi
  use iso_c_binding
  use m_vorts
  implicit none

  public

contains

  subroutine run( &
    nv, x, y, G, &
    dt, nt, &
    imethod, &
    xout, yout &
  ) bind(c, name='asdf')
    !
    integer(c_int), intent(in), value :: nv
    real(c_double), intent(in) :: x(nv), y(nv), G(nv)
    real(c_double), intent(in), value :: dt
    integer(c_int), intent(in), value :: nt
    integer(c_int), intent(in), value :: imethod
    !
    real(c_double), intent(inout) :: xout(nt+1), yout(nt+1)
    !
  
    print *, nv
    print *, x, y, G
    print *, dt, nt
    print *, imethod
    print *, xout, yout

  end subroutine run

end module
