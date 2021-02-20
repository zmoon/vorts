!> Parameters used by other modules
!> Kinds approach based on https://github.com/wavebitscientific/wavy/blob/master/src/lib/mod_precision.f90
module m_params
  use iso_fortran_env, only: int16, int32, int64, real32, real64, real128
  implicit none

  private
  public :: intkind, realkind, pi

  integer, parameter :: realkind = real64
  integer, parameter :: intkind = int32

  real(kind=realkind), parameter :: pi = 4*atan(1.0_realkind)

end module m_params
