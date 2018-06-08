!> Contains common definitions for the data and datatypes used by MONC
module datadefn_mod
  use mpi
  use iso_fortran_env
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, public, parameter :: STRING_LENGTH=150 !< Default length of strings
  integer, public, parameter :: LONG_STRING_LENGTH=STRING_LENGTH + 50!< Length of longer strings

  integer, public, parameter :: SINGLE_PRECISION = selected_real_kind(6,30)   !< Single precision (32 bit) kind
  integer, public, parameter :: DOUBLE_PRECISION = selected_real_kind(15,307) !< Double precision (64 bit) kind

  integer, public, parameter :: LONG_INTEGER = int64
  integer, public, parameter :: SHORT_INTEGER = int32

  integer, public, parameter :: PARCEL_INTEGER = LONG_INTEGER

  !< Default precision which is used for prognostic data and calculations
  integer, public, parameter :: DEFAULT_PRECISION = DOUBLE_PRECISION
  !< MPI communication type which we use for the prognostic and calculation data
  integer, public :: PRECISION_TYPE
  integer, public :: MPI_PARCEL_INT

  public init_data_defn

contains

  !> Will initialise the data definitions. This should be called as soon as MONC starts up
  subroutine init_data_defn()
    if (DEFAULT_PRECISION .eq. DOUBLE_PRECISION) then
      PRECISION_TYPE = MPI_DOUBLE_PRECISION
    else
      PRECISION_TYPE = MPI_REAL
    endif

    if (PARCEL_INTEGER .eq. LONG_INTEGER) then
      MPI_PARCEL_INT = MPI_INTEGER8
    else
      MPI_PARCEL_INT = MPI_INTEGER
    endif
  end subroutine init_data_defn
end module datadefn_mod
