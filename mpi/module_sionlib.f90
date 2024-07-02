module hydro_sionlib
#ifdef SIONLIB
  use hydro_utils,only:title
  implicit none

contains

subroutine output_sionlib
  use hydro_commons
  use hydro_mpi,only:rang
  implicit none

  if (rang==0) print *,'WARNING: fileformat_out ',trim(fileformat_out),' not yet implemented'
end subroutine output_sionlib


subroutine input_sionlib
  use hydro_commons
  use hydro_mpi,only:rang
  implicit none

  if (rang==0) print *,'WARNING: fileformat_in ',trim(fileformat_out),' not yet implemented'
end subroutine input_sionlib

#endif
end module hydro_sionlib
