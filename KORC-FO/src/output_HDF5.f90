module output_HDF5
use HDF5
implicit none

contains

subroutine initialize_HDF5()
use HDF5
implicit none
	INTEGER :: error  ! Error flag
	call h5open_f(error)
end subroutine initialize_HDF5

subroutine finalize_HDF5()
implicit none
	INTEGER :: error  ! Error flag
!	call h5close_f(error)
end subroutine finalize_HDF5

end module output_HDF5
