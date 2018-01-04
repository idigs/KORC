module korc_HDF5
	use korc_hpc
	use korc_types
	use korc_constants
	use HDF5

	implicit none

	INTEGER(HID_T), PRIVATE :: KORC_HDF5_REAL ! Real precision used in HDF5

	INTEGER(SIZE_T), PRIVATE :: rp_hdf5 ! Size of real precision used in HDF5

	INTERFACE load_from_hdf5
	  module procedure iload_from_hdf5, rload_from_hdf5
	END INTERFACE

	INTERFACE load_array_from_hdf5
	  module procedure rload_1d_array_from_hdf5, rload_3d_array_from_hdf5, rload_2d_array_from_hdf5
	END INTERFACE

	INTERFACE save_to_hdf5
	  module procedure i1save_to_hdf5,i2save_to_hdf5,i4save_to_hdf5,i8save_to_hdf5,rsave_to_hdf5
	END INTERFACE

	INTERFACE save_1d_array_to_hdf5
	  module procedure isave_1d_array_to_hdf5,rsave_1d_array_to_hdf5
	END INTERFACE

	INTERFACE save_2d_array_to_hdf5
	  module procedure rsave_2d_array_to_hdf5
	END INTERFACE

	INTERFACE save_3d_array_to_hdf5
	  module procedure rsave_3d_array_to_hdf5
	END INTERFACE

	INTERFACE save_array_to_hdf5
	  module procedure isave_1d_array_to_hdf5,rsave_1d_array_to_hdf5,rsave_2d_array_to_hdf5,rsave_3d_array_to_hdf5
	END INTERFACE

	PRIVATE :: rsave_to_hdf5,&
				isave_1d_array_to_hdf5,&
				rsave_1d_array_to_hdf5,&
				rsave_2d_array_to_hdf5,&
				iload_from_hdf5,&
				rload_from_hdf5,&
				rload_1d_array_from_hdf5,&
				rload_3d_array_from_hdf5,&
				rload_2d_array_from_hdf5,&
				i1save_to_hdf5,&
				i2save_to_hdf5,&
				i4save_to_hdf5,&
				i8save_to_hdf5
				
	PUBLIC :: initialize_HDF5,&
				finalize_HDF5,&
				save_simulation_parameters,&
				save_to_hdf5,&
				save_1d_array_to_hdf5,&
				save_2d_array_to_hdf5,&
				load_from_hdf5,&
				load_array_from_hdf5,&
				save_string_parameter,&
				get_last_iteration,&
				save_restart_variables,&
				load_particles_ic

contains


subroutine initialize_HDF5()
	INTEGER :: h5error  ! Error flag
	call h5open_f(h5error)
	
#ifdef HDF5_DOUBLE_PRESICION
	call h5tcopy_f(H5T_NATIVE_DOUBLE, KORC_HDF5_REAL, h5error)
#elif HDF5_SINGLE_PRESICION
	call h5tcopy_f(H5T_NATIVE_REAL, KORC_HDF5_REAL, h5error)
#endif

	call h5tget_size_f(KORC_HDF5_REAL, rp_hdf5, h5error)
end subroutine initialize_HDF5


subroutine finalize_HDF5()
	INTEGER :: h5error  ! Error flag
	call h5close_f(h5error)
end subroutine finalize_HDF5


subroutine iload_from_hdf5(h5file_id,dset,rdatum,attr)
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	INTEGER, INTENT(OUT) :: rdatum
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) :: attr
	CHARACTER(4) :: aname = "Info"
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: atype_id
	INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/1/)
	INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
	INTEGER :: rank = 1
	INTEGER :: arank = 1
	INTEGER(SIZE_T) :: tmplen
	INTEGER(SIZE_T) :: attrlen
	INTEGER :: h5error

	! * * * Read datum from file * * *

	call h5dopen_f(h5file_id, TRIM(dset), dset_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: iload_from_hdf5 --> h5dopen_f")')
	end if

	call h5dread_f(dset_id, H5T_NATIVE_INTEGER, rdatum, dims, h5error)

	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: iload_from_hdf5 --> h5dread_f")')
	end if

	call h5dclose_f(dset_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: iload_from_hdf5 --> h5dclose_f")')
	end if

	if (PRESENT(attr)) then
		! * * * Read attribute from file * * *

		! * * * Read attribute from file * * *
	end if

	! * * * Read datum from file * * *
end subroutine iload_from_hdf5


subroutine rload_from_hdf5(h5file_id,dset,rdatum,attr)
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	REAL(rp), INTENT(OUT) :: rdatum
	REAL :: raw_datum
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) :: attr
	CHARACTER(4) :: aname = "Info"
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: atype_id
	INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/1/)
	INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
	INTEGER :: rank = 1
	INTEGER :: arank = 1
	INTEGER(SIZE_T) :: tmplen
	INTEGER(SIZE_T) :: attrlen
	INTEGER :: h5error

	! * * * Read datum from file * * *

	call h5dopen_f(h5file_id, TRIM(dset), dset_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dopen_f")')
	end if

	call h5dread_f(dset_id, H5T_NATIVE_REAL, raw_datum, dims, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dread_f")')
	end if
	rdatum = REAL(raw_datum,rp)

	call h5dclose_f(dset_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dclose_f")')
	end if

	if (PRESENT(attr)) then
		! * * * Read attribute from file * * *

		! * * * Read attribute from file * * *
	end if

	! * * * Read datum from file * * *
end subroutine rload_from_hdf5


subroutine rload_1d_array_from_hdf5(h5file_id,dset,rdata,attr)
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: rdata
	REAL, DIMENSION(:), ALLOCATABLE :: raw_data
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: attr
	CHARACTER(MAX_STRING_LENGTH) :: aname
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: atype_id
	INTEGER(HSIZE_T), DIMENSION(1) :: dims
	INTEGER(HSIZE_T), DIMENSION(1) :: adims
	INTEGER(SIZE_T) :: tmplen
	INTEGER(SIZE_T) :: attrlen
	INTEGER :: h5error

	dims = (/ shape(rdata) /)

	ALLOCATE( raw_data(dims(1)) )

	! * * * Read data from file * * *

	call h5dopen_f(h5file_id, TRIM(dset), dset_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dopen_f")')
	end if

	call h5dread_f(dset_id, H5T_NATIVE_REAL, raw_data, dims, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dread_f")')
	end if
	rdata = REAL(raw_data,rp)

	call h5dclose_f(dset_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dclose_f")')
	end if

	DEALLOCATE( raw_data )

	if (PRESENT(attr)) then
		! * * * Read data attribute(s) from file * * *
	end if

	! * * * Read data from file * * *
end subroutine rload_1d_array_from_hdf5


subroutine rload_2d_array_from_hdf5(h5file_id,dset,rdata,attr)
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: rdata
	REAL, DIMENSION(:,:), ALLOCATABLE :: raw_data
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: attr
	CHARACTER(MAX_STRING_LENGTH) :: aname
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: atype_id
	INTEGER(HSIZE_T), DIMENSION(2) :: dims
	INTEGER(HSIZE_T), DIMENSION(2) :: adims
	INTEGER(SIZE_T) :: tmplen
	INTEGER(SIZE_T) :: attrlen
	INTEGER :: h5error

	dims = shape(rdata)

	ALLOCATE( raw_data(dims(1),dims(2)) )

	! * * * Read data from file * * *

	call h5dopen_f(h5file_id, TRIM(dset), dset_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dopen_f")')
	end if

	call h5dread_f(dset_id, H5T_NATIVE_REAL, raw_data, dims, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dread_f")')
	end if
	rdata = REAL(raw_data,rp)

	call h5dclose_f(dset_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dclose_f")')
	end if

	DEALLOCATE( raw_data )

	if (PRESENT(attr)) then
		! * * * Read data attribute(s) from file * * *
	end if

	! * * * Read data from file * * *
end subroutine rload_2d_array_from_hdf5


subroutine rload_3d_array_from_hdf5(h5file_id,dset,rdata,attr)
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: rdata
	REAL, DIMENSION(:,:,:), ALLOCATABLE :: raw_data
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: attr
	CHARACTER(MAX_STRING_LENGTH) :: aname
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: atype_id
	INTEGER(HSIZE_T), DIMENSION(3) :: dims
	INTEGER(HSIZE_T), DIMENSION(3) :: adims
	INTEGER(SIZE_T) :: tmplen
	INTEGER(SIZE_T) :: attrlen
	INTEGER :: h5error

	dims = shape(rdata)

	ALLOCATE( raw_data(dims(1),dims(2),dims(3)) )

	! * * * Read data from file * * *

	call h5dopen_f(h5file_id, TRIM(dset), dset_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dopen_f")')
	end if

	call h5dread_f(dset_id, H5T_NATIVE_REAL, raw_data, dims, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dread_f")')
	end if
	rdata = REAL(raw_data,rp)

	call h5dclose_f(dset_id, h5error)
	if (h5error .EQ. -1) then
		write(6,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dclose_f")')
	end if

	DEALLOCATE( raw_data )

	if (PRESENT(attr)) then
		! * * * Read data attribute(s) from file * * *
	end if

	! * * * Read data from file * * *
end subroutine rload_3d_array_from_hdf5


subroutine i1save_to_hdf5(h5file_id,dset,idata,attr)
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	INTEGER(KIND=1), INTENT(IN) :: idata
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) :: attr
	CHARACTER(4) :: aname = "Info"
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: atype_id
	INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/1/)
	INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
	INTEGER :: rank = 1
	INTEGER :: arank = 1
	INTEGER(SIZE_T) :: tmplen
	INTEGER(SIZE_T) :: attrlen
	INTEGER :: h5error

	! * * * Write data to file * * *
	call h5screate_simple_f(rank,dims,dspace_id,h5error)
	call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_INTEGER, dspace_id, dset_id, h5error)
	call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INT(idata,idef), dims, h5error)

	if (PRESENT(attr)) then
		! * * * Write attribute of data to file * * *
		attrlen = LEN_TRIM(attr)
		call h5screate_simple_f(arank,adims,aspace_id,h5error)
		call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
		call h5tset_size_f(atype_id, attrlen, h5error)
		call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
		call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

		call h5aclose_f(attr_id, h5error)
		call h5sclose_f(aspace_id, h5error)
		! * * * Write attribute of data to file * * *
	end if

	call h5sclose_f(dspace_id, h5error)
	call h5dclose_f(dset_id, h5error)
	! * * * Write data to file * * *
end subroutine i1save_to_hdf5


subroutine i2save_to_hdf5(h5file_id,dset,idata,attr)
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	INTEGER(KIND=2), INTENT(IN) :: idata
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) :: attr
	CHARACTER(4) :: aname = "Info"
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: atype_id
	INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/1/)
	INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
	INTEGER :: rank = 1
	INTEGER :: arank = 1
	INTEGER(SIZE_T) :: tmplen
	INTEGER(SIZE_T) :: attrlen
	INTEGER :: h5error

	! * * * Write data to file * * *
	call h5screate_simple_f(rank,dims,dspace_id,h5error)
	call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_INTEGER, dspace_id, dset_id, h5error)
	call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INT(idata,idef), dims, h5error)

	if (PRESENT(attr)) then
		! * * * Write attribute of data to file * * *
		attrlen = LEN_TRIM(attr)
		call h5screate_simple_f(arank,adims,aspace_id,h5error)
		call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
		call h5tset_size_f(atype_id, attrlen, h5error)
		call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
		call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

		call h5aclose_f(attr_id, h5error)
		call h5sclose_f(aspace_id, h5error)
		! * * * Write attribute of data to file * * *
	end if

	call h5sclose_f(dspace_id, h5error)
	call h5dclose_f(dset_id, h5error)
	! * * * Write data to file * * *
end subroutine i2save_to_hdf5


subroutine i4save_to_hdf5(h5file_id,dset,idata,attr)
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	INTEGER(KIND=4), INTENT(IN) :: idata
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) :: attr
	CHARACTER(4) :: aname = "Info"
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: atype_id
	INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/1/)
	INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
	INTEGER :: rank = 1
	INTEGER :: arank = 1
	INTEGER(SIZE_T) :: tmplen
	INTEGER(SIZE_T) :: attrlen
	INTEGER :: h5error

	! * * * Write data to file * * *
	call h5screate_simple_f(rank,dims,dspace_id,h5error)
	call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_INTEGER, dspace_id, dset_id, h5error)
	call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INT(idata,idef), dims, h5error)

	if (PRESENT(attr)) then
		! * * * Write attribute of data to file * * *
		attrlen = LEN_TRIM(attr)
		call h5screate_simple_f(arank,adims,aspace_id,h5error)
		call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
		call h5tset_size_f(atype_id, attrlen, h5error)
		call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
		call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

		call h5aclose_f(attr_id, h5error)
		call h5sclose_f(aspace_id, h5error)
		! * * * Write attribute of data to file * * *
	end if

	call h5sclose_f(dspace_id, h5error)
	call h5dclose_f(dset_id, h5error)
	! * * * Write data to file * * *
end subroutine i4save_to_hdf5


subroutine i8save_to_hdf5(h5file_id,dset,idata,attr)
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	INTEGER(KIND=8), INTENT(IN) :: idata
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) :: attr
	CHARACTER(4) :: aname = "Info"
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: atype_id
	INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/1/)
	INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
	INTEGER :: rank = 1
	INTEGER :: arank = 1
	INTEGER(SIZE_T) :: tmplen
	INTEGER(SIZE_T) :: attrlen
	INTEGER :: h5error

	! * * * Write data to file * * *
	call h5screate_simple_f(rank,dims,dspace_id,h5error)
	call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_DOUBLE, dspace_id, dset_id, h5error)
	call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, REAL(idata,8), dims, h5error)

	if (PRESENT(attr)) then
		! * * * Write attribute of data to file * * *
		attrlen = LEN_TRIM(attr)
		call h5screate_simple_f(arank,adims,aspace_id,h5error)
		call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
		call h5tset_size_f(atype_id, attrlen, h5error)
		call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
		call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

		call h5aclose_f(attr_id, h5error)
		call h5sclose_f(aspace_id, h5error)
		! * * * Write attribute of data to file * * *
	end if

	call h5sclose_f(dspace_id, h5error)
	call h5dclose_f(dset_id, h5error)
	! * * * Write data to file * * *
end subroutine i8save_to_hdf5


subroutine isave_1d_array_to_hdf5(h5file_id,dset,idata,attr)
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	INTEGER, DIMENSION(:), INTENT(IN) :: idata
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: attr
	CHARACTER(4) :: aname = "Info"
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: atype_id
	INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: dims
	INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: adims
	INTEGER :: rank
	INTEGER :: arank
	INTEGER(SIZE_T) :: tmplen
	INTEGER(SIZE_T) :: attrlen
	INTEGER :: h5error
	INTEGER :: rr,dd ! Iterators

	rank = size(shape(idata))
	ALLOCATE(dims(rank))
	dims = shape(idata)

	! * * * Write data to file * * *
	call h5screate_simple_f(rank,dims,dspace_id,h5error)
	call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_INTEGER, dspace_id, dset_id, h5error)
	call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, idata, dims, h5error)

	if (PRESENT(attr)) then
		arank = size(shape(attr))
		ALLOCATE(adims(arank))
		adims = shape(attr)

		! * * * Write attribute of data to file * * *
		tmplen = 0
		attrlen = 0
		do rr=1_idef,arank
			do dd=1_idef,adims(rr)
				tmplen = LEN_TRIM(attr(dd))
				if ( tmplen .GT. attrlen) then
					attrlen = tmplen
				end if
			end do
		end do

		call h5screate_simple_f(arank,adims,aspace_id,h5error)
		call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
		call h5tset_size_f(atype_id, attrlen, h5error)
		call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
		call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

		call h5aclose_f(attr_id, h5error)
		call h5sclose_f(aspace_id, h5error)
		! * * * Write attribute of data to file * * *

		DEALLOCATE(adims)
	end if

	call h5sclose_f(dspace_id, h5error)
	call h5dclose_f(dset_id, h5error)
	! * * * Write data to file * * *

	DEALLOCATE(dims)
end subroutine isave_1d_array_to_hdf5


subroutine rsave_to_hdf5(h5file_id,dset,rdata,attr)
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	REAL(rp), INTENT(IN) :: rdata
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) :: attr
	CHARACTER(4) :: aname = "Info"
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: atype_id
	INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/1/)
	INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
	INTEGER :: rank = 1
	INTEGER :: arank = 1
	INTEGER(SIZE_T) :: tmplen
	INTEGER(SIZE_T) :: attrlen
	INTEGER :: h5error

	! * * * Write data to file * * *

	call h5screate_simple_f(rank,dims,dspace_id,h5error)
	call h5dcreate_f(h5file_id, TRIM(dset), KORC_HDF5_REAL, dspace_id, dset_id, h5error)

	if (rp .EQ. INT(rp_hdf5)) then
		call h5dwrite_f(dset_id, KORC_HDF5_REAL, rdata, dims, h5error)
	else
		call h5dwrite_f(dset_id, KORC_HDF5_REAL, REAL(rdata,4), dims, h5error)
	end if

	if (PRESENT(attr)) then
		! * * * Write attribute of data to file * * *
		attrlen = LEN_TRIM(attr)
		call h5screate_simple_f(arank,adims,aspace_id,h5error)
		call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
		call h5tset_size_f(atype_id, attrlen, h5error)
		call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
		call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

		call h5aclose_f(attr_id, h5error)
		call h5sclose_f(aspace_id, h5error)
		! * * * Write attribute of data to file * * *
	end if

	call h5sclose_f(dspace_id, h5error)
	call h5dclose_f(dset_id, h5error)
	! * * * Write data to file * * *
end subroutine rsave_to_hdf5


subroutine rsave_1d_array_to_hdf5(h5file_id,dset,rdata,attr)
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	REAL(rp), DIMENSION(:), INTENT(IN) :: rdata
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: attr
	CHARACTER(4) :: aname = "Info"
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: atype_id
	INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: dims
	INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: adims
	INTEGER :: rank
	INTEGER :: arank
	INTEGER(SIZE_T) :: tmplen
	INTEGER(SIZE_T) :: attrlen
	INTEGER :: h5error
	INTEGER :: rr,dd ! Iterators

	rank = size(shape(rdata))
	ALLOCATE(dims(rank))
	dims = shape(rdata)

	! * * * Write data to file * * *

	call h5screate_simple_f(rank,dims,dspace_id,h5error)
	call h5dcreate_f(h5file_id, TRIM(dset), KORC_HDF5_REAL, dspace_id, dset_id, h5error)

	if (rp .EQ. INT(rp_hdf5)) then
		call h5dwrite_f(dset_id, KORC_HDF5_REAL, rdata, dims, h5error)
	else
		call h5dwrite_f(dset_id, KORC_HDF5_REAL, REAL(rdata,4), dims, h5error)
	end if

	if (PRESENT(attr)) then
		arank = size(shape(attr))
		ALLOCATE(adims(arank))
		adims = shape(attr)

		! * * * Write attribute of data to file * * *
		tmplen = 0
		attrlen = 0
		do rr=1_idef,arank
			do dd=1_idef,adims(rr)
				tmplen = LEN_TRIM(attr(dd))
				if ( tmplen .GT. attrlen) then
					attrlen = tmplen
				end if
			end do
		end do

		call h5screate_simple_f(arank,adims,aspace_id,h5error)
		call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
		call h5tset_size_f(atype_id, attrlen, h5error)
		call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
		call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

		call h5aclose_f(attr_id, h5error)
		call h5sclose_f(aspace_id, h5error)
		! * * * Write attribute of data to file * * *

		DEALLOCATE(adims)
	end if

	call h5sclose_f(dspace_id, h5error)
	call h5dclose_f(dset_id, h5error)
	! * * * Write data to file * * *

	DEALLOCATE(dims)
end subroutine rsave_1d_array_to_hdf5


subroutine rsave_2d_array_to_hdf5(h5file_id,dset,rdata,attr)
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	REAL(rp), DIMENSION(:,:), INTENT(IN) :: rdata
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: attr
	CHARACTER(4) :: aname = "Info"
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: atype_id
	INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: dims
	INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: adims
	INTEGER :: rank
	INTEGER :: arank
	INTEGER(SIZE_T) :: tmplen
	INTEGER(SIZE_T) :: attrlen
	INTEGER :: h5error
	INTEGER :: rr,dd ! Iterators

	rank = size(shape(rdata))
	ALLOCATE(dims(rank))
	dims = shape(rdata)

	! * * * Write data to file * * *

	call h5screate_simple_f(rank,dims,dspace_id,h5error)
	call h5dcreate_f(h5file_id, TRIM(dset), KORC_HDF5_REAL, dspace_id, dset_id, h5error)

	if (rp .EQ. INT(rp_hdf5)) then
		call h5dwrite_f(dset_id, KORC_HDF5_REAL, rdata, dims, h5error)
	else
		call h5dwrite_f(dset_id, KORC_HDF5_REAL, REAL(rdata,4), dims, h5error)
	end if

	if (PRESENT(attr)) then
		! * * * Write attribute of data to file * * *
	end if

	call h5sclose_f(dspace_id, h5error)
	call h5dclose_f(dset_id, h5error)
	! * * * Write data to file * * *

	DEALLOCATE(dims)
end subroutine rsave_2d_array_to_hdf5


subroutine rsave_3d_array_to_hdf5(h5file_id,dset,rdata,attr)
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: rdata
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: attr
	CHARACTER(4) :: aname = "Info"
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HID_T) :: aspace_id
	INTEGER(HID_T) :: attr_id
	INTEGER(HID_T) :: atype_id
	INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: dims
	INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: adims
	INTEGER :: rank
	INTEGER :: arank
	INTEGER(SIZE_T) :: tmplen
	INTEGER(SIZE_T) :: attrlen
	INTEGER :: h5error
	INTEGER :: rr,dd ! Iterators

	rank = size(shape(rdata))
	ALLOCATE(dims(rank))
	dims = shape(rdata)

	! * * * Write data to file * * *

	call h5screate_simple_f(rank,dims,dspace_id,h5error)
	call h5dcreate_f(h5file_id, TRIM(dset), KORC_HDF5_REAL, dspace_id, dset_id, h5error)

	if (rp .EQ. INT(rp_hdf5)) then
		call h5dwrite_f(dset_id, KORC_HDF5_REAL, rdata, dims, h5error)
	else
		call h5dwrite_f(dset_id, KORC_HDF5_REAL, REAL(rdata,4), dims, h5error)
	end if

	if (PRESENT(attr)) then
		! * * * Write attribute of data to file * * *
	end if

	call h5sclose_f(dspace_id, h5error)
	call h5dclose_f(dset_id, h5error)
	! * * * Write data to file * * *

	DEALLOCATE(dims)
end subroutine rsave_3d_array_to_hdf5


subroutine save_string_parameter(h5file_id,dset,string_array)
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), INTENT(IN) :: string_array
	INTEGER(HID_T) :: dset_id
	INTEGER(HID_T) :: dspace_id
	INTEGER(HSIZE_T), DIMENSION(1) :: dims
	INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
	INTEGER(SIZE_T), DIMENSION(:), ALLOCATABLE :: str_len
	INTEGER(HID_T) :: string_type
	INTEGER :: h5error
	
	ALLOCATE(str_len(SIZE(string_array)))

	dims = (/SIZE(string_array)/)
	data_dims = (/MAX_STRING_LENGTH,SIZE(string_array)/)
	str_len = (/LEN_TRIM(string_array)/)

	call h5tcopy_f(H5T_STRING,string_type,h5error)
	call h5tset_strpad_f(string_type,H5T_STR_SPACEPAD_F,h5error)

	call h5screate_simple_f(1,dims,dspace_id,h5error)

	call h5dcreate_f(h5file_id,TRIM(dset),string_type,dspace_id,dset_id,h5error)

	call h5dwrite_vl_f(dset_id,string_type,string_array,data_dims,str_len,h5error,dspace_id)

	call h5sclose_f(dspace_id, h5error)
	call h5dclose_f(dset_id, h5error)

	DEALLOCATE(str_len)
end subroutine save_string_parameter


subroutine save_simulation_parameters(params,spp,F,P)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp
	TYPE(FIELDS), INTENT(IN) :: F
	TYPE(PROFILES), INTENT(IN) :: P
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH) :: dset
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: dims
	REAL(rp), DIMENSION(:), ALLOCATABLE :: rdata
	INTEGER, DIMENSION(:), ALLOCATABLE :: idata
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
	CHARACTER(MAX_STRING_LENGTH) :: attr
	INTEGER :: h5error
	CHARACTER(19) :: tmp_str
	REAL(rp) :: units

	! * * * Error handling * * * !
	call h5eset_auto_f(params%HDF5_error_handling, h5error) ! Turn off: 0_idef. Turn on: 1_idef

	if (.NOT.params%restart) then

	if (SIZE(params%outputs_list).GT.1_idef) then
		write(tmp_str,'(I18)') params%mpi_params%rank
		filename = TRIM(params%path_to_outputs) // "file_" // TRIM(ADJUSTL(tmp_str)) // ".h5"
		call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)
		call h5fclose_f(h5file_id, h5error)
	end if

	if (params%mpi_params%rank .EQ. 0) then
		filename = TRIM(params%path_to_outputs) // "simulation_parameters.h5"

		call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

		! Simulation parameters group
		gname = "simulation"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

		ALLOCATE(attr_array(1))		
		ALLOCATE(idata(1))

		dset = TRIM(gname) // "/plasma_model"
		call save_string_parameter(h5file_id,dset,(/params%plasma_model/))

		dset = TRIM(gname) // "/simulation_time"
        attr = "Total aimed simulation time in seconds"
		call save_to_hdf5(h5file_id,dset,params%simulation_time*params%cpp%time,attr)

		dset = TRIM(gname) // "/snapshot_frequency"
        attr = "Time between snapshots in seconds"
		call save_to_hdf5(h5file_id,dset,params%snapshot_frequency*params%cpp%time,attr)

		dset = TRIM(gname) // "/dt"
        attr = "Time step in secs"
		call save_to_hdf5(h5file_id,dset,params%dt*params%cpp%time,attr)

		dset = TRIM(gname) // "/t_steps"
		attr_array(1) = "Number of time steps"
		idata = params%t_steps
		call save_1d_array_to_hdf5(h5file_id,dset,idata,attr_array)

		dset = TRIM(gname) // "/num_omp_threads"
		attr = "Number of omp threads"
		call save_to_hdf5(h5file_id,dset, params%num_omp_threads,attr)

		dset = TRIM(gname) // "/output_cadence"
		attr_array(1) = "Cadence of output files"
		idata = params%output_cadence
		call save_1d_array_to_hdf5(h5file_id,dset,idata,attr_array)

		dset = TRIM(gname) // "/HDF5_error_handling"
		attr_array(1) = "Error handling option: 0=OFF, 1=ON"
		idata = params%HDF5_error_handling
		call save_1d_array_to_hdf5(h5file_id,dset,idata,attr_array)

		dset = TRIM(gname) // "/restart_output_cadence"
		attr_array(1) = "Cadence of output files"
		idata = params%restart_output_cadence
		call save_1d_array_to_hdf5(h5file_id,dset,idata,attr_array)

		dset = TRIM(gname) // "/num_snapshots"
		attr_array(1) = "Number of outputs for each variable"
		idata = params%num_snapshots
		call save_1d_array_to_hdf5(h5file_id,dset,idata,attr_array)

		dset = TRIM(gname) // "/num_species"
		attr = "Number of particle species"
		call save_to_hdf5(h5file_id,dset,params%num_species,attr)

		dset = TRIM(gname) // "/nmpi"
		attr = "Number of mpi processes"
		call save_to_hdf5(h5file_id,dset,params%mpi_params%nmpi,attr)

		dset = TRIM(gname) // "/minimum_particle_energy"
        attr = "Minimum energy of simulated particles in eV"
		call save_to_hdf5(h5file_id,dset,params%minimum_particle_energy*params%cpp%energy/C_E,attr)

		dset = TRIM(gname) // "/radiation"
		attr = "Radiation losses included in simulation"
		if(params%radiation) then
			call save_to_hdf5(h5file_id,dset,1_idef,attr)
		else
			call save_to_hdf5(h5file_id,dset,0_idef,attr)
		end if

		dset = TRIM(gname) // "/collisions"
		attr = "Radiation losses included in simulation"
		if(params%collisions) then
			call save_to_hdf5(h5file_id,dset,1_idef,attr)
		else
			call save_to_hdf5(h5file_id,dset,0_idef,attr)
		end if

		dset = TRIM(gname) // "/outputs_list"
		call save_string_parameter(h5file_id,dset,params%outputs_list)

		DEALLOCATE(idata)
		DEALLOCATE(attr_array)

		call h5gclose_f(group_id, h5error)


		! Plasma species group
		gname = "species"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

		ALLOCATE(attr_array(params%num_species))

		dset = TRIM(gname) // "/spatial_distribution"
		call save_string_parameter(h5file_id,dset,spp%spatial_distribution)

		dset = TRIM(gname) // "/energy_distribution"
		call save_string_parameter(h5file_id,dset,spp%energy_distribution)

		dset = TRIM(gname) // "/pitch_distribution"
		call save_string_parameter(h5file_id,dset,spp%pitch_distribution)

		dset = TRIM(gname) // "/ppp"
		attr_array(1) = "Particles per (mpi) process"
		call save_1d_array_to_hdf5(h5file_id,dset,spp%ppp,attr_array)

		dset = TRIM(gname) // "/q"
		attr_array(1) = "Electric charge"
		call save_1d_array_to_hdf5(h5file_id,dset,spp%q*params%cpp%charge,attr_array)

		dset = TRIM(gname) // "/m"
		attr_array(1) = "Species mass in kg"
		call save_1d_array_to_hdf5(h5file_id,dset,spp%m*params%cpp%mass,attr_array)

		dset = TRIM(gname) // "/Eo"
		attr_array(1) = "Initial (average) energy in eV"
		call save_1d_array_to_hdf5(h5file_id,dset,spp%Eo*params%cpp%energy/C_E,attr_array)

		dset = TRIM(gname) // "/go"
		attr_array(1) = "Initial relativistic g factor."
		call save_1d_array_to_hdf5(h5file_id,dset,spp%go,attr_array)

		dset = TRIM(gname) // "/etao"
		attr_array(1) = "Initial pitch angle in degrees"
		call save_1d_array_to_hdf5(h5file_id,dset,spp%etao,attr_array)

		dset = TRIM(gname) // "/wc"
		attr_array(1) = "Average relativistic cyclotron frequency in Hz"
		call save_1d_array_to_hdf5(h5file_id,dset,spp%wc_r/params%cpp%time,attr_array)

		dset = TRIM(gname) // "/Ro"
		attr_array(1) = "Initial radial position of population"
		call save_1d_array_to_hdf5(h5file_id,dset,spp%Ro*params%cpp%length,attr_array)

		dset = TRIM(gname) // "/PHIo"
		attr_array(1) = "Azimuthal angle in degrees."
		call save_1d_array_to_hdf5(h5file_id,dset,spp%PHIo*180.0_rp/C_PI,attr_array)

		dset = TRIM(gname) // "/Zo"
		attr_array(1) = "Initial Z position of population"
		call save_1d_array_to_hdf5(h5file_id,dset,spp%Zo*params%cpp%length,attr_array)

		dset = TRIM(gname) // "/ri"
		attr_array(1) = "Inner radius of initial spatial distribution"
		call save_1d_array_to_hdf5(h5file_id,dset,spp%r_inner*params%cpp%length,attr_array)

		dset = TRIM(gname) // "/ro"
		attr_array(1) = "Outter radius of initial spatial distribution"
		call save_1d_array_to_hdf5(h5file_id,dset,spp%r_outter*params%cpp%length,attr_array)

		dset = TRIM(gname) // "/falloff_rate"
		attr_array(1) = "Falloff of gaussian or exponential radial profile in m"
		call save_1d_array_to_hdf5(h5file_id,dset,spp%falloff_rate*params%cpp%length,attr_array)

		call h5gclose_f(group_id, h5error)

		DEALLOCATE(attr_array)


		! Plasma profiles group

		gname = "profiles"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

		dset = TRIM(gname) // "/density_profile"
		call save_string_parameter(h5file_id,dset,(/P%ne_profile/))

		dset = TRIM(gname) // "/temperature_profile"
		call save_string_parameter(h5file_id,dset,(/P%Te_profile/))

		dset = TRIM(gname) // "/Zeff_profile"
		call save_string_parameter(h5file_id,dset,(/P%Zeff_profile/))

		dset = TRIM(gname) // "/neo"
		attr = "Density at the magnetic axis (m^-3)"
		call save_to_hdf5(h5file_id,dset,P%neo*params%cpp%density,attr)

		dset = TRIM(gname) // "/Teo"
		attr = "Temperature at the magnetic axis (eV)"
		call save_to_hdf5(h5file_id,dset,P%Teo*params%cpp%temperature/C_E,attr)

		dset = TRIM(gname) // "/Zeffo"
		attr = "Zeff at the magnetic axis"
		call save_to_hdf5(h5file_id,dset,P%Zeffo,attr)

		if (TRIM(params%plasma_model) .EQ. 'ANALYTICAL') then
			dset = TRIM(gname) // "/n_ne"
			attr = "Exponent of tanh(x)^n for density profile"
			call save_to_hdf5(h5file_id,dset,P%n_ne,attr)

			dset = TRIM(gname) // "/a_ne"
			attr = "Coefficients f=ao+a1*r+a2*r^2+a3*r^3. a_ne=[a0,a1,a2,a3]"
			call save_1d_array_to_hdf5(h5file_id,dset,P%a_ne)

			dset = TRIM(gname) // "/n_Te"
			attr = "Exponent of tanh(x)^n for density profile"
			call save_to_hdf5(h5file_id,dset,P%n_Te,attr)

			dset = TRIM(gname) // "/a_Te"
			attr = "Coefficients f=ao+a1*r+a2*r^2+a3*r^3. a_Te=[a0,a1,a2,a3]"
			call save_1d_array_to_hdf5(h5file_id,dset,P%a_Te)

			dset = TRIM(gname) // "/n_Zeff"
			attr = "Exponent of tanh(x)^n for Zeff profile"
			call save_to_hdf5(h5file_id,dset,P%n_Zeff,attr)

			dset = TRIM(gname) // "/a_Zeff"
			attr = "Coefficients f=ao+a1*r+a2*r^2+a3*r^3. a_Zeff=[a0,a1,a2,a3]"
			call save_1d_array_to_hdf5(h5file_id,dset,P%a_Zeff)
		else if (params%plasma_model .EQ. 'EXTERNAL') then
			ALLOCATE(attr_array(1))
			dset = TRIM(gname) // "/dims"
			attr_array(1) = "Mesh dimension of the profiles (NR,NPHI,NZ)"
			call save_1d_array_to_hdf5(h5file_id,dset,P%dims,attr_array)		

			dset = TRIM(gname) // "/R"
			attr_array(1) = "Grid nodes of profiles along the radial position"
			call save_1d_array_to_hdf5(h5file_id,dset,P%X%R*params%cpp%length,attr_array)

			if (ALLOCATED(F%X%PHI)) then
				dset = TRIM(gname) // "/PHI"
			attr_array(1) = "Grid nodes of profiles along the azimuthal position"
				call save_1d_array_to_hdf5(h5file_id,dset,P%X%PHI,attr_array)
			end if

			dset = TRIM(gname) // "/Z"
			attr_array(1) = "Grid nodes of profiles along the Z position"
			call save_1d_array_to_hdf5(h5file_id,dset,P%X%Z*params%cpp%length,attr_array)

			DEALLOCATE(attr_array)
		else if (params%plasma_model .EQ. 'UNIFORM') then
			! Something
		end if

		call h5gclose_f(group_id, h5error)


		! Electromagnetic fields group

		gname = "fields"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

		if (TRIM(params%plasma_model) .EQ. 'ANALYTICAL') then
			dset = TRIM(gname) // "/Bo"
			attr = "Toroidal field at the magnetic axis in T"
			call save_to_hdf5(h5file_id,dset,F%Bo*params%cpp%Bo,attr)

			dset = TRIM(gname) // "/current_direction"
			call save_string_parameter(h5file_id,dset,(/F%AB%current_direction/))

			dset = TRIM(gname) // "/a"
			attr = "Minor radius in m"
        	call save_to_hdf5(h5file_id,dset,F%AB%a*params%cpp%length,attr)

			dset = TRIM(gname) // "/Ro"
			attr = "Major radius in m"
        	call save_to_hdf5(h5file_id,dset,F%Ro*params%cpp%length,attr)

			dset = TRIM(gname) // "/qa"
			attr = "Safety factor at minor radius"
        	call save_to_hdf5(h5file_id,dset,F%AB%qa,attr)

			dset = TRIM(gname) // "/qo"
			attr = "Safety factor at the magnetic axis"
        	call save_to_hdf5(h5file_id,dset,F%AB%qo,attr)

			dset = TRIM(gname) // "/lambda"
			attr = "Parameter lamda in m"
        	call save_to_hdf5(h5file_id,dset,F%AB%lambda*params%cpp%length,attr)

			dset = TRIM(gname) // "/Bpo"
			attr = "Poloidal magnetic field in T"
        	call save_to_hdf5(h5file_id,dset,F%AB%Bpo*params%cpp%Bo,attr)

			dset = TRIM(gname) // "/Eo"
			attr = "Electric field at the magnetic axis in V/m"
			call save_to_hdf5(h5file_id,dset,F%Eo*params%cpp%Eo,attr)
		else if (params%plasma_model .EQ. 'EXTERNAL') then
			ALLOCATE(attr_array(1))
			dset = TRIM(gname) // "/dims"
			attr_array(1) = "Mesh dimension of the magnetic field (NR,NPHI,NZ)"
			call save_1d_array_to_hdf5(h5file_id,dset,F%dims,attr_array)		

			dset = TRIM(gname) // "/R"
			attr_array(1) = "Radial position of the magnetic field grid nodes"
			call save_1d_array_to_hdf5(h5file_id,dset,F%X%R*params%cpp%length,attr_array)

			if (ALLOCATED(F%X%PHI)) then
				dset = TRIM(gname) // "/PHI"
				attr_array(1) = "Azimuthal angle of the magnetic field grid nodes"
				call save_1d_array_to_hdf5(h5file_id,dset,F%X%PHI,attr_array)
			end if

			dset = TRIM(gname) // "/Z"
			attr_array(1) = "Z position of the magnetic field grid nodes"
			call save_1d_array_to_hdf5(h5file_id,dset,F%X%Z*params%cpp%length,attr_array)

			dset = TRIM(gname) // "/Bo"
			attr = "Toroidal field at the magnetic axis in T"
			call save_to_hdf5(h5file_id,dset,F%Bo*params%cpp%Bo,attr)

			dset = TRIM(gname) // "/Eo"
			attr = "Electric field at the magnetic axis in V/m"
			call save_to_hdf5(h5file_id,dset,F%Eo*params%cpp%Eo,attr)

			dset = TRIM(gname) // "/Ro"
			attr = "Radial position of magnetic axis"
			call save_to_hdf5(h5file_id,dset,F%Ro*params%cpp%length,attr)

			dset = TRIM(gname) // "/Zo"
			attr = "Radial position of magnetic axis"
			call save_to_hdf5(h5file_id,dset,F%Zo*params%cpp%length,attr)

			DEALLOCATE(attr_array)
		else if (params%plasma_model .EQ. 'UNIFORM') then
			dset = TRIM(gname) // "/Bo"
			attr = "Magnetic field in T"
			call save_to_hdf5(h5file_id,dset,F%Bo*params%cpp%Bo,attr)

			dset = TRIM(gname) // "/Eo"
			attr = "Electric field in V/m"
			call save_to_hdf5(h5file_id,dset,F%Eo*params%cpp%Eo,attr)
		end if

		call h5gclose_f(group_id, h5error)


		! Characteristic scales
		gname = "scales"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

		dset = TRIM(gname) // "/t"
		attr = "Characteristic time in secs"
		call save_to_hdf5(h5file_id,dset,params%cpp%time,attr)

		dset = TRIM(gname) // "/m"
		attr = "Characteristic mass in kg"
		call save_to_hdf5(h5file_id,dset,params%cpp%mass,attr)

		dset = TRIM(gname) // "/q"
		attr = "Characteristic charge in Coulombs"
		call save_to_hdf5(h5file_id,dset,params%cpp%charge,attr)

		dset = TRIM(gname) // "/l"
		attr = "Characteristic length in m"
		call save_to_hdf5(h5file_id,dset,params%cpp%length,attr)

		dset = TRIM(gname) // "/v"
		attr = "Characteristic velocity in m"
		call save_to_hdf5(h5file_id,dset,params%cpp%velocity,attr)

		dset = TRIM(gname) // "/K"
		attr = "Characteristic kinetic energy in J"
		call save_to_hdf5(h5file_id,dset,params%cpp%energy,attr)

		dset = TRIM(gname) // "/n"
		attr = "Characteristic plasma density in m^-3"
		call save_to_hdf5(h5file_id,dset,params%cpp%density,attr)

		dset = TRIM(gname) // "/E"
		attr = "Characteristic electric field in V/m"
		call save_to_hdf5(h5file_id,dset,params%cpp%Eo,attr)

		dset = TRIM(gname) // "/B"
		attr = "Characteristic magnetic field in T"
		call save_to_hdf5(h5file_id,dset,params%cpp%Bo,attr)

		dset = TRIM(gname) // "/P"
		attr = "Characteristic pressure in Pa"
		call save_to_hdf5(h5file_id,dset,params%cpp%pressure,attr)

		dset = TRIM(gname) // "/T"
		attr = "Characteristic plasma temperature in J"
		call save_to_hdf5(h5file_id,dset,params%cpp%temperature,attr)

		call h5gclose_f(group_id, h5error)

		call h5fclose_f(h5file_id, h5error)
	end if

	end if
end subroutine save_simulation_parameters


subroutine save_simulation_outputs(params,spp,F)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp
	TYPE(FIELDS), INTENT(IN) :: F
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH) :: subgname
	CHARACTER(MAX_STRING_LENGTH) :: dset
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	INTEGER(HID_T) :: subgroup_id
	INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: dims
	REAL(rp), DIMENSION(:), ALLOCATABLE :: rdata
	INTEGER, DIMENSION(:), ALLOCATABLE :: idata
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
	CHARACTER(MAX_STRING_LENGTH) :: attr
	INTEGER :: h5error
	CHARACTER(19) :: tmp_str
	REAL(rp) :: units
    INTEGER :: ss,jj
	LOGICAL :: object_exists

	if (params%mpi_params%rank .EQ. 0) then
		write(6,'("Saving snapshot: ",I15)') params%it/params%output_cadence
	end if

	if (SIZE(params%outputs_list).GT.1_idef) then
		write(tmp_str,'(I18)') params%mpi_params%rank
		filename = TRIM(params%path_to_outputs) // "file_" // TRIM(ADJUSTL(tmp_str)) // ".h5"
		call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, h5file_id, h5error)

		! Create group 'it'
		write(tmp_str,'(I18)') params%it
		gname = TRIM(ADJUSTL(tmp_str))
		call h5lexists_f(h5file_id,TRIM(gname),object_exists,h5error)

		if (.NOT.object_exists) then ! Check if group does exist.
			call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)
	
			dset = TRIM(gname) // "/time"
			attr = "Simulation time in secs"
			call save_to_hdf5(h5file_id,dset,REAL(params%it,rp)*params%dt*params%cpp%time,attr)

			do ss=1_idef,params%num_species
				write(tmp_str,'(I18)') ss
				subgname = "spp_" // TRIM(ADJUSTL(tmp_str))
				call h5gcreate_f(group_id, TRIM(subgname), subgroup_id, h5error)

				do jj=1_idef,SIZE(params%outputs_list)
					SELECT CASE (TRIM(params%outputs_list(jj)))
						CASE ('X')
							dset = "X"
							units = params%cpp%length
							call rsave_2d_array_to_hdf5(subgroup_id, dset, units*spp(ss)%vars%X)
						CASE('V')
							dset = "V"
							units = params%cpp%velocity
							call rsave_2d_array_to_hdf5(subgroup_id, dset, units*spp(ss)%vars%V)
						CASE('Rgc')
							dset = "Rgc"
							units = params%cpp%length
							call rsave_2d_array_to_hdf5(subgroup_id, dset, units*spp(ss)%vars%Rgc)
						CASE('g')
							dset = "g"
							call save_1d_array_to_hdf5(subgroup_id, dset, spp(ss)%vars%g)
						CASE('eta')
							dset = "eta"
							call save_1d_array_to_hdf5(subgroup_id, dset, spp(ss)%vars%eta)
						CASE('mu')
							dset = "mu"
							units = params%cpp%mass*params%cpp%velocity**2/params%cpp%Bo
							call save_1d_array_to_hdf5(subgroup_id, dset, units*spp(ss)%vars%mu)
						CASE('Prad')
							dset = "Prad"
							units = params%cpp%mass*(params%cpp%velocity**3)/params%cpp%length
							call save_1d_array_to_hdf5(subgroup_id, dset, units*spp(ss)%vars%Prad)
						CASE('Pin')
							dset = "Pin"
							units = params%cpp%mass*(params%cpp%velocity**3)/params%cpp%length
							call save_1d_array_to_hdf5(subgroup_id, dset, units*spp(ss)%vars%Pin)
						CASE('flag')
							dset = "flag"
							call save_1d_array_to_hdf5(subgroup_id,dset, INT(spp(ss)%vars%flag,idef))
						CASE('B')
							dset = "B"
							units = params%cpp%Bo
							call rsave_2d_array_to_hdf5(subgroup_id, dset, units*spp(ss)%vars%B)
						CASE('E')
							dset = "E"
							units = params%cpp%Eo
							call rsave_2d_array_to_hdf5(subgroup_id, dset, units*spp(ss)%vars%E)
						CASE('AUX')
							dset = "AUX"
							call save_1d_array_to_hdf5(subgroup_id, dset, spp(ss)%vars%AUX)
						CASE ('ne')
							dset = "ne"
							units = params%cpp%density
							call save_1d_array_to_hdf5(subgroup_id, dset, units*spp(ss)%vars%ne)
						CASE ('Te')
							dset = "Te"
							units = params%cpp%temperature
							call save_1d_array_to_hdf5(subgroup_id, dset, units*spp(ss)%vars%Te/C_E)
						CASE ('Zeff')
							dset = "Zeff"
							call save_1d_array_to_hdf5(subgroup_id, dset, spp(ss)%vars%Zeff)
						CASE DEFAULT
				
					END SELECT
				end do 

				call h5gclose_f(subgroup_id, h5error)
			end do

			call h5gclose_f(group_id, h5error)
		end if ! Check if group does exist.

		call h5fclose_f(h5file_id, h5error)
	end if
end subroutine save_simulation_outputs


subroutine save_restart_variables(params,spp,F)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp
	TYPE(FIELDS), INTENT(IN) :: F
    REAL(rp), DIMENSION(:), ALLOCATABLE :: send_buffer_rp, receive_buffer_rp
    INTEGER(is), DIMENSION(:), ALLOCATABLE :: send_buffer_is, receive_buffer_is
    REAL(rp), DIMENSION(:,:), ALLOCATABLE :: X
    REAL(rp), DIMENSION(:,:), ALLOCATABLE :: V
	INTEGER(is), DIMENSION(:), ALLOCATABLE :: flag
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH) :: subgname
	CHARACTER(MAX_STRING_LENGTH) :: dset
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	INTEGER(HID_T) :: subgroup_id
	INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: dims
	REAL(rp), DIMENSION(:), ALLOCATABLE :: rdata
	INTEGER, DIMENSION(:), ALLOCATABLE :: idata
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
	CHARACTER(MAX_STRING_LENGTH) :: attr
	INTEGER :: h5error
	CHARACTER(19) :: tmp_str
	REAL(rp) :: units
    INTEGER :: ss,jj
	INTEGER :: mpierr
	INTEGER :: numel_send, numel_receive


	if ( MODULO(params%it,params%restart_output_cadence) .EQ. 0_ip ) then
		if (params%mpi_params%rank.EQ.0_idef) then
			filename = TRIM(params%path_to_outputs) // "restart_file.h5"
			call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

			dset = "it"
			attr = "Iteration"
			call save_to_hdf5(h5file_id,dset,params%it,attr)
		
			dset = "time"
			attr = "Simulation time in secs"
			call save_to_hdf5(h5file_id,dset,REAL(params%it,rp)*params%dt*params%cpp%time,attr)
		end if

		do ss=1_idef,params%num_species
			numel_send = 3_idef*spp(ss)%ppp
			numel_receive = 3_idef*spp(ss)%ppp*params%mpi_params%nmpi

			if (params%mpi_params%rank.EQ.0_idef) then
				ALLOCATE(X(3,spp(ss)%ppp*params%mpi_params%nmpi))
				ALLOCATE(V(3,spp(ss)%ppp*params%mpi_params%nmpi))
				ALLOCATE(flag(spp(ss)%ppp*params%mpi_params%nmpi))
			end if

			ALLOCATE(send_buffer_rp(numel_send))
			ALLOCATE(receive_buffer_rp(numel_receive))

			send_buffer_rp = RESHAPE(spp(ss)%vars%X,(/numel_send/))
			receive_buffer_rp = 0.0_rp
			CALL MPI_GATHER(send_buffer_rp,numel_send,MPI_REAL8,receive_buffer_rp,numel_send,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
			if (params%mpi_params%rank.EQ.0_idef) then
				X = RESHAPE(receive_buffer_rp,(/3,spp(ss)%ppp*params%mpi_params%nmpi/))
			end if

			send_buffer_rp = RESHAPE(spp(ss)%vars%V,(/numel_send/))
			receive_buffer_rp = 0.0_rp
			CALL MPI_GATHER(send_buffer_rp,numel_send,MPI_REAL8,receive_buffer_rp,numel_send,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
			if (params%mpi_params%rank.EQ.0_idef) then
				V = RESHAPE(receive_buffer_rp,(/3,spp(ss)%ppp*params%mpi_params%nmpi/))
			end if

			DEALLOCATE(send_buffer_rp)
			DEALLOCATE(receive_buffer_rp)

			numel_send = spp(ss)%ppp
			numel_receive = spp(ss)%ppp*params%mpi_params%nmpi

			ALLOCATE(send_buffer_is(numel_send))
			ALLOCATE(receive_buffer_is(numel_receive))

			send_buffer_is = spp(ss)%vars%flag
			receive_buffer_is = 0_is
			CALL MPI_GATHER(send_buffer_is,numel_send,MPI_INTEGER1,receive_buffer_is,numel_send,&
							MPI_INTEGER1,0,MPI_COMM_WORLD,mpierr)
			if (params%mpi_params%rank.EQ.0_idef) then
				flag = receive_buffer_is
			end if

			DEALLOCATE(send_buffer_is)
			DEALLOCATE(receive_buffer_is)

			if (params%mpi_params%rank.EQ.0_idef) then
				write(tmp_str,'(I18)') ss
				subgname = "spp_" // TRIM(ADJUSTL(tmp_str))
				call h5gcreate_f(h5file_id, TRIM(subgname), group_id, h5error)

				dset = "X"
				call rsave_2d_array_to_hdf5(group_id, dset, X)

				dset = "V"
				call rsave_2d_array_to_hdf5(group_id, dset, V)

				dset = "flag"
				call save_1d_array_to_hdf5(group_id,dset, INT(flag,idef))

				call h5gclose_f(group_id, h5error)
			end if

			if (params%mpi_params%rank.EQ.0_idef) then
				DEALLOCATE(X)
				DEALLOCATE(V)
				DEALLOCATE(flag)
			end if
		end do

		if (params%mpi_params%rank.EQ.0_idef) then
			call h5fclose_f(h5file_id, h5error)
		end if

	end if
end subroutine save_restart_variables

! * * * * * * * * * * * * * * * * * * * * * * * * * !
! * * * SUBROUTINES FOR RESTARTING SIMULATION * * * !
! * * * * * * * * * * * * * * * * * * * * * * * * * !

subroutine get_last_iteration(params)
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: dset
	INTEGER(HID_T) :: h5file_id
	REAL(KIND=8) :: it
	CHARACTER(19) :: tmp_str
	INTEGER :: h5error
	INTEGER :: mpierr
	INTEGER :: ss

	if (params%mpi_params%rank.EQ.0_idef) then
		filename = TRIM(params%path_to_outputs) // "restart_file.h5"
		call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
		if (h5error .EQ. -1) then
			write(6,'("KORC ERROR: Something went wrong in: load_particles_ic --> h5fopen_f")')
		end if

		dset = "/it"
		call load_from_hdf5(h5file_id,dset,it)

		params%ito = INT(it,ip) + 1_ip

		call h5fclose_f(h5file_id, h5error)
	end if


	CALL MPI_BCAST(params%ito,1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
end subroutine get_last_iteration


subroutine load_particles_ic(params,spp)
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
    REAL(rp), DIMENSION(:), ALLOCATABLE :: X_send_buffer
    REAL(rp), DIMENSION(:), ALLOCATABLE :: X_receive_buffer
    REAL(rp), DIMENSION(:), ALLOCATABLE :: V_send_buffer
    REAL(rp), DIMENSION(:), ALLOCATABLE :: V_receive_buffer
    REAL(rp), DIMENSION(:), ALLOCATABLE :: AUX_send_buffer
    REAL(rp), DIMENSION(:), ALLOCATABLE :: AUX_receive_buffer
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: dset
	INTEGER(HID_T) :: h5file_id
	CHARACTER(19) :: tmp_str
	INTEGER :: h5error
	INTEGER :: mpierr
	INTEGER :: ss

	do ss=1_idef,params%num_species
		ALLOCATE(X_send_buffer(3*spp(ss)%ppp*params%mpi_params%nmpi))
		ALLOCATE(X_receive_buffer(3*spp(ss)%ppp))

		ALLOCATE(V_send_buffer(3*spp(ss)%ppp*params%mpi_params%nmpi))
		ALLOCATE(V_receive_buffer(3*spp(ss)%ppp))

		ALLOCATE(AUX_send_buffer(spp(ss)%ppp*params%mpi_params%nmpi))
		ALLOCATE(AUX_receive_buffer(spp(ss)%ppp))

		if (params%mpi_params%rank.EQ.0_idef) then
			filename = TRIM(params%path_to_outputs) // "restart_file.h5"
			call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
			if (h5error .EQ. -1) then
				write(6,'("KORC ERROR: Something went wrong in: load_particles_ic --> h5fopen_f")')
				call KORC_ABORT()
			end if

			write(tmp_str,'(I18)') ss

			dset = "/spp_" // TRIM(ADJUSTL(tmp_str)) // "/X"
			call load_array_from_hdf5(h5file_id,dset,X_send_buffer)

			call h5fclose_f(h5file_id, h5error)
		end if

		X_receive_buffer = 0.0_rp
		CALL MPI_SCATTER(X_send_buffer,3*spp(ss)%ppp,MPI_REAL8,X_receive_buffer,3*spp(ss)%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
		spp(ss)%vars%X = RESHAPE(X_receive_buffer,(/3,spp(ss)%ppp/))

		if (params%mpi_params%rank.EQ.0_idef) then
			filename = TRIM(params%path_to_outputs) // "restart_file.h5"
			call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
			if (h5error .EQ. -1) then
				write(6,'("KORC ERROR: Something went wrong in: load_particles_ic --> h5fopen_f")')
				call KORC_ABORT()
			end if

			write(tmp_str,'(I18)') ss

			dset = "/spp_" // TRIM(ADJUSTL(tmp_str)) // "/V"
			call load_array_from_hdf5(h5file_id,dset,V_send_buffer)

			call h5fclose_f(h5file_id, h5error)
		end if

		V_receive_buffer = 0.0_rp
		CALL MPI_SCATTER(V_send_buffer,3*spp(ss)%ppp,MPI_REAL8,V_receive_buffer,3*spp(ss)%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
		spp(ss)%vars%V = RESHAPE(V_receive_buffer,(/3,spp(ss)%ppp/))

		if (params%mpi_params%rank.EQ.0_idef) then
			filename = TRIM(params%path_to_outputs) // "restart_file.h5"
			call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
			if (h5error .EQ. -1) then
				write(6,'("KORC ERROR: Something went wrong in: load_particles_ic --> h5fopen_f")')
				call KORC_ABORT()
			end if

			write(tmp_str,'(I18)') ss

			dset = "/spp_" // TRIM(ADJUSTL(tmp_str)) // "/flag"
			call load_array_from_hdf5(h5file_id,dset,AUX_send_buffer)

			call h5fclose_f(h5file_id, h5error)
		end if

		AUX_receive_buffer = 0.0_rp
		CALL MPI_SCATTER(AUX_send_buffer,spp(ss)%ppp,MPI_REAL8,AUX_receive_buffer,spp(ss)%ppp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
		spp(ss)%vars%flag = INT(AUX_receive_buffer,is)

		DEALLOCATE(X_send_buffer)
		DEALLOCATE(X_receive_buffer)

		DEALLOCATE(V_send_buffer)
		DEALLOCATE(V_receive_buffer)

		DEALLOCATE(AUX_send_buffer)
		DEALLOCATE(AUX_receive_buffer)
	end do
end subroutine load_particles_ic

end module korc_HDF5
