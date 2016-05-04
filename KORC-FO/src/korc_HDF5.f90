module korc_HDF5

use korc_types
use HDF5

implicit none
	PRIVATE :: isave_to_hdf5, rsave_to_hdf5
	PUBLIC :: initialize_HDF5, finalize_HDF5, save_simulation_parameters
contains

subroutine initialize_HDF5()
implicit none
	INTEGER :: h5error  ! Error flag
	call h5open_f(h5error)
end subroutine initialize_HDF5


subroutine finalize_HDF5()
implicit none
	INTEGER :: h5error  ! Error flag
	call h5close_f(h5error)
end subroutine finalize_HDF5


subroutine isave_to_hdf5(h5file_id,dset,idata,attr)
implicit none
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: idata
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: attr
	CHARACTER(4) :: aname = "Type"
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

	arank = size(shape(attr))
	ALLOCATE(adims(arank))
	adims = shape(attr)

	! * * * Write data to file * * *
	call h5screate_simple_f(rank,dims,dspace_id,h5error)
	call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_INTEGER, dspace_id, dset_id, h5error)
	call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, idata, dims, h5error)

	! * * * Write attribute of data to file * * *
	tmplen = 0
	attrlen = 0
	do rr=1,arank
		do dd=1,adims(rr)
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
	call h5awrite_f(attr_id, atype_id, attr, dims, h5error)

	call h5aclose_f(attr_id, h5error)
	call h5sclose_f(aspace_id, h5error)
	! * * * Write attribute of data to file * * *

	call h5sclose_f(dspace_id, h5error)
	call h5dclose_f(dset_id, h5error)
	! * * * Write data to file * * *

	DEALLOCATE(dims)
	DEALLOCATE(adims)

end subroutine isave_to_hdf5


subroutine rsave_to_hdf5(h5file_id,dset,rdata,attr)
implicit none
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: rdata
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: attr
	CHARACTER(4) :: aname = "Type"
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

	arank = size(shape(attr))
	ALLOCATE(adims(arank))
	adims = shape(attr)

	! * * * Write data to file * * *
	call h5screate_simple_f(rank,dims,dspace_id,h5error)
	call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_DOUBLE, dspace_id, dset_id, h5error)
	call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, rdata, dims, h5error)

	! * * * Write attribute of data to file * * *
	tmplen = 0
	attrlen = 0
	do rr=1,arank
		do dd=1,adims(rr)
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
	call h5awrite_f(attr_id, atype_id, attr, dims, h5error)

	call h5aclose_f(attr_id, h5error)
	call h5sclose_f(aspace_id, h5error)
	! * * * Write attribute of data to file * * *

	call h5sclose_f(dspace_id, h5error)
	call h5dclose_f(dset_id, h5error)
	! * * * Write data to file * * *

	DEALLOCATE(dims)
	DEALLOCATE(adims)
end subroutine rsave_to_hdf5


subroutine save_simulation_parameters(params,spp,EB,cpp)
implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp
	TYPE(CHARCS_PARAMS), INTENT(IN) :: cpp
	TYPE(FIELDS), INTENT(IN) :: EB
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH) :: dset
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: dims
	REAL(rp), DIMENSION(:), ALLOCATABLE :: rdata
	INTEGER, DIMENSION(:), ALLOCATABLE :: idata
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_data
	INTEGER :: h5error
	INTEGER :: mpierror

	if (params%mpi_params%rank_topo .EQ. 0) then
		filename = TRIM(params%path_to_outputs) // "simulation_parameters.h5"

		call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

		! Simulation parameters group
		gname = "parameters"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)
		
		ALLOCATE(rdata(1))
		ALLOCATE(idata(1))
		ALLOCATE(attr_data(1))

		dset = TRIM(gname) // "/dt"
		attr_data(1) = "Time step in secs"
		rdata = params%dt*cpp%time
		call rsave_to_hdf5(h5file_id,dset,rdata,attr_data)

		dset = TRIM(gname) // "/t_steps"
		attr_data(1) = "Number of time steps"
		idata = params%t_steps
		call isave_to_hdf5(h5file_id,dset,idata,attr_data)

		dset = TRIM(gname) // "/num_omp_threads"
		attr_data(1) = "Number of omp threads"
		idata = params%num_omp_threads
		call isave_to_hdf5(h5file_id,dset,idata,attr_data)

		dset = TRIM(gname) // "/output_cadence"
		attr_data(1) = "Cadence of output files"
		idata = params%output_cadence
		call isave_to_hdf5(h5file_id,dset,idata,attr_data)

		dset = TRIM(gname) // "/num_snapshots"
		attr_data(1) = "Number of outputs for each variable"
		idata = params%num_snapshots
		call isave_to_hdf5(h5file_id,dset,idata,attr_data)

		dset = TRIM(gname) // "/num_species"
		attr_data(1) = "Number of particle species"
		idata = params%num_species
		call isave_to_hdf5(h5file_id,dset,idata,attr_data)

		dset = TRIM(gname) // "/nmpi"
		attr_data(1) = "Number of mpi processes"
		idata = params%mpi_params%nmpi
		call isave_to_hdf5(h5file_id,dset,idata,attr_data)

		DEALLOCATE(rdata)
		DEALLOCATE(idata)
		DEALLOCATE(attr_data)

		call h5gclose_f(group_id, h5error)


		! Plasma species group
		gname = "plasmaSpecies"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)
		call h5gclose_f(group_id, h5error)


		! Electromagnetic fields group
		gname = "emf"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)
		call h5gclose_f(group_id, h5error)

		call h5fclose_f(h5file_id, h5error)
	end if

	call MPI_BARRIER(params%mpi_params%rank_topo, mpierror)

end subroutine save_simulation_parameters

end module korc_HDF5
