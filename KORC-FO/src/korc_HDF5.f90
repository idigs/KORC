module korc_HDF5

use korc_types
use HDF5

implicit none

	INTEGER(HID_T), PRIVATE :: KORC_HDF5_REAL ! Real precision used in HDF5
	INTEGER(SIZE_T), PRIVATE :: rp_hdf5 ! Size of real precision used in HDF5

	PRIVATE :: isave_to_hdf5, rsave_to_hdf5
	PUBLIC :: initialize_HDF5, finalize_HDF5, save_simulation_parameters

contains


subroutine initialize_HDF5()
implicit none
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
implicit none
	INTEGER :: h5error  ! Error flag
	call h5close_f(h5error)
end subroutine finalize_HDF5


subroutine isave_to_hdf5(h5file_id,dset,idata,attr)
implicit none
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: idata
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: attr
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

		DEALLOCATE(adims)
	end if

	call h5sclose_f(dspace_id, h5error)
	call h5dclose_f(dset_id, h5error)
	! * * * Write data to file * * *

	DEALLOCATE(dims)

end subroutine isave_to_hdf5


subroutine rsave_to_hdf5(h5file_id,dset,rdata,attr)
implicit none
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: rdata
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: attr
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

		DEALLOCATE(adims)
	end if

	call h5sclose_f(dspace_id, h5error)
	call h5dclose_f(dset_id, h5error)
	! * * * Write data to file * * *

	DEALLOCATE(dims)
end subroutine rsave_to_hdf5


subroutine save_simulation_parameters(params,cpp,spp,EB)
implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(CHARCS_PARAMS), INTENT(IN) :: cpp
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp
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
	CHARACTER(19) :: tmp_str

	write(tmp_str,'(I18)') params%mpi_params%rank
	filename = TRIM(params%path_to_outputs) // "file_" // TRIM(ADJUSTL(tmp_str)) // ".h5"
	call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)
	call h5fclose_f(h5file_id, h5error)

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
!		call rsave_to_hdf5(h5file_id,dset,rdata,attr_data)
		call rsave_to_hdf5(h5file_id,dset,rdata)

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

		ALLOCATE(rdata(params%num_species))
		ALLOCATE(idata(params%num_species))
		ALLOCATE(attr_data(params%num_species))

		dset = TRIM(gname) // "/ppp"
		attr_data(1) = "Particles per (mpi) process"
		idata = spp%ppp
		call isave_to_hdf5(h5file_id,dset,idata,attr_data)

		dset = TRIM(gname) // "/q"
		attr_data(1) = "Electric charge"
		rdata = spp%q*cpp%charge
		call rsave_to_hdf5(h5file_id,dset,rdata,attr_data)

		dset = TRIM(gname) // "/m"
		attr_data(1) = "Species mass in kg"
		rdata = spp%m*cpp%mass
		call rsave_to_hdf5(h5file_id,dset,rdata,attr_data)

		dset = TRIM(gname) // "/Eo"
		attr_data(1) = "Initial (average) energy in eV"
		rdata = spp%Eo*cpp%energy
		call rsave_to_hdf5(h5file_id,dset,rdata,attr_data)

		dset = TRIM(gname) // "/wc"
		attr_data(1) = "Average cyclotron frequency in Hz"
		rdata = spp%wc/cpp%time
		call rsave_to_hdf5(h5file_id,dset,rdata,attr_data)

		call h5gclose_f(group_id, h5error)

		DEALLOCATE(rdata)
		DEALLOCATE(idata)
		DEALLOCATE(attr_data)

		! Electromagnetic fields group
		gname = "electromagneticFields"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

		ALLOCATE(rdata(1))
		ALLOCATE(idata(1))
		ALLOCATE(attr_data(1))

		dset = TRIM(gname) // "/Bo"
		attr_data(1) = "Characteristic (toroidal) field in T"
		rdata = EB%Bo*cpp%magnetic_field
		call rsave_to_hdf5(h5file_id,dset,rdata,attr_data)

		if (params%magnetic_field_model .EQ. 'ANALYTICAL') then
			dset = TRIM(gname) // "/a"
			attr_data(1) = "Minor radius in m"
			rdata = EB%AB%a*cpp%length
			call rsave_to_hdf5(h5file_id,dset,rdata,attr_data)

			dset = TRIM(gname) // "/Ro"
			attr_data(1) = "Major radius in m"
			rdata = EB%AB%Ro*cpp%length
			call rsave_to_hdf5(h5file_id,dset,rdata,attr_data)

			dset = TRIM(gname) // "/qa"
			attr_data(1) = "Safety factor at minor radius"
			rdata = EB%AB%qa
			call rsave_to_hdf5(h5file_id,dset,rdata,attr_data)

			dset = TRIM(gname) // "/lambda"
			attr_data(1) = "Parameter lamda in m"
			rdata = EB%AB%lambda*cpp%length
			call rsave_to_hdf5(h5file_id,dset,rdata,attr_data)

			dset = TRIM(gname) // "/Bpo"
			attr_data(1) = "Poloidal magnetic field in T"
			rdata = EB%AB%Bpo*cpp%magnetic_field
			call rsave_to_hdf5(h5file_id,dset,rdata,attr_data)
		end if

		DEALLOCATE(rdata)
		DEALLOCATE(idata)
		DEALLOCATE(attr_data)

		call h5gclose_f(group_id, h5error)

		call h5fclose_f(h5file_id, h5error)
	end if

end subroutine save_simulation_parameters



subroutine rsave_2d_array_to_hdf5(h5file_id,dset,rdata,attr)
implicit none
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: rdata
	CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: attr
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

		DEALLOCATE(adims)
	end if

	call h5sclose_f(dspace_id, h5error)
	call h5dclose_f(dset_id, h5error)
	! * * * Write data to file * * *

	DEALLOCATE(dims)
end subroutine rsave_2d_array_to_hdf5


subroutine save_simulation_outputs(params,cpp,spp,EB,it)
implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(CHARCS_PARAMS), INTENT(IN) :: cpp
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp
	TYPE(FIELDS), INTENT(IN) :: EB
	INTEGER(ip), INTENT(IN) :: it
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
	CHARACTER(19) :: tmp_str

	write(tmp_str,'(I18)') params%mpi_params%rank
!	write(6,*) TRIM(ADJUSTL(tmp_str))

	filename = TRIM(params%path_to_outputs) // "file_" // TRIM(ADJUSTL(tmp_str)) // ".h5"
!	call h5fcreate_f(TRIM(filename), H5F_ACC_EXCL_F, h5file_id, h5error)
	call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, h5file_id, h5error)

	write(tmp_str,'(I18)') it
	gname = TRIM(ADJUSTL(tmp_str))
	call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)
	call h5gclose_f(group_id, h5error)

	call h5fclose_f(h5file_id, h5error)
end subroutine save_simulation_outputs

end module korc_HDF5
