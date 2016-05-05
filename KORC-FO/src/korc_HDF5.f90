module korc_HDF5

use korc_types
use HDF5

implicit none

	INTEGER(HID_T), PRIVATE :: KORC_HDF5_REAL ! Real precision used in HDF5
	INTEGER(SIZE_T), PRIVATE :: rp_hdf5 ! Size of real precision used in HDF5

	PRIVATE :: isave_1d_array_to_hdf5, rsave_1d_array_to_hdf5, rsave_2d_array_to_hdf5,&
                isave_to_hdf5, rsave_to_hdf5
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
	INTEGER, INTENT(IN) :: idata
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
	call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, idata, dims, h5error)

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
end subroutine isave_to_hdf5


subroutine isave_1d_array_to_hdf5(h5file_id,dset,idata,attr)
implicit none
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: idata
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
implicit none
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
implicit none
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: rdata
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
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
	CHARACTER(MAX_STRING_LENGTH) :: attr
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

		ALLOCATE(attr_array(1))		
		ALLOCATE(idata(1))

		dset = TRIM(gname) // "/dt"
        attr = "Time step in secs"
		call rsave_to_hdf5(h5file_id,dset,params%dt*cpp%time,attr)

		dset = TRIM(gname) // "/t_steps"
		attr_array(1) = "Number of time steps"
		idata = params%t_steps
		call isave_1d_array_to_hdf5(h5file_id,dset,idata,attr_array)

		dset = TRIM(gname) // "/num_omp_threads"
		attr = "Number of omp threads"
		call isave_to_hdf5(h5file_id,dset, params%num_omp_threads,attr)

		dset = TRIM(gname) // "/output_cadence"
		attr_array(1) = "Cadence of output files"
		idata = params%output_cadence
		call isave_1d_array_to_hdf5(h5file_id,dset,idata,attr_array)

		dset = TRIM(gname) // "/num_snapshots"
		attr_array(1) = "Number of outputs for each variable"
		idata = params%num_snapshots
		call isave_1d_array_to_hdf5(h5file_id,dset,idata,attr_array)

		dset = TRIM(gname) // "/num_species"
		attr = "Number of particle species"
		call isave_to_hdf5(h5file_id,dset,params%num_species,attr)

		dset = TRIM(gname) // "/nmpi"
		attr = "Number of mpi processes"
		call isave_to_hdf5(h5file_id,dset,params%mpi_params%nmpi,attr)

		DEALLOCATE(idata)
		DEALLOCATE(attr_array)

		call h5gclose_f(group_id, h5error)


		! Plasma species group
		gname = "plasmaSpecies"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

		ALLOCATE(rdata(params%num_species))
		ALLOCATE(idata(params%num_species))
		ALLOCATE(attr_array(params%num_species))

		dset = TRIM(gname) // "/ppp"
		attr_array(1) = "Particles per (mpi) process"
		attr_array(2) = "A second attribute"
		idata = spp%ppp
		call isave_1d_array_to_hdf5(h5file_id,dset,idata,attr_array)

		dset = TRIM(gname) // "/q"
		attr_array(1) = "Electric charge"
		rdata = spp%q*cpp%charge
		call rsave_1d_array_to_hdf5(h5file_id,dset,rdata,attr_array)

		dset = TRIM(gname) // "/m"
		attr_array(1) = "Species mass in kg"
		rdata = spp%m*cpp%mass
		call rsave_1d_array_to_hdf5(h5file_id,dset,rdata,attr_array)

		dset = TRIM(gname) // "/Eo"
		attr_array(1) = "Initial (average) energy in eV"
		rdata = spp%Eo*cpp%energy
		call rsave_1d_array_to_hdf5(h5file_id,dset,rdata,attr_array)

		dset = TRIM(gname) // "/wc"
		attr_array(1) = "Average cyclotron frequency in Hz"
		rdata = spp%wc/cpp%time
		call rsave_1d_array_to_hdf5(h5file_id,dset,rdata,attr_array)

		call h5gclose_f(group_id, h5error)

		DEALLOCATE(rdata)
		DEALLOCATE(idata)
		DEALLOCATE(attr_array)

		! Electromagnetic fields group
		gname = "electromagneticFields"
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

		dset = TRIM(gname) // "/Bo"
		attr = "Characteristic (toroidal) field in T"
		call rsave_to_hdf5(h5file_id,dset,EB%Bo*cpp%magnetic_field,attr)

		if (params%magnetic_field_model .EQ. 'ANALYTICAL') then
			dset = TRIM(gname) // "/a"
			attr = "Minor radius in m"
        		call rsave_to_hdf5(h5file_id,dset,EB%AB%a*cpp%length,attr)

			dset = TRIM(gname) // "/Ro"
			attr = "Major radius in m"
        		call rsave_to_hdf5(h5file_id,dset,EB%AB%Ro*cpp%length,attr)

			dset = TRIM(gname) // "/qa"
			attr = "Safety factor at minor radius"
        		call rsave_to_hdf5(h5file_id,dset,EB%AB%Ro*cpp%length,attr)

			dset = TRIM(gname) // "/lambda"
			attr = "Parameter lamda in m"
        		call rsave_to_hdf5(h5file_id,dset,EB%AB%lambda*cpp%length,attr)

			dset = TRIM(gname) // "/Bpo"
			attr = "Poloidal magnetic field in T"
        		call rsave_to_hdf5(h5file_id,dset,EB%AB%Bpo*cpp%magnetic_field,attr)
		end if

		call h5gclose_f(group_id, h5error)

		call h5fclose_f(h5file_id, h5error)
	end if

end subroutine save_simulation_parameters



subroutine rsave_2d_array_to_hdf5(h5file_id,dset,rdata,attr)
implicit none
	INTEGER(HID_T), INTENT(IN) :: h5file_id
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: dset
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: rdata
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


subroutine save_simulation_outputs(params,cpp,spp,EB,it)
implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(CHARCS_PARAMS), INTENT(IN) :: cpp
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp
	TYPE(FIELDS), INTENT(IN) :: EB
	INTEGER(ip), INTENT(IN) :: it
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
	INTEGER :: mpierror
	CHARACTER(19) :: tmp_str
    INTEGER :: ii

	write(tmp_str,'(I18)') params%mpi_params%rank
	filename = TRIM(params%path_to_outputs) // "file_" // TRIM(ADJUSTL(tmp_str)) // ".h5"
	call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, h5file_id, h5error)

    ! Create group 'it'
	write(tmp_str,'(I18)') it
	gname = TRIM(ADJUSTL(tmp_str))
	call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)
    
	dset = TRIM(gname) // "/time"
	attr = "Simulation time in secs"
	call rsave_to_hdf5(h5file_id,dset,REAL(it,ip)*params%dt*cpp%time,attr)

    do ii=1,params%num_species
        write(tmp_str,'(I18)') ii
	    subgname = "spp_" // TRIM(ADJUSTL(tmp_str))
	    call h5gcreate_f(group_id, TRIM(subgname), subgroup_id, h5error)

	    dset = "X"
	    call rsave_2d_array_to_hdf5(subgroup_id, dset, spp(ii)%vars%X)

	    dset = "V"
	    call rsave_2d_array_to_hdf5(subgroup_id, dset, spp(ii)%vars%V)

	    dset = "Rgc"
	    call rsave_2d_array_to_hdf5(subgroup_id, dset, spp(ii)%vars%Rgc)

	    dset = "gamma"
	    call rsave_1d_array_to_hdf5(subgroup_id, dset, spp(ii)%vars%gamma)

	    dset = "eta"
	    call rsave_1d_array_to_hdf5(subgroup_id, dset, spp(ii)%vars%eta)

	    dset = "mu"
	    call rsave_1d_array_to_hdf5(subgroup_id, dset, spp(ii)%vars%mu)

	    dset = "kappa"
	    call rsave_1d_array_to_hdf5(subgroup_id, dset, spp(ii)%vars%kappa)

	    dset = "tau"
	    call rsave_1d_array_to_hdf5(subgroup_id, dset, spp(ii)%vars%tau)

        call h5gclose_f(subgroup_id, h5error)
    end do

	call h5gclose_f(group_id, h5error)

	call h5fclose_f(h5file_id, h5error)
end subroutine save_simulation_outputs

end module korc_HDF5
