MODULE korc_binning_diagnostic
	USE korc_types
	USE korc_constants
	USE korc_HDF5
    USE korc_hpc

	IMPLICIT NONE

	TYPE, PRIVATE :: BINNING
		LOGICAL :: diagnostic_on
		REAL(rp) :: start_at ! In seconds
		INTEGER, DIMENSION(2) :: num_bins ! Number of bins
		REAL(rp), DIMENSION(2) :: rlim ! in meters
		REAL(rp), DIMENSION(2) :: zlim ! in meters
		REAL(rp) :: rmin
		REAL(rp) :: rmax
		REAL(rp) :: zmin
		REAL(rp) :: zmax
		REAL(rp) :: dr
		REAL(rp) :: dz
		REAL(rp), DIMENSION(:), ALLOCATABLE :: rnodes
		REAL(rp), DIMENSION(:), ALLOCATABLE :: znodes
		LOGICAL :: toroidal_sections
		INTEGER(idef) :: ntor_sections
	END TYPE BINNING

	TYPE(BINNING), PRIVATE :: binning_params

	INTERFACE save_snapshot_var
	  module procedure save_snapshot_var_1d,save_snapshot_var_2d,save_snapshot_var_3d,save_snapshot_var_4d
	END INTERFACE

	PRIVATE :: clockwise_rotation,&
				anticlockwise_rotation,&
				cross,&
				save_binning_diagnostic_params,&
				save_snapshot_var_1d,&
				save_snapshot_var_2d,&
				save_snapshot_var_3d,&
				save_snapshot_var_4d
	PUBLIC :: initialize_binning_diagnostic,&
				binning_diagnostic

	CONTAINS

SUBROUTINE initialize_binning_diagnostic(params)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	LOGICAL :: diagnostic_on
	REAL(rp) :: start_at ! In seconds
	INTEGER, DIMENSION(2) :: num_bins ! Number of bins
	REAL(rp), DIMENSION(2) :: rlim ! in meters
	REAL(rp), DIMENSION(2) :: zlim ! in meters
	LOGICAL :: toroidal_sections
	INTEGER(idef) :: ntor_sections
	INTEGER :: ii

	NAMELIST /BinningDiagnostic/ diagnostic_on,start_at,num_bins,rlim,zlim,toroidal_sections,ntor_sections

	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=BinningDiagnostic)
	close(default_unit_open)

	
	binning_params%diagnostic_on = diagnostic_on

	if (binning_params%diagnostic_on) then
		if (params%mpi_params%rank .EQ. 0) then
			write(6,'(/,"* * * * * * * * * * * * * * * * * *")')
			write(6,'("* Initializing Binning diagnostic *")')
		end if
	
		binning_params%start_at = start_at
		binning_params%num_bins = num_bins
		binning_params%rlim = rlim
		binning_params%zlim = zlim
		binning_params%toroidal_sections = toroidal_sections
		if (binning_params%toroidal_sections) then
			binning_params%ntor_sections = ntor_sections
		else
			binning_params%ntor_sections = 1_idef
		end if	

		binning_params%dr = (binning_params%rlim(2) - binning_params%rlim(1))/REAL(binning_params%num_bins(1),rp)
		binning_params%dz = (binning_params%zlim(2) - binning_params%zlim(1))/REAL(binning_params%num_bins(2),rp)

		binning_params%rmin = MINVAL(binning_params%rlim)
		binning_params%rmax = MAXVAL(binning_params%rlim)

		binning_params%zmin = MINVAL(binning_params%zlim)
		binning_params%zmax = MAXVAL(binning_params%zlim)
	
		ALLOCATE(binning_params%rnodes(binning_params%num_bins(1)))
		ALLOCATE(binning_params%znodes(binning_params%num_bins(2)))

		do ii=1_idef,binning_params%num_bins(1)
			binning_params%rnodes(ii) = binning_params%rmin + (REAL(ii-1_idef,rp)+0.5_rp)*binning_params%dr
		end do

		do ii=1_idef,binning_params%num_bins(2)
			binning_params%znodes(ii) = binning_params%zmin + (REAL(ii-1_idef,rp)+0.5_rp)*binning_params%dr
		end do

		call save_binning_diagnostic_params(params)

		if (params%mpi_params%rank .EQ. 0) then
			write(6,'("*    Binning diagnostic ready!    *")')
			write(6,'("* * * * * * * * * * * * * * * * * *")')
		end if
	end if
END SUBROUTINE initialize_binning_diagnostic


! * * * * * * * * * * * * * * * !
! * * * * * FUNCTIONS * * * * * !
! * * * * * * * * * * * * * * * !

FUNCTION cross(a,b)
	REAL(rp), DIMENSION(3), INTENT(IN) :: a
	REAL(rp), DIMENSION(3), INTENT(IN) :: b
	REAL(rp), DIMENSION(3) :: cross

	cross(1) = a(2)*b(3) - a(3)*b(2)
	cross(2) = a(3)*b(1) - a(1)*b(3)
	cross(3) = a(1)*b(2) - a(2)*b(1)
END FUNCTION cross


FUNCTION clockwise_rotation(x,t)
	IMPLICIT NONE
	REAL(rp), DIMENSION(2), INTENT(IN) :: x
	REAL(rp), INTENT(IN) :: t ! Angle in radians
	REAL(rp), DIMENSION(2) :: clockwise_rotation

	clockwise_rotation(1) = x(1)*COS(t) + x(2)*SIN(t)
	clockwise_rotation(2) = -x(1)*SIN(t) + x(2)*COS(t)
END FUNCTION clockwise_rotation


FUNCTION anticlockwise_rotation(x,t)
	IMPLICIT NONE
	REAL(rp), DIMENSION(2), INTENT(IN) :: x
	REAL(rp), INTENT(IN) :: t ! Angle in radians
	REAL(rp), DIMENSION(2) :: anticlockwise_rotation

	anticlockwise_rotation(1) = x(1)*COS(t) - x(2)*SIN(t)
	anticlockwise_rotation(2) = x(1)*SIN(t) + x(2)*COS(t)
END FUNCTION anticlockwise_rotation

! * * * * * * * * * * * * * * * !
! * * * * * FUNCTIONS * * * * * !
! * * * * * * * * * * * * * * * !


SUBROUTINE bin_variables(params,spp)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp
	CHARACTER(MAX_STRING_LENGTH) :: var_name
	REAL(rp), DIMENSION(3) :: X
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: eta
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: g
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE :: N
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: array3D
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: array2D
	REAL(rp), DIMENSION(:), ALLOCATABLE :: array1D
	REAL(rp) :: R,Z,phi
	REAL(rp) :: Dtor
	REAL(rp) :: q,m
	INTEGER :: ii,jj,kk,ss,pp
    REAL(rp), DIMENSION(:), ALLOCATABLE :: send_buffer, receive_buffer
    INTEGER :: numel, mpierr
	REAL(rp) :: units

	ALLOCATE(eta(binning_params%num_bins(1),binning_params%num_bins(2),binning_params%ntor_sections,params%num_species))
	ALLOCATE(g(binning_params%num_bins(1),binning_params%num_bins(2),binning_params%ntor_sections,params%num_species))
	ALLOCATE(N(binning_params%num_bins(1),binning_params%num_bins(2),binning_params%ntor_sections,params%num_species))

	eta = 0.0_rp
	g = 0.0_rp
	N = 0.0_rp

	Dtor = 2.0_rp*C_PI/REAL(binning_params%ntor_sections,rp)

	do ss=1_idef,params%num_species
		q = ABS(spp(ss)%q)*params%cpp%charge
		m = spp(ss)%m*params%cpp%mass
!$OMP PARALLEL DO FIRSTPRIVATE(q,m,Dtor) PRIVATE(X,ii,jj,kk,pp)&
!$OMP& SHARED(params,spp,ss,eta,g,N)
		do pp=1_idef,spp(ss)%ppp
			if ( spp(ss)%vars%flag(pp) .EQ. 1_idef ) then
				X = spp(ss)%vars%X(:,pp)*params%cpp%length

				R = SQRT(SUM(X(1:2)**2))
				Z = X(3)
				
				ii = FLOOR((R - binning_params%rmin)/binning_params%dr) + 1_idef
				jj = FLOOR((Z + ABS(binning_params%zmin))/binning_params%dz) + 1_idef

				phi = ATAN2(X(2),X(1))
				if (phi.LT.0.0_rp) phi = 2.0_rp*C_PI + phi
				kk = floor(phi/Dtor) + 1_idef

				eta(ii,jj,kk,ss) = eta(ii,jj,kk,ss) + spp(ss)%vars%eta(pp)
				g(ii,jj,kk,ss) = g(ii,jj,kk,ss) + spp(ss)%vars%g(pp)
				N(ii,jj,kk,ss) = N(ii,jj,kk,ss) + 1.0_rp
			end if ! if confined
		end do ! particles
!$OMP END PARALLEL DO
	end do ! species

	if (params%mpi_params%rank.EQ.0_idef) then
		if (.NOT.binning_params%toroidal_sections) then
			ALLOCATE(array3D(binning_params%num_bins(1),binning_params%num_bins(2),params%num_species))
		end if
	end if

	if (params%mpi_params%nmpi.GT.1_idef) then
		numel = binning_params%num_bins(1)*binning_params%num_bins(2)*binning_params%ntor_sections*params%num_species

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(eta,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    eta = RESHAPE(receive_buffer,(/binning_params%num_bins(1),binning_params%num_bins(2),&
							binning_params%ntor_sections,params%num_species/))

			var_name = 'eta'
			if (binning_params%toroidal_sections) then
				call save_snapshot_var(params,eta,var_name)
			else
				array3D = SUM(eta,3)
				call save_snapshot_var(params,array3D,var_name)
			end if
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)

		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(g,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    g = RESHAPE(receive_buffer,(/binning_params%num_bins(1),binning_params%num_bins(2),&
							binning_params%ntor_sections,params%num_species/))

			var_name = 'g'
			if (binning_params%toroidal_sections) then
				call save_snapshot_var(params,g,var_name)
			else
				array3D = SUM(g,3)
				call save_snapshot_var(params,array3D,var_name)
			end if
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)
		
		ALLOCATE(send_buffer(numel))
		ALLOCATE(receive_buffer(numel))

		send_buffer = RESHAPE(N,(/numel/))
		receive_buffer = 0.0_rp
		CALL MPI_REDUCE(send_buffer,receive_buffer,numel,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
		if (params%mpi_params%rank.EQ.0_idef) then
		    N = RESHAPE(receive_buffer,(/binning_params%num_bins(1),binning_params%num_bins(2),&
							binning_params%ntor_sections,params%num_species/))

			var_name = 'N'
			if (binning_params%toroidal_sections) then
				call save_snapshot_var(params,N,var_name)
			else
				array3D = SUM(N,3)
				call save_snapshot_var(params,array3D,var_name)
			end if
		end if

		DEALLOCATE(send_buffer)
		DEALLOCATE(receive_buffer)
	else
		var_name = 'eta'
		if (binning_params%toroidal_sections) then
			call save_snapshot_var(params,eta,var_name)
		else
			array3D = SUM(eta,3)
			call save_snapshot_var(params,array3D,var_name)
		end if

		var_name = 'g'
		if (binning_params%toroidal_sections) then
			call save_snapshot_var(params,g,var_name)
		else
			array3D = SUM(g,3)
			call save_snapshot_var(params,array3D,var_name)
		end if
		
		var_name = 'N'
		if (binning_params%toroidal_sections) then
			call save_snapshot_var(params,N,var_name)
		else
			array3D = SUM(N,3)
			call save_snapshot_var(params,array3D,var_name)
		end if
	end if

	DEALLOCATE(eta)
	DEALLOCATE(g)
	DEALLOCATE(N)
	
	if (ALLOCATED(array3D)) DEALLOCATE(array3D)
	if (ALLOCATED(array2D)) DEALLOCATE(array2D)
	if (ALLOCATED(array1D)) DEALLOCATE(array1D)
END SUBROUTINE bin_variables


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * SUBROUTINES TO GENERATE OUTPUTS OF THE SYNTHETIC CAMERA * * * * 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


SUBROUTINE save_binning_diagnostic_params(params)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
	CHARACTER(MAX_STRING_LENGTH) :: dset
	CHARACTER(MAX_STRING_LENGTH) :: attr
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	CHARACTER(19) :: tmp_str
	INTEGER :: h5error
	REAL(rp) :: units
	
	if (.NOT.params%restart) then

		if (params%mpi_params%rank .EQ. 0) then
			filename = TRIM(params%path_to_outputs) // "binning_diagnostic.h5"
			call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

			gname = "binning_diagnostic_params"
			call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

			dset = TRIM(gname) // "/start_at"
			attr = "Time at which camera starts working (s)"
			call save_to_hdf5(h5file_id,dset,binning_params%start_at,attr)

			dset = TRIM(gname) // "/num_bins"
			call save_1d_array_to_hdf5(h5file_id,dset,binning_params%num_bins)

			dset = TRIM(gname) // "/rlim"
			call save_1d_array_to_hdf5(h5file_id,dset,binning_params%rlim)

			dset = TRIM(gname) // "/zlim"
			call save_1d_array_to_hdf5(h5file_id,dset,binning_params%zlim)

			dset = TRIM(gname) // "/dr"
			attr = "Size of bin along the radial direction (m)"
			call save_to_hdf5(h5file_id,dset,binning_params%dr,attr)

			dset = TRIM(gname) // "/dz"
			attr = "Size of bin along the radial direction (m)"
			call save_to_hdf5(h5file_id,dset,binning_params%dr,attr)

			dset = TRIM(gname) // "/rnodes"
			call save_1d_array_to_hdf5(h5file_id,dset,binning_params%rnodes)

			dset = TRIM(gname) // "/znodes"
			call save_1d_array_to_hdf5(h5file_id,dset,binning_params%znodes)

			dset = TRIM(gname) // "/toroidal_sections"
			attr = "Logical variable: 1=decomposed in toroidal sections, 0=no toroidal decomposition"
			if (binning_params%toroidal_sections) then
				call save_to_hdf5(h5file_id,dset,1_idef,attr)

				dset = TRIM(gname) // "/ntor_sections"
				attr = "Number of toroidal sections"
				call save_to_hdf5(h5file_id,dset,binning_params%ntor_sections,attr)
			else
				call save_to_hdf5(h5file_id,dset,0_idef,attr)
			end if

			call h5gclose_f(group_id, h5error)

			call h5fclose_f(h5file_id, h5error)
		end if

		if (params%mpi_params%rank.EQ.0_idef) then
			filename = TRIM(params%path_to_outputs) //"binning_diagnostic_snapshots.h5"
			call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)
			call h5fclose_f(h5file_id, h5error)
		end if
    
    end if
END SUBROUTINE save_binning_diagnostic_params


SUBROUTINE save_snapshot_var_1d(params,var,var_name)
IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: var
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: var_name
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH) :: subgname
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
	CHARACTER(MAX_STRING_LENGTH) :: dset
	CHARACTER(MAX_STRING_LENGTH) :: attr
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	INTEGER(HID_T) :: subgroup_id
	CHARACTER(19) :: tmp_str
	INTEGER :: h5error
	INTEGER :: ss
	LOGICAL :: object_exists

	filename = TRIM(params%path_to_outputs) //"binning_diagnostic_snapshots.h5"
	call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, h5file_id, h5error)

    ! Create group 'it' if it doesn't exist
	write(tmp_str,'(I18)') params%it
	gname = TRIM(ADJUSTL(tmp_str))
	call h5lexists_f(h5file_id,TRIM(gname),object_exists,h5error)

	if (.NOT.object_exists) then
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)
		
		dset = TRIM(gname) // "/time"
		attr = "Simulation time in secs"
		call save_to_hdf5(h5file_id,dset,REAL(params%it,rp)*params%dt*params%cpp%time,attr)
	else
		call h5gopen_f(h5file_id, TRIM(gname), group_id, h5error)
	end if

	do ss=1_idef,params%num_species
		write(tmp_str,'(I18)') ss
		subgname = "spp_" // TRIM(ADJUSTL(tmp_str))
		call h5lexists_f(group_id,TRIM(subgname),object_exists,h5error)

		if (.NOT.object_exists) then
			call h5gcreate_f(group_id, TRIM(subgname), subgroup_id, h5error)
		else
			call h5gopen_f(group_id, TRIM(subgname), subgroup_id, h5error)
		end if

		dset = TRIM(var_name)
		call h5lexists_f(subgroup_id,TRIM(dset),object_exists,h5error)
		if (.NOT.object_exists) then
			call save_to_hdf5(subgroup_id,dset,var(ss))
		end if

		call h5gclose_f(subgroup_id, h5error)
	end do	

	call h5gclose_f(group_id, h5error)

	call h5fclose_f(h5file_id, h5error)
END SUBROUTINE save_snapshot_var_1d


SUBROUTINE save_snapshot_var_2d(params,var,var_name)
IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: var
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: var_name
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH) :: subgname
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
	CHARACTER(MAX_STRING_LENGTH) :: dset
	CHARACTER(MAX_STRING_LENGTH) :: attr
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	INTEGER(HID_T) :: subgroup_id
	CHARACTER(19) :: tmp_str
	INTEGER :: h5error
	INTEGER :: ss
	LOGICAL :: object_exists

	filename = TRIM(params%path_to_outputs) //"binning_diagnostic_snapshots.h5"
	call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, h5file_id, h5error)

    ! Create group 'it' if it doesn't exist
	write(tmp_str,'(I18)') params%it
	gname = TRIM(ADJUSTL(tmp_str))
	call h5lexists_f(h5file_id,TRIM(gname),object_exists,h5error)

	if (.NOT.object_exists) then
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)
		
		dset = TRIM(gname) // "/time"
		attr = "Simulation time in secs"
		call save_to_hdf5(h5file_id,dset,REAL(params%it,rp)*params%dt*params%cpp%time,attr)
	else
		call h5gopen_f(h5file_id, TRIM(gname), group_id, h5error)
	end if

	do ss=1_idef,params%num_species
		write(tmp_str,'(I18)') ss
		subgname = "spp_" // TRIM(ADJUSTL(tmp_str))
		call h5lexists_f(group_id,TRIM(subgname),object_exists,h5error)

		if (.NOT.object_exists) then
			call h5gcreate_f(group_id, TRIM(subgname), subgroup_id, h5error)
		else
			call h5gopen_f(group_id, TRIM(subgname), subgroup_id, h5error)
		end if

		dset = TRIM(var_name)
		call h5lexists_f(subgroup_id,TRIM(dset),object_exists,h5error)
		if (.NOT.object_exists) then
			call save_array_to_hdf5(subgroup_id, dset, var(:,ss))
		end if

		call h5gclose_f(subgroup_id, h5error)
	end do	

	call h5gclose_f(group_id, h5error)

	call h5fclose_f(h5file_id, h5error)
END SUBROUTINE save_snapshot_var_2d


SUBROUTINE save_snapshot_var_3d(params,var,var_name)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(IN) :: var
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: var_name
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH) :: subgname
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
	CHARACTER(MAX_STRING_LENGTH) :: dset
	CHARACTER(MAX_STRING_LENGTH) :: attr
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	INTEGER(HID_T) :: subgroup_id
	CHARACTER(19) :: tmp_str
	INTEGER :: h5error
	INTEGER :: ss
	LOGICAL :: object_exists

	filename = TRIM(params%path_to_outputs) //"binning_diagnostic_snapshots.h5"
	call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, h5file_id, h5error)

    ! Create group 'it' if it doesn't exist
	write(tmp_str,'(I18)') params%it
	gname = TRIM(ADJUSTL(tmp_str))
	call h5lexists_f(h5file_id,TRIM(gname),object_exists,h5error)

	if (.NOT.object_exists) then
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)
		
		dset = TRIM(gname) // "/time"
		attr = "Simulation time in secs"
		call save_to_hdf5(h5file_id,dset,REAL(params%it,rp)*params%dt*params%cpp%time,attr)
	else
		call h5gopen_f(h5file_id, TRIM(gname), group_id, h5error)
	end if

	do ss=1_idef,params%num_species
		write(tmp_str,'(I18)') ss
		subgname = "spp_" // TRIM(ADJUSTL(tmp_str))
		call h5lexists_f(group_id,TRIM(subgname),object_exists,h5error)

		if (.NOT.object_exists) then
			call h5gcreate_f(group_id, TRIM(subgname), subgroup_id, h5error)
		else
			call h5gopen_f(group_id, TRIM(subgname), subgroup_id, h5error)
		end if

		dset = TRIM(var_name)
		call h5lexists_f(subgroup_id,TRIM(dset),object_exists,h5error)
		if (.NOT.object_exists) then
			call save_array_to_hdf5(subgroup_id, dset, var(:,:,ss))
		end if

		call h5gclose_f(subgroup_id, h5error)
	end do	

	call h5gclose_f(group_id, h5error)

	call h5fclose_f(h5file_id, h5error)
END SUBROUTINE save_snapshot_var_3d


SUBROUTINE save_snapshot_var_4d(params,var,var_name)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	REAL(rp), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(IN) :: var
	CHARACTER(MAX_STRING_LENGTH), INTENT(IN) :: var_name
	CHARACTER(MAX_STRING_LENGTH) :: filename
	CHARACTER(MAX_STRING_LENGTH) :: gname
	CHARACTER(MAX_STRING_LENGTH) :: subgname
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
	CHARACTER(MAX_STRING_LENGTH) :: dset
	CHARACTER(MAX_STRING_LENGTH) :: attr
	INTEGER(HID_T) :: h5file_id
	INTEGER(HID_T) :: group_id
	INTEGER(HID_T) :: subgroup_id
	CHARACTER(19) :: tmp_str
	INTEGER :: h5error
	INTEGER :: ss
	LOGICAL :: object_exists

	filename = TRIM(params%path_to_outputs) //"binning_diagnostic_snapshots.h5"
	call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, h5file_id, h5error)

    ! Create group 'it' if it doesn't exist
	write(tmp_str,'(I18)') params%it
	gname = TRIM(ADJUSTL(tmp_str))
	call h5lexists_f(h5file_id,TRIM(gname),object_exists,h5error)

	if (.NOT.object_exists) then
		call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)
		
		dset = TRIM(gname) // "/time"
		attr = "Simulation time in secs"
		call save_to_hdf5(h5file_id,dset,REAL(params%it,rp)*params%dt*params%cpp%time,attr)
	else
		call h5gopen_f(h5file_id, TRIM(gname), group_id, h5error)
	end if

	do ss=1_idef,params%num_species
		write(tmp_str,'(I18)') ss
		subgname = "spp_" // TRIM(ADJUSTL(tmp_str))
		call h5lexists_f(group_id,TRIM(subgname),object_exists,h5error)

		if (.NOT.object_exists) then
			call h5gcreate_f(group_id, TRIM(subgname), subgroup_id, h5error)
		else
			call h5gopen_f(group_id, TRIM(subgname), subgroup_id, h5error)
		end if

		dset = TRIM(var_name)
		call h5lexists_f(subgroup_id,TRIM(dset),object_exists,h5error)
		if (.NOT.object_exists) then
			call save_array_to_hdf5(subgroup_id, dset, var(:,:,:,ss))
		end if

		call h5gclose_f(subgroup_id, h5error)
	end do	

	call h5gclose_f(group_id, h5error)

	call h5fclose_f(h5file_id, h5error)
END SUBROUTINE save_snapshot_var_4d


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * MAIN CALL TO SYNTHETIC CAMERA SUBROUTINES * * * * 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


SUBROUTINE binning_diagnostic(params,spp)
	IMPLICIT NONE
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: spp

	if (binning_params%diagnostic_on.AND.(params%time*params%cpp%time >= binning_params%start_at)) then
		if (params%mpi_params%rank .EQ. 0) then
			write(6,'("Binning diagnostic: ON")',advance="no")
		end if

		call bin_variables(params,spp)

		if (params%mpi_params%rank .EQ. 0) then
			write(6,'(" ---> OFF")')
		end if
	end if
END SUBROUTINE binning_diagnostic

END MODULE korc_binning_diagnostic
