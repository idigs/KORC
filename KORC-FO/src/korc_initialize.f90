module korc_initialize

    use korc_types
    use korc_constants
    use korc_hpc
    use korc_HDF5
    use korc_fields
    use rnd_numbers

	use korc_avalanche ! external module

    implicit none
	

	PRIVATE :: set_paths,load_korc_params,initialization_sanity_check,unitVectors,random_norm,fth_1V,fth_3V,iso_thermal_distribution
	PUBLIC :: initialize_korc_parameters,initialize_particles,initialize_fields

    contains

! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! ** SUBROUTINES FOR INITIALIZING KORC PARAMETERS ** !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !

subroutine set_paths(params)
	implicit none
	INTEGER :: argn
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params

	argn = command_argument_count()
	call get_command_argument(1,params%path_to_inputs)
	call get_command_argument(2,params%path_to_outputs)

	if (params%mpi_params%rank .EQ. 0) then
		write(6,'(/,"* * * * * PATHS * * * * *")')
		write(6,'("The input file is:",A50)') TRIM(params%path_to_inputs)
		write(6,'("The output folder is:",A50,/)') TRIM(params%path_to_outputs)
	end if
end subroutine set_paths


subroutine load_korc_params(params)
	implicit none
	TYPE (KORC_PARAMS), INTENT(INOUT) :: params
	LOGICAL :: restart ! Not used, yet.
	INTEGER(ip) :: t_steps
	REAL(rp) :: dt
	LOGICAL :: radiation
	LOGICAL :: collisions
	CHARACTER(MAX_STRING_LENGTH) :: collisions_model
	CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_model
	LOGICAL :: poloidal_flux
	CHARACTER(MAX_STRING_LENGTH) :: magnetic_field_filename
	CHARACTER(MAX_STRING_LENGTH) :: outputs_list
	INTEGER(ip) :: output_cadence
	INTEGER :: num_species
	INTEGER :: num_impurity_species
	INTEGER :: imax,imin,ii,jj,num_outputs
	INTEGER, DIMENSION(2) :: indices

	NAMELIST /input_parameters/ magnetic_field_model,poloidal_flux,magnetic_field_filename,t_steps,dt,&
            output_cadence,num_species,radiation,collisions,collisions_model,outputs_list
	
	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=input_parameters)
	close(default_unit_open)

	! params%restart = restart
	params%t_steps = t_steps
	params%output_cadence = output_cadence
	params%num_snapshots = t_steps/output_cadence
	params%dt = dt
	params%num_species = num_species
	params%magnetic_field_model = TRIM(magnetic_field_model)
	params%poloidal_flux = poloidal_flux
	params%magnetic_field_filename = TRIM(magnetic_field_filename)
	params%radiation = radiation
	params%collisions = collisions
	params%collisions_model = TRIM(collisions_model)

	! Loading list of output parameters
	imin = SCAN(outputs_list,'{')
	imax = SCAN(outputs_list,'}')

	ii = 1_idef
	jj = 1_idef
	num_outputs = 1_idef
	do while (ii.NE.0)
		ii = SCAN(outputs_list(jj:),",")
		if (ii.NE.0) then
			jj = jj + ii
			num_outputs = num_outputs + 1_idef
		end if
	end do

	ALLOCATE(params%outputs_list(num_outputs))

	if (num_outputs.GT.1_idef) then
		indices = 0_idef
		indices(2) = SCAN(outputs_list,",")
		params%outputs_list(1) = TRIM(outputs_list(imin+1_idef:indices(2)-1_idef))
		indices(1) = indices(1) + indices(2) + 1_idef
		do ii=2_idef,num_outputs
			indices(2) = SCAN(outputs_list(indices(1):),",")
			if (indices(2).EQ.0_idef) then
				params%outputs_list(ii) = TRIM(outputs_list(indices(1):imax-1_idef))
			else
				params%outputs_list(ii) = TRIM(outputs_list(indices(1):indices(1)+indices(2)-2_idef))
				indices(1) = indices(1) + indices(2)
			end if
		end do
	else
		params%outputs_list(1) = TRIM(outputs_list(imin+1_idef:imax-1_idef))
	end if

!	do ii=1_idef,SIZE(params%outputs_list)
!		write(6,*) TRIM(params%outputs_list(ii))
!	end do
	! Loading list of output parameters

	if (params%mpi_params%rank .EQ. 0) then
		write(6,'(/,"* * * * * SIMULATION PARAMETERS * * * * *")')
		write(6,'("Number of time steps: ",I16)') params%t_steps
		write(6,'("Output cadence: ",I16)') params%output_cadence
		write(6,'("Number of outputs: ",I16)') params%num_snapshots
		write(6,'("Time step in fraction of relativistic gyro-period: ",F15.10)') params%dt
		write(6,'("Number of electron populations: ",I16)') params%num_species
		write(6,'("Magnetic field model: ",A50)') TRIM(params%magnetic_field_model)
		write(6,'("USINg (JFIT) poloidal flux: ", L1)') params%poloidal_flux
		write(6,'("Magnetic field model: ",A100)') TRIM(params%magnetic_field_filename)
		write(6,'("Radiation losses included: ",L1)') params%radiation
		write(6,'("collisions losses included: ",L1)') params%collisions
		write(6,'("collisions model: ",A50)') TRIM(params%collisions_model)
		write(6,'(/)')
	end if	
end subroutine load_korc_params


subroutine initialize_korc_parameters(params)
	use korc_types
	implicit none
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params
	INTEGER :: mpierr


	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

	call set_paths(params)

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

	call load_korc_params(params)

	call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
end subroutine initialize_korc_parameters


! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! * * * SUBROUTINES FOR INITIALIZING PARTICLES * * * !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !

subroutine initialize_particles(params,F,spp) 
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: F
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: spp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: ppp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: q
	REAL(rp), DIMENSION(:), ALLOCATABLE :: m
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Eo
	REAL(rp), DIMENSION(:), ALLOCATABLE :: etao
	LOGICAL, DIMENSION(:), ALLOCATABLE :: runaway
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: energy_distribution
	CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: pitch_distribution
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Ro
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Zo
	REAL(rp), DIMENSION(:), ALLOCATABLE :: r
	INTEGER :: ii,jj ! Iterator

	NAMELIST /plasma_species/ ppp, q, m, Eo, etao, runaway, energy_distribution, pitch_distribution, Ro, Zo, r

	! Allocate array containing variables of particles for each species
	ALLOCATE(spp(params%num_species))

	ALLOCATE(ppp(params%num_species))
	ALLOCATE(q(params%num_species))
	ALLOCATE(m(params%num_species))
	ALLOCATE(Eo(params%num_species))
	ALLOCATE(etao(params%num_species))
	ALLOCATE(runaway(params%num_species))
	ALLOCATE(energy_distribution(params%num_species))
	ALLOCATE(pitch_distribution(params%num_species))
	ALLOCATE(Ro(params%num_species))
	ALLOCATE(Zo(params%num_species))
	ALLOCATE(r(params%num_species))

	open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
	read(default_unit_open,nml=plasma_species)
	close(default_unit_open)

	do ii=1_idef,params%num_species
		spp(ii)%Eo = Eo(ii)*C_E
		spp(ii)%etao = etao(ii)
		spp(ii)%runaway = runaway(ii)
		spp(ii)%energy_distribution = TRIM(energy_distribution(ii))
		spp(ii)%pitch_distribution = TRIM(pitch_distribution(ii))
		spp(ii)%q = q(ii)*C_E
		spp(ii)%m = m(ii)*C_ME
		spp(ii)%ppp = ppp(ii)

		spp(ii)%Ro = Ro(ii)
		spp(ii)%Zo = Zo(ii)
		spp(ii)%r = r(ii)

		spp(ii)%go = (spp(ii)%Eo + spp(ii)%m*C_C**2)/(spp(ii)%m*C_C**2)

		ALLOCATE( spp(ii)%vars%X(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%V(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%Rgc(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%Y(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%E(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%B(3,spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%g(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%eta(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%mu(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%Prad(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%Pin(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%flag(spp(ii)%ppp) )
		ALLOCATE( spp(ii)%vars%AUX(spp(ii)%ppp) )

		SELECT CASE (TRIM(spp(ii)%energy_distribution))
			CASE ('MONOENERGETIC')
				spp(ii)%vars%g = spp(ii)%go ! Monoenergetic
			CASE ('AVALANCHE')
				call get_avalanche_PDF_params(params,spp(ii)%vars%g,spp(ii)%vars%eta)
				spp(ii)%go = SUM(spp(ii)%vars%g)/SIZE(spp(ii)%vars%g)
				spp(ii)%Eo = spp(ii)%m*C_C**2*spp(ii)%go - spp(ii)%m*C_C**2
			CASE DEFAULT
				if (params%mpi_params%rank .EQ. 0) then
					write(6,'(/,"* * * * * * * * * * * * * * * * * * * * * * * * *")')
					write(6,'("Energy distribution of species ",I4 " is: ",A50)') ii,TRIM(spp(ii)%energy_distribution)
					write(6,'("* * * * * * * * * * * * * * * * * * * * * * * * *",/)')
				end if
		END SELECT

		SELECT CASE (TRIM(spp(ii)%pitch_distribution))
			CASE ('MONOPITCH')
				spp(ii)%vars%eta = spp(ii)%etao ! Mono-pitch-angle
			CASE ('AVALANCHE')
				spp(ii)%etao = SUM(spp(ii)%vars%eta)/SIZE(spp(ii)%vars%eta)
			CASE DEFAULT
				if (params%mpi_params%rank .EQ. 0) then
					write(6,'(/,"* * * * * * * * * * * * * * * * * * * * * * * * *")')
					write(6,'("The energy distribution will be: ",A50)') TRIM(spp(ii)%pitch_distribution)
					write(6,'("* * * * * * * * * * * * * * * * * * * * * * * * *",/)')
				end if
		END SELECT
	
		! Initialize to zero
		spp(ii)%vars%X = 0.0_rp
		spp(ii)%vars%V = 0.0_rp
		spp(ii)%vars%Rgc = 0.0_rp
		spp(ii)%vars%Y = 0.0_rp
		spp(ii)%vars%E = 0.0_rp
		spp(ii)%vars%B = 0.0_rp
		spp(ii)%vars%mu = 0.0_rp
		spp(ii)%vars%Prad = 0.0_rp
		spp(ii)%vars%Pin = 0.0_rp
		spp(ii)%vars%flag = 1_idef
		spp(ii)%vars%AUX = 0.0_rp
	end do

	DEALLOCATE(ppp)
	DEALLOCATE(q)
	DEALLOCATE(m)
	DEALLOCATE(Eo)
	DEALLOCATE(etao)
	DEALLOCATE(runaway)
	DEALLOCATE(Ro)
	DEALLOCATE(Zo)
	DEALLOCATE(r)
end subroutine initialize_particles


FUNCTION fth_3V(Vth,V)
	IMPLICIT NONE
	REAL(rp), DIMENSION(3), INTENT(IN) :: V
    REAL(rp), INTENT(IN) :: Vth
	REAL(rp) :: fth_3V

    fth_3V = EXP(-0.5_rp*DOT_PRODUCT(V,V)/Vth**2.0_rp)
END FUNCTION fth_3V


FUNCTION fth_1V(Vth,V)
	IMPLICIT NONE
	REAL(rp), INTENT(IN) :: V
    REAL(rp), INTENT(IN) :: Vth
	REAL(rp) :: fth_1V

    fth_1V = EXP(-0.5_rp*(V/Vth)**2)
END FUNCTION fth_1V


FUNCTION random_norm(mean,sigma)
	REAL(rp), INTENT(IN) :: mean
	REAL(rp), INTENT(IN) :: sigma
	REAL(rp) :: random_norm
	REAL(rp) :: rand1, rand2

	call RANDOM_NUMBER(rand1)
	call RANDOM_NUMBER(rand2)

	random_norm = SQRT(-2.0_rp*LOG(1.0_rp-rand1))*COS(2.0_rp*C_PI*rand2);
END FUNCTION random_norm


subroutine iso_thermal_distribution(params,spp)
    implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
    TYPE(SPECIES), INTENT(INOUT) :: spp
    REAL(rp) :: Vmax,Vth, sv
    REAL(rp) :: ratio, rand_unif
    REAL(rp), DIMENSION(3) :: V, U
    REAL(rp), DIMENSION(3) :: b = (/1.0_rp,0.0_rp,0.0_rp/)
	INTEGER :: ii,ppp

	Vmax = 0.8_rp
    Vth = SQRT(spp%Eo*ABS(spp%q)/spp%m)
    ppp = spp%ppp

    V = (/0.0_rp,0.0_rp,0.0_rp/)
    sv = Vth/10.0_rp

    ii=2_idef
	do while (ii .LE. 1000_idef)
		U(1) = V(1) + random_norm(0.0_rp,sv)
		do while (ABS(U(1)) .GT. Vmax)
			U(1) = V(1) + random_norm(0.0_rp,sv)
		end do
		U(2) = V(2) + random_norm(0.0_rp,sv)
		do while (ABS(U(2)) .GT. Vmax)
			U(2) = V(2) + random_norm(0.0_rp,sv)
		end do
		U(3) = V(3) + random_norm(0.0_rp,sv)
		do while (ABS(U(3)) .GT. Vmax)
			U(3) = V(3) + random_norm(0.0_rp,sv)
		end do

		ratio = fth_3V(Vth,U)/fth_3V(Vth,V)

		if (ratio .GE. 1.0_rp) then
			V = U
			ii = ii + 1_idef
		else 
			call RANDOM_NUMBER(rand_unif)
			if (ratio .GT. rand_unif) then
				V = U
				ii = ii + 1_idef
			end if
		end if
	end do	

    spp%vars%V(:,1) = V
    ii=2_idef
	do while (ii .LE. ppp)
		U(1) = spp%vars%V(1,ii-1) + random_norm(0.0_rp,sv)
		do while (ABS(U(1)) .GT. Vmax)
			U(1) = spp%vars%V(1,ii-1) + random_norm(0.0_rp,sv)
		end do
		U(2) = spp%vars%V(2,ii-1) + random_norm(0.0_rp,sv)
		do while (ABS(U(2)) .GT. Vmax)
			U(2) = spp%vars%V(2,ii-1) + random_norm(0.0_rp,sv)
		end do
		U(3) = spp%vars%V(3,ii-1) + random_norm(0.0_rp,sv)
		do while (ABS(U(3)) .GT. Vmax)
			U(3) = spp%vars%V(3,ii-1) + random_norm(0.0_rp,sv)
		end do

		ratio = fth_3V(Vth,U)/fth_3V(Vth,spp%vars%V(:,ii-1))

		if (ratio .GE. 1.0_rp) then
			spp%vars%V(:,ii) = U
			ii = ii + 1_idef
		else 
			call RANDOM_NUMBER(rand_unif)
			if (ratio .GT. rand_unif) then
				spp%vars%V(:,ii) = U
				ii = ii + 1_idef
			end if
		end if
	end do

    do ii=1_idef,ppp
        spp%vars%g(ii) = 1.0_rp/SQRT(1.0_rp - SUM(spp%vars%V(:,ii)**2,1))
        spp%vars%eta(ii) = ACOS(DOT_PRODUCT(b,spp%vars%V(:,ii)/SQRT(SUM(spp%vars%V(:,ii)**2,1))))
    end do
end subroutine iso_thermal_distribution


subroutine set_up_particles_ic(params,F,spp) 
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(IN) :: F
	TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
	REAL(rp), DIMENSION(:), ALLOCATABLE :: Vo
	REAL(rp), DIMENSION(:), ALLOCATABLE :: V1
	REAL(rp), DIMENSION(:), ALLOCATABLE :: V2
	REAL(rp), DIMENSION(:), ALLOCATABLE :: V3
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: b1, b2, b3
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Xo
	REAL(rp), DIMENSION(:), ALLOCATABLE :: theta, zeta, R ! temporary vars
	REAL(rp), DIMENSION(3) :: x = (/1.0_rp,0.0_rp,0.0_rp/)
	REAL(rp), DIMENSION(3) :: y = (/0.0_rp,1.0_rp,0.0_rp/)
	REAL(rp), DIMENSION(3) :: z = (/0.0_rp,0.0_rp,1.0_rp/)
	INTEGER :: ii,jj ! Iterator

	do ii=1_idef,params%num_species
		ALLOCATE( Xo(3,spp(ii)%ppp) )
		ALLOCATE( Vo(spp(ii)%ppp) )
		ALLOCATE( V1(spp(ii)%ppp) )
		ALLOCATE( V2(spp(ii)%ppp) )
		ALLOCATE( V3(spp(ii)%ppp) )
		ALLOCATE( b1(3,spp(ii)%ppp) )
		ALLOCATE( b2(3,spp(ii)%ppp) )
		ALLOCATE( b3(3,spp(ii)%ppp) )

		
		ALLOCATE( theta(spp(ii)%ppp) )
		ALLOCATE( zeta(spp(ii)%ppp) )
		ALLOCATE( R(spp(ii)%ppp) )

		! * * * * INITIALIZE POSITION * * * * 
		if (params%magnetic_field_model .EQ. 'UNIFORM') then
			spp(ii)%vars%X = 0.0_rp
		else
			! Initial condition of uniformly distributed particles on a disk in the xz-plane
			! A unique velocity direction
			call init_u_random(10986546_8)

			call init_random_seed()
			call RANDOM_NUMBER(theta)
			theta = 2.0_rp*C_PI*theta

			call init_random_seed()
			call RANDOM_NUMBER(zeta)
			zeta = 2.0_rp*C_PI*zeta

			! Uniform distribution on a disk at a fixed azimuthal theta		
			call init_random_seed()
			call RANDOM_NUMBER(R)
		
			Xo(1,:) = ( spp(ii)%Ro + spp(ii)%r*SQRT(R)*COS(theta) )*SIN(zeta)
			Xo(2,:) = ( spp(ii)%Ro + spp(ii)%r*SQRT(R)*COS(theta) )*COS(zeta)
			Xo(3,:) = spp(ii)%Zo + spp(ii)%r*SQRT(R)*SIN(theta)

			spp(ii)%vars%X(1,:) = Xo(1,:)
			spp(ii)%vars%X(2,:) = Xo(2,:)
			spp(ii)%vars%X(3,:) = Xo(3,:)
		end if


		! * * * * INITIALIZE VELOCITY * * * * 
		if ((TRIM(spp(ii)%energy_distribution)).EQ.'THERMAL') then
            call iso_thermal_distribution(params,spp(ii))
		else
			call init_random_seed()
			call RANDOM_NUMBER(theta)
			theta = 2.0_rp*C_PI*theta

			Vo = SQRT( 1.0_rp - 1.0_rp/(spp(ii)%vars%g(:)**2) )
		    V1 = Vo*COS(C_PI*spp(ii)%vars%eta/180.0_rp)
		    V2 = Vo*SIN(C_PI*spp(ii)%vars%eta/180.0_rp)*COS(theta)
		    V3 = Vo*SIN(C_PI*spp(ii)%vars%eta/180.0_rp)*SIN(theta)

		    call unitVectors(params,Xo,F,b1,b2,b3,spp(ii)%vars%flag)

			do jj=1_idef,spp(ii)%ppp
				if ( spp(ii)%vars%flag(jj) .EQ. 1_idef ) then
					spp(ii)%vars%V(1,jj) = V1(jj)*DOT_PRODUCT(b1(:,jj),x) + &
				                            V2(jj)*DOT_PRODUCT(b2(:,jj),x) + &
				                            V3(jj)*DOT_PRODUCT(b3(:,jj),x)

					spp(ii)%vars%V(2,jj) = V1(jj)*DOT_PRODUCT(b1(:,jj),y) + &
				                            V2(jj)*DOT_PRODUCT(b2(:,jj),y) + &
				                            V3(jj)*DOT_PRODUCT(b3(:,jj),y)

					spp(ii)%vars%V(3,jj) = V1(jj)*DOT_PRODUCT(b1(:,jj),z) + &
				                            V2(jj)*DOT_PRODUCT(b2(:,jj),z) + &
				                            V3(jj)*DOT_PRODUCT(b3(:,jj),z)
				end if
			end do
		end if

		DEALLOCATE(theta)
		DEALLOCATE(zeta)
		DEALLOCATE(R)
		DEALLOCATE(Xo)
		DEALLOCATE(Vo)
		DEALLOCATE(V1)
		DEALLOCATE(V2)
		DEALLOCATE(V3)
		DEALLOCATE(b1)
		DEALLOCATE(b2)
		DEALLOCATE(b3)
	end do
end subroutine set_up_particles_ic


! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! ** SUBROUTINES FOR INITIALIZING COMMUNICATIONS  ** !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !

subroutine initialize_communications(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(INOUT) :: params

	call initialize_mpi(params)

!$OMP PARALLEL SHARED(params)
        params%num_omp_threads = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL

	call initialization_sanity_check(params) 
end subroutine initialize_communications


subroutine initialization_sanity_check(params)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	CHARACTER(MAX_STRING_LENGTH) :: env_variable
	INTEGER :: ierr, mpierr
	LOGICAL :: flag = .FALSE.

!	call GET_ENVIRONMENT_VARIABLE("OMP_PLACES",env_variable)
!	call GET_ENVIRONMENT_VARIABLE("GOMP_CPU_AFFINITY",env_variable)
!	write(6,*) TRIM(env_variable)

	call MPI_INITIALIZED(flag, ierr)

!$OMP PARALLEL SHARED(params) FIRSTPRIVATE(ierr,flag)
	!$OMP CRITICAL
	write(6,'("MPI: ",I3," OMP/of: ",I3," / ",I3," Procs: ",I3," Init: ",l1)') &
	params%mpi_params%rank,OMP_GET_THREAD_NUM(),OMP_GET_NUM_THREADS(),OMP_GET_NUM_PROCS(),flag
	!$OMP END CRITICAL
!$OMP END PARALLEL
end subroutine initialization_sanity_check


! * * * * * * * * * * * *  * * * * * * * * * * * * * !
! * * *  SUBROUTINES FOR INITIALIZING FIELDS   * * * !
! * * * * * * * * * * * *  * * * * * * * * * * * * * !

subroutine initialize_fields(params,F)
	implicit none
	TYPE(KORC_PARAMS), INTENT(IN) :: params
	TYPE(FIELDS), INTENT(OUT) :: F
	TYPE(KORC_STRING) :: field
	REAL(rp) :: Bo
	REAL(rp) :: minor_radius
	REAL(rp) :: major_radius
	REAL(rp) :: qa
	REAL(rp) :: qo
    CHARACTER(MAX_STRING_LENGTH) :: electric_field_mode
	REAL(rp) :: Eo
    REAL(rp) :: pulse_maximum
    REAL(rp) :: pulse_duration

	NAMELIST /analytic_mag_field_params/ Bo,minor_radius,major_radius,&
			qa,qo,electric_field_mode,Eo,pulse_maximum,pulse_duration

	SELECT CASE (TRIM(params%magnetic_field_model))
		CASE('ANALYTICAL')
			! Load the parameters of the analytical magnetic field
			open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
			read(default_unit_open,nml=analytic_mag_field_params)
			close(default_unit_open)

			F%AB%Bo = Bo
			F%AB%a = minor_radius
			F%AB%Ro = major_radius
			F%Ro = major_radius
			F%AB%qa = qa
			F%AB%qo = qo
			F%AB%lambda = F%AB%a/SQRT(qa/qo - 1.0_rp)
			F%AB%Bpo = F%AB%lambda*F%AB%Bo/(F%AB%qo*F%AB%Ro)

			F%Eo = Eo
			F%Bo = F%AB%Bo

		    F%electric_field_mode = TRIM(electric_field_mode)
			F%to = pulse_maximum
			F%sig = pulse_duration
		CASE('EXTERNAL')
			! Load the magnetic field from an external HDF5 file
		    call load_dim_data_from_hdf5(params,F%dims)

		   	call ALLOCATE_FIELDS_ARRAYS(F,params%poloidal_flux)

		    call load_field_data_from_hdf5(params,F)

			if (.NOT. params%poloidal_flux) then
				field%str = 'B'
				call mean_F_field(F,F%Bo,field)
			end if
		CASE('UNIFORM')
			! Load the parameters of the analytical magnetic field
			open(unit=default_unit_open,file=TRIM(params%path_to_inputs),status='OLD',form='formatted')
			read(default_unit_open,nml=analytic_mag_field_params)
			close(default_unit_open)

			F%Eo = Eo
			F%Bo = Bo
		CASE DEFAULT
	END SELECT
end subroutine initialize_fields

end module korc_initialize
