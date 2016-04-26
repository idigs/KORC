program main

    use constants
    use korc_types
    use main_mpi
    use initialize
    use omp_lib

    implicit none
	TYPE (KORC_PARAMS) :: params

    call initialize_korc(params)

end program main
