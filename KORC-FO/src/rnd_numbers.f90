module rnd_numbers

#ifdef INTEL
    use ifport
#endif

    use korc_types

    implicit none

    contains

subroutine init_random_seed()
!use iso_fortran_env, only: int64
    implicit none
	INTEGER, allocatable :: seed(:)
	INTEGER(8), DIMENSION(8) :: dt
	INTEGER(8) :: i, istat, pid
	INTEGER(4) :: n
	INTEGER(8) :: t
		      
	call random_seed(size = n)
	allocate(seed(n))
	! First try if the OS provides a random number generator
	open(default_unit_open, file="/dev/urandom", access="stream", &
	form="unformatted", action="read", status="old", iostat=istat)
	if (istat == 0) then
		read(default_unit_open) seed
		close(default_unit_open)
	else
		! Fallback to XOR:ing the current time and pid. The PID is
		! useful in case one launches multiple instances of the same
		! program in parallel.
		call system_clock(t)
	if (t == 0) then
		call date_and_time(values=dt)
		t = (dt(1) - 1970_8) * 365_8 * 24_8 * 60_8 * 60_8 * 1000_8 &
			+ dt(2) * 31_8 * 24_8 * 60_8 * 60_8 * 1000_8 &
			+ dt(3) * 24_8 * 60_8 * 60_8 * 1000_8 &
			+ dt(5) * 60_8 * 60_8 * 1000_8 &
			+ dt(6) * 60_8 * 1000_8 &
			+ dt(7) * 1000_8 &
			+ dt(8)
	end if
		pid = getpid()
		write(6,'("PID: ",I15)') pid
		t = ieor(t, int(pid, kind(t)))
		do i = 1, n
		seed(i) = lcg(t)
		end do
	end if
	call random_seed(put=seed)
	contains

	! This simple PRNG might not be good enough for real work, but is
	! sufficient for seeding a better PRNG.
	function lcg(s)
	INTEGER :: lcg
	INTEGER(8) :: s
	if (s == 0) then
		s = 104729_8
	else
		s = mod(s, 4294967296_8)
	end if
	s = mod(s * 279470273_8, 4294967291_8)
	lcg = int(mod(s, int(huge(0), 8)), kind(0))
	end function lcg

end subroutine init_random_seed

end module rnd_numbers
