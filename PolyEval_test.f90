program PolyEval_test
	use mpi
	use modPolyEval
	implicit none

	integer, parameter :: MAXITER = 1!000
	real(kind=prec), allocatable, dimension(:) :: rndArray
	integer :: i, numTerms
	integer, parameter :: numTestCases = 8
	integer, dimension(numTestCases) :: testCases
	type(polynomial) :: poly
	
	!MPI variables
	integer :: ierr, comm, rank, size

	!Initialise MPI
	comm = MPI_COMM_WORLD

	call MPI_Init(ierr)

    call MPI_Comm_Rank(comm, rank, ierr)

    call MPI_Comm_Size(comm, size, ierr)

	!Number of processes must either be 1 or a power of two
	if (size .ne. 1) then
	    if ((iand(size, (size - 1))) .ne. 0) then
	    	if (rank==0) then
	    		write(*,*) 'Number of processes must be a power of two'
	    	end if
	    	call MPI_Finalize(ierr)
	    	stop
	    end if
    end if

	testCases = (/5,15,30,100,500,1000,2000,4000/)

	do i = 1, numTestCases

		if (i > 1) then
			deallocate(poly%f)
			deallocate(rndArray)
		end if

		numTerms = testCases(i)

		poly%x = 1.1
		poly%n = numTerms
		allocate(poly%f(poly%n))
		allocate(rndArray(poly%n))

		if (rank == 0) then

			write(*,*) 'Test Case ', i
			write(*,*) 'Number of terms: ', numTerms
			write(*,*)

			call generateRandomArray(poly%n, rndArray)

		end if

		call MPI_Bcast(rndArray, numTerms, MPI_DOUBLE_PRECISION, 0, comm, ierr)

		poly%f(:) = rndArray(:)

		call evaluationMethods(poly, rank)

		if (rank == 0) then
			write(*,*)
			write(*,*) '==============='
			write(*,*)
		end if

	end do
	
	call MPI_Finalize(ierr)

contains

	subroutine generateRandomArray(size, rndArray)
		integer, intent(IN) :: size
		real(kind=prec), allocatable, intent(OUT), dimension(:) :: rndArray
		integer :: seedSize,date(8)
		integer, allocatable  :: seed(:)

		call date_and_time(values=date)
	 	call random_seed(SIZE=seedSize)
	 	allocate( seed(seedSize) )
		call random_seed(GET=seed)
	 	seed = seed * date(8)      ! date(8) is milliseconds
	 	call random_seed(PUT=seed)

		allocate(rndArray(size))

		call random_number(rndArray)

		rndArray(:) = rndArray(:) * 10

	end subroutine generateRandomArray

	subroutine evaluationMethods(poly, rank)

		type(polynomial), intent(IN) :: poly
		integer, parameter :: numMethods =  4
		real (kind=prec), dimension(numMethods) :: results, averageIterTime
		real (kind=prec) :: startTime, endTime, timeDiff
		integer :: i, loop, rank

		if (rank == 0) then

			!Brute force

			startTime = MPI_Wtime()
			do loop = 1, MAXITER
				results(1) = Eval(poly)
			end do
			endTime = MPI_Wtime()

			timeDiff = endTime - startTime
			averageIterTime(1) = timeDiff/MAXITER

		    !Brute force (optimised)

			startTime = MPI_Wtime()
			do loop = 1, MAXITER
				results(2) = EvalOpt(poly)
			end do
			endTime = MPI_Wtime()

			timeDiff = endTime - startTime
			averageIterTime(2) = timeDiff/MAXITER

			!Horners Form

			startTime = MPI_Wtime()
			do loop = 1, MAXITER
				results(3) = EvalHorner(poly)
			end do
			endTime = MPI_Wtime()

			timeDiff = endTime - startTime
			averageIterTime(3) = timeDiff/MAXITER

		end if

		!Estrins Method

		startTime = MPI_Wtime()
		do loop = 1, MAXITER
			results(4) = EvalEstrin(poly, rank, size, comm)
		end do
		endTime = MPI_Wtime()

		timeDiff = endTime - startTime
		averageIterTime(4) = timeDiff/MAXITER

		if (rank == 0) then
			write(*,*) 'Brute Force               |  Brute Force - Opt        |  Horners Form             |  Estrins Method'
			write(*,*) 'Results:'
			write(*,*) (results(i),',', i=1,numMethods)
			write(*,*) 'Average Iteration Times:'
			write(*,*) (averageIterTime(i),',', i=1,numMethods)
		end if

	end subroutine evaluationMethods

end program PolyEval_test
