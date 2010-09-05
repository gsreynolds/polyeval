program PolyEval_test
	use omp_lib
	use modPolyEval
	implicit none

	integer, parameter :: MAXITER = 1!000
	real(kind=prec), allocatable, dimension(:) :: rndArray
	integer :: i, numTerms
	integer, parameter :: numTestCases = 8
	integer, dimension(numTestCases) :: testCases
	type(polynomial) :: poly

	testCases = (/5,15,30,100,500,1000,2000,4000/)

	do i = 1, numTestCases

		if (i > 1) then
			deallocate(poly%f)
		end if

		numTerms = testCases(i)

		write(*,*) 'Test Case ', i
		write(*,*) 'Number of terms: ', numTerms
		write(*,*)

		poly%x = 1.1
		poly%n = numTerms
		allocate(poly%f(poly%n))

		call generateRandomArray(poly%n, rndArray)

		poly%f(:) = rndArray(:)
		deallocate(rndArray)

		call evaluationMethods(poly)

		write(*,*)
		write(*,*) '==============='
		write(*,*)

	end do

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

	subroutine evaluationMethods(poly)

		type(polynomial), intent(IN) :: poly
		integer, parameter :: numMethods =  4
		real (kind=prec), dimension(numMethods) :: results, averageIterTime
		real (kind=prec) :: startTime, endTime, timeDiff
		integer :: i, loop

		!Brute force

		startTime = omp_get_wtime()
		do loop = 1, MAXITER
			results(1) = Eval(poly)
		end do
		endTime = omp_get_wtime()

		timeDiff = endTime - startTime
		averageIterTime(1) = timeDiff/MAXITER

	    !Brute force (optimised)

		startTime = omp_get_wtime()
		do loop = 1, MAXITER
			results(2) = EvalOpt(poly)
		end do
		endTime = omp_get_wtime()

		timeDiff = endTime - startTime
		averageIterTime(2) = timeDiff/MAXITER

		!Horners Form

		startTime = omp_get_wtime()
		do loop = 1, MAXITER
			results(3) = EvalHorner(poly)
		end do
		endTime = omp_get_wtime()

		timeDiff = endTime - startTime
		averageIterTime(3) = timeDiff/MAXITER

		!Estrins Method

		startTime = omp_get_wtime()
		do loop = 1, MAXITER
			results(4) = EvalEstrin(poly)
		end do
		endTime = omp_get_wtime()

		timeDiff = endTime - startTime
		averageIterTime(4) = timeDiff/MAXITER

		write(*,*) 'Brute Force               |  Brute Force - Opt        |  Horners Form             |  Estrins Method'
		write(*,*) 'Results:'
		write(*,*) (results(i),',', i=1,numMethods)
		write(*,*) 'Average Iteration Times:'
		write(*,*) (averageIterTime(i),',', i=1,numMethods)

	end subroutine evaluationMethods

end program PolyEval_test
