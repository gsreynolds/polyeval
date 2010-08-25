program PolyEval_test
	use mpi
	use modPolyEval
	implicit none

	integer, parameter :: MAXITER = 1!000000000
	real (kind=prec), parameter :: a = 3.068474, b=0.0d0, c=20.857847, d=0.0d0, e=0.757463, f=0.0d0, g=8.673527
	real (kind=prec), parameter :: h=765.638467, i=0.0d0, j=-20.889708, k=67.786429, l=-0.754380, m= 1120.000000
	real (kind=prec) :: x, estrin, horner, brute, brute_multi, bruteopt, bruteopt_multi
	real (kind=prec) :: startTime, endTime
	real (kind=prec) :: timeDiff, averageIterTime
	integer :: loop
	type(polynomial) :: poly
	type(polynomial_multi) :: poly_multi

	!MPI variables
	integer, dimension(MPI_STATUS_SIZE) :: status
	integer :: ierr, comm, rank, size

	!Initialise MPI
	comm = MPI_COMM_WORLD

	call MPI_Init(ierr)

    call MPI_Comm_Rank(comm, rank, ierr)

    call MPI_Comm_Size(comm, size, ierr)

	poly%x = 2.0
	poly%n = 7
	allocate(poly%f(poly%n))
	poly%f = (/a,a,a,a,a,a,a/)

	! Only executing Brute force and Horners form on rank 0

	if(rank==0) then

		write(*,*) 'Brute force x'

		startTime = MPI_Wtime()
		do loop = 1, MAXITER
			brute = Eval(poly)
		end do
		endTime = MPI_Wtime()

		timeDiff = endTime - startTime
	    write(*,*) 'result =', brute
	    write(*,*) 'Completed in', timeDiff
	    write(*,*)

	    write(*,*) 'Brute force (optimised) x'

		startTime = MPI_Wtime()
		do loop = 1, MAXITER
			bruteopt = EvalOpt(poly)
		end do
		endTime = MPI_Wtime()

		timeDiff = endTime - startTime
	    write(*,*) 'result =', bruteopt
	    write(*,*) 'Completed in', timeDiff
	    write(*,*)

	    write(*,*) 'Brute force multi'

	    poly_multi%m = 2
		allocate(poly_multi%vars(poly_multi%m))
	    poly_multi%vars(1)=2.0
	    poly_multi%vars(2)=3.0
	    poly_multi%n=poly%n
	    allocate(poly_multi%f(poly_multi%n))
	    allocate(poly_multi%powers(poly_multi%m,poly_multi%n-1))
	    poly_multi%f = poly%f
		poly_multi%powers(1,:) = (/2,2,1,1,1,0/)
		poly_multi%powers(2,:) = (/2,1,2,1,0,1/)

		startTime = MPI_Wtime()
		do loop = 1, MAXITER
			brute_multi = Eval_multi(poly_multi)
		end do
		endTime = MPI_Wtime()

		timeDiff = endTime - startTime
	    write(*,*) 'result =', brute_multi
	    write(*,*) 'Completed in', timeDiff
	    write(*,*)

	    write(*,*) 'Brute force multi (optimised)'
		startTime = MPI_Wtime()
		do loop = 1, MAXITER
			bruteopt_multi = EvalOpt_multi(poly_multi)
		end do
		endTime = MPI_Wtime()

		timeDiff = endTime - startTime
	    write(*,*) 'result =', bruteopt_multi
	    write(*,*) 'Completed in', timeDiff
	    write(*,*)

		write(*,*) '==============='
		write(*,*)
		write(*,*) 'Horners Form'

		startTime = MPI_Wtime()
		do loop = 1, MAXITER
			horner = EvalHorner(poly)
		end do
		endTime = MPI_Wtime()

		timeDiff = endTime - startTime
	    write(*,*) 'result =', horner
	    write(*,*) 'Completed in', timeDiff

    end if

	write(*,*)
	write(*,*) 'Estrins Method'

	startTime = MPI_Wtime()
	do loop = 1, MAXITER
		estrin = EvalEstrin(poly, rank)
	end do
	endTime = MPI_Wtime()

	timeDiff = endTime - startTime
    write(*,*) 'result =', estrin
	write(*,*) 'Completed in', timeDiff

!	write(*,*)
!	write(*,*) 'Horner - Estrin'
!	write(*,*) 'diff   =', horner-estrin

	call MPI_Finalize(ierr)

end program PolyEval_test
