program PolyEval_test
	use omp_lib
	use modPolyEval
	implicit none

	real (kind=prec) :: x, y, z, estrin, horner
	real (kind=prec) :: startTime, endTime
	real (kind=prec) :: timeDiff, averageIterTime
	integer :: loop, numcoeff
	type(polynomial) :: poly
	type(polynomial2) :: poly2
	type(polynomial3) :: poly3
	real (kind=prec), allocatable, dimension(:) :: coeff

	numcoeff = 6
	x = 2.0
	y = 3.0
	z = 4.0
	allocate(coeff(numcoeff))
	coeff = (/4,3,3,2,1,-80/)

	poly%x = x
	poly%n = numcoeff
	allocate(poly%f(poly%n))
	poly%f = coeff(:)

	write(*,*) 'Horners Form (1 variable)'

	startTime = omp_get_wtime()
	horner = EvalHorner(poly)
	endTime = omp_get_wtime()

	timeDiff = endTime - startTime
    write(*,*) 'result =', horner
    write(*,*) 'Completed in', timeDiff

	write(*,*)
    write(*,*) 'Horners Form (2 variables)'

    poly2%x = x
    poly2%y = y
	poly2%n = numcoeff

	allocate(poly2%f(poly2%n))
	allocate(poly2%powers(2,poly2%n-1))

	poly2%f = coeff(:)
	poly2%powers(1,:) = (/3,2,2/)
	poly2%powers(2,:) = (/1,0,1/)

	startTime = omp_get_wtime()
	horner = EvalHorner(poly2)
	endTime = omp_get_wtime()

	timeDiff = endTime - startTime
    write(*,*) 'result =', horner
    write(*,*) 'Completed in', timeDiff

    write(*,*)
    write(*,*) 'Horners Form (3 variables)'

    poly3%x = x
    poly3%y = y
    poly3%z = z
	poly3%n = numcoeff

	allocate(poly3%f(poly3%n))
	allocate(poly3%powers(3,poly3%n-1))

	poly3%f = coeff(:)
	poly3%powers(1,:) = (/5,4,3,1,0/)
	poly3%powers(2,:) = (/0,4,1,1,0/)
	poly3%powers(3,:) = (/0,4,3,0,1/)

	startTime = omp_get_wtime()
	horner = EvalHorner(poly3)
	endTime = omp_get_wtime()

	timeDiff = endTime - startTime
    write(*,*) 'result =', horner
    write(*,*) 'Completed in', timeDiff

	write(*,*)
	write(*,*) 'Estrins Method (1 variable)'

	startTime = omp_get_wtime()
	estrin = EvalEstrin(poly)
	endTime = omp_get_wtime()

	timeDiff = endTime - startTime
    write(*,*) 'result =', estrin
	write(*,*) 'Completed in', timeDiff

	!write(*,*)
	!write(*,*) 'diff   =', horner-estrin

end program PolyEval_test
