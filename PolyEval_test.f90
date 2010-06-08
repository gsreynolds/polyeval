program PolyEval_test
	use omp_lib
	use modPolyEval
	implicit none

	real (kind=prec), parameter :: a = 3.068474, b=0.0d0, c=20.857847, d=0.0d0, e=0.757463, f=0.0d0, g=8.673527
	real (kind=prec), parameter :: h=765.638467, i=0.0d0, j=-20.889708, k=67.786429, l=-0.754380, m= 1120.000000
	real (kind=prec) :: x, estrin, horner
	real (kind=prec) :: startTime, endTime
	real (kind=prec) :: timeDiff, averageIterTime
	integer :: loop
	type(polynomial) poly

	poly%x = 2.0
	poly%n = 128
	allocate(poly%f(poly%n))
	poly%f = (/0.0d0, 0.0d0, 0.0d0,m,l,k,j, i, h, g, f, e, d, c, b, a, &
		0.0d0, 0.0d0, 0.0d0,m,l,k,j, i, h, g, f, e, d, c, b, a, &
		0.0d0, 0.0d0, 0.0d0,m,l,k,j, i, h, g, f, e, d, c, b, a, &
		0.0d0, 0.0d0, 0.0d0,m,l,k,j, i, h, g, f, e, d, c, b, a, &
		0.0d0, 0.0d0, 0.0d0,m,l,k,j, i, h, g, f, e, d, c, b, a, &
		0.0d0, 0.0d0, 0.0d0,m,l,k,j, i, h, g, f, e, d, c, b, a, &
		0.0d0, 0.0d0, 0.0d0,m,l,k,j, i, h, g, f, e, d, c, b, a, &
		0.0d0, 0.0d0, 0.0d0,m,l,k,j, i, h, g, f, e, d, c, b, a/)

	write(*,*) 'Horners Form'

	startTime = omp_get_wtime()
	horner = EvalHorner(poly)
	endTime = omp_get_wtime()

	timeDiff = endTime - startTime
    write(*,*) 'result =', horner
    write(*,*) 'Completed in', timeDiff

	write(*,*)
	write(*,*) 'Estrins Method'

	startTime = omp_get_wtime()
	estrin = EvalEstrin(poly)
	endTime = omp_get_wtime()

	timeDiff = endTime - startTime
    write(*,*) 'result =', estrin
	write(*,*) 'Completed in', timeDiff

	write(*,*)
	write(*,*) 'diff   =', horner-estrin

end program PolyEval_test
