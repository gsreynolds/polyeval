program PolyEval_test
	use omp_lib
	use modPolyEval
	implicit none

	real (kind=prec), parameter :: a = 3.068474, b=0.0d0, c=20.857847, d=0.0d0, e=0.757463, f=0.0d0, g=8.673527
	real (kind=prec), parameter :: h=765.638467, i=0.0d0, j=-20.889708, k=67.786429, l=-0.754380, m= 1120.000000
	real (kind=prec) :: x, estrin, horner, brutex, brutexy, brutexyz, hornerxy, hornerxyz
	real (kind=prec) :: startTime, endTime
	real (kind=prec) :: timeDiff, averageIterTime
	integer :: loop
	type(polynomial) :: poly
	type(polynomial2) :: poly2
	type(polynomial3) :: poly3

	poly%x = 2.0
	poly%n = 7
	allocate(poly%f(poly%n))
	poly%f = (/g, f, e, d, c, b, a/)

	write(*,*) 'Brute force x'

	startTime = omp_get_wtime()
	brutex = Eval(poly)
	endTime = omp_get_wtime()

	timeDiff = endTime - startTime
    write(*,*) 'result =', brutex
    write(*,*) 'Completed in', timeDiff
    write(*,*)

    write(*,*) 'Brute force xy'

    poly2%x=2.0
    poly2%y=3.0
    poly2%n=poly%n
    allocate(poly2%f(poly2%n))
    allocate(poly2%powers(2,poly2%n-1))
    poly2%f = poly%f
	poly2%powers(1,:) = (/2,2,1,1,1,0/)
	poly2%powers(2,:) = (/2,1,2,1,0,1/)

	startTime = omp_get_wtime()
	brutexy = Eval(poly2)
	endTime = omp_get_wtime()

	timeDiff = endTime - startTime
    write(*,*) 'result =', brutexy
    write(*,*) 'Completed in', timeDiff
    write(*,*)
    
 	write(*,*) 'Horner xy'

	startTime = omp_get_wtime()
	hornerxy = EvalHorner(poly2)
	endTime = omp_get_wtime()

	timeDiff = endTime - startTime
    write(*,*) 'result =', hornerxy
    write(*,*) 'Completed in', timeDiff
    write(*,*)

    write(*,*) 'Brute force xyz'

    poly3%x=2.0
    poly3%y=3.0
    poly3%z=4.0
    poly3%n=poly%n
    allocate(poly3%f(poly3%n))
    allocate(poly3%powers(3,poly3%n-1))
    poly3%f = poly%f
	poly3%powers(1,:) = (/2,2,1,1,1,0/)
	poly3%powers(2,:) = (/2,1,2,1,0,1/)
	poly3%powers(3,:) = (/0,0,0,0,0,0/)

	startTime = omp_get_wtime()
	brutexyz = Eval(poly3)
	endTime = omp_get_wtime()

	timeDiff = endTime - startTime
    write(*,*) 'result =', brutexyz
    write(*,*) 'Completed in', timeDiff
    write(*,*)

    write(*,*) 'Horner xyz'

	startTime = omp_get_wtime()
	hornerxyz = EvalHorner(poly3)
	endTime = omp_get_wtime()

	timeDiff = endTime - startTime
    write(*,*) 'result =', hornerxyz
    write(*,*) 'Completed in', timeDiff
    write(*,*)

	write(*,*) '==============='
	write(*,*)
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
	write(*,*) 'Horner - Estrin'
	write(*,*) 'diff   =', horner-estrin

end program PolyEval_test
