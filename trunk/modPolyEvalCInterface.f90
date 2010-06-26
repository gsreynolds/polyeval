module modPolyEvalCInterface
	use modPolyEval
    implicit none

    contains

    double precision function CtoFtest(order, eps)
    	integer, intent(IN) :: order !C int
    	real(kind=prec), intent(IN) :: eps !C double

    	write(*,*) 'Fortran function'

    	write(*,*) 'Poly Order', order
    	write(*,*) 'eps', eps

    	CtoFtest = order*eps
    	write(*,*) 'result = ', CtoFtest
    end function CtoFtest

	!Test C interface to Evalx, as opposed to passing a C struct to function
	!with Fortran derived datatype as argument.
    double precision function CEvalx(n, f, x, type)
    	integer, intent(IN) :: n, type
    	real(kind=prec), dimension(n) :: f
    	real (kind=prec), intent(IN) :: x
    	type(polynomial) :: poly

		poly%x = x
		poly%n = n
		allocate(poly%f(poly%n))
		poly%f(:) = f(:)

		select case (type)

			case (0)
				write(*,*) 'Evalx'
				CEvalx = Evalx(poly)

			case (1)
				write(*,*) 'EvalOptx'
				CEvalx = EvalOptx(poly)

			case (2)
				write(*,*) 'EvalHorner'
				CEvalx = EvalHorner(poly)

			case (3)
				write(*,*) 'EvalEstrins'
				CEvalx = EvalEstrin(poly)

			case default
				write(*,*) 'Invalid type identifier, defaulting to Evalx'
				CEvalx = Evalx(poly)

		end select

    end function CEvalx

end module modPolyEvalCInterface
