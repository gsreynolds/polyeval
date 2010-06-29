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

    double precision function CEvalx(n, f, x, type)
    	integer, intent(IN) :: n, type
    	real(kind=prec), intent(IN), dimension(n) :: f
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

   	double precision function CEvalxy(n, f, x, y, powers1, powers2, type)
    	integer, intent(IN) :: n, type
    	real(kind=prec), intent(IN), dimension(n) :: f
    	real(kind=prec), intent(IN) :: x, y
    	real(kind=prec), intent(IN), dimension(n-1) :: powers1, powers2
    	type(polynomial2) :: poly

		poly%x = x
		poly%y = y
		poly%n = n
		allocate(poly%f(poly%n))
		poly%f(:) = f(:)
		allocate(poly%powers(2,poly%n-1))
		poly%powers(1,:) = powers1(:)
		poly%powers(2,:) = powers2(:)

		select case (type)

			case (0)
				write(*,*) 'Evalxy'
				CEvalxy = Evalxy(poly)

			case (1)
				write(*,*) 'EvalOptxy'
				CEvalxy = EvalOptxy(poly)

!			case (2)
!				write(*,*) 'EvalHorner'
!				CEvalxy = EvalHornerxy(poly)

			case default
				write(*,*) 'Invalid type identifier, defaulting to Evalxy'
				CEvalxy = Evalxy(poly)

		end select

    end function CEvalxy

end module modPolyEvalCInterface
