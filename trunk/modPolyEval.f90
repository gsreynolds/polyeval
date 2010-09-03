module modPolyEval
    implicit none

	integer, parameter :: prec = kind(1.0d0)
	integer, parameter :: long = selected_int_kind(12)

   	type polynomial
   		integer :: n
   		real (kind=prec), allocatable, dimension(:) :: f !Coefficients of polynomial
		real (kind=prec) :: x !Value of independent variable
	end type

	type polynomial_multi
		integer :: n
		integer :: m
		real (kind=prec), allocatable, dimension(:) :: f !Coefficients of polynomial
		real (kind=prec), allocatable, dimension(:) :: vars !Independent variables
		integer (kind=long), allocatable, dimension(:,:) :: powers
	end type

	!Evaluate by brute force
	interface Eval
		module procedure Eval, Eval_multi
	end interface

	!Evaluate by optimised brute force method
	interface EvalOpt
		module procedure EvalOpt, EvalOpt_multi
	end interface

    contains

    !Function to evaluate a univariate polynomial by brute force
   	double precision function Eval(poly)
   		type(polynomial), intent(IN) :: poly
   		integer :: i
   		real (kind=prec), dimension(poly%n) :: monomial

   		do i = 1, poly%n-1
			!Build up the monomial value by taking the coefficient
   			!and then multiplying by the variable values raised to the respective power
			monomial(i) = poly%f(i)*poly%x**(poly%n-i)
   		end do

   		monomial(poly%n) = poly%f(poly%n)

   		Eval = sum(monomial(:))

   	end function Eval

   	 !Function to evaluate a univariate polynomial by brute force, with optimisations
   	double precision function EvalOpt(poly)
   		type(polynomial), intent(IN) :: poly
   		integer :: i, j, numSteps
   		real (kind=prec), dimension(poly%n) :: monomial
   		real (kind=prec) :: x2

   		x2 = poly%x*poly%x

   		do i = 1, poly%n-1
			monomial(i) = poly%f(i)
			!Build up the power required from x^2 and x
			!i.e. x^7 = x^2 * x^2 * x^2 * x
			numSteps = (poly%n-i)/2
			do j = 1, numSteps
				monomial(i) = monomial(i) * x2
			end do
			if (mod(poly%n-i, 2) .ne. 0) then
				monomial(i) = monomial(i) * poly%x
   			end if
   		end do

   		monomial(poly%n) = poly%f(poly%n)

   		EvalOpt = sum(monomial(:))

   	end function EvalOpt

	!Function to evaluate a multivariate polynomial by brute force
   	double precision function Eval_multi(poly)
   		type(polynomial_multi), intent(IN) :: poly
   		integer :: i, j
   		real (kind=prec), dimension(poly%n) :: monomial

   		do i = 1, poly%n-1
   			!Build up the monomial value by taking the coefficient
   			!and then multiplying by the variable values raised to the respective power
			monomial(i) = poly%f(i)
			do j = 1, poly%m
				monomial(i) = monomial(i)*poly%vars(j)**poly%powers(j,i)
			end do
   		end do

   		monomial(poly%n) = poly%f(poly%n)

   		Eval_multi = sum(monomial(:))

   	end function Eval_multi

	!Function to evaluate a multivariate polynomial by brute force, with optimisations
   	double precision function EvalOpt_multi(poly)
   		type(polynomial_multi), intent(IN) :: poly
   		integer :: i, j, k, numSteps
   		real (kind=prec), dimension(poly%n) :: monomial
   		real (kind=prec), dimension(poly%m) :: vars2 !variables, squared

		!Pre-calculate all the x^2, y^2 etc variables
   		do i = 1, poly%m
   			vars2(i) = poly%vars(i)*poly%vars(i)
   		end do

   		do i = 1, poly%n-1
			monomial(i) = poly%f(i)

			do j = 1, poly%m
				!Build up the power required from x^2 and x
				!i.e. x^7 = x^2 * x^2 * x^2 * x
				numSteps = (poly%powers(j,i))/2
				do k = 1, numSteps
					monomial(i) = monomial(i) * vars2(j)
				end do
				if (mod(poly%powers(j,i), 2) .ne. 0) then
					monomial(i) = monomial(i) * poly%vars(j)
	   			end if
   			end do
   		end do

   		monomial(poly%n) = poly%f(poly%n)

   		EvalOpt_multi = sum(monomial(:))

   	end function EvalOpt_multi

	!Function to evaluate a polynomial using Horner's form
    double precision function EvalHorner(poly)

		type(polynomial), intent(IN) :: poly
		integer :: i

		!Evaluate using Horner's Method
		EvalHorner = poly%f(1)*poly%x
		do i = 2, poly%n-1
			EvalHorner = (EvalHorner + poly%f(i))*poly%x
		end do
		EvalHorner = EvalHorner + poly%f(poly%n)

	end function EvalHorner

	!Function to evaluate a polynomial using Estrin's method
	double precision function EvalEstrin(poly)
		use omp_lib

		type(polynomial), intent(IN) :: poly
		real (kind=prec), allocatable, dimension(:) :: powers
		real (kind=prec), allocatable, dimension(:,:) :: coeff
		integer :: i,j,numsteps, shift, nearestpoweroftwo, npow2

		!Round up the number of coefficients to the nearest power of two
		nearestpoweroftwo = 2**ceiling(log(real(poly%n))/log(2.0d0))
		if (mod(poly%n, nearestpoweroftwo) .ne. 0) then
			shift = nearestpoweroftwo-mod(poly%n, nearestpoweroftwo)
			npow2 = poly%n + shift
		else
			shift = 0
			npow2 = poly%n
		end if

		allocate(coeff(npow2, 0:npow2))
		coeff(:,:) = 0

		!Calculate the number of calculation steps required to reach the result
		numsteps = log(real(npow2))/log(2.0d0)
		allocate(powers(numsteps))

		!Build the powers array
		powers(1) = poly%x
		do i = 2, numsteps
			powers(i) = powers(i-1)**2
		end do

		!Populate the first (0) row of the coeff array with the polynomial coefficients
		coeff(1+shift:npow2, 0) = poly%f(:)

		!Evaluate using Estrin's Method
		do i = 1, numsteps

			!$omp parallel do default(none) shared(coeff, powers, npow2, i) private(j)
			do j = 1, npow2/(2**i)
				!write(*,*) 'Thread', OMP_GET_THREAD_NUM(), 'i', i, 'j', j
				coeff(j, i) = coeff(2*j-1, i-1)*powers(i)+coeff(2*j, i-1)
			end do
			!$omp end parallel do

		end do

		EvalEstrin = coeff(1, numsteps)

	end function EvalEstrin

end module modPolyEval
