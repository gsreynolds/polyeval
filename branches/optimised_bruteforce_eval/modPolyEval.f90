module modPolyEval
    implicit none

	integer, parameter :: prec = kind(1.0d0)

   	type polynomial
   		integer :: n
   		real (kind=prec), allocatable, dimension(:) :: f !Coefficients of polynomial
		real (kind=prec) :: x !Value of independent variable
	end type

	type polynomial2
   		integer :: n
   		real (kind=prec), allocatable, dimension(:) :: f !Coefficients of polynomial
		real (kind=prec) :: x !Value of independent variable
		real (kind=prec) :: y !Value of independent variable
		real (kind=prec), allocatable, dimension(:,:) :: powers
	end type

	type polynomial3
   		integer :: n
   		real (kind=prec), allocatable, dimension(:) :: f !Coefficients of polynomial
		real (kind=prec) :: x !Value of independent variable
		real (kind=prec) :: y !Value of independent variable
		real (kind=prec) :: z !Value of independent variable
		real (kind=prec), allocatable, dimension(:,:) :: powers
	end type

	!Evaluate by brute force
	interface Eval
		module procedure Evalx, Evalxy, Evalxyz
	end interface

	!Evaluate by optimised brute force method
	interface EvalOpt
		module procedure EvalOptx!, EvalOptxy, EvalOptxyz
	end interface

    contains

    !Function to evaluate a univariate polynomial by brute force
   	double precision function Evalx(poly)
   		type(polynomial), intent(IN) :: poly
   		integer :: i
   		real (kind=prec), dimension(poly%n) :: monomial

   		do i = 1, poly%n-1
			monomial(i) = poly%f(i)*poly%x**(poly%n-i)
   		end do

   		monomial(poly%n) = poly%f(poly%n)

   		Evalx = sum(monomial(:))

   	end function Evalx

   	 !Function to evaluate a univariate polynomial by brute force, with optimisations
   	double precision function EvalOptx(poly)
   		type(polynomial), intent(IN) :: poly
   		integer :: i, j, numSteps
   		real (kind=prec), dimension(poly%n) :: monomial
   		real (kind=prec) :: x2

   		x2 = poly%x*poly%x

   		do i = 1, poly%n-1
			monomial(i) = poly%f(i)
			numSteps = (poly%n-i)/2
			do j = 1, numSteps
				monomial(i) = monomial(i) * x2
			end do
			if (mod(poly%n-i, 2) .ne. 0) then
				monomial(i) = monomial(i) * poly%x
   			end if
   		end do

   		monomial(poly%n) = poly%f(poly%n)

   		EvalOptx = sum(monomial(:))

   	end function EvalOptx

   	 !Function to evaluate a bivariate polynomial by brute force
   	double precision function Evalxy(poly)
   		type(polynomial2), intent(IN) :: poly
   		integer :: i
   		real (kind=prec), dimension(poly%n) :: monomial

   		do i = 1, poly%n-1
			monomial(i) = poly%f(i)*poly%x**poly%powers(1,i)*poly%y**poly%powers(2,i)
   		end do

   		monomial(poly%n) = poly%f(poly%n)

   		Evalxy = sum(monomial(:))

   	end function Evalxy

   	 !Function to evaluate a trivariate polynomial by brute force
   	double precision function Evalxyz(poly)
   		type(polynomial3), intent(IN) :: poly
   		integer :: i
   		real (kind=prec), dimension(poly%n) :: monomial

   		do i = 1, poly%n-1
			monomial(i) = poly%f(i)*poly%x**poly%powers(1,i)*poly%y**poly%powers(2,i)*poly%z**poly%powers(3,i)
   		end do

   		monomial(poly%n) = poly%f(poly%n)

   		Evalxyz = sum(monomial(:))

   	end function Evalxyz

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
		real (kind=prec), allocatable, dimension(:) :: func, powers
		real (kind=prec), allocatable, dimension(:,:) :: coeff
		integer :: i,j,numsteps, shift, nearestpoweroftwo, npow2

		nearestpoweroftwo = 2**ceiling(log(real(poly%n))/log(2.0d0))
		!write(*,*) 'nearestpoweroftwo', nearestpoweroftwo
		if (mod(poly%n, nearestpoweroftwo) .ne. 0) then
			!write(*,*) 'ne 0'
			shift = nearestpoweroftwo-mod(poly%n, nearestpoweroftwo)
			!write(*,*) 'shift', shift
			npow2 = poly%n + shift
			!write(*,*) 'npow2', npow2
		else
			shift = 0
			npow2 = poly%n
		end if

		allocate(func(npow2))
		func(:) = 0
		allocate(coeff(npow2, 0:npow2))
		coeff(:,:) = 0

		numsteps = log(real(npow2))/log(2.0d0)
		!write(*,*) 'numsteps', numsteps
		allocate(powers(numsteps))

		func(1+shift:npow2) = poly%f(:)
		!write(*,*) 'func'
		!write(*,*) func

		!Evaluate using Estrin's Method
		powers(1) = poly%x
		do i = 2, numsteps
			powers(i) = powers(i-1)**2
		end do

		coeff(:, 0) = func

		do i = 1, numsteps

			!$omp parallel do default(none) shared(coeff, powers, npow2, i) private(j)
			do j = 1, npow2/(2**i)
				!write(*,*) 'Thread', OMP_GET_THREAD_NUM(), 'i', i, 'j', j
				coeff(j, i) = coeff(2*j-1, i-1)*powers(i)+coeff(2*j, i-1)
			end do
			!$omp end parallel do

			!write(*,*)
			!write(*,*) coeff(:, i)

		end do

		EvalEstrin = coeff(1, numsteps)

	end function EvalEstrin

end module modPolyEval
