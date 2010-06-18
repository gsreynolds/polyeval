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

	interface Eval
		module procedure Evalx, Evalxy, Evalxyz
	end interface

	interface EvalHorner
		module procedure EvalHornerx, EvalHornerxy, EvalHornerxyz
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

	!Function to evaluate a univariate polynomial using Horner's form
    double precision function EvalHornerx(poly)

		type(polynomial), intent(IN) :: poly
		integer :: i

		!Evaluate using Horner's Method
		EvalHornerx = poly%f(1)*poly%x
		do i = 2, poly%n-1
			EvalHornerx = (EvalHornerx + poly%f(i))*poly%x
		end do
		EvalHornerx = EvalHornerx + poly%f(poly%n)

	end function EvalHornerx

	!Function to evaluate a bivariate polynomial using Horner's form
    double precision function EvalHornerxy(poly)

		type(polynomial2), intent(IN) :: poly
		real(kind=prec), dimension(2) :: variables
		real(kind=prec), dimension(poly%n) :: b
		real(kind=prec) :: temp
		integer :: i, j

		variables(1) = poly%x
		variables(2) = poly%y

		b(1) = poly%f(1)
		do i = 2,poly%n
			if(poly%powers(1,i) .eq. 1) then
				temp = variables(1)
			else
				temp = 1
			end if
			do j = 2,2
				if(poly%powers(j,i) .eq. 1) then
			    	temp = variables(j)*temp
			    end if
			end do
			b(i) = poly%f(i) + temp*b(i-1)
		end do

		EvalHornerxy = sum(b)

	end function EvalHornerxy

	!Function to evaluate a trivariate polynomial using Horner's form
    double precision function EvalHornerxyz(poly)

		type(polynomial3), intent(IN) :: poly
		real(kind=prec), dimension(3) :: variables
		real(kind=prec), dimension(poly%n) :: b
		real(kind=prec) :: temp
		integer :: i, j

		variables(1) = poly%x
		variables(2) = poly%y
		variables(3) = poly%z

		b(1) = poly%f(1)
		do i = 2,poly%n
			if(poly%powers(1,i) .eq. 1) then
				temp = variables(1)
			else
				temp = 1
			end if
			do j = 2,3
				if(poly%powers(j,i) .eq. 1) then
			    	temp = variables(j)*temp
			    end if
			end do
			b(i) = poly%f(i) + temp*b(i-1)
		end do

		EvalHornerxyz = sum(b)

	end function EvalHornerxyz

	!Function to evaluate a univariate polynomial using Estrin's method
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
