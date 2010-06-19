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

	!Function to evaluate a polynomial using Horner's form
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

    double precision function EvalHornerxy(poly)

		type(polynomial2), intent(IN) :: poly
		EvalHornerxy = 0

	end function EvalHornerxy

    double precision function EvalHornerxyz(poly)

		type(polynomial3), intent(IN) :: poly

		integer :: di !effective number of variables
		real(kind=prec), dimension(0:3) :: r !registers for summation
		integer, dimension(3) :: i !counters for linear functions
		integer, dimension(3) :: tempjarr
		integer :: k, n, var, M, j, l
		real(kind=prec) :: tempsum

		di = 3
		r(:) = 0
		r(0) = poly%f(poly%n)
		!do k = 1, di
		i(:) = 0
		!end do

		M = poly%n - 1


		do n = 2, M
			!determine k = max(1 <= j <= di : a(n)j != a(n-1)j)
			tempjarr(:) = 0
			do j = 1,3
				if (poly%powers(j, n) .ne. poly%powers(j, n-1)) then
					tempjarr(j) = j
				end if
			end do
			!write(*,*) tempjarr
			k = maxval(tempjarr)
			write(*,*) 'k', k

			!set ik = ik-1
			i(k) = i(k) - 1
			write(*,*) 'i(k)', i(k)

			!set rk = lk,ik(x)(r0 + r1 + ... + rk)
			tempsum = 0
			do l = 0,k
				tempsum = tempsum + r(l)
			end do
			r(k) = tempsum!*lk,ik(x)????
			write(*,*) 'r(k)', r(k)

			!set r0=ca(n), r1...rk-1=0
			r(0) = poly%f(n)
			r(1:k-1) = 0

			!set i1=a(n)1,..., ik-1=a(n)k-1
			do l=1, k-1
				i(l) = poly%powers(l,n)
				write(*,*) 'i(l)', i(l)
			end do

			write(*,*) 'i', i

		end do

		write(*,*) 'r', r

		EvalHornerxyz = sum(r)

	end function EvalHornerxyz

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
