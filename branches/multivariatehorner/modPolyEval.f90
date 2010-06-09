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
	end type

	type polynomial3
   		integer :: n
   		real (kind=prec), allocatable, dimension(:) :: f !Coefficients of polynomial
		real (kind=prec) :: x !Value of independent variable
		real (kind=prec) :: y !Value of independent variable
		real (kind=prec) :: z !Value of independent variable
	end type

	interface EvalHorner
		module procedure EvalHornerx, EvalHornerxy, EvalHornerxyz
	end interface

    contains

	!Function to evaluate a polynomial using Horner's form
    double precision function EvalHornerx(poly)

		type(polynomial), intent(IN) :: poly
		integer :: i

		write(*,*) 'Hornerx'

		!Evaluate using Horner's Method
		EvalHornerx = poly%f(1)*poly%x
		do i = 2, poly%n-1
			EvalHornerx = (EvalHornerx + poly%f(i))*poly%x
		end do
		EvalHornerx = EvalHornerx + poly%f(poly%n)

	end function EvalHornerx

	double precision function EvalHornerxy(poly)

		type(polynomial2), intent(IN) :: poly
		integer :: i

		write(*,*) 'Hornerxy'

		!Evaluate using Horner's Method
		EvalHornerxy = poly%f(1)*poly%x
		do i = 2, poly%n-1
			EvalHornerxy = (EvalHornerxy + poly%f(i))*poly%x
		end do
		EvalHornerxy = EvalHornerxy + poly%f(poly%n)

	end function EvalHornerxy

	double precision function EvalHornerxyz(poly)

		type(polynomial3), intent(IN) :: poly
		integer :: i

		write(*,*) 'Hornerxyz'

		!Evaluate using Horner's Method
		EvalHornerxyz = poly%f(1)*poly%x
		do i = 2, poly%n-1
			EvalHornerxyz = (EvalHornerxyz + poly%f(i))*poly%x
		end do
		EvalHornerxyz = EvalHornerxyz + poly%f(poly%n)

	end function EvalHornerxyz

	!Function to evaluate a polynomial using Estrin's method
	double precision function EvalEstrin(poly)
		use omp_lib

		type(polynomial), intent(IN) :: poly
		real (kind=prec), allocatable, dimension(:) :: func, powers
		real (kind=prec), allocatable, dimension(:,:) :: coeff
		integer :: i,j,numsteps, shift, nearestpoweroftwo, npow2

		nearestpoweroftwo = 2**ceiling(log(real(poly%n))/log(2.0d0))
		if (mod(poly%n, nearestpoweroftwo) .ne. 0) then
			shift = nearestpoweroftwo-mod(poly%n, nearestpoweroftwo)
			npow2 = poly%n + shift
		else
			shift = 0
			npow2 = poly%n
		end if

		allocate(func(npow2))
		allocate(coeff(npow2, 0:npow2))

		numsteps = log(real(npow2))/log(2.0d0)
		allocate(powers(numsteps))

		func(1+shift:npow2) = poly%f(:)

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
