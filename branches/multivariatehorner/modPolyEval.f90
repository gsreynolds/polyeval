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
		integer, dimension(2) :: N

		write(*,*) 'Hornerxy'

		EvalHornerxy = 0

	end function EvalHornerxy

	double precision function EvalHornerxyz(poly)

		type(polynomial3), intent(IN) :: poly
		integer :: var, i, j, k, numSteps
		integer, dimension(3) :: numFactors
		real(kind=prec), dimension(3, poly%n-1) ::powers, factorisation
		real(kind=prec), dimension(3) :: monomial, variables
		real(kind=prec), allocatable, dimension(:) :: factors

		write(*,*) 'Hornerxyz'

		powers = poly%powers

		write(*,*) 'xpowers', powers(1, :)
		write(*,*) 'ypowers', powers(2, :)
		write(*,*) 'zpowers', powers(3, :)
		write(*,*)
		write(*,*) 'x', poly%x
		write(*,*) 'y', poly%y
		write(*,*) 'z', poly%z
		write(*,*)

		variables(1) = poly%x
		variables(2) = poly%y
		variables(3) = poly%z

		do var = 1,3

			write(*,*) 'var', var, 'value', variables(var)
			i = 0
			numFactors(var) = 0
			do while (sum(powers(var,:)) >= 1)

				i = i + 1

				!write(*,*)
				write(*,*) 'i', i

				factorisation(var,:) = powers(var,:)

				where (powers(var,:) > 0) powers(var,:) = powers(var,:) - 1

				write(*,*) 'var', var, powers(var, :)
			end do

			write(*,*)
			if (i .eq. 0) then
				numFactors(var) = 0
			else
				numFactors(var) = i-1
			end if
		end do

		write(*,*) 'numFactors', numFactors

		monomial(:) = 0
		!by the monomials
		do i = 1, poly%n-1
			write(*,*) 'i', i
			write(*,*) 'poly%f(i)', poly%f(i)
			monomial(i) = poly%f(i)

			write(*,*)
			write(*,*) 'factorisation(1, i)', factorisation(1, i)
			write(*,*) 'factorisation(2, i)', factorisation(2, i)
			write(*,*) 'factorisation(3, i)', factorisation(3, i)

			do j = 1, 3
				write(*,*) 'j', j
				write(*,*) variables(j), factorisation(j, i), variables(j)**factorisation(j, i)
				monomial(i) = monomial(i) * variables(j)**factorisation(j, i)
				write(*,*) monomial(i)
			end do
			write(*,*) 'monomial', i, monomial(i)
			write(*,*)

		end do

		numSteps = maxval(numFactors)
		allocate(factors(numSteps))
		factors(:) = 1

		write(*,*)
		write(*,*) numFactors

		do i = 1, numSteps
			!write(*,*) 'i', i
			do var = 1, 3

					!write(*,*) 'var', var

					if(numFactors(var) .ne. 0) then
						factors(i) = factors(i)*variables(var)
					end if
					if (numFactors(var) > 0) then
						numFactors(var) = numFactors(var) - 1
						!write(*,*) 'var', var, 'numFactors(var)', numFactors(var)
					end if

			end do
		end do

		write(*,*) 'factors', factors
		EvalHornerxyz = monomial(1)+monomial(2)

		do i = 1, numSteps
			if(2+i > poly%n-1) then
				EvalHornerxyz = factors(i)*EvalHornerxyz
			else
				EvalHornerxyz = factors(i)*EvalHornerxyz+monomial(3)
			end if
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
