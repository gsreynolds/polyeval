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
		real(kind=prec), dimension(3,poly%n-1) :: powers
		real(kind=prec), dimension(3) :: factorValue
		real(kind=prec), dimension(100, 3, poly%n-1) :: factorisation, unfactorised
		real(kind=prec) :: xfactpower, yfactpower, zfactpower, xunfactpower, yunfactpower, zunfactpower, unfactorisedSum

		integer, dimension(3) :: N
		integer :: maxNloc, minAloc, i, j, k, l, unfactorisableCount, V, W, numsteps, numUnfactorisedMonomials
		integer, dimension(100) :: factor, M, constantCoeff
		integer, dimension(poly%n-1) :: A, B
		logical, dimension(3) :: selectMask
		logical, dimension(poly%n-1) :: factoredOutMask

		write(*,*) 'Hornerxyz'

		powers = poly%powers
		i = 0
		write(*,*) 'xpowers', powers(1, :)
		write(*,*) 'ypowers', powers(2, :)
		write(*,*) 'zpowers', powers(3, :)

		factorisation(:,:,:) = 0
		unfactorised(:,:,:) = 0
		constantCoeff(:) = 0
		A(:) = 0
		factoredOutMask(:) = .true.
		!do while (sum(powers) > poly%n-1)
		do while (i < 10)
			i = i + 1
			write(*,*) '============'
			write(*,*)
			write(*,*) 'i', i
			write(*,*) sum(powers)
			write(*,*) 'factored out mask', factoredOutMask(:)

			write(*,*)
			write(*,*) 'powers'
			write(*,*) powers(1,:)
			write(*,*) powers(2,:)
			write(*,*) powers(3,:)

			!by the variables
			N(1) = count(powers(1,:) > 0)
			N(2) = count(powers(2,:) > 0)
			N(3) = count(powers(3,:) > 0)
			write(*,*) 'N', N

			maxNloc = maxloc(N,1)
			write(*,*) 'maxNloc', maxNloc, 'i.e. the factor, x=1, y=2, z=3'
			factor(i) = maxNloc

			do j = 1, 3
				!write(*,*) j, maxNloc
				if (j .eq. maxNloc) then
					selectMask(j) = .false.
				else
					selectMask(j) = .true.
				end if
			end do


			!(*,*) powers(maxNloc, :)
			unfactorisableCount = 0
			do j = 1, poly%n-1

				if (factoredOutMask(j) .eqv. .true.) then
					if (powers(maxNloc,j) .eq. 0) then
						unfactorisableCount = unfactorisableCount + 1
					end if
				end if

			end do
!			unfactorisableCount = count(powers(maxNloc, :) .eq. 0)
			!count is the number of monomials that do not contain the current variable being considered, x,y or z.
			write(*,*) 'unfactorisable count', unfactorisableCount

			if (unfactorisableCount > 0) then
				do j=1, unfactorisableCount
						write(*,*) '-----'
					!by the monomials
					A(:) = 0

					write(*,*) powers(maxNloc, :)
					write(*,*) 'factored out mask', factoredOutMask(:)
					minAloc = minloc(powers(maxNloc, :), 1, factoredOutMask)

					write(*,*) 'minAloc', minAloc, 'i.e. monomial', minAloc, 'is not factorisable'

					do k = 1, 3

						write(*,*) 'k', k
						write(*,*) '!!!!'
						write(*,*) powers(k, minAloc)
						unfactorised(i, k, minAloc) = powers(k, minAloc)
						powers(k, minAloc) = 0

						if (sum(powers(:, minAloc)) .eq. 0) then
							!write(*,*) '!!!!!!!!!!!', minAloc
							factoredOutMask(minAloc) = .false.
							constantCoeff(i) = minAloc
							!write(*,*) 'i', i
						else
							factoredOutMask(minAloc) = .true.
						end if

					end do

					write(*,*) 'factored out mask', factoredOutMask(:)

				end do
			end if

			where (powers(maxNloc,:) > 0) powers(maxNloc,:) = powers(maxNloc,:) - 1

			factorisation(i,:,:) = powers(:,:)
			write(*,*) 'factorised'
			write(*,*) 'x', factorisation(i, 1, :)
			write(*,*) 'y', factorisation(i, 2, :)
			write(*,*) 'z', factorisation(i, 3, :)

			write(*,*) 'unfactorised'
			write(*,*) 'x', unfactorised(i, 1, :)
			write(*,*) 'y', unfactorised(i, 2, :)
			write(*,*) 'z', unfactorised(i, 3, :)

		end do

		numsteps = i
		write(*,*)

		write(*,*) 'constantCoeff', constantCoeff(1:numsteps)

		EvalHornerxyz = 0

		factorValue(1) = poly%x
		factorValue(2) = poly%y
		factorValue(3) = poly%z

		!write(*,*) 'numsteps', numsteps
		numUnfactorisedMonomials = count(factoredOutMask .eqv. .true.)
		write(*,*) 'numUnfactorisedMonomials', numUnfactorisedMonomials

		do j = 1, numUnfactorisedMonomials
			!j by the monomials

			!write(*,*) ';;;', j

			xfactpower = factorisation(numsteps, 1, j)
			yfactpower = factorisation(numsteps, 2, j)
			zfactpower = factorisation(numsteps, 3, j)
			A(:) = factorisation(numsteps, :,j)

			xunfactpower = unfactorised(numsteps, 1, j)
			yunfactpower = unfactorised(numsteps, 2, j)
			zunfactpower = unfactorised(numsteps, 3, j)
			B(:) = unfactorised(numsteps, :,j)

			if (sum(A(:)) .eq. 0) then
				V = 0
			else
				V = 1
			end if

			if (sum(B(:)) .eq. 0) then
				!do i = 1, poly%n-1
				W = 0
				!end do
			else
				W = 1
			end if

			!write(*,*) 'V', V, 'W', W
			write(*,*) 'factor', factor(1)
			write(*,*) 'factor value',factorValue(factor(1))
			write(*,*) 'V', V
			write(*,*) poly%f(j),poly%x,xfactpower,poly%y,yfactpower,poly%z,zfactpower
			write(*,*) 'W', W
			write(*,*) poly%f(j),poly%x,xunfactpower,poly%y,yunfactpower,poly%z,zunfactpower

			EvalHornerxyz = EvalHornerxyz + factorValue(factor(numsteps)) *&
							(V*poly%f(j)*poly%x**xfactpower*poly%y**yfactpower*poly%z**zfactpower) +&
 							(W*poly%f(j)*poly%x**xunfactpower*poly%y**yunfactpower*poly%z**zunfactpower)
			write(*,*) 'result', EvalHornerxyz

		end do


		write(*,*)
		write(*,*)

		write(*,*) 'factor', factor(1:numsteps)
		!numFactorisedMonomials = poly%n-1-numUnfactorisedMonomials

		do k=numsteps-1, 1, -1

			write(*,*) 'k', k
			unfactorisedSum = 0

			do j = numUnfactorisedMonomials+1, poly%n-1
				write(*,*) 'j', j
				xunfactpower = unfactorised(k, 1, j)
				yunfactpower = unfactorised(k, 2, j)
				zunfactpower = unfactorised(k, 3, j)

				B(:) = unfactorised(k, :,j)

				write(*,*) 'unfactpow'
				write(*,*) B(:)


				if (sum(B(:)) .eq. 0) then
					if(j .eq. constantCoeff(k)) then
						unfactorisedSum = unfactorisedSum + (poly%f(j))
						write(*,*) 'j .eq. constantCoeff(k)', j, constantCoeff(k)
					end if
				else
					unfactorisedSum = unfactorisedSum+ &
								(poly%f(j) * poly%x ** xunfactpower * poly%y ** yunfactpower * poly%z ** zunfactpower)
				end if


				write(*,*) 'unfactorisedSum', unfactorisedSum

			end do

			write(*,*) 'numsteps-k', numsteps-k+1
			write(*,*) 'factor', factor(numsteps-k+1)
			write(*,*) 'factor value',factorValue(factor(numsteps-k+1))
			EvalHornerxyz = factorValue(factor(numsteps-k+1))*(EvalHornerxyz) + unfactorisedSum
			write(*,*) EvalHornerxyz
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
