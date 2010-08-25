module modPolyEval
    implicit none

	integer, parameter :: prec = kind(1.0d0)

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
		real (kind=prec), allocatable, dimension(:,:) :: powers
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

   		do i = 1, poly%m
   			vars2(i) = poly%vars(i)*poly%vars(i)
   		end do

   		do i = 1, poly%n-1
			monomial(i) = poly%f(i)

			do j = 1, poly%m
				numSteps = (poly%powers(j,i))/2
				do k = 1, numSteps
					monomial(i) = monomial(i) * vars2(j)
				end do
				if (mod(poly%powers(j,i), 2.0d0) .ne. 0) then
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
	double precision function EvalEstrin(poly, rank, size, comm)
		use mpi

		type(polynomial), intent(IN) :: poly
		integer, intent(IN) :: rank, size, comm
		real (kind=prec), allocatable, dimension(:) :: func, powers
		real (kind=prec), allocatable, dimension(:,:) :: coeff
		integer :: i,j,numsteps, shift, nearestpoweroftwo, npow2, ll, ul, source, sourcell, sourceul, sizeofchunk

		!MPI Variables
		integer, dimension(MPI_STATUS_SIZE) :: status
		integer :: ierr

		nearestpoweroftwo = 2**ceiling(log(real(poly%n))/log(2.0d0))

		if (size > nearestpoweroftwo) then
			if (rank==0) then
	    		write(*,*) 'Number of processes is greater than the width of the array.'
	    		write(*,*) 'Num procs: ', size, ', numcoeffs: ', poly%n, ' nearestpoweroftwo: ', nearestpoweroftwo
	    	end if
	    	call MPI_Finalize(ierr)
	    	stop
	    end if

		if (mod(poly%n, nearestpoweroftwo) .ne. 0) then
			shift = nearestpoweroftwo-mod(poly%n, nearestpoweroftwo)
			npow2 = poly%n + shift
		else
			shift = 0
			npow2 = poly%n
		end if

		allocate(func(npow2))
		func(:) = 0

		numsteps = log(real(npow2))/log(2.0d0)
		allocate(coeff(npow2, 0:numsteps))
		coeff(:,:) = 0
		allocate(powers(numsteps))

		func(1+shift:npow2) = poly%f(:)

		!Evaluate using Estrin's Method
		powers(1) = poly%x
		do i = 2, numsteps
			powers(i) = powers(i-1)**2
		end do

		coeff(:, 0) = func

!		write(*,*) 'rank ', rank
!		write(*,*) 'of ', size
!		write(*,*) 'npow2', npow2

		do i = 1, numsteps

!			write(*,*) 'i', i

			sizeofchunk = npow2/size

!			write(*,*) 'sizeofchunk', sizeofchunk

			if (i .ne. numsteps) then
				ll = rank*sizeofchunk+1
				ul = ll+sizeofchunk-1
			else
				ll=1
				ul=1
			end if

!			write(*,*) 'll', ll, 'ul', ul

			do j = ll, ul
				coeff(j, i) = coeff(2*j-1, i-1)*powers(i)+coeff(2*j, i-1)
			end do

!			write(*,*)
!			write(*,*) 'rank ', rank, ' data'
!			write(*,*) coeff(:, i)


			if (i .ne. numsteps) then

				if (rank .ne. 0) then

!					write(*,*) 'rank', rank, 'sending to root'
!					write(*,*) 'data: ', coeff(ll:ul, i)
					call MPI_SSend(coeff(ll:ul, i), ul-ll+1, MPI_DOUBLE_PRECISION, 0, 0, comm, ierr)

				else
					do source = 1, size-1
!						write(*,*) 'root receiving from rank ', source
						sourcell = source*sizeofchunk+1
						sourceul = sourcell+sizeofchunk-1
!						write(*,*) 'root inserting data from rank ', source, ' into coeff(', sourcell,':',sourceul, ',', i, ')'
						call MPI_Recv(coeff(sourcell:sourceul, i), ul-ll+1, MPI_DOUBLE_PRECISION, source, 0, comm, status, ierr)
					end do

!					write(*,*)
!					write(*,*) 'updated data on rank 0'
!					write(*,*) coeff(:, i)
				end if

				!Broadcasting updated data back to processes other than root
				call MPI_Bcast(coeff(:, i), npow2, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!				if (rank .ne. 0) then
!					write(*,*) 'rank ', rank, ' recieved updated coeff array slice'
!				end if

			end if

!			write(*,*)

		end do

!		if (rank == 0) then
!
!			do i = 0, numsteps
!
!				write(*,*) coeff(:, i)
!
!			end do
!
!		end if

		EvalEstrin = coeff(1, numsteps)

	end function EvalEstrin

end module modPolyEval
