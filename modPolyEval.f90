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
	double precision function EvalEstrin(poly, rank, size, comm)
		use mpi

		type(polynomial), intent(IN) :: poly
		integer, intent(IN) :: rank, size, comm
		real (kind=prec), allocatable, dimension(:) :: powers
		real (kind=prec), allocatable, dimension(:,:) :: coeff
		integer :: i,j,numsteps, shift, nearestpoweroftwo, npow2, ll, ul, source, sourcell, sourceul, sizeofchunk, arraywidthinuse
		integer :: effectivesize

		!MPI Variables
		integer, dimension(MPI_STATUS_SIZE) :: status
		integer :: ierr

		!Round up the number of coefficients to the nearest power of two
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

		!Calculate the number of calculation steps required to reach the result
		numsteps = log(real(npow2))/log(2.0d0)
		allocate(coeff(npow2, 0:numsteps))
		coeff(:,:) = 0
		allocate(powers(numsteps))

		!Build the powers array
		powers(1) = poly%x
		do i = 2, numsteps
			powers(i) = powers(i-1)**2
		end do

		!Populate the first (0) row of the coeff array with the polynomial coefficients
		coeff(1+shift:npow2, 0) = poly%f(:)

		effectivesize = size

		!Evaluate using Estrin's Method
		do i = 1, numsteps
			!Work out the number of elements in the current array row that are "active"/in use.
			!i.e. for a 512 coefficient poly, at step 1 there are 512 coefficients but at step two there are only 256 coefficients
			!therefore, elements > 256 are unused i.e. zero
			arraywidthinuse = npow2/(2**(i-1))

			if (arraywidthinuse < effectivesize) then
				!If the array width in use at this step is greater than the active/effective number of MPI processes at this step
				!then reduce the effective number of MPI processes.
				!i.e. if at this step there is only 8 coefficients (i.e. array width of 8) and there are 16 MPI processes,
				!then reduce the effective number of MPI processes by a factor of two.
				effectivesize = effectivesize/2
			end if

			!Size of chunk is the portion of the array row at this step assigned to each MPI process 
			sizeofchunk = arraywidthinuse/effectivesize

			if (i .ne. numsteps) then
				ll = rank*sizeofchunk+1
				ul = ll+sizeofchunk-1
			else
				!If this is the last step, then lower limit and upper limit are 1, no matter what the rank
				ll=1
				ul=1
			end if

			if (rank <= effectivesize - 1) then
				!If the rank is less than or equal to the active/effective number of MPI processes at this step then do the calc
				do j = ll, ul
					coeff(j, i) = coeff(2*j-1, i-1)*powers(i)+coeff(2*j, i-1)
				end do
			end if

			if (i .ne. numsteps) then

				if (rank .ne. 0) then
					if (rank <= effectivesize - 1) then
						!Don't send data outside of the width of the array in use at this step.
						!Data outside the effective width of the array in use at this step will only consist of zeros.

						call MPI_SSend(coeff(ll:ul, i), ul-ll+1, MPI_DOUBLE_PRECISION, 0, 0, comm, ierr)
					end if
				else
					do source = 1, effectivesize-1
						!Only expect to recieve from processes with data inside the effective width of the array.
						sourcell = source*sizeofchunk+1
						sourceul = sourcell+sizeofchunk-1

						!Root recieves and inserts data from other processes into coeff array
						call MPI_Recv(coeff(sourcell:sourceul, i), ul-ll+1, MPI_DOUBLE_PRECISION, source, 0, comm, status, ierr)

					end do
				end if

				!Broadcasting updated data back to processes other than root
				call MPI_Bcast(coeff(:, i), npow2, MPI_DOUBLE_PRECISION, 0, comm, ierr)

			end if
		end do

		EvalEstrin = coeff(1, numsteps)

	end function EvalEstrin

end module modPolyEval
