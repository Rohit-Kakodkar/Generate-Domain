MODULE essentials

USE mp_module

IMPLICIT NONE

CONTAINS

	!----------------------------------------------------------------------------------
	FUNCTION mldivide(A, b, n)
	!----------------------------------------------------------------------------------
		USE lapack95
		USE kinds
		IMPLICIT NONE 
		INCLUDE 'mkl.fi'
		! solve linear equation Ax=b
		INTEGER :: n
		REAL(RP), DIMENSION(n,n), INTENT(IN) :: A
		REAL(RP), DIMENSION(n,n) :: Atemp
		REAL(RP), DIMENSION(n) :: b
		REAL(RP), DIMENSION(n) :: mldivide
		mldivide=b
		Atemp = A
		CALL GESV( Atemp, mldivide)
		RETURN
	END FUNCTION
	
	!------------------------------------------------------------------------------------
	SUBROUTINE Supercell(	atsc, at, tau, ityp, nat, tausc, itypsc, natsc, z, zsc, &
							r_cell, na_vec, ib_vec	)
	!------------------------------------------------------------------------------------
	! generates a supercell defined by atsc
		USE lapack95
		USE kinds
		IMPLICIT NONE
		INCLUDE 'mkl.fi'
		INTEGER :: nat, natsc
		REAL(KIND = RP) :: atsc(3,3), at(3,3), V, Vsc, r(3), tau_temp(3), r_temp(3)
		REAL(KIND = RP), INTENT (IN) :: tau(3,nat)
		INTEGER, INTENT(IN) :: ityp(nat)
		INTEGER, ALLOCATABLE :: itypsc(:)
		REAL(KIND = RP), ALLOCATABLE :: tausc(:,:)
		INTEGER :: counter, counter1, nr1, nr2, nr3, nsc, i1, i2, i3, na, &
					i, nb, j, ia
		COMPLEX(KIND = CP), optional :: z(3*nat)
		COMPLEX(KIND = CP), optional, ALLOCATABLE :: zsc(:)
		REAL(KIND = RP), optional, ALLOCATABLE :: r_cell(:,:)
		INTEGER, optional, ALLOCATABLE :: ib_vec(:,:,:,:,:), na_vec(:)
		LOGICAL :: logic
		INTEGER :: sta, check
		
		CALL cell_volume(at, 1.D0, V)
		
		CALL cell_volume(atsc, 1.D0, Vsc)
		
		IF ((Vsc/V-NINT(Vsc/V)).gt.1.0E-4) THEN
			WRITE(*, *) 'ERROR : Volume of super cell ~= n * Volume of primitive lattice'
			STOP
		ENDIF
		
		nsc = NINT(Vsc/V)
		
		natsc = nat*nsc
		
		ALLOCATE(tausc(3, natsc), itypsc(natsc))
		
		
		nr1 = 8
		nr2 = 8
		nr3 = 8
		
		IF (present(na_vec)) ALLOCATE(na_vec(natsc))
		IF (present(r_cell)) ALLOCATE(r_cell(3, natsc))
		IF (present(zsc)) ALLOCATE(zsc(3*natsc))

		counter= 0
		counter1 = 0
		
		IF (present(ib_vec)) ib_vec(:,:,:,:,:) = -1
					
		DO i1 = -2*nr1,2*nr1
			DO i2 = -2*nr2,2*nr2
				DO i3 = -2*nr3,2*nr3
					DO na = 1, nat
						DO i = 1, 3
							r(i) = i1*at(i,1) + i2*at(i,2) + i3*at(i,3) + tau(i,na)
						ENDDO
						r = mldivide(atsc, r, 3)
						IF ((r(1).ge.0) .and. (r(2).ge.0) .and. (r(3).ge.0) &
						   &.and. (r(1).lt.1) .and. (r(2).lt.1) .and. (r(3).lt.1)) THEN
								counter = counter +1
								tausc(:,counter) = r(1)*atsc(:,1) + r(2)*atsc(:,2) + r(3)*atsc(:,3)
								itypsc(counter) = ityp(na)
								IF (present(na_vec))na_vec(counter) = na
								IF (present(zsc)) THEN
									DO i = 1, 3
										counter1 = counter1+1
										zsc(counter1) = z(3*na-(3-i))
									ENDDO
								ENDIF
								IF (present(r_cell)) THEN
									DO i = 1,3
										r_cell (i, counter) = i1*at(i,1) + i2*at(i,2) + i3*at(i,3)
									ENDDO
								ENDIF
						ENDIF
					ENDDO
				ENDDO
			ENDDO
		ENDDO
		
		IF (counter.ne.natsc) WRITE(*,*) 'WARNING: conventional unit cell might not be correct'
		
		IF (counter.lt.natsc) WRITE(*,*) 'WARNING: Might need to increase the value of nr1 in this code (Supercell in essentials.f90)'
			
		nr1 = 2
		nr2 = 2
		nr3 = 2
		
		IF (present(ib_vec)) ALLOCATE(ib_vec(nat, -2*nr3:2*nr3, -2*nr2:2*nr2, -2*nr1:2*nr1, natsc))

		IF (present(ib_vec)) THEN
			DO ia = 1, natsc
				DO i1 = -2*nr1, 2*nr1
					DO i2 = -2*nr2, 2*nr2
						DO i3 = -2*nr3, 2*nr3
							DO nb = 1, nat
								DO i = 1, 3
									r(i) = i1*at(i,1) + i2*at(i,2) + i3*at(i,3) + tau(i,nb) + r_cell(i,ia)
								ENDDO
!								r = mldivide(atsc, r, 3)
								tau_temp = r - floor(r)
								check = 0
								DO na = 1, natsc
									logic = .true.
									DO i = 1, 3
										IF (abs(tau_temp(i) - tausc(i, na)).le.1E-3) THEN
											logic = (logic .and. (.true.))
										ELSE
											logic = .false.
										ENDIF
									ENDDO
									IF (logic) THEN
										ib_vec(nb,i3,i2,i1,ia) = na
										check = check + 1
									ENDIF
								ENDDO
								IF (check.ne.1) THEN
									WRITE (stdout, '(a)') 'WARNING : Supercell is not calculated correctly'
								ENDIF
							ENDDO	
						ENDDO
					ENDDO
				ENDDO
			ENDDO
		ENDIF
	END SUBROUTINE Supercell
	
	!====================================================================!
	function ind2sub(iG,nSub) result(iSub)
	!! Compute the indices in each dimension from the global index
	!====================================================================!
	USE kinds
	
	IMPLICIT NONE
	
	integer(KIND = IP), intent(in) :: iG !! Index into a global vector
	integer(KIND = IP), intent(in) :: nSub(:) !! Size in each dimension
	integer(KIND = IP) :: iSub(size(nSub)) !! Indices in each dimension to return
	integer(KIND = IP) :: i,iGtmp,iTmp
	integer(KIND = IP) :: nDims
	integer(KIND = IP) :: prod
	
	nDims=size(nSub)
	if (nDims == 1) then
	iSub(1) = iG
	return
	end if
	
	prod = product(nSub)
	iGtmp = iG
	do i = nDims, 1, -1
		prod = prod / nSub(i)
		iTmp = mod(iGtmp-1,prod)+1
		iSub(i) = (iGtmp - iTmp)/prod + 1
		iGtmp = iTmp
	end do
	
	end function
	!====================================================================!
	!====================================================================!
	function sub2ind(iSub,nSub) result(iG)
	!! Given component indices, get the global vector location.
	!====================================================================!
	USE kinds
	
	IMPLICIT NONE
	
	integer(KIND = IP), intent(in) :: iSub(:) !! Indices in each dimension. The first entry in iL is the left most index
	integer(KIND = IP), intent(in) :: nSub(:) !! Size in each dimension
	integer(KIND = IP) :: iG !! Index in the global vector
	integer(KIND = IP) :: i
	integer(KIND = IP) :: nDims
	integer(KIND = IP) :: prod
	nDims=size(iSub)
	
	prod = 1
	iG = 1
	do i = 1,nDims
		iG = iG + (iSub(i)-1)*prod
		prod = prod * nSub(i)
	end do
	
	end function
	!====================================================================!
	
	!--------------------------------------------
	SUBROUTINE cell_volume(at, alat, omega)
	!--------------------------------------------
		USE kinds
		IMPLICIT NONE
		!
		REAL(RP) :: at(3,3), bg_blk(3,3), alat, omega
		INTEGER :: i,j
		
		omega =alat**3 * ABS(at(1,1)*(at(2,2)*at(3,3)-at(3,2)*at(2,3))- &
					   at(1,2)*(at(2,1)*at(3,3)-at(2,3)*at(3,1))+ &
					   at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)))
		RETURN					 
	END SUBROUTINE 
	
	!-------------------------
	SUBROUTINE recips(at, bg)
	!-------------------------
	! Calculate the reciprocal lattice vectors
		USE kinds
		IMPLICIT NONE
		REAL(KIND = RP) :: at(3,3), bg(3,3)
		REAL(KIND = RP) :: a1(3), a2(3), a3(3)
		
		REAL(KIND = RP) :: den, s
		INTEGER :: iperm, i, j, k, l, ipol
		! counter on the permutations
		!\
		!  Auxiliary variables
		!/
		!
		! Counter on the polarizations
		!
		!    first we compute the denominator
		!
		a1 = at(:,1)
		a2 = at(:,2)
		a3 = at(:,3)
		den = 0
		i = 1
		j = 2
		k = 3
		s = 1.d0
		100 do iperm = 1, 3
			den = den + s * a1 (i) * a2 (j) * a3 (k)
			l = i
			i = j
			j = k
			k = l
		enddo
		i = 2
		j = 1
		k = 3
		s = - s
		if (s.lt.0.d0) goto 100
		!
		!    here we compute the reciprocal vectors
		!
		i = 1
		j = 2
		k = 3
		DO ipol = 1, 3
			bg (ipol , 1) = (a2 (j) * a3 (k) - a2 (k) * a3 (j) ) / den
			bg (ipol , 2) = (a3 (j) * a1 (k) - a3 (k) * a1 (j) ) / den
			bg (ipol , 3) = (a1 (j) * a2 (k) - a1 (k) * a2 (j) ) / den
			l = i
			i = j
			j = k
			k = l
		ENDDO
		RETURN
	END SUBROUTINE 
	
	SUBROUTINE heapsort(I, J, A, nnz)
	
!	--------------------------------------------------
!	SUBROUTINE DISCRIPTION
!	--------------------------------------------------
!	Sort A-matrix in COO format based on index I and J
!	For example, Before sorting : 	I = {1, 2, 1, 3, 2} 
!									J = {2, 1, 3, 2, 3} 
!									A = {3, 4, 6, 1, 2}	
!				 After sorting : 	I = {1, 1, 2, 2, 3} 
!									J = {2, 3, 1, 3, 2} 
!									A = {3, 6, 4, 2, 2}
!	Check algorithm for heapsort and pseudocode at ||https://rosettacode.org/wiki/Sorting_algorithms/Heapsort#Fortran|| 

!	--------------------------------------------------
! 	INPUT DATA DISCRIPTION:
!	--------------------------------------------------
! 		I - Array of I indices
!		J - Array of J indices
!		A - Array of values

!	-------------------------------------------------
!	OUTPUT DATA DISCRIPTION
!	-------------------------------------------------
!		I, J, A - Sorted version of I, J, A

		USE kinds
		IMPLICIT NONE

		INTEGER(KIND=IP) :: start, nnz, bottom, n
		COMPLEX(KIND=CP), DIMENSION(nnz) :: A
		INTEGER(KIND=IP), DIMENSION(nnz) :: I, J
		COMPLEX(KIND=CP) :: tempA
		INTEGER(KIND=IP) :: tempI, tempJ
		
		n = nnz
		start = n/2
		DO WHILE (start.ge.1)
			CALL siftdown(I, J, A, start, nnz, nnz)
			start = start - 1
		END DO
		bottom = nnz
		DO WHILE (bottom.gt.1)
			tempI = I(1)
			I(1) = I(bottom)
			I(bottom) = tempI
			tempJ = J(1)
			J(1)=J(bottom)
			J(bottom) = tempJ
			tempA = A(1)
			A(1)=A(bottom)
			A(bottom) = tempA
			bottom = bottom - 1
			CALL siftdown(I, J, A, 1, bottom, nnz)
		END DO
		
	END SUBROUTINE heapsort

	SUBROUTINE siftdown(I, J, A, start, bottom, n)
	
		USE kinds
		IMPLICIT NONE
		INTEGER(KIND=IP) :: child, root, counter, n
		COMPLEX(KIND=CP), DIMENSION(n) :: A
		INTEGER(KIND=IP), DIMENSION(n) :: I, J
		INTEGER(KIND=IP), INTENT(IN) :: start, bottom
		COMPLEX(KIND=CP) :: tempA
		INTEGER(KIND=IP) :: tempI, tempJ
		root = start
		
		counter = 0
		DO WHILE((root*2) .le. bottom)
			child = root * 2
			IF ((child + 1) .le. bottom) THEN
				IF (((I(child)) .lt. I(child+1)) .OR. ((I(child) .lt. I(child+1)) .OR. ((I(child).eq.I(child+1)) .AND. (J(child).lt.J(child+1))))) THEN 
					child = child + 1
				ENDIF
			ENDIF
 
			IF (I(root).lt.I(child)) THEN
				tempI = I(child)
				I(child) = I(root)
				I(root) = tempI
				tempJ = J(child)
				J(child) = J(root)
				J(root) = tempJ
				tempA = A(child)
				A(child) = A(root)
				A(root) = tempA
				root = child
			ELSEIF ((I(root).eq.I(child)) .and. (J(root).lt.J(child))) THEN
				tempJ = J(child)
				J(child) = J(root)
				J(root) = tempJ
				tempA = A(child)
				A(child) = A(root)
				A(root) = tempA
				root = child
			ELSE
				RETURN
			ENDIF 
		END DO   
	END SUBROUTINE siftdown
	
END MODULE essentials
