MODULE io_module

	! All variables read from file that need dynamical allocation
	!
	USE kinds
	
	IMPLICIT NONE
	
	CONTAINS
	
	!-----------------------------------------------------------------------
	SUBROUTINE readfc ( flfrc, frc, tau, zeu, m_loc, ityp, nr1, nr2, nr3, &
					    epsil, nat, ibrav, alat, at, ntyp, amass, omega, has_zstar )
	!-----------------------------------------------------------------------
	  
	  USE mp_module
	  USE constants,  ONLY : amu_ry
	  USE essentials
	  !
	  IMPLICIT NONE
	  ! I/O variable
	  INCLUDE 'mpif.h'
	  CHARACTER(LEN=256) :: flfrc
	  INTEGER :: ibrav, nr1,nr2,nr3,nat, ntyp
	  REAL(KIND = RP) :: alat, at(3,3), epsil(3,3)
	  LOGICAL :: has_zstar
	  ! local variables
	  INTEGER :: i, j, na, nb, m1,m2,m3
	  INTEGER :: ibid, jbid, nabid, nbbid, m1bid,m2bid,m3bid
	  REAL(KIND = RP), ALLOCATABLE :: amass(:)
	  REAL(KIND = RP) :: amass_from_file, omega
	  INTEGER :: nt
	  REAL(KIND = RP) :: celldm(6)
	  CHARACTER(LEN=3) :: atm
	  REAL(KIND = RP), ALLOCATABLE :: frc(:,:,:,:,:,:,:), tau(:,:),  zeu(:,:,:), m_loc(:,:)
	  ! frc : interatomic force constants in real space
	  ! tau_blk : atomic positions for the original cell
	  ! zeu : effective charges for the original cell
	  ! m_loc: the magnetic moments of each atom
	  INTEGER, ALLOCATABLE  :: ityp(:)
	  ! ityp_blk : atomic types for each atom of the original cell
	  !
	  
	  ! Default values
	  celldm(:) = 0.D0
	  !
	  !
	  IF (io_node) OPEN (unit=1,file=flfrc,status='old',form='formatted')
	  !
	  !  read cell data
	  !
	  IF (io_node)THEN
		 READ(1,*) ntyp,nat,ibrav,(celldm(i),i=1,6)
		 if (ibrav==0) then
			read(1,*) ((at(i,j),i=1,3),j=1,3)
		 end if
	  ENDIF
	  CALL MPI_BCAST(ntyp, 1, MPI_INT, root_process, comm, ierr)
	  CALL MPI_BCAST(nat, 1, MPI_INT, root_process, comm, ierr)
	  CALL MPI_BCAST(ibrav, 1, MPI_INT, root_process, comm, ierr)
	  CALL MPI_BCAST(celldm, 6, mp_real, root_process, comm, ierr)

	  IF (ibrav==0) THEN
		 CALL MPI_BCAST(at, 9, mp_real, root_process, comm, ierr)
	  ENDIF
	  
	  ALLOCATE(amass(ntyp))
	  
	  amass(:)=0.0_RP
	  !
	  CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
	  alat = celldm(1)
	  at = at / alat !  bring at in units of alat
	  CALL cell_volume(at, alat, omega)
	  !
	  !  read atomic types, positions and masses
	  !
	  DO nt = 1,ntyp
		 IF (io_node) READ(1,*) i,atm,amass_from_file

		 CALL MPI_BCAST(i, 1, MPI_INT, root_process, comm, ierr)
		 CALL MPI_BCAST(atm, ntyp, MPI_CHARACTER, root_process, comm, ierr)
		 CALL MPI_BCAST(amass_from_file, 1, mp_real, root_process, comm, ierr)
		 IF (i.NE.nt) THEN
			write(*,*) ' ERROR'
			write(*,*) ' Error reading mass from file '
			write(*,*) ' ABORTING....'
			STOP
		 ENDIF
		 IF (amass(nt).EQ.0.d0) THEN
			amass(nt) = amass_from_file/amu_ry
		 ELSE
			WRITE(*,*) 'for atomic type',nt,' mass from file not used'
		 END IF
	  END DO
	  !
	  ALLOCATE (tau(3,nat), ityp(nat), zeu(3,3,nat))
	  !
	  DO na=1,nat
		 IF (io_node) READ(1,*) i,ityp(na),(tau(j,na),j=1,3)
		 CALL MPI_BCAST(i, 1, MPI_INT, root_process, comm, ierr)
		 IF (i.NE.na) THEN
			write(*,*) ' ERROR'
			write(*,*) ' Error reading ityp from file'
			write(*,*) ' ABORTING....'
			STOP
		 ENDIF
	  END DO
	  CALL MPI_BCAST(ityp, nat, MPI_INT, root_process, comm, ierr)
	  CALL MPI_BCAST(tau, 3*nat, mp_real, root_process, comm, ierr)
	  !
	  !  read macroscopic variable
	  !
	  IF (io_node) READ (1,*) has_zstar
	  CALL MPI_BCAST(has_zstar, 1, mp_logical, root_process, comm, ierr)
	  IF (has_zstar) THEN
		 IF (io_node) READ(1,*) ((epsil(i,j),j=1,3),i=1,3)
		 CALL MPI_BCAST(epsil, 9, mp_real, root_process, comm, ierr)
		 IF (io_node) THEN
			DO na=1,nat
			   READ(1,*)
			   READ(1,*) ((zeu(i,j,na),j=1,3),i=1,3)
			END DO
		 ENDIF
		  CALL MPI_BCAST(zeu, 9, mp_real, root_process, comm, ierr)
	  ELSE
		 zeu  (:,:,:) = 0.d0
		 epsil(:,:) = 0.d0
	  END IF
	  !
	  IF (io_node) READ (1,*) nr1,nr2,nr3
	  CALL MPI_BCAST(nr1, 1, MPI_INT, root_process, comm, ierr)
      CALL MPI_BCAST(nr2, 1, MPI_INT, root_process, comm, ierr)
	  CALL MPI_BCAST(nr3, 1, MPI_INT, root_process, comm, ierr)
	  !
	  !  read real-space interatomic force constants
	  !
	  ALLOCATE ( frc(nr1,nr2,nr3,3,3,nat,nat) )
	  frc(:,:,:,:,:,:,:) = 0.d0
	  DO i=1,3
		 DO j=1,3
			DO na=1,nat
			   DO nb=1,nat
				  IF (io_node) READ (1,*) ibid, jbid, nabid, nbbid
				  CALL MPI_BCAST(ibid, 1, MPI_INT, root_process, comm, ierr)
				  CALL MPI_BCAST(jbid, 1, MPI_INT, root_process, comm, ierr)
				  CALL MPI_BCAST(nabid, 1, MPI_INT, root_process, comm, ierr)
				  CALL MPI_BCAST(nbbid, 1, MPI_INT, root_process, comm, ierr)
				  IF(i .NE.ibid  .OR. j .NE.jbid .OR.                   &
					 na.NE.nabid .OR. nb.NE.nbbid) THEN
					 write(*,*) ' ERROR'
					 write(*,*) ' Error in reading force constants from frc file'
					 write(*,*) ' ABORTING....'
					 STOP
				  ENDIF
				  IF (io_node) READ (1,*) (((m1bid, m2bid, m3bid,        &
							  frc(m1,m2,m3,i,j,na,nb),                  &
							   m1=1,nr1),m2=1,nr2),m3=1,nr3)
				   
				  CALL MPI_BCAST(frc(:,:,:,i,j,na,nb), (nr1*nr2*nr3), mp_real, root_process, comm, ierr)
			   END DO
			END DO
		 END DO
	  END DO
	  !
	  IF (io_node) CLOSE(unit=1)
	  !
	  RETURN
	END SUBROUTINE readfc

	!-------------------------------------------------------------------------
	subroutine latgen(ibrav,celldm,a1,a2,a3,omega)
	!-----------------------------------------------------------------------
	  !     sets up the crystallographic vectors a1, a2, and a3.
	  !
	  !     ibrav is the structure index:
	  !       1  cubic P (sc)                8  orthorhombic P
	  !       2  cubic F (fcc)               9  1-face (C) centered orthorhombic
	  !       3  cubic I (bcc)              10  all face centered orthorhombic
	  !       4  hexagonal and trigonal P   11  body centered orthorhombic
	  !       5  trigonal R, 3-fold axis c  12  monoclinic P (unique axis: c)
	  !       6  tetragonal P (st)          13  one face (base) centered monoclinic
	  !       7  tetragonal I (bct)         14  triclinic P
	  !     Also accepted:
	  !       0  "free" structure          -12  monoclinic P (unique axis: b)
	  !      -5  trigonal R, threefold axis along (111) 
	  !      -9  alternate description for base centered orthorhombic
	  !     -13  one face (base) centered monoclinic (unique axis: b)
	  !      91  1-face (A) centered orthorombic
	  !
	  !     celldm are parameters which fix the shape of the unit cell
	  !     omega is the unit-cell volume
	  !
	  !     NOTA BENE: all axis sets are right-handed
	  !     Boxes for US PPs do not work properly with left-handed axis
	  !
	  implicit none
	  integer, intent(in) :: ibrav
	  real(RP), intent(inout) :: celldm(6)
	  real(RP), intent(inout) :: a1(3), a2(3), a3(3)
	  real(RP), intent(out) :: omega
	  !
	  real(RP), parameter:: sr2 = 1.414213562373d0, &
							sr3 = 1.732050807569d0
	  integer :: i,j,k,l,iperm,ir
	  real(RP) :: term, cbya, s, term1, term2, singam, sen
	  !
	  !  user-supplied lattice vectors
	  !
	  if (ibrav == 0) then
		 if (SQRT( a1(1)**2 + a1(2)**2 + a1(3)**2 ) == 0 )  THEN
			 write(*,*) ' Error in input data at lattice generation. Check input.'
			 STOP
		 endif
		 if (SQRT( a2(1)**2 + a2(2)**2 + a2(3)**2 ) == 0 )  then
			 write(*,*) ' Error in input data at lattice generation. Check input.'
			 stop
		 endif
		 if (SQRT( a3(1)**2 + a3(2)**2 + a3(3)**2 ) == 0 )  then
			 write(*,*) ' Error in input data at lattice generation. Check input.'
			 stop
		 endif

		 if ( celldm(1) /= 0.D0 ) then
		 !
		 ! ... input at are in units of alat => convert them to a.u.
		 !
			 a1(:) = a1(:) * celldm(1)
			 a2(:) = a2(:) * celldm(1)
			 a3(:) = a3(:) * celldm(1)
		 else
		 !
		 ! ... input at are in atomic units: define celldm(1) from a1
		 !
			 celldm(1) = SQRT( a1(1)**2 + a1(2)**2 + a1(3)**2 )
		 end if
		 !
	  else
		 a1(:) = 0.d0
		 a2(:) = 0.d0
		 a3(:) = 0.d0
	  end if
	  !
	  if (celldm (1) <= 0.d0) then
		 write(*,*) ' Error in input data at lattice generation. Check input.'
		 stop
	  endif
	  !
	  !  index of bravais lattice supplied
	  !
	  if (ibrav == 1) then
		 !
		 !     simple cubic lattice
		 !
		 a1(1)=celldm(1)
		 a2(2)=celldm(1)
		 a3(3)=celldm(1)
		 !
	  else if (ibrav == 2) then
		 !
		 !     fcc lattice
		 !
		 term=celldm(1)/2.d0
		 a1(1)=-term
		 a1(3)=term
		 a2(2)=term
		 a2(3)=term
		 a3(1)=-term
		 a3(2)=term
		 !
	  else if (ibrav == 3) then
		 !
		 !     bcc lattice
		 !
		 term=celldm(1)/2.d0
		 do ir=1,3
			a1(ir)=term
			a2(ir)=term
			a3(ir)=term
		 end do
		 a2(1)=-term
		 a3(1)=-term
		 a3(2)=-term
		 !
	  else if (ibrav == 4) then
		 !
		 !     hexagonal lattice
		 !
		 if (celldm (3) <= 0.d0) then 
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 !
		 cbya=celldm(3)
		 a1(1)=celldm(1)
		 a2(1)=-celldm(1)/2.d0
		 a2(2)=celldm(1)*sr3/2.d0
		 a3(3)=celldm(1)*cbya
		 !
	  else if (ABS(ibrav) == 5) then
		 !
		 !     trigonal lattice
		 !
		 if (celldm (4) <= -0.5_RP .or. celldm (4) >= 1.0_RP) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 !
		 term1=sqrt(1.0_RP + 2.0_RP*celldm(4))
		 term2=sqrt(1.0_RP - celldm(4))
		 !
		 IF ( ibrav == 5) THEN
			!     threefold axis along c (001)
			a2(2)=sr2*celldm(1)*term2/sr3
			a2(3)=celldm(1)*term1/sr3
			a1(1)=celldm(1)*term2/sr2
			a1(2)=-a1(1)/sr3
			a1(3)= a2(3)
			a3(1)=-a1(1)
			a3(2)= a1(2)
			a3(3)= a2(3)
		 ELSE IF ( ibrav == -5) THEN
			!     threefold axis along (111)
			! Notice that in the cubic limit (alpha=90, celldm(4)=0, term1=term2=1)
			! does not yield the x,y,z axis, but an equivalent rotated triplet:
			!    a/3 (-1,2,2), a/3 (2,-1,2), a/3 (2,2,-1)
			! If you prefer the x,y,z axis as cubic limit, you should modify the
			! definitions of a1(1) and a1(2) as follows:'
			!    a1(1) = celldm(1)*(term1+2.0_dp*term2)/3.0_dp
			!    a1(2) = celldm(1)*(term1-term2)/3.0_dp
			! (info by G. Pizzi and A. Cepellotti)
			!
			a1(1) = celldm(1)*(term1-2.0_RP*term2)/3.0_RP
			a1(2) = celldm(1)*(term1+term2)/3.0_RP
			a1(3) = a1(2)
			a2(1) = a1(3)
			a2(2) = a1(1)
			a2(3) = a1(2)
			a3(1) = a1(2)
			a3(2) = a1(3)
			a3(3) = a1(1)
		 END IF
	  else if (ibrav == 6) then
		 !
		 !     tetragonal lattice
		 !
		 if (celldm (3) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 !
		 cbya=celldm(3)
		 a1(1)=celldm(1)
		 a2(2)=celldm(1)
		 a3(3)=celldm(1)*cbya
		 !
	  else if (ibrav == 7) then
		 !
		 !     body centered tetragonal lattice
		 !
		 if (celldm (3) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 !
		 cbya=celldm(3)
		 a2(1)=celldm(1)/2.d0
		 a2(2)=a2(1)
		 a2(3)=cbya*celldm(1)/2.d0
		 a1(1)= a2(1)
		 a1(2)=-a2(1)
		 a1(3)= a2(3)
		 a3(1)=-a2(1)
		 a3(2)=-a2(1)
		 a3(3)= a2(3)
		 !
	  else if (ibrav == 8) then
		 !
		 !     Simple orthorhombic lattice
		 !
		 if (celldm (2) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (celldm (3) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 !
		 a1(1)=celldm(1)
		 a2(2)=celldm(1)*celldm(2)
		 a3(3)=celldm(1)*celldm(3)
		 !
	  else if ( ABS(ibrav) == 9) then
		 !
		 !     One face (base) centered orthorhombic lattice  (C type)
		 !
		 if (celldm (2) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (celldm (3) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 !
		 IF ( ibrav == 9 ) THEN
			!   old PWscf description
			a1(1) = 0.5d0 * celldm(1)
			a1(2) = a1(1) * celldm(2)
			a2(1) = - a1(1)
			a2(2) = a1(2)
		 ELSE
			!   alternate description
			a1(1) = 0.5d0 * celldm(1)
			a1(2) =-a1(1) * celldm(2)
			a2(1) = a1(1)
			a2(2) =-a1(2)
		 END IF
		 a3(3) = celldm(1) * celldm(3)
		 !
	  else if ( ibrav == 91 ) then
		 !
		 !     One face (base) centered orthorhombic lattice  (A type)
		 !
		 if (celldm (2) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (celldm (3) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 !
		 a1(1) = celldm(1)
		 a2(2) = celldm(1) * celldm(2) * 0.5_RP
		 a2(3) = - celldm(1) * celldm(3) * 0.5_RP
		 a3(2) = a2(2)
		 a3(3) = - a2(3)
		 !
	  else if (ibrav == 10) then
		 !
		 !     All face centered orthorhombic lattice
		 !
		 if (celldm (2) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (celldm (3) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 !
		 a2(1) = 0.5d0 * celldm(1)
		 a2(2) = a2(1) * celldm(2)
		 a1(1) = a2(1)
		 a1(3) = a2(1) * celldm(3)
		 a3(2) = a2(1) * celldm(2)
		 a3(3) = a1(3)
		 !
	  else if (ibrav == 11) then
		 !
		 !     Body centered orthorhombic lattice
		 !
		 if (celldm (2) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (celldm (3) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 !
		 a1(1) = 0.5d0 * celldm(1)
		 a1(2) = a1(1) * celldm(2)
		 a1(3) = a1(1) * celldm(3)
		 a2(1) = - a1(1)
		 a2(2) = a1(2)
		 a2(3) = a1(3)
		 a3(1) = - a1(1)
		 a3(2) = - a1(2)
		 a3(3) = a1(3)
		 !
	  else if (ibrav == 12) then
		 !
		 !     Simple monoclinic lattice, unique (i.e. orthogonal to a) axis: c
		 !
		 if (celldm (2) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (celldm (3) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (abs(celldm(4))>=1.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 !
		 sen=sqrt(1.d0-celldm(4)**2)
		 a1(1)=celldm(1)
		 a2(1)=celldm(1)*celldm(2)*celldm(4)
		 a2(2)=celldm(1)*celldm(2)*sen
		 a3(3)=celldm(1)*celldm(3)
		 !
	  else if (ibrav ==-12) then
		 !
		 !     Simple monoclinic lattice, unique axis: b (more common)
		 !
		 if (celldm (2) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (celldm (3) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (abs(celldm(5))>=1.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 !
		 sen=sqrt(1.d0-celldm(5)**2)
		 a1(1)=celldm(1)
		 a2(2)=celldm(1)*celldm(2)
		 a3(1)=celldm(1)*celldm(3)*celldm(5)
		 a3(3)=celldm(1)*celldm(3)*sen
		 !
	  else if (ibrav == 13) then
		 !
		 !     One face centered monoclinic lattice unique axis c
		 !
		 if (celldm (2) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (celldm (3) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (abs(celldm(4))>=1.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 !
		 sen = sqrt( 1.d0 - celldm(4) ** 2 )
		 a1(1) = 0.5d0 * celldm(1) 
		 a1(3) =-a1(1) * celldm(3)
		 a2(1) = celldm(1) * celldm(2) * celldm(4)
		 a2(2) = celldm(1) * celldm(2) * sen
		 a3(1) = a1(1)
		 a3(3) =-a1(3)
	  else if (ibrav == -13) then
		 !
		 !     One face centered monoclinic lattice unique axis b
		 !
		 if (celldm (2) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (celldm (3) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (abs(celldm(5))>=1.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 !
		 sen = sqrt( 1.d0 - celldm(5) ** 2 )
		 a1(1) = 0.5d0 * celldm(1) 
		 a1(2) =-a1(1) * celldm(2)
		 a2(1) = a1(1)
		 a2(2) =-a1(2)
		 a3(1) = celldm(1) * celldm(3) * celldm(5)
		 a3(3) = celldm(1) * celldm(3) * sen
		 !
	  else if (ibrav == 14) then
		 !
		 !     Triclinic lattice
		 !
		 if (celldm (2) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (celldm (3) <= 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (abs(celldm(4))>=1.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (abs(celldm(5))>=1.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 if (abs(celldm(6))>=1.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 !
		 singam=sqrt(1.d0-celldm(6)**2)
		 term= (1.d0+2.d0*celldm(4)*celldm(5)*celldm(6)             &
			  -celldm(4)**2-celldm(5)**2-celldm(6)**2)
		 if (term < 0.d0) then
			write(*,*) ' Error in input data at lattice generation. Check input.'
			stop
		 endif
		 term= sqrt(term/(1.d0-celldm(6)**2))
		 a1(1)=celldm(1)
		 a2(1)=celldm(1)*celldm(2)*celldm(6)
		 a2(2)=celldm(1)*celldm(2)*singam
		 a3(1)=celldm(1)*celldm(3)*celldm(5)
		 a3(2)=celldm(1)*celldm(3)*(celldm(4)-celldm(5)*celldm(6))/singam
		 a3(3)=celldm(1)*celldm(3)*term
		 !
	  else
		 !
		 write(*,*) ' Error lattice type non-existent'
		 stop
		 !
	  end if
	  !
	  !  calculate unit-cell volume omega
	  !
	  omega=0.d0
	  s=1.d0
	  i=1
	  j=2
	  k=3
	  !
	101 do iperm=1,3
		 omega=omega+s*a1(i)*a2(j)*a3(k)
		 l=i
		 i=j
		 j=k
		 k=l
	  end do
	!
	  i=2
	  j=1
	  k=3
	  s=-s
	  if(s < 0.d0) go to 101
	  omega=abs(omega)
	  return
	!
	end subroutine latgen
	
	!-----------------------------------------------------------------------
	subroutine write_dyn_on_file (xq, phi, nat)
	!-----------------------------------------------------------------------
	  USE kinds, only : RP
	  implicit none
	  ! input variables
	  integer :: nat
	  ! unit number
	  ! number of atom in the unit cell
	  complex(RP) :: phi (3, 3, nat, nat)
	  !  the dynamical matrix
	  real(RP) :: xq (3)
	  ! the q vector
	  ! local variables
	
	  integer :: na, nb, icar, jcar
	  ! counters on atoms
	  ! cartesian coordinate counters
	  write (*, 9000) (xq (icar), icar = 1, 3)
	  do na = 1, nat
	     do nb = 1, nat
	        write (*, '(2i5)') na, nb
	        do icar = 1, 3
	!           write (iudyn, '(3e24.12)') (phi(icar,jcar,na,nb), jcar=1,3)
	           write (*, '(3(2f12.8,2x))') (phi(icar,jcar,na,nb), jcar=1,3)
	        enddo
	     enddo
	  enddo
	
	  return
	9000 format(/,5x,'Dynamical  Matrix in cartesian axes', &
	       &       //,5x,'q = ( ',3f14.9,' ) ',/)
	end subroutine write_dyn_on_file

END MODULE io_module
