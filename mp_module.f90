MODULE mp_module

USE kinds
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
											  stdout=>output_unit, &
											  stderr=>error_unit

IMPLICIT NONE
! Global variables-------------
INTEGER(KIND = RP), PARAMETER :: root_process = 0 ! If you wanna mess with root processor. root_process.lt.world_size
INTEGER(KIND = RP), PARAMETER :: io_process = 0 ! If you wanna mess with I/O processor. Set it .ge. world_size if you do not want any printed outputs
INTEGER :: my_id, world_size = 0, comm ! mpi variables
LOGICAL :: root_node, io_node 
INTEGER :: mp_real, mp_complex, mp_int, mp_logical, mp_sumr, mp_sumc, mp_sumi
INTEGER :: ierr
!------------------------------

CONTAINS
	
	SUBROUTINE mp_init ( )
!	======================================================================
!	You shouldn't be required to change anything in this routine
!	======================================================================
!	
!	Initializes mpi framework
!
!	----------------------------------------------------------------------
!	OUTPUT DATA DISCRIPTION
!	----------------------------------------------------------------------
!	world_size - Total number of processors required
!	my_id - my processor tag
!	comm - MPI world communicator handle
!	ierr - error handle 
!	root_node - logical defining if I am root process
!	io_node - logical defining if I am I/O process
!	mp_logical - MPI logical datatype
!	mp_int - MPI integer datatype of KIND = IP
!	mp_real - MPI real datatype of KIND = RP
!	mp_complex - MPI complax datatype of KIND = CP
!	mp_sumr - summation operation for KIND = mp_real
!	mp_sumc - summation operation for KIND = mp_complex
!	mp_sumi - summation operation for KIND = mp_int


		INCLUDE 'mpif.h' ! MPI header file
		
		CALL MPI_INIT(ierr)
		comm = MPI_COMM_WORLD
		CALL MPI_COMM_RANK (comm, my_id, ierr)	! call the rank of every processor
		CALL MPI_COMM_SIZE (comm, world_size, ierr) 		! Total number of processors
		
		IF (root_process.ge.world_size) THEN
			WRITE (stdout,'(a)') 'root_process > world_size'
			STOP
		ENDIF
		
		IF (io_process.ge.world_size) THEN
			WRITE (stdout,'(a)') 'No printing output'
			STOP
		ENDIF
		
		root_node = .false.
		io_node = .false.
		IF (my_id.eq.root_process) THEN
			root_node = .true.
		ENDIF
		
		IF (my_id.eq.io_process) THEN
			io_node = .true.
		ENDIF
		
		! define MPI datatypes here, These are called in MPI routines
		! Make sure these are consistent with the precision-datatypes defined in kinds
		mp_logical = MPI_LOGICAL
		CALL MPI_TYPE_CONTIGUOUS(IP/4, MPI_INTEGER4, mp_int, ierr)
		CALL MPI_TYPE_COMMIT(mp_int, ierr)
		CALL MPI_TYPE_CONTIGUOUS(RP/4, MPI_REAL, mp_real, ierr)
		CALL MPI_TYPE_COMMIT(mp_real, ierr)
		CALL MPI_TYPE_CONTIGUOUS(CP/4, MPI_COMPLEX, mp_complex, ierr)
		CALL MPI_TYPE_COMMIT(mp_complex, ierr)
		
		CALL MPI_OP_CREATE(mp_sum_real, .true., mp_sumr, ierr)
		CALL MPI_OP_CREATE(mp_sum_cmplx, .true., mp_sumc, ierr)
		CALL MPI_OP_CREATE(mp_sum_int, .true., mp_sumi, ierr)
		
		IF (root_node) THEN
			write (stdout, *) '	'
			write (stdout, *) '====================================================='
			write (stdout, *) '	MPI INITIALIZED '
			write (stdout, '(a, I3)') '	World Size =', world_size
			write (stdout, *) '====================================================='
			write (stdout, *) '	'
		ENDIF
		
	END SUBROUTINE
	
	SUBROUTINE check_mpinit( )
		! Check if the mpi environment has been initialized
		! If not then stop
		IF (world_size.eq.0) THEN
			WRITE (stdout,'(a)') 'MPI has not been initialized'
			WRITE (stdout,'(a)') 'Call mp_init( )'
			STOP
		ENDIF
	
	END SUBROUTINE

	SUBROUTINE calculate_displs(counts, displs)
		INTEGER(KIND = RP), DIMENSION(world_size) :: counts, displs
		INTEGER :: i
		displs(1)=0
		DO i=2,world_size
			displs(i) = displs(i-1) + counts(i-1)
		ENDDO
	END SUBROUTINE

	SUBROUTINE mp_finalize ( )
!	Shuts of the MPI process
		INCLUDE 'mpif.h' ! MPI header file
		INTEGER ::  ierr
		
		IF (io_node) write(stdout, '(a)') '	'
		IF (io_node) write(stdout, '(a)') '==================================='
		IF (io_node) write(stdout, '(a)') 'Closing MPI session'
		
		CALL MPI_FINALIZE(ierr)
		
		IF (io_node) write(stdout, '(a)') 'DONE'
		IF (io_node) write(stdout, '(a)') '==================================='


		
	END SUBROUTINE
	
!	**********************************************************************************
!	USER defined MPI functions

	SUBROUTINE mp_sum_real( invec, inoutvec, datatype )
!	mp_real sum operation (To be used in ALLREDUCE or REDUCE MPI operations) 
		INTEGER :: datatype 
		REAL(KIND = RP) :: invec(world_size), inoutvec(world_size) 
		INTEGER :: i 
		
		DO i= 1, world_size 
			inoutvec(i) = invec(i) + inoutvec(i) 
		return 
		ENDDO 
	END SUBROUTINE mp_sum_real
	
	SUBROUTINE mp_sum_cmplx( invec, inoutvec, datatype )
!	mp_complex sum operation (To be used in ALLREDUCE or REDUCE MPI operations) 
		INTEGER :: datatype 
		COMPLEX(KIND = CP) :: invec(world_size), inoutvec(world_size) 
		INTEGER :: i 
		
		DO i= 1, world_size 
			inoutvec(i) = invec(i) + inoutvec(i) 
		return 
		ENDDO 
	END SUBROUTINE mp_sum_cmplx
	
	SUBROUTINE mp_sum_int( invec, inoutvec, datatype )
!	mp_int sum operation (To be used in ALLREDUCE or REDUCE MPI operations) 
		INTEGER :: datatype 
		COMPLEX(KIND = IP) :: invec(world_size), inoutvec(world_size) 
		INTEGER :: i 
		
		DO i= 1, world_size 
			inoutvec(i) = invec(i) + inoutvec(i) 
		return 
		ENDDO 
	END SUBROUTINE mp_sum_int
	
	

END MODULE mp_module
