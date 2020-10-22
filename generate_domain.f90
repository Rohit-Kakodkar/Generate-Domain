PROGRAM generate_domain
	
	USE mp_module
	USE io_module
	USE kinds
	USE lapack95
	USE blas95
	USE f95_precision
	USE essentials
	use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
											  stdout=>output_unit, &
											  stderr=>error_unit
											  
	IMPLICIT NONE
	
	INCLUDE 'mpif.h' ! MPI header file
	INCLUDE 'mkl.fi' ! MKL header file
	
	CHARACTER(LEN = 256) 			:: flfrc1, flfrc2, mass_prefix, domain_prefix
	REAL							:: PD(3), center(3), zmid, r(3), rn
	LOGICAL 						:: nanoparticle ! if nanoparticle is false interface = .true.
	REAL(KIND = RP) 				:: Radius
	REAL(KIND = RP), ALLOCATABLE 	:: my_amass_PD(:,:,:,:)
	INTEGER, ALLOCATABLE 			:: my_ityp_PD(:,:,:,:)
	INTEGER 						:: i, n1,n2, n3, na
	INTEGER 						:: nr1, nr2, nr3, nat(2), ibrav(2), ntyp(2), &
									   natsc, my_nr3, nr3_start
	REAL(KIND = RP) 				:: epsil(3,3,2), at1(3,3), at2(3,3), at_conv(3,3)
	REAL(KIND = RP) 				:: alat1, alat2, omega1, omega2, q(3), distance
	REAL(KIND = RP), ALLOCATABLE 	:: amass1(:), amass2(:)
	LOGICAL 						:: has_zstar(2), fd, na_ifc
	LOGICAL 						:: crystal_coordinates, mass_input
	REAL(KIND = RP), ALLOCATABLE 	:: frc1(:,:,:,:,:,:,:), tau1(:,:),  zeu1(:,:,:), m_loc(:,:)
	INTEGER, ALLOCATABLE  			:: ityp1(:)
	REAL(KIND = RP), ALLOCATABLE 	:: frc2(:,:,:,:,:,:,:), tau2(:,:),  zeu2(:,:,:)
	INTEGER, ALLOCATABLE  			:: ityp2(:)
	INTEGER, ALLOCATABLE 			:: itypsc(:)
	REAL(KIND = RP), ALLOCATABLE 	:: tausc(:,:)
	LOGICAL							:: THERE
	INTEGER							:: type_size
	INTEGER(KIND = MPI_OFFSET_KIND)	:: offset
	INTEGER							:: status(MPI_STATUS_SIZE), fh
	
	
	
	CALL mp_init() ! Intialize MPI
	
	IF (io_node) THEN
		write (stdout, *) 'Generating Domain for Frequency Domain Perfectly Matched Layer'
		write (stdout, *) '	'
	ENDIF
	
	NAMELIST /system/ PD, crystal_coordinates, domain_prefix, mass_prefix, mass_input, flfrc1, flfrc2, &
			 /domain/ nanoparticle, Radius
			 
	
	! Initalize default values
	nanoparticle = .false.
	Radius = 0.D0
	
	IF (root_node) THEN
		READ(stdin, system)
		READ(stdin, domain)
	ENDIF
	
	IF (io_node) THEN
		IF (.not. crystal_coordinates) THEN
			WRITE (stdout, *) 'Generating domain in convetional co-ordinate system (cartesian)'
		ELSE
			WRITE (stdout, *) 'Generating domain in crystal co-ordinate system'
		ENDIF
	ENDIF
	
	CALL MPI_BCAST(PD, 3, MPI_REAL, root_process, comm, ierr)
	CALL MPI_BCAST(crystal_coordinates, 1, mp_logical, root_process, comm, ierr)
	CALL MPI_BCAST(domain_prefix, 256, MPI_CHAR, root_process, comm, ierr)
	CALL MPI_BCAST(mass_prefix, 256, MPI_CHAR, root_process, comm, ierr)
	CALL MPI_BCAST(mass_input, 1, mp_logical, root_process, comm, ierr)
	CALL MPI_BCAST(flfrc1, 256, MPI_CHAR, root_process, comm, ierr)
	CALL MPI_BCAST(flfrc2, 256, MPI_CHAR, root_process, comm, ierr)
	CALL MPI_BCAST(nanoparticle, 1, mp_logical, root_process, comm, ierr)
	CALL MPI_BCAST(Radius, 1, mp_real, root_process, comm, ierr)
	
	IF (nanoparticle .and. (Radius.eq.0)) THEN
		WRITE(stdout, '(a)') ' Radius == 0'
		CALL MPI_ABORT(comm, ierr)
	ENDIF
	
	CALL readfc( flfrc1, frc1, tau1, zeu1, m_loc, ityp1, nr1, nr2, nr3, epsil(:,:,1),&
				 nat(1), ibrav(1), alat1, at1, ntyp(1), amass1, omega1, has_zstar(1) )
	CALL readfc( flfrc2, frc2, tau2, zeu2, m_loc, ityp2, nr1, nr2, nr3, epsil(:,:,2),&
				 nat(2), ibrav(2), alat2, at2, ntyp(2), amass2, omega2, has_zstar(2) )


	nr1 = PD(1)
	nr2 = PD(2)
	nr3 = PD(3)
	
	IF (nanoparticle) THEN
		center = PD/2 + (/1.0, 1.0, 1.0/)
	ELSE
		zmid = nr3/2
	ENDIF
	
!	Check if I want to work in conventional or crystal coordinates
	
	IF (.not. crystal_coordinates) THEN
		at_conv(:,:) = 0.D0
		DO i = 1, 3
			at_conv(i,i) = 1.D0
		ENDDO
		CALL Supercell(at_conv, at1, tau1, ityp1, nat(1), tausc, itypsc, natsc)
		DEALLOCATE (tau1, ityp1, ityp2)
		ALLOCATE(tau1(3, natsc), ityp1(natsc), ityp2(natsc))
		at1 = at_conv
		tau1 = tausc
		ityp1 = itypsc
		ityp2 = itypsc
	ENDIF
	
!	Split up the work between different processors

	CALL get_nr3(nr3, my_nr3, nr3_start)
	
	ALLOCATE(my_ityp_PD(natsc,nr1,nr2,my_nr3), my_amass_PD(natsc,nr1,nr2,my_nr3))
	
	DO n3 = 1, my_nr3
		DO n2 = 1, nr2
			DO n1 = 1, nr1
				DO na = 1, natsc
					r = n1*at_conv(:,1) + n2*at_conv(:,2) + (n3 + nr3_start)*at_conv(:,3) + tausc(:,na)
					IF (nanoparticle) THEN
						distance = nrm2((center-r(:)))
						IF (distance.le.Radius) THEN 
							my_ityp_PD(na,n1,n2,n3) = 2
							my_amass_PD(na,n1,n2,n3) = amass2(itypsc(na))
						ELSE
							my_ityp_PD(na,n1,n2,n3) = 1
							my_amass_PD(na,n1,n2,n3)  = amass1(itypsc(na))
						ENDIF
					ELSE
						IF (r(3).gt.zmid) THEN
							my_ityp_PD(na,n1,n2,n3) = 2
							my_amass_PD(na,n1,n2,n3) = amass2(itypsc(na))
						ELSE
							my_ityp_PD(na,n1,n2,n3) = 1
							my_amass_PD(na,n1,n2,n3)  = amass1(itypsc(na))
						ENDIF
					ENDIF
				ENDDO
			ENDDO
		ENDDO
	ENDDO
	
!	If mass and domain files already exist then delete them
	
	IF (io_node) THEN
		INQUIRE(FILE = domain_prefix, EXIST = THERE)
		IF (THERE) CALL MPI_FILE_DELETE(domain_prefix, MPI_INFO_NULL, ierr)
	END IF 
	
	IF (io_node) THEN
		INQUIRE(FILE = mass_prefix, EXIST = THERE)
		IF (THERE) CALL MPI_FILE_DELETE(mass_prefix, MPI_INFO_NULL, ierr)
	END IF
	
!	Now start writing the files 

	IF (io_node) WRITE (stdout, fmt = '(a)') 'Storing primary domain in ', domain_prefix
	IF (io_node) THEN
		IF (mass_input) WRITE (stdout, fmt = '(a)') 'Storing mass of every atom in domain in ', mass_prefix
	END IF
	
	CALL MPI_BARRIER(comm, ierr)
	
!	These are MPI-IO calls to write to same file at the same time
	CALL MPI_File_open(comm, domain_prefix, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)
	CALL MPI_Type_size(MPI_INT, type_size, ierr)
	offset = (nr1*nr2*natsc*nr3_start)*type_size
	CALL MPI_File_seek(fh, offset, MPI_SEEK_SET, ierr)
	CALL MPI_File_write_all(fh, my_ityp_PD, size(my_ityp_PD), MPI_INT, status, ierr)
	CAlL MPI_File_close(fh, ierr)
	
	IF (mass_input) THEN
		CALL MPI_File_open(comm, mass_prefix, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)
		CALL MPI_TYPE_size(mp_real, type_size, ierr)
		offset = (nr1*nr2*natsc*nr3_start)*type_size
		CALL MPI_File_seek(fh, offset, MPI_SEEK_SET, ierr)
		CALL MPI_File_write_all(fh, my_ityp_PD, size(my_amass_PD), mp_real, status, ierr)
		CAlL MPI_File_close(fh, ierr)
	END IF
	
	IF (io_node) WRITE(stdout, fmt = '(a)') 'DONE'
	
	CALL mp_finalize()

END PROGRAM

SUBROUTINE get_nr3(nr3, my_nr3, nr3_start)

	USE mp_module
	IMPLICIT NONE
	INCLUDE 'mpif.h' ! MPI header file
	
	INTEGER				:: nr3, my_nr3, everyones_nr3(world_size), nr3_start, rem, i, cumulative
	

	my_nr3 = nr3/world_size
	rem = MOD(nr3,world_size)
		
	IF (my_id.gt.(world_size-rem-1)) THEN
		my_nr3 = my_nr3+1
	ENDIF
	
	
	CALL MPI_ALLGATHER(my_nr3, 1, MPI_INT, everyones_nr3, 1, MPI_INT, comm, ierr)
	
	cumulative = 0
	DO i = 0, world_size-1
		IF (i .eq. my_id) THEN
			nr3_start = cumulative
		END IF
		cumulative = cumulative + everyones_nr3(i+1)
	END DO
	
	IF (cumulative .ne. nr3) THEN
		IF (io_node) WRITE(stdout, '(a)') 'Problems with splitting up the domain'
		CALL MPI_ABORT(comm, ierr)
	END IF
	
END SUBROUTINE
