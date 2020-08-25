!
! Copyright (C) 2002-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE constants
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : RP
  !
  ! ... The constants needed everywhere
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! ... Mathematical constants
  ! 
  REAL(RP), PARAMETER :: pi     = 3.14159265358979323846_RP 
  REAL(RP), PARAMETER :: tpi    = 2.0_RP * pi
  REAL(RP), PARAMETER :: fpi    = 4.0_RP * pi
  REAL(RP), PARAMETER :: sqrtpi = 1.77245385090551602729_RP 
  REAL(RP), PARAMETER :: sqrtpm1= 1.0_RP / sqrtpi
  REAL(RP), PARAMETER :: sqrt2  = 1.41421356237309504880_RP
  !
  ! ... Physical constants, SI (NIST CODATA 2006), Web Version 5.1
  !     http://physics.nist.gov/constants
  REAL(RP), PARAMETER :: H_PLANCK_SI      = 6.62606896E-34_RP   ! J s
  REAL(RP), PARAMETER :: K_BOLTZMANN_SI   = 1.3806504E-23_RP    ! J K^-1 
  REAL(RP), PARAMETER :: ELECTRON_SI      = 1.602176487E-19_RP  ! C
  REAL(RP), PARAMETER :: ELECTRONVOLT_SI  = 1.602176487E-19_RP  ! J  
  REAL(RP), PARAMETER :: ELECTRONMASS_SI  = 9.10938215E-31_RP   ! Kg
  REAL(RP), PARAMETER :: HARTREE_SI       = 4.35974394E-18_RP   ! J
  REAL(RP), PARAMETER :: RYDBERG_SI       = HARTREE_SI/2.0_RP   ! J
  REAL(RP), PARAMETER :: BOHR_RADIUS_SI   = 0.52917720859E-10_RP ! m
  REAL(RP), PARAMETER :: AMU_SI           = 1.660538782E-27_RP  ! Kg
  REAL(RP), PARAMETER :: C_SI             = 2.99792458E+8_RP    ! m sec^-1
  REAL(RP), PARAMETER :: MUNOUGHT_SI      = fpi*1.0E-7_RP       ! N A^-2
  REAL(RP), PARAMETER :: EPSNOUGHT_SI     = 1.0_RP / (MUNOUGHT_SI * &
                                                       C_SI**2) ! F m^-1
  !
  ! ... Physical constants, atomic units:
  ! ... AU for "Hartree" atomic units (e = m = hbar = 1)
  ! ... RY for "Rydberg" atomic units (e^2=2, m=1/2, hbar=1)
  !
  REAL(RP), PARAMETER :: K_BOLTZMANN_AU   = K_BOLTZMANN_SI / HARTREE_SI
  REAL(RP), PARAMETER :: K_BOLTZMANN_RY   = K_BOLTZMANN_SI / RYDBERG_SI
  !
  ! ... Unit conversion factors: energy and masses
  !
  REAL(RP), PARAMETER :: AUTOEV           = HARTREE_SI / ELECTRONVOLT_SI
  REAL(RP), PARAMETER :: RYTOEV           = AUTOEV / 2.0_RP
  REAL(RP), PARAMETER :: AMU_AU           = AMU_SI / ELECTRONMASS_SI
  REAL(RP), PARAMETER :: AMU_RY           = AMU_AU / 2.0_RP
  !
  ! ... Unit conversion factors: atomic unit of time, in s and ps
  !
  REAL(RP), PARAMETER :: AU_SEC           = H_PLANCK_SI/tpi/HARTREE_SI
  REAL(RP), PARAMETER :: AU_PS            = AU_SEC * 1.0E+12_RP
  !
  ! ... Unit conversion factors: pressure (1 Pa = 1 J/m^3, 1GPa = 10 Kbar )
  !
  REAL(RP), PARAMETER :: AU_GPA           = HARTREE_SI / BOHR_RADIUS_SI ** 3 &
                                            / 1.0E+9_RP 
  REAL(RP), PARAMETER :: RY_KBAR          = 10.0_RP * AU_GPA / 2.0_RP
  !
  ! ... Unit conversion factors: 1 debye = 10^-18 esu*cm 
  ! ...                                  = 3.3356409519*10^-30 C*m 
  ! ...                                  = 0.208194346 e*A
  ! ... ( 1 esu = (0.1/c) Am, c=299792458 m/s)
  !
  REAL(RP), PARAMETER :: DEBYE_SI         = 3.3356409519_RP * 1.0E-30_RP ! C*m 
  REAL(RP), PARAMETER :: AU_DEBYE         = ELECTRON_SI * BOHR_RADIUS_SI / &
                                            DEBYE_SI
  !
  REAL(RP), PARAMETER :: eV_to_kelvin = ELECTRONVOLT_SI / K_BOLTZMANN_SI
  REAL(RP), PARAMETER :: ry_to_kelvin = RYDBERG_SI / K_BOLTZMANN_SI
  !
  ! .. Unit conversion factors: Energy to wavelength
  !
  REAL(RP), PARAMETER :: EVTONM = 1E+9_RP * H_PLANCK_SI * C_SI / &
                                  &ELECTRONVOLT_SI
  REAL(RP), PARAMETER :: RYTONM = 1E+9_RP * H_PLANCK_SI * C_SI / RYDBERG_SI
  !
  !  Speed of light in atomic units
  !
  REAL(RP), PARAMETER :: C_AU             = C_SI / BOHR_RADIUS_SI * AU_SEC
  !
  ! ... zero up to a given accuracy
  !
  REAL(RP), PARAMETER :: eps4  = 1.0E-4_RP
  REAL(RP), PARAMETER :: eps6  = 1.0E-6_RP
  REAL(RP), PARAMETER :: eps8  = 1.0E-8_RP
  REAL(RP), PARAMETER :: eps12 = 1.0E-12_RP
  REAL(RP), PARAMETER :: eps14 = 1.0E-14_RP
  REAL(RP), PARAMETER :: eps16 = 1.0E-16_RP
  REAL(RP), PARAMETER :: eps24 = 1.0E-24_RP
  REAL(RP), PARAMETER :: eps32 = 1.0E-32_RP
  !
  REAL(RP), PARAMETER :: gsmall = 1.0E-12_RP
  !
  REAL(RP), PARAMETER :: e2 = 2.0_RP      ! the square of the electron charge
  REAL(RP), PARAMETER :: degspin = 2.0_RP ! the number of spins per level
  !
  !!!!!! COMPATIBIILITY
  !
  REAL(RP), PARAMETER :: BOHR_RADIUS_CM = BOHR_RADIUS_SI * 100.0_RP
  REAL(RP), PARAMETER :: BOHR_RADIUS_ANGS = BOHR_RADIUS_CM * 1.0E8_RP
  REAL(RP), PARAMETER :: ANGSTROM_AU = 1.0_RP/BOHR_RADIUS_ANGS
  REAL(RP), PARAMETER :: DIP_DEBYE = AU_DEBYE
  REAL(RP), PARAMETER :: AU_TERAHERTZ  = AU_PS
  REAL(RP), PARAMETER :: AU_TO_OHMCMM1 = 46000.0_RP ! (ohm cm)^-1
  REAL(RP), PARAMETER :: RY_TO_THZ = 1.0_RP / AU_TERAHERTZ / FPI
  REAL(RP), PARAMETER :: RY_TO_GHZ = RY_TO_THZ*1000.0_RP
  REAL(RP), PARAMETER :: RY_TO_CMM1 = 1.E+10_RP * RY_TO_THZ / C_SI
  !

END MODULE constants
