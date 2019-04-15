module field
  
  !     .. variables
    use precision_definition

    implicit none
    
    real(dp), dimension(:), allocatable :: xpolAB  ! total volume fraction of polymer on sphere
    real(dp), dimension(:), allocatable :: xpolC   ! total volume fraction of polymer on sphere: hydrocarbon chain
    real(dp), dimension(:), allocatable :: rhopolA ! density A monomer of polymer on sphere
    real(dp), dimension(:), allocatable :: rhopolB ! density B monomer of polymer on sphere
    real(dp), dimension(:), allocatable :: rhopolC ! density C monomer of polymer on sphere
    real(dp), dimension(:), allocatable :: xsol    ! volume fraction solvent
    real(dp), dimension(:), allocatable :: psi     ! electrostatic potential 
    real(dp), dimension(:), allocatable :: xNa     ! volume fraction of positive Na+ ion
    real(dp), dimension(:), allocatable :: xK      ! volume fraction of positive K+ ion
    real(dp), dimension(:), allocatable :: xTB     ! volume fraction of psitive TB (tetra butyl ammonium) ion
    real(dp), dimension(:), allocatable :: xTM     ! volume fraction of psitive TB (tetra methyl ammonium) ion
    real(dp), dimension(:), allocatable :: xCa     ! volume fraction of positive Ca2+ ion
    real(dp), dimension(:), allocatable :: xNaCl   ! volume fraction of NaCl ion pair
    real(dp), dimension(:), allocatable :: xKCl    ! volume fraction of KCl  ion pair
    real(dp), dimension(:), allocatable :: xCl     ! volume fraction of Cl- ion
    real(dp), dimension(:), allocatable :: xNO3    ! volume fraction of NO3- ion, nitrate
    real(dp), dimension(:), allocatable :: xHplus  ! volume fraction of Hplus
    real(dp), dimension(:), allocatable :: xOHmin  ! volume fraction of OHmin 
    real(dp), dimension(:), allocatable :: rhoq    ! total charge density in units of vsol
    real(dp), dimension(:), allocatable :: qpol    ! charge density of polymer
    real(dp), dimension(:,:), allocatable :: fdisA ! degree of dissociation 
    real(dp), dimension(:,:), allocatable :: fdisB ! degree of dissociation
    real(dp), dimension(:,:), allocatable :: xpp   ! volume fraction pp ligand

    real(dp) :: qAB             ! normalization partion fnc polymer 
    real(dp) :: qC              ! normalization partion fnc polymer 

    real(dp), dimension(:), allocatable :: rhopolAL ! density A monomer of polymer on sphere
    real(dp), dimension(:), allocatable :: rhopolBL ! density B monomer of polymer on sphere
    real(dp), dimension(:), allocatable :: rhopolAR ! density A monomer of polymer on sphere
    real(dp), dimension(:), allocatable :: rhopolBR ! density B monomer of polymer on sphere

    real(dp) :: qABL,qABR

  
contains

    subroutine allocate_field(N)
        implicit none

        integer, intent(in) :: N

        allocate(xpolAB(N))
        allocate(xpolC(N))
        allocate(rhopolA(N))
        allocate(rhopolB(N))
        allocate(rhopolC(N))
        allocate(xsol(N))
        allocate(psi(N+1))
        allocate(xNa(N))
        allocate(xK(N))
        allocate(xTB(N))
        allocate(xTM(N))
        allocate(xCa(N))
        allocate(xNaCl(N)) 
        allocate(xKCl(N)) 
        allocate(xCl(N))
        allocate(xNO3(N))  
        allocate(xHplus(N))
        allocate(xOHmin(N))
        allocate(rhoq(N))
        allocate(qpol(N))
        allocate(fdisA(5,N))
        allocate(fdisB(5,N))
        allocate(rhopolAL(N))
        allocate(rhopolAR(N))
        allocate(rhopolBL(N))
        allocate(rhopolBR(N))
        allocate(xpp(N,6))
        
    end subroutine allocate_field


    subroutine deallocate_field()
        implicit none
        
        
        deallocate(xpolAB)
        deallocate(xpolC)
        deallocate(rhopolA)
        deallocate(rhopolB)
        deallocate(rhopolC)
        deallocate(xsol)
        deallocate(psi)
        deallocate(xNa)
        deallocate(xK)
        deallocate(xTB)
        deallocate(xTM)
        deallocate(xCa)
        deallocate(xNaCl) 
        deallocate(xKCl) 
        deallocate(xCl)
        deallocate(xNO3)  
        deallocate(xHplus)
        deallocate(xOHmin)
        deallocate(rhoq)
        deallocate(qpol)
        deallocate(fdisA)
        deallocate(fdisB)
        deallocate(rhopolAL)
        deallocate(rhopolAR)
        deallocate(rhopolBL)
        deallocate(rhopolBR)
        deallocate(xpp)
        
    end subroutine deallocate_field



    logical function myIsNaN(x)
        implicit none
        real(dp) :: x
        if (x /= x) then
            myIsNaN=.true.
        else
            myIsNaN=.false.
        endif 

    end function myIsNaN

  
end module field

