module field
  
  !     .. variables
    use precision_definition

    implicit none
   
    real(dp), dimension(:), allocatable :: xsol    ! volume fraction solvent
    real(dp), dimension(:), allocatable :: xpa     ! volume fraction pa-fiber
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
    real(dp), dimension(:), allocatable :: rhoqpa  ! pa charge density in units of vsol
    
    real(dp), dimension(:,:), allocatable :: xpp   ! volume fraction pp ligand

    !real(dp) :: qAB             ! normalization partion fnc polymer 
    !real(dp) :: qC              ! normalization partion fnc polymer 

    
  
contains

    subroutine allocate_field(N)
        implicit none

        integer, intent(in) :: N

       
        allocate(xsol(N))
        allocate(xpa(N))
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
        allocate(rhoqpa(N))
        
        allocate(xpp(N,6))
        
    end subroutine allocate_field


    subroutine deallocate_field()
        implicit none
        
        
        deallocate(xsol)
        deallocate(xpa)
        
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
        deallocate(rhoqpa)
        
        deallocate(xpp)
        
    end subroutine deallocate_field



    ! set volume fraction of pa 
    ! split pa fiber in tree core-shell region 
    ! r< R=radius is 

    subroutine init_xpa_elect_volume_dist
        

        use volume, only : nr

        integer ::i 

        ! init volumer fraction pa fiber 
        do i=1,nr
            xpa(i)=0.0_dp
        enddo    
        
        deltai=(radiuspahgr-radiuspacore)
        radiuspahgrend = 4.5_dp 
        rhohgr         = 5.0_dp 
        xpalinker      = 0.8_dp

    end subroutine
    
    


    subroutine init_rhoqpa_charge_dist
        
        use volume, only : nr

        integer ::i 

        ! .. init 
        do i=1,nr
            rhoqpa(i)=0.0_dp
        enddo    

            


    end subroutine

  
end module field

