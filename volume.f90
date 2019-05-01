!     makes volume elements 
!     for spherical coordinates 
module volume      

    use precision_definition
    implicit none


    !     .. variables

    real(dp) :: delta             ! delta  spacing of lattice site in z-direction
    integer :: nr                 ! nr number of lattice sites in z-direction  nr <= nsize

    real(dp)  :: radius           ! radius spherical/cylinderical/invcylinderical surface
    real(dp)  :: Asurf            ! area of spherical/cylinderical/invcylinderical surface  

    integer :: nrmax              ! nzmax  maximum number of lattice sites in z-direction
    integer :: nrmin              ! nzmin minimumal number of lattice sites in z-direction  
    integer :: nrstep             ! nzstep number of lattice sites stepped over or reduced 

    real(dp), dimension(:), allocatable :: rc ! z-coordinate
    real(dp), dimension(:), allocatable :: G ! geometrical factor 
    real(dp), dimension(:), allocatable :: deltaG ! geometrical factor
    real(dp), dimension(:), allocatable :: Fplus ! factor in Poisson Eq  
    real(dp), dimension(:), allocatable :: Fmin


    character(len=14) :: geometry

    real(dp), private, parameter :: voltol       = 0.0001_dp  ! tolerance of volume
    real(dp), private, parameter :: epsradius   = 1.0e-4_dp  ! tolerance for radius to be zero  

    logical :: isPACore  
   
    ! .. switch
    ! .false. if radial cooridate runs from zero 
    ! .true of radial coordiante runs from R=radius 


contains
  
subroutine allocate_geometry(N)
    implicit none
    integer, intent(in) :: N
    allocate(rc(N))
    allocate(G(N))
    allocate(deltaG(N+1))
    allocate(Fplus(N))
    allocate(Fmin(N))
    
end subroutine allocate_geometry
  
subroutine  make_geometry()

    use globals
    use myutils

    implicit none

    integer ::  i
    real(dp)  ::  vol
    real(dp)  ::  Vtest
    character(len=lenText) :: text

    vol=0.0_dp

    if(radius > epsradius) then  
        isPACore=.true.

        select case (geometry) 
        case ("spherical")
            do i=1,nr
                rc(i)= (i-0.5_dp) * delta + radius      ! radial coordinate 
                G(i) =  (rc(i) /radius)**2              ! geometrical factor 
                deltaG(i) = G(i) + (delta*delta/(12.0_dp*radius*radius)) ! delta G(i)= (1/delta) \int dr G(r) 
                Fplus(i)  = 1.0_dp+  delta/rc(i)
                Fmin(i)   = 2.0_dp - Fplus(i)               ! factors in Poisson Equation
                vol=vol+ deltaG(i)
            enddo
            Asurf=4.0_dp*pi*(radius**2)  
            Vtest=(4.0_dp/3.0_dp)*pi*((nr*delta+radius)**3-radius**3)/(4.0_dp*pi*(radius**2)*delta)
        case("cylindrical")
            do i=1,nr
                rc(i)= (i-0.5_dp) * delta + radius  ! radial coordinate 
                G(i) =  (rc(i) /radius)             ! geometrical factor 
                deltaG(i) = G(i)                    ! delta G(i)= (1/delta) \int dr G(r) 
                Fplus(i)=1.0_dp+ delta/(2.0_dp*rc(i))
                Fmin(i) =2.0_dp-Fplus(i)            ! factors in Poisson Equation
                vol=vol+ deltaG(i)
            enddo
            Asurf=2.0_dp*pi*(radius)  
            Vtest=pi*((nr*delta+radius)**2-radius**2)/(2.0_dp*pi*radius*delta)
        case("invcylindrical")
            do i=1,nr
                rc(i)= (i-0.5_dp) * delta           ! radial coordinate 
                G(i) =  (rc(i) /radius)             ! geometrical factor 
                deltaG(i) = G(i)                    ! delta G(i)= (1/delta) \int dr G(r) 
                Fplus(i)=1.0_dp+ delta/(2.0_dp*rc(i))
                Fmin(i) =2.0_dp-Fplus(i)            ! factors in Poisson Equation
                vol=vol+ deltaG(i)
            enddo
            Asurf=2.0_dp*pi*(radius) 
            Vtest=pi*((nr*delta)**2)/(2.0_dp*pi*radius*delta)
        case("planar")
            do i=1,nr
                rc(i)= (i-0.5_dp) * delta            ! coordinate 
                G(i) = 1.0_dp                        ! geometrical factor 
                deltaG(i) = 1.0_dp                   ! delta G(i)= (1/delta) \int dr G(r) 
                Fplus(i)=1.0_dp
                Fmin(i) =1.0_dp                      ! factors in Poisson Equation
                vol=vol+ deltaG(i) 
            enddo
            Asurf=1.0_dp
            Vtest=nr
        case default
            print*,"Error: geometry not PLANAR, CYLINDRICAL, INVCYLINDERICAL, or SPHERICAL"
            print*,"stopping program"
            stop
        end select
    else 
         isPACore=.false.  

        select case (geometry) 
        case ("spherical")
        
            vol=0.0d0
  
            do i=1,2*nr
                rc(i)= (i-0.5_dp) * delta        ! radial coordinate 
                G(i) =  4.0_dp*pi*rc(i)**2       ! geometrical factor 
                deltaG(i) = G(i) + (4.0_dp*pi*delta*delta/12.0_dp) ! delta G(i)= (1/delta) \int dr G(r) 
         
            Fplus(i)=1.0d0+ delta/rc(i)
            Fmin(i) = 2.0d0-Fplus(i)      ! factors in Poisson Equation
            vol=vol+ deltaG(i)
            enddo
            Asurf=1.0_dp
            vol=vol*delta
            Vtest=(4.0/3.0)*pi*((nr*delta)**3)
    
        case("cylindrical")
             do i=1,nr
                rc(i)= (i-0.5_dp) * delta           ! radial coordinate 
                G(i) =  2.0_dp*pi*rc(i)             ! geometrical factor 
                deltaG(i) = G(i)                    ! delta G(i)= (1/delta) \int dr G(r) 
                Fplus(i)=1.0_dp+ delta/(2.0_dp*rc(i))
                Fmin(i) = 2.0_dp-Fplus(i)            ! factors in Poisson Equation
                vol=vol+ deltaG(i)
            enddo 
            Asurf=1.0_dp
            
            vol=vol*delta
            Vtest=pi*(nr*delta)**2

        case("planar")
            do i=1,nr
                rc(i)= (i-0.5_dp) * delta            ! coordinate 
                G(i) = 1.0_dp                        ! geometrical factor 
                deltaG(i) = 1.0_dp                   ! delta G(i)= (1/delta) \int dr G(r) 
                Fplus(i)=1.0_dp
                Fmin(i) =1.0_dp                      ! factors in Poisson Equation
                vol=vol+ deltaG(i) 
            enddo
            Asurf=1.0_dp
            
            Vtest=nr
        case default
            print*,"Error: geometry not PLANAR, CYLINDRICAL,INVCYLINDERICAL, or SPHERICAL"
            print*,"stopping program"
            stop
        end select

    endif     
  
    write(text,'(A20)')geometry
    text="geometry= "//trim(text)
    call print_to_log(LogUnit,text)

    if(abs(vtest/vol-1.0_dp)>voltol) then 
        print*,"geometry=",geometry  
        print*,"Error: volume incorrect"
        print*,"vol=",vol,"vtest=",vtest
    endif

end subroutine make_geometry


    
end module volume
  
