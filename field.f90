module field
  
  !     .. variables
    use precision_definition

    implicit none
   
    real(dp), dimension(:), allocatable :: xsol    ! volume fraction solvent
    real(dp), dimension(:), allocatable :: xpa     ! volume fraction pa-fiber
    real(dp), dimension(:), allocatable :: psi     ! electrostatic potential 
    real(dp), dimension(:), allocatable :: xNa     ! volume fraction of positive Na+ ion
    real(dp), dimension(:), allocatable :: xRb     ! volume fraction of positive Rb+ ion
    real(dp), dimension(:), allocatable :: xIm     ! volume fraction of positive Imadozolium +ion

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
    real(dp), dimension(:), allocatable :: fdispa  ! fraction pa  EE charged 
    real(dp), dimension(:), allocatable :: rhoEpa  ! pa of E AA number density in units of vsol
    
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
        allocate(xRb(N))
        allocate(xIm(N))
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
        allocate(fdispa(N))
        allocate(rhoEpa(N))
        
        allocate(xpp(N,6))
        
    end subroutine allocate_field


    subroutine deallocate_field()
        implicit none
        
        
        deallocate(xsol)
        deallocate(xpa)
        
        deallocate(psi)
        deallocate(xNa)
        deallocate(xRb)
        deallocate(xIm)
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
        deallocate(fdispa)
        deallocate(rhoEpa)
        
        deallocate(xpp)
        
    end subroutine deallocate_field



    ! set volume fraction of pa 
    ! split pa fiber in tree core-shell region 
    ! r< R=radius is 

    subroutine init_xpa_elect_volume_dist
        

        use volume, only : nr

        integer ::i 

        ! init volume fraction pa fiber 
        do i=1,nr
            xpa(i)=0.0_dp
        enddo    
       
        call read_xpa_dist


    end subroutine
    
    


    subroutine init_rhoEpa_dist
        
        use volume, only : nr, deltaG, Asurf
        use parameters, only : totalEpa
        
        integer :: i 

        ! .. init 
        do i=1,nr
            rhoEpa(i)=0.0_dp
        enddo  

        call read_rhoEpa_dist

        totalEpa=0.0_dp
        do i=1,nr
            totalEpa=totalEpa+deltaG(i)*rhoEpa(i)
        enddo

        totalEpa=totalEpa*Asurf

    end subroutine


    ! read file rhoqpa.dat assumed to start at postion 
    ! first element at r=delta/2 second element r=delta3/2 etc

    subroutine read_rhoEpa_dist(info)
        
        use volume, only : nr, radius, isPACore, delta
        use parameters, only :zpa
        use myutils

        integer, intent(out), optional :: info

        character(len=9) :: fname
        integer :: ios, un_input  ! un = unit number    
        integer :: i, line , nradius
        real(dp) :: rcoor,rhoEpa_value

        !     .. reading in of variables from file
        write(fname,'(A9)')'rhoEpa.in'
        open(unit=newunit(un_input),file=fname,iostat=ios,status='old')
        if(ios >0 ) then
            print*, 'Error opening rhoqpa file : iostat =', ios
            if (present(info)) info = 1 ! myio_err_inputfile
            return
        endif

        nradius=int(radius/delta)

        if(isPACore) then 
          
            ios=0 
            line = 0
            i=1
            do while (ios == 0)
                read(un_input, * , iostat=ios) rcoor, rhoEpa_value
                line=line+1
                if(rcoor>radius) then
                    rhoEpa(i)=rhoEpa_value
                    i=i+1
                endif    
                if(line==(nr+nradius)) ios=1 ! do not read beyond line nr+nradius
            enddo

        else 
            ios=0 
            line = 0
            i=1
            do while (ios == 0)
                read(un_input, * , iostat=ios) rcoor, rhoEpa_value
                line=line+1
                rhoEpa(i)=rhoEpa_value
                i=i+1
                if(line==nr) ios=1 ! do not read beyond line nradius
            enddo

        endif    


    end subroutine



    ! read file rhoqpa.dat assumed to start at postion 
    ! first element at r=delta/2 second element r=delta3/2 etc

    subroutine read_xpa_dist(info)
        
        use volume, only : nr, radius, isPACore, delta
        use myutils
        use parameters, only : zpa

        integer, intent(out), optional :: info

        character(len=6) :: fname
        integer :: ios, un_input  ! un = unit number    
        integer :: i, line , nradius
        real(dp) :: rcoor,xpa_value

        !     .. reading in of variables from file
        write(fname,'(A6)')'xpa.in'
        open(unit=newunit(un_input),file=fname,iostat=ios,status='old')
        if(ios >0 ) then
            print*, 'Error opening xpa file : iostat =', ios
            if (present(info)) info = 1 ! myio_err_inputfile
            return
        endif

        nradius=int(radius/delta)

        if(isPACore) then 
          
            ios=0 
            line = 0
            i=1
            do while (ios == 0)
                read(un_input, * , iostat=ios) rcoor, xpa_value
                line=line+1
                if(rcoor>radius) then
                    xpa(i)=xpa_value
                    i=i+1
                endif    
                if(line==(nr+nradius)) ios=1 ! do not read beyond line nr+nradius
            enddo

        else 
            ios=0 
            line = 0
            i=1
            do while (ios == 0)
                read(un_input, * , iostat=ios) rcoor, xpa_value
                line=line+1
                xpa(i)=xpa_value
                if(line==nr) ios=1 ! do not read beyond line nradius
            enddo

        endif    
            
    end subroutine


    function total_charge(rhoq,sigmaqSurf) result(sumcharge)

        use mathconst 
        use volume, only : delta,deltaG, Asurf
        use parameters, only : vsol,lb

        implicit none

        real(dp), intent(in) :: rhoq(:)
        real(dp), intent(in) :: sigmaqSurf
        real(dp) :: sumcharge
        integer :: n, i

        n=size(rhoq)

        sumcharge=0.0_dp
        do i=1, n
            sumcharge=sumcharge+rhoq(i)*deltaG(i)
        enddo
        sumcharge = Asurf*( (delta/vsol)*sumcharge +sigmaqSurf/(4.0_dp*pi*lb*delta))

    end function total_charge

    function average_charge_pa() result(avfdispa)     ! .. post : return average charge of state of polymers

        use volume, only : deltaG, nr
        use globals, only : sysflag

        implicit none 

        real(dp) :: avfdispa

        integer :: i
        real(dp) :: sumpa

        if(sysflag/="pa-fiber") then !
            
            avfdispa=0.0_dp
            sumpa =0.0_dp

            do i=1,nr
                avfdispa=avfdispa+fdispa(i)*rhoEpa(i)*deltaG(i)
                sumpa =sumpa  + rhoEpa(i)*deltaG(i)
            enddo
                
            avfdispa=avfdispa/sumpa
                
        else
            avfdispa=0.0_dp
        endif    

    end function average_charge_pa

    
end module field

