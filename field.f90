module field
  
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
    real(dp), dimension(:,:), allocatable :: fdisA ! fraction pa  EE charged
    real(dp), dimension(:), allocatable :: rhoEpa  ! pa of E AA number density in units of vsol 
    real(dp), dimension(:,:), allocatable :: xpp   ! volume fraction pp ligand
    real(dp), dimension(:), allocatable :: epsfcn  ! dielectric constant 
    real(dp), dimension(:), allocatable :: Depsfcn ! derivative dielectric constant


contains

    subroutine allocate_field(N)
        
        integer, intent(in) :: N

        integer :: ier
       
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
        allocate(fdisA(N,6))
        allocate(xpp(N,6))
        allocate(epsfcn(N),stat=ier)    ! relative dielectric constant
        allocate(Depsfcn(N),stat=ier)   ! derivate relative dielectric constant

        if( ier/=0 ) then
            print*, 'Allocation error : stat =', ier
            stop
        endif

    end subroutine allocate_field


    subroutine deallocate_field()        
        
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
        deallocate(fdisA)
        deallocate(xpp)
        deallocate(epsfcn)
        deallocate(Depsfcn)
        
    end subroutine deallocate_field



    ! set volume fraction of pa 
    ! split pa fiber in tree core-shell region 
    ! r< R=radius is 

    subroutine init_xpa_volume_dist
        
        use volume, only : nr

        integer ::i 

        ! init volume fraction pa fiber 
        do i=1,nr
            xpa(i)=0.0_dp
        enddo    
       
        call read_xpa_dist

    end subroutine
    
    

    ! computes total line density of Glu residues

    subroutine init_rhoEpa_dist
        
        use volume, only : nr, deltaG, Asurf, delta
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

        totalEpa=totalEpa*Asurf*delta  ! Asurf = pi*r*L 

    end subroutine


    ! read file rhoEpa.dat assumed to start at postion
    ! first element at r=delta/2 second element r=delta3/2 etc

    subroutine read_rhoEpa_dist(info)
        
        use volume, only : nr, radius, isPACore, delta
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
                    rhoEpa(i)=1.0_dp*rhoEpa_value
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
                rhoEpa(i)=1.0_dp*rhoEpa_value
                i=i+1
                if(line==nr) ios=1 ! do not read beyond line nradius
            enddo

        endif    

        close(un_input)

    end subroutine



    ! read file rhoqpa.dat assumed to start at postion 
    ! first element at r=delta/2 second element r=delta3/2 etc

    subroutine read_xpa_dist(info)
        
        use volume, only : nr, radius, isPACore, delta
        use myutils

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
                    xpa(i)=1.0_dp*xpa_value
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
                xpa(i)=1.0_dp*xpa_value
                i=i+1
                if(line==nr) ios=1 ! do not read beyond line nradius
            enddo

        endif       
           
        close(un_input)    

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

        if(sysflag=="pafiber".or.sysflag=="pafiberIm".or.sysflag=="pafibervarelec".or.&
            sysflag=="pafiberborn".or.sysflag=="pafiberbornscf") then !
            
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

    
    function average_charge_pa_Ca() result(avfdispa)     ! .. post : return average charge of state of polymers

        use volume, only : deltaG, nr
        use globals, only : sysflag
        use parameters, only : isCabinding

        implicit none 

        real(dp) :: avfdispa(6)

        integer :: i, k
        real(dp) :: sumpa

    
        if(sysflag=="pafiber".or.sysflag=="pafiberIm".or.sysflag=="pafibervarelec".or.&
            sysflag=="pafiberborn".or.sysflag=="pafiberbornscf") then !
            
            if(.not.isCabinding) then
                sumpa =0.0_dp
                avfdispa=0.0_dp
                
                do i=1,nr
                    avfdispa(1)=avfdispa(1)+fdispa(i)*rhoEpa(i)*deltaG(i)
                    sumpa =sumpa  + rhoEpa(i)*deltaG(i)
                enddo
                    
                avfdispa(1)=avfdispa(1)/sumpa
            
            else
               
                do k=1,6
                    avfdispa(k)=0.0_dp
                    do i=1,nr
                        avfdispa(k)=avfdispa(k)+fdisA(i,k)*rhoEpa(i)*deltaG(i)
                    enddo
                enddo

                 sumpa =0.0_dp
                do i=1,nr     
                    sumpa =sumpa  + rhoEpa(i)*deltaG(i)
                enddo
                 
                do k=1,6    
                    avfdispa(k)=avfdispa(k)/sumpa
                enddo
            
            endif
                
        else
            avfdispa=0.0_dp
        endif    

    end function average_charge_pa_Ca

    ! computes amount of Rb number density  within maxDebye*LD distance for maximum of rhoE pa distruction
    ! and amount of Rb numbger denisty within disntance of maximum of rhoEpa such that (xRb(i)-xbulk%Rb)/xbulk%Rb > epsdeltaxRb

    subroutine  charge_pa_ratio_freeRb()     ! .. post : return average charge of state of polymers

        use volume, only : deltaG, nr, delta
        use globals, only : sysflag, nsize
        use parameters, only : cRbCl,xbulk, vRb, vsol,DebyeLength,lB
        use parameters, only : maxpalayer, maxdeltaRblayer, numlDs, epsdeltaxRb
        use parameters, only : ratio_free_Rb,ratio_free_Rb_Debye

        
        integer :: i,maxlayer,range
        real(dp) :: sumpa,deltaxRb,lD

        ! maxpalayer =location of maximum  rhoEpa in unit of delta is input 
        ! maxDeybe = number of Deybe length to 
        lD=DebyeLength(lB,cRbCl)
        maxlayer=int(numlDs*lD/delta)
        
        if(sysflag=="pafiber".or.sysflag=="pafiberIm".or.sysflag=="pafibervarelec".or.&
            sysflag=="pafiberborn".or.sysflag=="pafiberbornscf") then !
            
    
            ratio_free_Rb=0.0_dp
            ratio_free_Rb_Debye=0.0_dp

            do i=1,maxpalayer
                ratio_free_Rb_Debye=ratio_free_Rb_Debye+xRb(i)*deltaG(i)
            enddo
            ratio_free_Rb=ratio_free_Rb_Debye

            ! integrate out 4 Debye length from maxpalayer
            range=maxpalayer+1+maxlayer
            if(range>nsize) range=nsize
            do i=maxpalayer+1,range
                ratio_free_Rb_Debye=ratio_free_Rb_Debye+xRb(i)*deltaG(i)
            enddo

            deltaxRb=(xRb(maxpalayer+1)-xbulk%Rb)/xbulk%Rb

            i=maxpalayer        
            do while(abs(deltaxRb)>epsdeltaxRb.and.i<nr-1)
                i=i+1
                deltaxRb=(xRb(i+1)-xbulk%Rb)/xbulk%Rb
                ratio_free_Rb=ratio_free_Rb+xRb(i)*deltaG(i)
            enddo
            maxdeltaRblayer=i

            sumpa =0.0_dp
            do i=1,nr     
                sumpa =sumpa+rhoEpa(i)*deltaG(i)
            enddo
              
            ! normalize
              
            ratio_free_Rb = ratio_free_Rb/(sumpa*vRb*vsol)    
            ratio_free_Rb_Debye= ratio_free_Rb_Debye/(sumpa*vRb*vsol)    
            
                
        else
            maxdeltaRblayer=0
            ratio_free_Rb=0.0_dp
        
        endif    

    end subroutine charge_pa_ratio_freeRb
    
    
end module field

