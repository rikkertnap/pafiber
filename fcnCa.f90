! --------------------------------------------------------------|
! fcnCa.f90:                                                    |
! constructs the vector function  needed by the                 |
! routine solver, which solves the SCMFT eqs for weak poly-     |
! electrolytes onto a tethered planar surface                   |
! --------------------------------------------------------------|


module fcnpointer

    implicit none

    abstract interface
        subroutine fcn(x,f,n)
            use precision_definition
            implicit none
            integer(8), intent(in) :: n
            real(dp), dimension(n), intent(in) :: x
            real(dp), dimension(n), intent(out) :: f
        end subroutine fcn
    end interface

    procedure(fcn), pointer :: fcnptr => null()

end module fcnpointer


module listfcn

    implicit none

    contains

    ! pre:  vector x=xsol+psi+rhopolA+rhoolB+xpolC                   
    ! post: vector f=xpolAB+xpolC+xsol+Sum_i x_i -1,poisson equation,rhoA,rhoB,xpolC              

    subroutine fcnelectHC(x,f,nn)

    !     .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW
        use surface
        use vectornorm
 
        implicit none

        !     .. scalar arguments
        !     .. array arguments
        integer(8), intent(in) :: nn
        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        

        !     .. declare local variables

        real(dp) :: exppiA(nsize),exppiB(nsize),exppiC(nsize)    ! auxilairy variable for computing P(\alpha) 
        real(dp) :: rhopolAin(nsize),rhopolBin(nsize),xpolCin(nsize)
        real(dp) :: xA(3),xB(3),sumxA,sumxB
        real(dp) :: constA,constB
        real(dp) :: pro,rhopolAB0,rhopolC0
        integer :: n                 ! half of n
        integer :: i,j,k,c,s         ! dummy indices
        real(dp) :: tmp,expVdW 
        real(dp) :: norm
        integer :: conf              ! counts number of conformations
        real(dp) :: cn                 ! auxilary variable for Poisson Eq


        real(dp), parameter :: tolconst = 1.0e-9_dp  ! tolerance for constA and constB 


        !     .. executable statements 

        n=nr                      ! size vector neq=5*nz x=(pi,psi,rhopolA,rhopolB,xpolC)

        do i=1,n                  ! init x 
            xsol(i)= x(i)          ! solvent volume fraction 
            psi(i) = x(i+n)        ! potential
            rhopolAin(i)=x(i+2*n)
            rhopolBin(i)=x(i+3*n)
            xpolCin(i)=x(i+4*n)    ! volume fraction C-polymer
        enddo

        psiSurf = psi(1)          ! surface potentail

        do i=1,n                  ! init volume fractions 
            xpolAB(i)  = 0.0_dp     ! AB polymer volume fraction 
            xpolC(i)   = 0.0_dp     ! C polymer volume fraction 
            rhopolA(i) = 0.0_dp     ! A polymer density 
            rhopolB(i) = 0.0_dp
            rhopolC(i) = 0.0_dp
        
            xNa(i)   = expmu%Na*(xsol(i)**vNa)*dexp(-psi(i)*zNa) ! ion plus volume fraction
            xK(i)    = expmu%K*(xsol(i)**vK)*dexp(-psi(i)*zK)    ! ion plus volume fraction
            xCa(i)   = expmu%Ca*(xsol(i)**vCa)*dexp(-psi(i)*zCa) ! ion divalent pos volume fraction
            xNaCl(i) = expmu%NaCl*(xsol(i)**vNaCl)               ! ion pair  volume fraction
            xKCl(i)  = expmu%KCl*(xsol(i)**vKCl)                 ! ion pair  volume fraction
            xCl(i)   = expmu%Cl*(xsol(i)**vCl)*dexp(-psi(i)*zCl) ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*dexp(-psi(i))      ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*dexp(+psi(i))      ! OH-  volume fraction
       
            xA(1)= xHplus(i)/(K0a(1)*(xsol(i)**deltavA(1)))     ! AH/A-
            xA(2)= (xNa(i)/vNa)/(K0a(2)*(xsol(i)**deltavA(2)))  ! ANa/A-
            xA(3)= (xCa(i)/vCa)/(K0a(3)*(xsol(i)**deltavA(3)))  ! ACa+/A-
       
            sumxA=xA(1)+xA(2)+xA(3)
            constA=(2.0_dp*(rhopolAin(i)*vsol)*(xCa(i)/vCa))/(K0a(4)*(xsol(i)**deltavA(4))) ! A2Ca/(A-)^2
            if(constA<=tolconst) then 
                fdisA(1,i)=1.0_dp/(1.0_dp+sumxA)
                fdisA(5,i)=0.0_dp
            else
                fdisA(1,i)= (-1.0_dp+dsqrt(1.0_dp+4.0_dp*constA/((sumxA+1.0_dp)**2)))
                fdisA(1,i)= fdisA(1,i)*(sumxA+1.0_dp)/(2.0_dp*constA)
                fdisA(5,i)= (fdisA(1,i)**2)*constA
            endif    
       
            fdisA(2,i)  = fdisA(1,i)*xA(1)                      ! AH 
            fdisA(3,i)  = fdisA(1,i)*xA(2)                      ! ANa 
            fdisA(4,i)  = fdisA(1,i)*xA(3)                      ! ACa+ 
       
            xB(1)= xHplus(i)/(K0b(1)*(xsol(i) **deltavB(1)))    ! BH/B-
            xB(2)= (xNa(i)/vNa)/(K0b(2)*(xsol(i)**deltavB(2)))  ! BNa/B-
            xB(3)= (xCa(i)/vCa)/(K0b(3)*(xsol(i)**deltavB(3)))  ! BCa+/B-
       
       
            sumxB=xB(1)+xB(2)+xB(3)
            constB=(2.0_dp*(rhopolBin(i)*vsol)*(xCa(i)/vCa))/(K0b(4)*(xsol(i)**deltavB(4)))
            if(constB<=tolconst) then
                fdisB(1,i)=1.0_dp/(1.0_dp+sumxB)
                fdisB(5,i)=0.0_dp
            else
                fdisB(1,i)= (-1.0_dp+dsqrt(1.0_dp+4.0_dp*constB/((sumxB+1.0_dp)**2)))
                fdisB(1,i)= fdisB(1,i)*(sumxB+1.0_dp)/(2.0_dp*constB) !B^-
                fdisB(5,i)= (fdisB(1,i)**2)*constB                    ! B2Ca
            endif
       
            fdisB(2,i)  = fdisB(1,i)*xB(1)                      ! BH 
            fdisB(3,i)  = fdisB(1,i)*xB(2)                      ! BNa 
            fdisB(4,i)  = fdisB(1,i)*xB(3)                      ! BCa+ 
       

    !        exppiA(i)=(xsol(i)**vpolA(1))*dexp(-zpolA(1)*psi(i))/fdisA(1,i) ! auxiliary variable
    !        exppiB(i)=(xsol(i)**vpolB(1))*dexp(-zpolB(1)*psi(i))/fdisB(1,i) ! auxiliary variable

    !       exppiA(i)=(xsol(i)**vpolA(2))*dexp(-zpolA(2)*psi(i))/fdisA(2,i) ! auxiliary variable                                           
    !       exppiB(i)=(xsol(i)**vpolB(2))*dexp(-zpolB(2)*psi(i))/fdisB(2,i) ! auxiliary variable   
            ! Na condensed ANa reference state
            exppiA(i)=(xsol(i)**vpolA(3))*dexp(-zpolA(3)*psi(i))/fdisA(3,i) ! auxiliary variable
            exppiB(i)=(xsol(i)**vpolB(3))*dexp(-zpolB(3)*psi(i))/fdisB(3,i) ! auxiliary variable   

       
            !     .. VdW interaction   
            tmp = 0.0_dp
            if((i+VdWcutoffdelta)<=nsize) then 
                do j=minrange(i),i+VdWcutoffdelta
                    tmp = tmp + chis(i,j)*xpolCin(j)
                enddo
            endif
            expVdW=dexp(-VdWepsC*tmp)
            exppiC(i)=(xsol(i)**vpolC)*expVdW ! auxiliary variable
        enddo

        !   .. computation polymer volume fraction 

        qAB = 0.0_dp                 ! init q

        do c=1,cuantasAB             ! loop over cuantas
            pro=1.0_dp               ! initial weight conformation 
            do s=1,nsegAB            ! loop over segments 
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment 
                    pro = pro*exppiA(k)
                else
                    pro = pro*exppiB(k)
                endif
            enddo

            qAB = qAB+pro
            do s=1,nsegAB
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment  !        if(isAmonomer(s).eqv..TRUE.) then ! A segment 
                    rhopolA(k)=rhopolA(k)+pro
                else
                    rhopolB(k)=rhopolB(k)+pro
                endif
            enddo
        enddo

        qC = 0.0_dp                 ! init q   
        do c=1,cuantasC             ! loop over cuantas                                                      
            pro=1.0_dp              ! initial weight conformation                                                   
            do s=1,nsegC            ! loop over segments                
                k=indexchainC(c,s)
                pro = pro*exppiC(k)
            enddo
            qC = qC+pro
            do s=1,nsegC
                k=indexchainC(c,s)
                rhopolC(k)=rhopolC(k)+pro
            enddo
        enddo

        !   .. construction of fcn and volume fraction polymer        

        rhopolAB0=sigmaAB/qAB
        rhopolC0=sigmaC/qC

        do i=1,n
            rhopolA(i)= rhopolAB0*rhopolA(i)/deltaG(i)
            rhopolB(i)= rhopolAB0*rhopolB(i)/deltaG(i)
            rhopolC(i)= rhopolC0*rhopolC(i)/deltaG(i)
       
            do k=1,4               ! polymer volume fraction
                xpolAB(i)=xpolAB(i)+rhopolA(i)*fdisA(k,i)*vpolA(k)*vsol  & 
                    +rhopolB(i)*fdisB(k,i)*vpolB(k)*vsol
            enddo    
            xpolAB(i)=xpolAB(i)+rhopolA(i)*(fdisA(5,i)*vpolA(5)*vsol/2.0_dp)
            xpolAB(i)=xpolAB(i)+rhopolB(i)*(fdisB(5,i)*vpolB(5)*vsol/2.0_dp)
       
            xpolC(i)=rhopolC(i)*vpolC*vsol
       
            f(i)=xpolAB(i)+xpolC(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0_dp
       
            rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                zpolA(1)*fdisA(1,i)*rhopolA(i)*vsol+ zpolA(4)*fdisA(4,i)*rhopolA(i)*vsol+ &
                zpolB(1)*fdisB(1,i)*rhopolB(i)*vsol+ zpolB(4)*fdisB(4,i)*rhopolB(i)*vsol         
       
            !   ..  total charge density in units of vsol
        enddo  !  .. end computation polymer density and charge density  

        ! .. electrostatics 

        sigmaqSurf=0.0_dp ! charge regulating surface charge 
        psi(n+1)=0.0_dp   ! bulk potential

        !    .. Poisson Eq 

        f(n+1)= -0.5_dp*((psi(2)-psi(1)) + sigmaqSurf +rhoq(1)*constqW)      !     boundary

        do i=2,n
            f(n+i)= -0.5_dp*(psi(i+1)-2.0_dp*psi(i) + psi(i-1) +rhoq(i)*constqW)
        enddo

        do i=1,n
            f(2*n+i)=rhopolA(i)-rhopolAin(i)
            f(3*n+i)=rhopolB(i)-rhopolBin(i)
            f(4*n+i)=xpolC(i)-xpolCin(i)
        enddo

        iter=iter+1

    end subroutine fcnelectHC

    subroutine fcnelectNoPoly(x,f,nn)

    !     .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW
        use surface 
        use vectornorm

        implicit none

        !     .. scalar arguments
        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        integer(8), intent(in) :: nn


        !     .. declare local variables

        integer :: n                 ! half of n
        integer :: i,j,k,c,s         ! dummy indices
        real(dp) :: norm
        integer :: conf              ! counts number of conformations
!
        integer :: neq_bc           

        !     .. executable statements 
 
        n=nr                       ! size vector neq=5*nz x=(pi,psi,rhopolA,rhopolB,xpolC)

        do i=1,n                   ! init x 
            xsol(i)= x(i)          ! solvent volume fraction 
            psi(i) = x(i+n)        ! potential
        enddo
        

        if(bcflag/="cc") then
            neq_bc=1 
            psiSurf =x(2*n+neq_bc) ! surface potential
        endif 
    
        do i=1,n                  ! init volume fractions 
        
            xNa(i)    = expmu%Na  *(xsol(i)**vNa)*dexp(-psi(i)*zNa)  ! ion plus volume fraction
            xK(i)     = expmu%K   *(xsol(i)**vK) *dexp(-psi(i)*zK)   ! ion plus volume fraction
            xCa(i)    = expmu%Ca  *(xsol(i)**vCa)*dexp(-psi(i)*zCa)  ! ion divalent pos volume fraction
            xNaCl(i)  = expmu%NaCl*(xsol(i)**vNaCl)                  ! ion pair  volume fraction
            xKCl(i)   = expmu%KCl *(xsol(i)**vKCl)                   ! ion pair  volume fraction
            xCl(i)    = expmu%Cl  *(xsol(i)**vCl)*dexp(-psi(i)*zCl)  ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*dexp(-psi(i))          ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*dexp(+psi(i))          ! OH-  volume fraction
            xTB(i)    = expmu%TB  *(xsol(i)**vTB) *dexp(-psi(i))   ! ion plus volume fraction
    
        enddo

        !   .. construction of fcn 
        

        do i=1,n

            f(i)=xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)+xTB(i)-1.0_dp
       
            rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+xTB(i)/vTB
       
            !   ..  total charge density in units of vsol
        enddo 

        ! .. electrostatics 

        ! .. charge regulating surface charge 
        sigmaqSurf=surface_charge(bcflag,psiSurf)
        psi(n+1)= 0.0_dp
       
        ! .. Poisson Eq  
        ! .. Poisson Eq 
  
        f(n+1)= -0.5_dp*(Fplus(1)*(psi(2)-psi(1)) + Fmin(1)*sigmaqSurf +rhoq(1)*constqW)      !     boundary
  
        do i=2,n
            f(n+i)= -0.5_dp*(Fplus(i)*psi(i+1)-2.0_dp*psi(i) + Fmin(i)*psi(i-1) +rhoq(i)*constqW)
        enddo

        ! self consistent boundary conditions

        if(bcflag/='cc') then 
            f(2*n+neq_bc)=psi(1)-psisurf+sigmaqSurf/2.0_dp
        else    
            psisurf=psi(1)+sigmaqSurf/2.0_dp
        endif   
       
        iter=iter+1 

    end subroutine fcnelectNoPoly


    subroutine fcnelect(x,f,nn)

    !     .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW
        use surface 
        use vectornorm

        implicit none

        !     .. scalar arguments
        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        integer(8), intent(in) :: nn


        !     .. declare local variables

        real(dp) :: exppiA(nsize),exppiB(nsize),exppiC(nsize)    ! auxilairy variable for computing P(\alpha) 
        real(dp) :: rhopolAin(nsize),rhopolBin(nsize),xpolCin(nsize)
        real(dp) :: xA(3),xB(3),sumxA,sumxB
        real(dp) :: constA,constB
        real(dp) :: pro,rhopolAB0,rhopolC0
        integer :: n                 ! half of n
        integer :: i,j,k,c,s         ! dummy indices
        real(dp) :: tmp,expVdW 
        real(dp) :: norm
        integer :: conf              ! counts number of conformations
        real(dp) :: cn               ! auxilary variable for Poisson Eq
        integer :: neq_bc           

        real(dp), parameter :: tolconst = 1.0e-9_dp  ! tolerance for constA and constB 


        !     .. executable statements 
        
        n=nr                       ! size vector neq=5*nz x=(pi,psi,rhopolA,rhopolB,xpolC)
        
        do i=1,n                   ! init x 
            xsol(i)= x(i)          ! solvent volume fraction 
            psi(i) = x(i+n)        ! potential
            rhopolAin(i)=x(i+2*n)
            rhopolBin(i)=x(i+3*n)
        enddo
        

        if(bcflag/="cc") then
            neq_bc=1 
            psiSurf =x(4*n+neq_bc)          ! surface potentail
        endif   
  
        do i=1,n                    ! init volume fractions 
            xpolAB(i)  = 0.0_dp     ! AB polymer volume fraction 
            xpolC(i)   = 0.0_dp     ! C polymer volume fraction 
            rhopolA(i) = 0.0_dp     ! A polymer density 
            rhopolB(i) = 0.0_dp

            xNa(i)   = expmu%Na*(xsol(i)**vNa)*dexp(-psi(i)*zNa) ! ion plus volume fraction
            xK(i)    = expmu%K*(xsol(i)**vK)*dexp(-psi(i)*zK)    ! ion plus volume fraction
            xCa(i)   = expmu%Ca*(xsol(i)**vCa)*dexp(-psi(i)*zCa) ! ion divalent pos volume fraction
            xNaCl(i) = expmu%NaCl*(xsol(i)**vNaCl)               ! ion pair  volume fraction
            xKCl(i)  = expmu%KCl*(xsol(i)**vKCl)                 ! ion pair  volume fraction
            xCl(i)   = expmu%Cl*(xsol(i)**vCl)*dexp(-psi(i)*zCl) ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*dexp(-psi(i))      ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*dexp(+psi(i))      ! OH-  volume fraction
       
            xA(1)= xHplus(i)/(K0a(1)*(xsol(i)**deltavA(1)))     ! AH/A-
            xA(2)= (xNa(i)/vNa)/(K0a(2)*(xsol(i)**deltavA(2)))  ! ANa/A-
            xA(3)= (xCa(i)/vCa)/(K0a(3)*(xsol(i)**deltavA(3)))  ! ACa+/A-
       
            sumxA=xA(1)+xA(2)+xA(3)
            constA=(2.0_dp*(rhopolAin(i)*vsol)*(xCa(i)/vCa))/(K0a(4)*(xsol(i)**deltavA(4))) ! A2Ca/(A-)^2
            if(constA<=tolconst) then 
                fdisA(1,i)=1.0_dp/(1.0_dp+sumxA)
                fdisA(5,i)=0.0_dp
            else
                fdisA(1,i)= (-1.0_dp+dsqrt(1.0_dp+4.0_dp*constA/((sumxA+1.0_dp)**2)))
                fdisA(1,i)= fdisA(1,i)*(sumxA+1.0_dp)/(2.0_dp*constA)
                fdisA(5,i)= (fdisA(1,i)**2)*constA
            endif    
       
            fdisA(2,i)  = fdisA(1,i)*xA(1)                      ! AH 
            fdisA(3,i)  = fdisA(1,i)*xA(2)                      ! ANa 
            fdisA(4,i)  = fdisA(1,i)*xA(3)                      ! ACa+ 
       
            xB(1)= xHplus(i)/(K0b(1)*(xsol(i) **deltavB(1)))    ! BH/B-
            xB(2)= (xNa(i)/vNa)/(K0b(2)*(xsol(i)**deltavB(2)))  ! BNa/B-
            xB(3)= (xCa(i)/vCa)/(K0b(3)*(xsol(i)**deltavB(3)))  ! BCa+/B-
       
       
            sumxB=xB(1)+xB(2)+xB(3)
            constB=(2.0_dp*(rhopolBin(i)*vsol)*(xCa(i)/vCa))/(K0b(4)*(xsol(i)**deltavB(4)))
            if(constB<=tolconst) then
                fdisB(1,i)=1.0_dp/(1.0_dp+sumxB)
                fdisB(5,i)=0.0_dp
            else
                fdisB(1,i)= (-1.0_dp+dsqrt(1.0_dp+4.0_dp*constB/((sumxB+1.0_dp)**2)))
                fdisB(1,i)= fdisB(1,i)*(sumxB+1.0_dp)/(2.0_dp*constB) !B^-
                fdisB(5,i)= (fdisB(1,i)**2)*constB                    ! B2Ca
            endif
       
            fdisB(2,i)  = fdisB(1,i)*xB(1)                      ! BH 
            fdisB(3,i)  = fdisB(1,i)*xB(2)                      ! BNa 
            fdisB(4,i)  = fdisB(1,i)*xB(3)                      ! BCa+ 
       

            exppiA(i)=(xsol(i)**vpolA(1))*dexp(-zpolA(1)*psi(i))/fdisA(1,i) ! auxiliary variable
            exppiB(i)=(xsol(i)**vpolB(1))*dexp(-zpolB(1)*psi(i))/fdisB(1,i) ! auxiliary variable

    !       exppiA(i)=(xsol(i)**vpolA(2))*dexp(-zpolA(2)*psi(i))/fdisA(2,i) ! auxiliary variable                                           
    !       exppiB(i)=(xsol(i)**vpolB(2))*dexp(-zpolB(2)*psi(i))/fdisB(2,i) ! auxiliary variable   
    !        Na condensed ANa reference state
    !        exppiA(i)=(xsol(i)**vpolA(3))*dexp(-zpolA(3)*psi(i))/fdisA(3,i) ! auxiliary variable
    !        exppiB(i)=(xsol(i)**vpolB(3))*dexp(-zpolB(3)*psi(i))/fdisB(3,i) ! auxiliary variable   

       
        enddo

        
        !   .. computation polymer volume fraction 

        qAB = 0.0_dp                  ! init q

        do c=1,cuantasAB              ! loop over cuantas
            pro=1.0_dp                ! initial weight conformation 
            do s=1,nsegAB             ! loop over segments 
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment 
                    pro = pro*exppiA(k)
                else
                    pro = pro*exppiB(k)
                endif
            enddo

            qAB = qAB+pro
            do s=1,nsegAB
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment  !        if(isAmonomer(s).eqv..TRUE.) then ! A segment 
                    rhopolA(k)=rhopolA(k)+pro
                else
                    rhopolB(k)=rhopolB(k)+pro
                endif
            enddo
        enddo

      
        !   .. construction of fcn and volume fraction polymer        

        rhopolAB0=sigmaAB/qAB
       

        do i=1,n
            rhopolA(i)= rhopolAB0*rhopolA(i)/deltaG(i)
            rhopolB(i)= rhopolAB0*rhopolB(i)/deltaG(i)
           
       
            do k=1,4               ! polymer volume fraction
                xpolAB(i)=xpolAB(i)+rhopolA(i)*fdisA(k,i)*vpolA(k)*vsol  & 
                    +rhopolB(i)*fdisB(k,i)*vpolB(k)*vsol
            enddo    
            xpolAB(i)=xpolAB(i)+rhopolA(i)*(fdisA(5,i)*vpolA(5)*vsol/2.0_dp)
            xpolAB(i)=xpolAB(i)+rhopolB(i)*(fdisB(5,i)*vpolB(5)*vsol/2.0_dp)
       
            f(i)=xpolAB(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0_dp
       
            rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                zpolA(1)*fdisA(1,i)*rhopolA(i)*vsol+ zpolA(4)*fdisA(4,i)*rhopolA(i)*vsol+ &
                zpolB(1)*fdisB(1,i)*rhopolB(i)*vsol+ zpolB(4)*fdisB(4,i)*rhopolB(i)*vsol         
       
            !   ..  total charge density in units of vsol
        enddo  !  .. end computation polymer density and charge density  

        ! .. electrostatics 
        ! .. charge regulating surface charge 
        sigmaqSurf=surface_charge(bcflag,psiSurf)
        psi(n+1)=0.0_dp

        ! .. Poisson Eq 
  
        f(n+1)= -0.5_dp*(Fplus(1)*(psi(2)-psi(1)) + Fmin(1)*sigmaqSurf +rhoq(1)*constqW)      !     boundary
  
        do i=2,n
            f(n+i)= -0.5_dp*(Fplus(i)*psi(i+1)-2.0_dp*psi(i) + Fmin(i)*psi(i-1) +rhoq(i)*constqW)
        enddo


        do i=1,n
            f(2*n+i)=rhopolA(i)-rhopolAin(i)
            f(3*n+i)=rhopolB(i)-rhopolBin(i)
        enddo

        ! self consistent boundary conditions
       
        if(bcflag/='cc') then 
            f(4*n+neq_bc)=psi(1)-psisurf+sigmaqSurf/2.0_dp
        else    
            psisurf=psi(1)+sigmaqSurf/2.0_dp
        endif   

   
        
        iter=iter+1

        qABL=qAB   ! defined both qAB and qABL, communcicates value     

    end subroutine fcnelect

   
    subroutine fcnneutral(x,f,nn)

    !     .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW
        use vectornorm

        implicit none

        !     .. scalar arguments
        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        integer(8), intent(in) :: nn

        !     .. declare local variables

        real(dp) :: exppiA(nsize),exppiB(nsize),exppiC(nsize)    ! auxilairy variable for computing P(\alpha) 
        real(dp) :: xpolBin(nsize)
        real(dp) :: xA(3),xB(3),sumxA,sumxB
        real(dp) :: constA,constB
        real(dp) :: pro,rhopolAB0,rhopolC0
        integer :: n                 ! half of n
        integer :: i,j,k,c,s         ! dummy indices
        real(dp) :: tmp,expVdW 
        real(dp) :: norm
        integer :: conf              ! counts number of conformations
        real(dp) :: cn               ! auxilary variable for Poisson Eq


        real(dp), parameter :: tolconst = 1.0e-9_dp  ! tolerance for constA and constB 


        !     .. executable statements 

        n=nr                        ! size vector neq=5*nz x=(pi,psi,rhopolA,rhopolB,xpolC)

        do i=1,n                    ! init x 
            xsol(i)= x(i)           ! solvent volume fraction 
            xpolBin(i)=x(i+n)       ! volume fraction B-polymer 
        enddo

        do i=1,n                    ! init volume fractions 
            xpolAB(i)  = 0.0_dp     ! AB polymer volume fraction 
            xpolC(i)   = 0.0_dp     ! C polymer volume fraction 
            rhopolA(i) = 0.0_dp     ! A polymer density 
            rhopolB(i) = 0.0_dp
            rhopolC(i) = 0.0_dp
        
            exppiA(i)=(xsol(i)**vpolA(3)) !*dexp(-zpolA(3)*psi(i))/fdisA(3,i) ! auxiliary variable
            exppiB(i)=(xsol(i)**vpolB(3)) !*dexp(-zpolB(3)*psi(i))/fdisB(3,i) ! auxiliary variable   
            exppiC(i)=(xsol(i)**vpolC)
       
            !     .. VdW interaction   
            tmp = 0.0_dp
            if((i+VdWcutoffdelta)<=nsize) then 
                do j=minrange(i),i+VdWcutoffdelta
                    tmp = tmp + chis(i,j)*xpolBin(j)
                enddo
            endif
            expVdW=dexp(-VdWepsB*tmp)
            exppiB(i)=exppiB(i)*expVdW ! auxiliary variable
        enddo

        !   .. computation polymer volume fraction 

        qAB = 0.0_dp                 ! init q

        do c=1,cuantasAB            ! loop over cuantas
            pro=1.0_dp                ! initial weight conformation 
            do s=1,nsegAB            ! loop over segments 
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment 
                    pro = pro*exppiA(k)
                else
                    pro = pro*exppiB(k)
                endif
            enddo

            qAB = qAB+pro
            do s=1,nsegAB
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment  !        if(isAmonomer(s).eqv..TRUE.) then ! A segment 
                    rhopolA(k)=rhopolA(k)+pro
                else
                    rhopolB(k)=rhopolB(k)+pro
                endif
            enddo
        enddo

        qC = 0.0_dp                 ! init q   
        do c=1,cuantasC             ! loop over cuantas                                                      
            pro=1.0_dp              ! initial weight conformation                                                   
            do s=1,nsegC            ! loop over segments                
                k=indexchainC(c,s)
                pro = pro*exppiC(k)
            enddo
            qC = qC+pro
            do s=1,nsegC
                k=indexchainC(c,s)
                rhopolC(k)=rhopolC(k)+pro
            enddo
        enddo

        !   .. construction of fcn and volume fraction polymer        

        rhopolAB0=sigmaAB/qAB
        rhopolC0=sigmaC/qC

        do i=1,n
            rhopolA(i)= rhopolAB0*rhopolA(i)/deltaG(i)
            rhopolB(i)= rhopolAB0*rhopolB(i)/deltaG(i)
            rhopolC(i)= rhopolC0*rhopolC(i)/deltaG(i)
       
            ! polymer volume fraction A and B all in state 3 ANa and BNa
            xpolAB(i)=(rhopolA(i)*vpolA(3)+rhopolB(i)*vpolB(3))*vsol
            xpolC(i)=rhopolC(i)*vpolC*vsol
       
            f(i)=xpolAB(i)+xpolC(i)+xsol(i)-1.0_dp
            f(n+i)=rhopolB(i)*vpolB(3)*vsol-xpolBin(i)
        enddo  !  .. end computation polymer density 

    !    norm=l2norm(f,2*n)
        iter=iter+1

    !    print*,'iter=', iter ,'norm=',norm

    end subroutine fcnneutral


    !     .. function solves for bulk volume fraction 

    subroutine fcnbulk(x,f,nn)   

        !     .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use physconst
        use vectornorm

        implicit none

        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out):: f(neq)


        !     .. local variables

        real(dp) :: phiNaCl,phiNa,phiCl,phiK,phiKCl
        real(dp) :: deltavolNaCl,deltavolKCl,norm

        !     .. executable statements 


        deltavolNaCl=(vNaCl-vNa-vCl)
        deltavolKCl=(vKCl-vNa-vCl)

        phiNa   =x(1)
        phiCl   =x(2)
        phiNaCl =x(3)
        phiK    =x(4)
        phiKCl  =x(5)

        !     .. condensation equilibrium eq NaCl

        f(1)= phiNaCl-(phiNa*phiCl*K0ionNa*vsol*vNaCl/(vNa*vCl*vsol))*      & 
         ((1.0_dp -phiNaCl-phiNa-phiCl-phiK-phiKCl-xbulk%Hplus-xbulk%OHmin-xbulk%Ca)**deltavolNaCl)

        !     .. charge neutrality
        f(2)=phiNa/vNa-phiCl/vCl+2.0_dp*xbulk%Ca/vCa+xbulk%Hplus-xbulk%OHmin+phiK/vK

        !     .. conservation of number NaCl

        f(3)=phiNa/vNa+phiCl/vCl+2.0_dp*phiNaCl/vNaCl + phiKCL/vKCL &    
             -2.0_dp*vsol*(Na/1.0e24_dp)*cNaCl-dabs(xbulk%Hplus-xbulk%OHmin) &
            -2.0_dp*xbulk%Ca/vCa-vsol*(Na/1.0e24_dp)*cKCl

        !     .. condensation equilibrium eq KCl

        f(4)= phiKCl-(phiK*phiCl*K0ionK*vsol*vKCl/(vK*vCl*vsol))* &
             ((1.0_dp -phiNaCl-phiNa-phiCl-phiK-phiKCl-xbulk%Hplus-xbulk%OHmin-xbulk%Ca)**deltavolKCl)

        !     .. conservation of number K
        f(5)= phiK/vK+phiKCl/vKCl-vsol*(Na/1.0e24_dp)*cKCl

    !    norm=l2norm(f,5)
        iter=iter+1
     
    !    print*,'iter=', iter ,'norm=',norm

    end subroutine fcnbulk



 ! selects correct fcn function 

    subroutine set_fcn

        use globals
        use fcnpointer
        
        implicit none   

        select case (sysflag)
            case ("elect") 
                fcnptr => fcnelect
            case ("electnopoly") 
                fcnptr => fcnelectNoPoly 
            case ("electHC") 
                fcnptr => fcnelectHC
            case ("neutral") 
                fcnptr => fcnneutral
            case ("bulk water") 
                 fcnptr => fcnbulk
            case default
                print*,"Error in call to solver subroutine"    
                print*,"Wrong value sysflag : ", sysflag
                stop
        end select  
    
    end subroutine set_fcn





end module listfcn
