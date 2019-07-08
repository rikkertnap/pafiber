! --------------------------------------------------------------|
! fcn.f90:                                                    |
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

   
    !
    subroutine fcnelectNoPoly(x,f,nn)

       !     .. variables and constant declaractions 

        use globals
        use volume
        use field
        use parameters
        use surface 
        use vectornorm

        !     .. arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        integer(8), intent(in) :: nn


        !     .. declare local variables

        integer :: n                 ! n=nr
        integer :: i,j,k,c,s         ! dummy indices
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
        
            xNa(i)    = expmu%Na  *(xsol(i)**vNa)*exp(-psi(i)*zNa)  ! ion plus volume fraction
            xK(i)     = expmu%K   *(xsol(i)**vK) *exp(-psi(i)*zK)   ! ion plus volume fraction
            xCa(i)    = expmu%Ca  *(xsol(i)**vCa)*exp(-psi(i)*zCa)  ! ion divalent pos volume fraction
            xNaCl(i)  = expmu%NaCl*(xsol(i)**vNaCl)                  ! ion pair  volume fraction
            xKCl(i)   = expmu%KCl *(xsol(i)**vKCl)                   ! ion pair  volume fraction
            xCl(i)    = expmu%Cl  *(xsol(i)**vCl)*exp(-psi(i)*zCl)  ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*exp(-psi(i))          ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*exp(+psi(i))          ! OH-  volume fraction
            ! xTB(i)    = expmu%TB  *(xsol(i)**vTB) *exp(-psi(i)*zTB)   ! ion plus volume fraction
            ! .. warning fcn bulk solution see fcnbulk has not TBCl thus XTb=0
            xTB(i) = 0.0_dp     
        enddo

        !   .. construction of fcn 
        

        do i=1,n

            f(i)=xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)+xTB(i)-1.0_dp
       
            rhoq(i)= zNa*xNa(i)/vNa+zCa*xCa(i)/vCa+zK*xK(i)/vK+zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+zTB*xTB(i)/vTB
       
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


    ! ligand pp adsorbtion 
    subroutine fcnelectligand(x,f,nn)

        !     .. variables and constant declaractions 

        use globals
        use volume
        use field
        use parameters
        use surface 
        use vectornorm

        !     .. scalar arguments
        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        integer(8), intent(in) :: nn   ! nn=neq 

        !     .. declare local variables

        integer :: n                 ! n=nr 
        integer :: i,t               ! dummy indices
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
            xNa(i)    = expmu%Na  *(xsol(i)**vNa)*exp(-psi(i)*zNa)  ! ion plus volume fraction
            xK(i)     = expmu%K   *(xsol(i)**vK) *exp(-psi(i)*zK)   ! ion plus volume fraction
            xCa(i)    = expmu%Ca  *(xsol(i)**vCa)*exp(-psi(i)*zCa)  ! ion divalent pos volume fraction
            xNaCl(i)  = expmu%NaCl*(xsol(i)**vNaCl)                  ! ion pair  volume fraction
            xKCl(i)   = expmu%KCl *(xsol(i)**vKCl)                   ! ion pair  volume fraction
            xCl(i)    = expmu%Cl  *(xsol(i)**vCl)*exp(-psi(i)*zCl)  ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*exp(-psi(i))          !mo H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*exp(+psi(i))           ! OH-  volume fraction
            xTB(i)    = expmu%TB  *(xsol(i)**vTB) *exp(-psi(i)*zTB)   ! ion plus volume fraction
            xTM(i)    = expmu%TM  *(xsol(i)**vTM) *exp(-psi(i)*zTM)   ! ion plus volume fraction  
            xNO3(i)   = expmu%NO3  *(xsol(i)**vNO3)*exp(-psi(i)*zNO3) 
            
        enddo
            
        do t=1,5 ! loop ligand types 
            do i=1,n
                xpp(i,t)  = expmu%pp(t)  *(xsol(i)**vpp(t)) *dexp(-psi(i)*zpp(t))   
            end do   
        enddo

        !   .. construction of fcn 
        do i=1,n
            f(i) = 0.0_dp
            rhoq(i) = 0.0_dp
            do t=1,5
                f(i) = f(i)+xpp(i,t)
                rhoq(i) = rhoq(i) + zpp(t) * xpp(i,t)/vpp(t)
            enddo   
            f(i)=f(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)+xTB(i)+&
                xTM(i)+xNO3(i)-1.0_dp
            rhoq(i)=rhoq(i)+zNa*xNa(i)/vNa+zCa*xCa(i)/vCa +zK*xK(i)/vK +zCl*xCl(i)/vCl+xHplus(i)-xOHmin(i)+&
                zTB*xTB(i)/vTB+zTM*xTM(i)/vTM+zNO3*xNO3(i)/vNO3
            
            !   ..  total charge density in units of vsol
        enddo 

        ! .. electrostatics 

        ! .. charge regulating surface charge 
        sigmaqSurf=surface_charge(bcflag,psiSurf)
        
        if(runflag=="rangenr") then
            psi(n+1) = psi(n)
        else 
            psi(n+1) = 0.0_dp
        endif 
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

    end subroutine fcnelectligand


    subroutine fcnpafiber(x,f,nn)

        !     .. variables and constant declaractions 

        use globals
        use volume
        use field
        use parameters
        use surface 
        use vectornorm

        !     .. scalar arguments
        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        integer(8), intent(in) :: nn   ! nn=neq 

        !     .. declare local variables

        integer :: n                 ! n=nr 
        integer :: i,t               ! dummy indices
        integer :: neq_bc           
        real(dp) :: norm
        real(dp) :: xA

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
            xNa(i)    = expmu%Na  *(xsol(i)**vNa)*exp(-psi(i)*zNa)  ! ion Na+ volume fraction
            xK(i)     = expmu%K   *(xsol(i)**vK) *exp(-psi(i)*zK)   ! ion K+ volume fraction
            xRb(i)    = expmu%Rb  *(xsol(i)**vRb) *exp(-psi(i)*zRb) ! ion Rb+ volume fraction
            xIm(i)    = expmu%Im  *(xsol(i)**vIm) *exp(-psi(i)*zIm)
            xCa(i)    = expmu%Ca  *(xsol(i)**vCa)*exp(-psi(i)*zCa)   ! ion divalent pos volume fraction
            xNaCl(i)  = expmu%NaCl*(xsol(i)**vNaCl)                  ! ion pair  volume fraction
            xKCl(i)   = expmu%KCl *(xsol(i)**vKCl)                   ! ion pair  volume fraction
            xCl(i)    = expmu%Cl  *(xsol(i)**vCl)*exp(-psi(i)*zCl)   ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*exp(-psi(i))           ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*exp(+psi(i))           ! OH-  volume fraction    
        enddo   
        

        if(isChargeRegularization) then 
            do i=1,n
                xA = xHplus(i)/(K0a*xsol(i))     ! AH/A-                                                       
                fdispa(i)  = 1.0_dp/(1.0_dp+xA)               ! A-
            enddo
        else
            do i=1,n
                fdispa(i) = 1.0_dp 
            enddo   
        endif    

        
        !   .. construction of fcn 
        do i=1,n
            f(i)=xpa(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i) +&
                    xRb(i)+ xIm(i)/vIm  -1.0_dp
              
            rhoq(i)=zpa*fdispa(i)*rhoEpa(i)*vsol+zNa*xNa(i)/vNa+zCa*xCa(i)/vCa +zK*xK(i)/vK +& 
                zCl*xCl(i)/vCl+zRb*xRb(i)/vRb +zIm*xIm(i)/vIm  + xHplus(i)-xOHmin(i) 
            !   ..  total charge density in units of vsol
        enddo 

        ! .. electrostatics 

        ! .. charge regulating surface charge 
        sigmaqSurf=surface_charge(bcflag,psiSurf)
        
        if(runflag=="rangenr") then
            psi(n+1) = psi(n)
        else 
            psi(n+1) = 0.0_dp
        endif 
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

        n=neq
        norm=l2norm(f,n)
        print*,'iter=', iter ,'norm=',norm

    end subroutine fcnpafiber

    

    subroutine fcnpafiberInt(x,f,nn)

        !     .. variables and constant declaractions 

        use globals
        use volume
        use field
        use parameters
        use surface 
        use vectornorm

        !     .. scalar arguments
        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        integer(8), intent(in) :: nn   ! nn=neq 

        !     .. declare local variables
 
        real(dp) :: rhoIm(nsize)
        integer :: n                 ! n=nr 
        integer :: i,t               ! dummy indices
        integer :: neq_bc           
        real(dp) :: norm
        real(dp) :: xA

        !     .. executable statements 
 
        n=nr                       ! size vector neq=5*nz x=(pi,psi,rhopolA,rhopolB,xpolC)

        do i=1,n                   ! init x 
            xsol(i)= x(i)          ! solvent volume fraction 
            psi(i) = x(i+n)        ! potential
            rhoIm(i) = x(i+2*n)    ! imidazolium density
        enddo
        

        if(bcflag/="cc") then
            neq_bc=1 
            psiSurf =x(2*n+neq_bc) ! surface potential
        endif 
    
        do i=1,n                  ! init volume fractions 
            xNa(i)    = expmu%Na  *(xsol(i)**vNa)*exp(-psi(i)*zNa)  ! ion Na+ volume fraction
            xK(i)     = expmu%K   *(xsol(i)**vK) *exp(-psi(i)*zK)   ! ion K+ volume fraction
            xRb(i)    = expmu%Rb  *(xsol(i)**vRb) *exp(-psi(i)*zRb) ! ion Rb+ volume fraction
           
            xIm(i)    = expmu%Im  *(xsol(i)**vIm) *exp(-psi(i)*zIm)*exp(rhoIm(i)*epsIm)
           
            xCa(i)    = expmu%Ca  *(xsol(i)**vCa)*exp(-psi(i)*zCa)   ! ion divalent pos volume fraction
            xNaCl(i)  = expmu%NaCl*(xsol(i)**vNaCl)                  ! ion pair  volume fraction
            xKCl(i)   = expmu%KCl *(xsol(i)**vKCl)                   ! ion pair  volume fraction
            xCl(i)    = expmu%Cl  *(xsol(i)**vCl)*exp(-psi(i)*zCl)   ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*exp(-psi(i))           ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*exp(+psi(i))           ! OH-  volume fraction    
        enddo   
        

        if(isChargeRegularization) then 
            do i=1,n
                xA = xHplus(i)/(K0a*xsol(i))     ! AH/A-                                                       
                fdispa(i)  = 1.0_dp/(1.0_dp+xA)               ! A-
            enddo
        else
            do i=1,n
                fdispa(i) = 1.0_dp 
            enddo   
        endif    

        
        !   .. construction of fcn 
        do i=1,n
            f(i)=xpa(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i) +&
                    xRb(i)+ xIm(i)/vIm  -1.0_dp
              
            rhoq(i)=zpa*fdispa(i)*rhoEpa(i)*vsol+zNa*xNa(i)/vNa+zCa*xCa(i)/vCa +zK*xK(i)/vK +& 
                zCl*xCl(i)/vCl+zRb*xRb(i)/vRb +zIm*xIm(i)/vIm  + xHplus(i)-xOHmin(i) 
            !   ..  total charge density in units of vsol
        enddo 

        ! .. electrostatics 

        ! .. charge regulating surface charge 
        sigmaqSurf=surface_charge(bcflag,psiSurf)
        
        if(runflag=="rangenr") then
            psi(n+1) = psi(n)
        else 
            psi(n+1) = 0.0_dp
        endif 
        ! .. Poisson Eq 
  
        f(n+1)= -0.5_dp*(Fplus(1)*(psi(2)-psi(1)) + Fmin(1)*sigmaqSurf +rhoq(1)*constqW)      !     boundary
  
        do i=2,n
            f(n+i)= -0.5_dp*(Fplus(i)*psi(i+1)-2.0_dp*psi(i) + Fmin(i)*psi(i-1) +rhoq(i)*constqW)
            f(2*n+i) = xIm(i)-rhoIm(i)*vIm*vsol
        enddo




        ! self consistent boundary conditions

        if(bcflag/='cc') then 
            f(2*n+neq_bc)=psi(1)-psisurf+sigmaqSurf/2.0_dp
        else    
            psisurf=psi(1)+sigmaqSurf/2.0_dp
        endif   
       
        iter=iter+1 

        n=neq
        norm=l2norm(f,n)
        print*,'iter=', iter ,'norm=',norm

    end subroutine fcnpafiberInt


    !     .. function solves for bulk volume fraction 

    subroutine fcnbulk(x,f,nn)   

        !     .. variables and constant declaractions 

        use globals
        use volume
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



    !     .. function solves for bulk volume fraction 

    subroutine fcnbulkligand(x,f,nn)   

        !     .. variables and constant declaractions 

        use globals
        use volume
        use field
        use parameters
        use physconst
        use vectornorm
        use molecules
        use vectornorm

        implicit none

        !     .. scalar arguments
        integer(8), intent(in) :: nn

        !     .. array arguments
        real(dp), intent(in) :: x(neq)
        real(dp), intent(out):: f(neq)

        !     .. local variables

        real(dp) :: fppin(5),fppout(5)
        real(dp) :: cppbulk, psisol
        real(dp) :: xppbulkin, rhoqppbulkin
        real(dp) :: xppbulkout, rhoqppbulkout
        real(dp) :: phisol,phiKin,phiClin,phiClout,phiKout
        real(dp) :: xA,xB,xBprime,xE,xF
        real(dp) :: sumx, sumfpp, deltacharge
        integer :: t
        real(dp) :: norm

        !     .. executable statements 
        !     .. input vector 

        sumfpp=0.0_dp
        do t=1,4
            fppin(t)=x(t)
            sumfpp=sumfpp+fppin(t)
        enddo    
        fppin(5)=1.0_dp-sumfpp
        phiClin   = x(5)
        phiKin    = x(6)


        cppbulk=cpp*(Na/1.0e24_dp) !  .. bulk concentration  cpp in mol/liter cppbulk in ligands/nm^3
        xppbulkin=0.0_dp             
        rhoqppbulkin=0.0_dp
        do t=1,5
            xppbulkin=xppbulkin+fppin(t)*vpp(t)
            rhoqppbulkin=rhoqppbulkin+fppin(t)*zpp(t)
        enddo
        xppbulkin = xppbulkin*cppbulk*vsol      !  .. bulk ligand volume fraction 
        rhoqppbulkin = rhoqppbulkin*cppbulk     !  .. charge density ligand 

        phisol=1.0_dp-phiClin-phiKin-xbulk%TB-xbulk%Hplus-xbulk%OHmin-xppbulkin
      
        phisol=phisol-xbulk%TM-xbulk%NO3! -xbulk%Na this needs to be checked

    !    print*,"phisol=",phisol


        ! pKpp(1)  = 2.26_dp  ! POH2COOH <=> POHCOOH- + H+ : A<=> B
        ! pKpp(2) =  4.6_dp   ! POHCOOH- <=> POHCOO2- + H+ : B<=> E
        ! pKpp(3) =  5.4_dp   ! POHCOOH- <=> POCOOH2- + H+ : B<=> C
        ! pKpp(4) =  6.9_dp   ! POCOOH2- <=> POCOO3- + H+  : C<=> F
        ! pKpp(5) =  7.8_dp   ! POHCOO2- <=> POCOO3- + H+  : E<=> F

        !  .. equilibrium eq AH2BH <=> AHBH^- +H^+  A<=> B

        xA = K0pp(1)*(vpp(AHBH)/vpp(AH2BH))*(phisol)/xbulk%Hplus

        !  .. equilibrium eq AHBH^- <=> ABH^2- +H^+    B<=>C 

        xB = K0pp(3)*(vpp(ABH)/vpp(AHBH))*(phisol)/xbulk%Hplus
       
        !   .. equilibrium eq AHBH^- <=> AHB^2- +H^+    B<=>E 

        xBprime = K0pp(2)*(vpp(AHB)/vpp(AHBH))*(phisol)/xbulk%Hplus

        !   .. equilibrium eq AHB^2- <=> AB^3- +H^+    E<=>F

        xE = K0pp(5)*(vpp(AB)/vpp(AHB))*(phisol)/xbulk%Hplus


        sumx=xA + xA*xB + xA*xBprime + xA*xBprime*xE

        fppout(AH2BH) = 1.0_dp/(1.0_dp+sumx)
        fppout(AHBH)  = fppout(AH2BH) * xA
        fppout(ABH)   = fppout(AHBH)  * xB
        fppout(AHB)   = fppout(AHBH)  * xBprime 
        fppout(AB)    = fppout(AHB)   * xE

        rhoqppbulkout=0.0_dp
        do t=1,5
            rhoqppbulkout=rhoqppbulkout+fppout(t)*zpp(t)
        enddo 

        rhoqppbulkout= rhoqppbulkout*cppbulk 
    
        deltacharge = xbulk%Hplus/vsol-xbulk%OHmin/vsol+rhoqppbulkout  ! number density   

        if(deltacharge<0) then
            ! delta [Cl^-]=0
            phiClout = xbulk%Cl    ! NaCl and TBCl no extra Cl added
            phiKout  = abs(deltacharge)*vK*vsol   ! added KOH  
        else if(deltacharge>0) then 
            ! delta [K^+]=0
            phiClout = deltacharge*vCl*vsol +xbulk%Cl  ! added HCl 
            phiKout  = 0.0_dp     ! no added KCl   
        else  !deltacharge==0
            phiClout = xbulk%Cl   ! no HCL
            phiKout  = 0.0_dp     ! no KOH 
        endif    


        do t=1,4
            f(t)=fppout(t)-fppin(t)
        enddo
        f(5)=phiClout-phiClin
        f(6)=phiKout -phiKin
    
       ! norm=l2norm(f,6)
        iter=iter+1
     
       ! print*,'iter=', iter ,'norm=',norm

    end subroutine fcnbulkligand


        !     .. function solves for bulk volume fraction 

    subroutine fcnbulkligandHCL(x,f,nn)   

        !     .. variables and constant declaractions 

        use globals
        use volume
        !use chains
        use field
        use parameters
        use physconst
        use vectornorm
        use molecules
        use vectornorm

        implicit none

        !     .. scalar arguments
        integer(8), intent(in) :: nn

        !     .. array arguments
        real(dp), intent(in) :: x(neq)
        real(dp), intent(out):: f(neq)

        !     .. local variables

        real(dp) :: fppin(5),fppout(5)
        real(dp) :: cppbulk, psisol
        real(dp) :: xppbulkin, rhoqppbulkin
        real(dp) :: xppbulkout, rhoqppbulkout
        real(dp) :: phisol,phiKin,phiClin,phiClout,phiKout
        real(dp) :: phiNO3in,phiNO3out
        real(dp) :: xA,xB,xBprime,xE,xF
        real(dp) :: sumx, sumfpp, deltacharge
        integer :: t
        real(dp) :: norm

        !     .. executable statements 
        !     .. input vector 

        sumfpp=0.0_dp
        do t=1,4
            fppin(t)=x(t)
            sumfpp=sumfpp+fppin(t)
        enddo    
        fppin(5)=1.0_dp-sumfpp

        if(isbulkHCl) then 
            phiClin   = x(5)
        else
            phiNO3in   = x(5)
        endif    
        phiKin    = x(6)


        cppbulk=cpp*(Na/1.0e24_dp) !  .. bulk concentration  cpp in mol/liter cppbulk in ligands/nm^3
        xppbulkin=0.0_dp             
        rhoqppbulkin=0.0_dp
        do t=1,5
            xppbulkin=xppbulkin+fppin(t)*vpp(t)
            rhoqppbulkin=rhoqppbulkin+fppin(t)*zpp(t)
        enddo
        xppbulkin = xppbulkin*cppbulk*vsol      !  .. bulk ligand volume fraction 
        rhoqppbulkin = rhoqppbulkin*cppbulk     !  .. charge density ligand 

        if(isbulkHCl) then 
            phisol=1.0_dp-phiClin-phiKin-xbulk%TB-xbulk%Hplus-xbulk%OHmin-xppbulkin
            phisol=phisol-xbulk%TM-xbulk%NO3! -xbulk%Na this needs to be checked assumes no NaCl added !!!
        else
            phisol=1.0_dp-phiNO3in-phiKin-xbulk%TB-xbulk%Hplus-xbulk%OHmin-xppbulkin
            phisol=phisol-xbulk%TM-xbulk%Cl!
        endif    
    !    print*,"phisol=",phisol


        ! pKpp(1)  = 2.26_dp  ! POH2COOH <=> POHCOOH- + H+ : A<=> B
        ! pKpp(2) =  4.6_dp   ! POHCOOH- <=> POHCOO2- + H+ : B<=> E
        ! pKpp(3) =  5.4_dp   ! POHCOOH- <=> POCOOH2- + H+ : B<=> C
        ! pKpp(4) =  6.9_dp   ! POCOOH2- <=> POCOO3- + H+  : C<=> F
        ! pKpp(5) =  7.8_dp   ! POHCOO2- <=> POCOO3- + H+  : E<=> F

        !  .. equilibrium eq AH2BH <=> AHBH^- +H^+  A<=> B

        xA = K0pp(1)*(vpp(AHBH)/vpp(AH2BH))*(phisol)/xbulk%Hplus

        !  .. equilibrium eq AHBH^- <=> ABH^2- +H^+    B<=>C 

        xB = K0pp(3)*(vpp(ABH)/vpp(AHBH))*(phisol)/xbulk%Hplus
       
        !   .. equilibrium eq AHBH^- <=> AHB^2- +H^+    B<=>E 

        xBprime = K0pp(2)*(vpp(AHB)/vpp(AHBH))*(phisol)/xbulk%Hplus

        !   .. equilibrium eq AHB^2- <=> AB^3- +H^+    E<=>F

        xE = K0pp(5)*(vpp(AB)/vpp(AHB))*(phisol)/xbulk%Hplus


        sumx=xA + xA*xB + xA*xBprime + xA*xBprime*xE

        fppout(AH2BH) = 1.0_dp/(1.0_dp+sumx)
        fppout(AHBH)  = fppout(AH2BH) * xA
        fppout(ABH)   = fppout(AHBH)  * xB
        fppout(AHB)   = fppout(AHBH)  * xBprime 
        fppout(AB)    = fppout(AHB)   * xE

        rhoqppbulkout=0.0_dp
        do t=1,5
            rhoqppbulkout=rhoqppbulkout+fppout(t)*zpp(t)
        enddo 

        rhoqppbulkout= rhoqppbulkout*cppbulk 
    
        deltacharge = xbulk%Hplus/vsol-xbulk%OHmin/vsol+rhoqppbulkout  ! number density   

        if(deltacharge<0) then

            ! delta [Cl^-]=0
            if(isbulkHCl) then
                phiClout = xbulk%Cl    ! NaCl and TBCl no extra Cl added
            else
                phiNO3out = xbulk%NO3
            endif    
            phiKout  = abs(deltacharge)*vK*vsol   ! added KOH  
        
        else if(deltacharge>0) then 
        
            ! delta [K^+]=0
            if(isbulkHCl) then
                phiClout = deltacharge*vCl*vsol +xbulk%Cl  ! added HCl 
            else
                phiNO3out = deltacharge*vNO3*vsol +xbulk%NO3  ! added HCl     
            endif
            phiKout  = 0.0_dp     ! no added KCl   
        
        else  !deltacharge==0
        
            if(isbulkHCl) then
                phiClout = xbulk%Cl   ! no HCL
            else
                phiNO3out= xbulk%NO3   ! no HNO3
            endif    
            phiKout  = 0.0_dp     ! no KOH 
        
        endif    


        do t=1,4
            f(t)=fppout(t)-fppin(t)
        enddo

        if(isbulkHCl) then
            f(5)=phiClout-phiClin
        else
            f(5)=phiNO3out-phiNO3in
        endif    
        f(6)=phiKout -phiKin
    
       ! norm=l2norm(f,6)
        iter=iter+1
     
       ! print*,'iter=', iter ,'norm=',norm

    end subroutine fcnbulkligandHCl


    !     set constrains on vector x  depending on systype value

    subroutine set_contraints(constr)

        use precision_definition
        use globals, only : sysflag, neq, nsize ,bcflag 

        implicit none
            
        real(dp), intent(inout):: constr(:)

        integer :: i, neqint ,neq_bc 

        neqint=int(neq,kind(neqint))     ! explict conversion from integer(8) to integer
          
        neq_bc=0 
        if(bcflag/="cc") neq_bc=neq_bc+1
        
        select case (sysflag)
        case ("electnopoly")     
            do i=1,nsize
                constr(i)=1.0_dp
                constr(i+nsize)=0.0_dp   !  electrostatic potential
            enddo
            do i=1,neq_bc                  ! surface electrostatic potential if bcflag/=cc
                constr(i+2*nsize)=0.0_dp
            enddo    
        case ("electligand")     
            do i=1,nsize
                constr(i)=1.0_dp
                constr(i+nsize)=0.0_dp   !  electrostatic potential
            enddo   
            do i=1,neq_bc                  ! surface electrostatic potential if bcflag/=cc
                constr(i+2*nsize)=0.0_dp
            enddo 
      
        case ("pafiberIm")     
            do i=1,nsize
                constr(i)=1.0_dp
                constr(i+nsize)=0.0_dp   !  electrostatic potential
                constr(i+2*nsize)=0.0_dp   !  electrostatic potential
                
            enddo  
            do i=1,neq_bc                  ! surface electrostatic potential if bcflag/=cc
                constr(i+3*nsize)=0.0_dp
            enddo  
        case ("bulk water")             
            do i=1,nsize
                constr(i)=1.0_dp     
            enddo     
        case ("bulk ligand")                 !  neutral polymers
            do i=1,neqint
                constr(i)=1.0_dp
            enddo 
        case default
            do i=1,neqint
                constr(i)=1.0_dp
            enddo 
        end select  

    end subroutine set_contraints


    subroutine set_fcn

        use globals
        use fcnpointer
        
        implicit none   

        select case (sysflag)
            case ("electnopoly") 
                fcnptr => fcnelectNoPoly 
            case ("electligand") 
                fcnptr => fcnelectligand 
            case ("bulk water") 
                fcnptr => fcnbulk
            case ("bulk ligand") 
                !fcnptr => fcnbulkligand
                fcnptr => fcnbulkligandHCl
            case ("pafiber")
                fcnptr => fcnpafiber
            case ("pafiberIm")
                fcnptr => fcnpafiber
            case default
                print*,"Error in call to solver subroutine"    
                print*,"Wrong value sysflag : ", sysflag
                stop
        end select  
    
    end subroutine set_fcn





end module listfcn
