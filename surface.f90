module surface 
   
    use globals, only :    AH2BH, AHBH, AHB, ABH, AB, SuOH, SuCl, Su, SuNO3
    use mathconst
   
    implicit none
    
    !   different surface states

    real(dp) :: fdisS(9)            ! fraction of different surface states
    real(dp) :: gdisS(9)            ! fraction of different surface states, only used for bc=pc
    
    real(dp) :: KS(7)               ! experimemtal equilibruim constant 
    real(dp) :: pKS(7)              ! experimental equilibruim constant pKS= -log[KS]   
    real(dp) :: K0S(7)              ! intrinsic equilibruim constant
    real(dp) :: qS(9)               ! charge  real because possitbilito offractional charge 
    real(dp) :: cap                 ! capacitance
 
    real(dp) :: sigmaSurf           ! surface density of acid on surface in nm^2
    real(dp) :: sigmaqSurf          ! surface charge density on surface in nm^2
    real(dp) :: fdisR               ! fraction of surface ligand density sites not bound 

    real(dp) :: psiSurf             ! surface potential     
    
   ! taurine
    real(dp) :: fdisTaL(4),fdisTaR(4) ! fraction of different surface states
    real(dp) :: KTa(3)               ! experimemtal equilibruim constant 
    real(dp) :: pKTa(3)              ! experimental equilibruim constant pKS= -log[KS]   
    real(dp) :: K0Ta(3)              ! intrinsic equilibruim constant
    real(dp) :: qTa(4)               ! charge 

contains
 
    subroutine init_surface(bc)

        character(len=2) :: bc
    
        select case (bc)
            case ("qu")  
                call  init_surface_quartz()
            case ("cl")
                call  init_surface_clay()
            case ("ca")  
                call init_surface_calcite()
            case ("ta") 
                call init_surface_taurine()
            case ("cc")
                call init_surface_constcharge()
            case ("pp")
                call init_surface_pp(sigmaSurf)  ! pp == phoshonatepropionate
            case ("pd")
                call init_surface_pp_dynamic(sigmaSurf)   
            case ("pc")
                call init_surface_pp_dynamic_cond(sigmaSurf)      
            case default
                print*,"bc does not match qu, cl, ca, ta, cc, pp, pd, or pc"
        end select 
            
    end subroutine init_surface

    function surface_charge(bc,psiSurf) result(sigma_surface_charge)
             
        real(dp), intent(in)  :: psiSurf
        character(len=2), intent(in) :: bc

        real(dp) :: sigma_surface_charge
        
        select case (bc)
            case ("qu") 
                sigma_surface_charge = surface_charge_quartz(psiSurf)
            case ("cl") 
                sigma_surface_charge = surface_charge_clay(psiSurf)
            case ("ca") 
                sigma_surface_charge = surface_charge_calcite(psiSurf)
            case ("ta")
                sigma_surface_charge = surface_charge_taurine(psiSurf)
            case ("cc")
                sigma_surface_charge = sigmaSurf
            case ("pp")
                sigma_surface_charge = surface_charge_pp(psisurf)
            case ("pd")
                sigma_surface_charge = surface_charge_pp_dynamic(psiSurf)
            case ("pc")
                sigma_surface_charge = surface_charge_pp_dynamic_cond(psiSurf)       
            case default
                print*,"bc does not match qu, cl, ca, ta, cc, pp, pd, or pc"    
                sigma_surface_charge = 0.0_dp
        end select


    end function surface_charge
 
    subroutine init_surface_quartz()
   
        use mathconst 
        use physconst, only : Na
        use parameters,  only : vsol,delta,lb

        implicit none

        integer :: i
    
        cap=1.0e-18_dp ! in F/nm^2 1F= 1C/V
 
        pKS(1)=   8.1_dp  !  8.1_dp !  >SOH <=> >SO- + H+ ! see Luetzenkirchen book chapter 14
        pKS(2)=  -4.4_dp  ! -4.4_dp !  >SOH2+ <=> >SOH + H+
        pKS(3)=  -1.5_dp  ! >SONa <=> >SO- + Na+  
        pKS(4)=  -5.20_dp ! >SOCa+ <=> >SO- + Ca2+  Ber. Bunsenges. Phys. Chem., 98, pp 1062-1067, 1994 
        pKS(5)=  -1.5_dp  ! >SOH2Cl <=>>SOH2+ Cl- 

        do i=1,5
            KS(i)  = 10.0_dp**(-pKS(i))       ! experimental equilibruim constant surface acid
            K0S(i) = (KS(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
        enddo

        ! charges surface states
        qS(1)=-1.0_dp ! SO^-
        qS(2)=0.0_dp  ! SOH
        qS(3)=1.0_dp  ! SOH_2^+
        qS(4)=0.0_dp  ! SONa
        qS(5)=1.0_dp  ! SOCa
        qS(6)=0.0_dp  ! SOH2Cl
    
        ! sites density
        sigmaSurf = 8.0_dp
        sigmaSurf = sigmaSurf * (4.0_dp*pi*lb)*delta ! dimensionless surface charge     

    end subroutine init_surface_quartz

    subroutine init_surface_calcite()
        
        use mathconst
        use physconst, only : Na
        use parameters,  only : vsol,delta,lb

        implicit none

        integer :: i

    !    cap=1.0d-18 ! in F/nm^2 1F= 1C/V                                                                                                                              
        !   1=CO3-, 2=CO3H, 3=CO3Ca+, 4=CaO^-, 5=CaOH, 6= CaOH2^+ 
        
        pKS(1)=   4.9_dp   ! >CO3H         <=> >CO3^- + H+ 
        pKS(2)=   2.8_dp   ! >CO3H + Ca^2+ <=> >CO3Ca^+ + H+ 
        pKS(3)=   12.2_dp  ! >CaOH2+       <=> >CaOH + H+     
        pKS(4)=   17.0_dp  ! >CaOH         <=> >CaO- + H+                                                                                                                     

        do i=1,4
           KS(i)  = 10.0_dp**(-pKS(i))       ! experimental equilibruim constant surface acid                                                                           
           K0S(i) = (KS(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant                                                                                       
        enddo
        ! i=2 is different
        K0S(2)=KS(2)

        ! charges surface states 
        qS(1)=-1.0_dp ! >CO3^-
        qS(2)=0.0_dp  ! >CO3H                                                
        qS(3)=1.0_dp  ! >CO3Ca^+

        qS(4)=-1.0_dp ! >CaO^-                              
        qS(5)=0.0_dp  ! >CaOH 
        qS(6)=1.0_dp  ! >CaOH2^+

        ! site density
        sigmaSurf = 5.0_dp
        sigmaSurf = sigmaSurf * (4.0_dp*pi*lb)*delta ! dimensionless surface charge     

    end subroutine init_surface_calcite

   
    subroutine init_surface_clay()

        use mathconst 
        use physconst, only : Na
        use parameters,  only : vsol,delta,lb

        integer :: i
        
        ! see Luetzenkirchen book chapter 7 pages 204-206                                     

        cap=0.90e-18_dp ! in F/nm^2 1F= 1C/V  

        pKS(1)=   10.0_dp  ! 10.0  >SOH^0.5   <=> >SO^-0.5 + H+ 
        pKS(2)=   5.3_dp   !  5.3  >SOH_2^0.5 <=> >SOH_1.5^0 + 1/2H+ 
        pKS(3)=   0.20_dp  !       >SONa^0.5  <=> >SO^-0.5 + Na+ 
        pKS(4)=  -5.20_dp  !       >SOCa^1.5  <=> >SO^-0.5 + Ca2+             
        pKS(5)=  -0.40_dp !        >SOH_2Cl^-0.5 <=> >SOH_2^0.5+Cl-         
       
        do i=1,5
           KS(i)  = 10.0_dp**(-pKS(i))       ! experimental equilibruim constant surface acid  
           K0S(i) = (KS(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
        enddo
        ! i=2
        !K0S(2) = KS(2)*((vsol*Na/1.0d24)**0.5) ! intrinstic equilibruim constant                                                  

        ! charges surface states                                                                                                 
        qS(1)=-0.5_dp  ! SOH^-0.5                                                                            
        qS(2)= 0.5_dp  ! SOH_2^0.5   
        qS(3)= 0.0_dp  ! SOH_1.5^0  
        qS(4)= 0.5_dp  ! SOHNa^0.5                                  
        qS(5)= 1.5_dp  ! SOHCa^1.5
        qS(6)=-0.5_dp  ! SOH_2Cl^-0.5                                                                                                        

        ! site density 
        sigmaSurf = 8.5_dp
        sigmaSurf = sigmaSurf * (4.0_dp*pi*lb)*delta ! dimensionless surface charge     

    end subroutine init_surface_clay

    subroutine init_surface_taurine()

        use mathconst
        use physconst, only : Na
        use parameters,  only : vsol,delta,lb

        implicit none      

        integer :: i

        pKTa(1)=   -2.0_dp    !  >SOH <=> >SO^- + H+ 
        pKTa(2)=   -0.42_dp   !  >SONa <=> >SO^- + Na+
        pKTa(3)=   -0.72243_dp  !  >SOCa^+ <=> >SO^- +Ca2+
       
        do i=1,3
            KTa(i)  = 10.0_dp**(-pKTa(i))       ! experimental equilibruim constant surface acid  
            K0Ta(i) = (KTa(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
        enddo

        ! charges surface states                                                                                                 
        qTa(1)=-1.0_dp ! SO^-                                                                            
        qTa(2)=0.0_dp  ! SOH  
        qTa(3)=0.0_dp  ! SONa 
        qTa(4)=1.0_dp  ! SOHCa^+                                  
                                                                                                            
        sigmaSurf = sigmaSurf * 4.0_dp*pi*lb*delta ! dimensionless surface charge     
       
        ! site density
        !if(side==RIGHT) sigmaSurfR = sigmaSurfR * 4.0_dp*pi*lb*delta ! dimensionless surface charge     
        !if(side==LEFT)  sigmaSurfL = sigmaSurfL * 4.0_dp*pi*lb*delta ! dimensionless surface charge     

    end subroutine init_surface_taurine


    
    subroutine init_surface_pp_dynamic(sigmaSurf)

        real(dp) ,intent(inout) :: sigmaSurf 
        call init_surface_pp(sigmaSurf)

    end subroutine init_surface_pp_dynamic


    subroutine init_surface_pp(sigmaSurf)
   
        use mathconst 
        use physconst, only : Na
        use parameters,  only : vsol,delta,lb

        real(dp) ,intent(inout) :: sigmaSurf 

        integer :: i
        real(dp) :: pKSa,pKSb

        pKS(1)=   4.6_dp  !  >SPOHCOOH <=> >SPOHCOO- + H+ ! carboxlic group of pp ligand
        pKS(2)=   5.4_dp  !  >SPOHCOOH <=> >SPOCOOH- + H+ ! second phosphonate state of pp ligand
        pKSa  =   6.9_dp  !  >SPOCOOH- <=> >SPOCOO2- + H+ ! carboxylic group of pp ligand 
        pKSb  =   7.8_dp  !  >SPOHCOO- <=> >SPOCOO2- + H+ ! both 
        pKS(3)=   pKS(1)+pKSb !  >SPOHCOOH <=> >SPOCOO2- + H+ 
        

        
        do i=1,4
            KS(i)  = 10.0_dp**(-pKS(i))       ! experimental equilibruim constant surface acid
            K0S(i) = (KS(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
        enddo

        ! charges surface states
        qS(AHBH)= 0.0_dp    !>SPOHCOOH
        qS(AHB) =-1.0_dp   !>SPOHCOO-
        qS(ABH) =-1.0_dp   !>SPOCOOH-
        qS(AB)  =-2.0_dp   !>SPOCOO2-
        qS(AH2BH)= 0.0_dp
       
    
        ! sites density
        sigmaSurf = sigmaSurf * (4.0_dp*pi*lb)*delta ! dimensionless surface charge     

    end subroutine init_surface_pp

    subroutine init_surface_pp_dynamic_cond(sigmaSurf)
    
        use mathconst 
        use physconst, only : Na
        use parameters,  only : vsol,delta,lb

        real(dp) ,intent(inout) :: sigmaSurf 

        integer :: i
        real(dp) :: pKSa,pKSb

        pKS(1)=   4.6_dp  !  >SPOHCOOH <=> >SPOHCOO- + H+ ! carboxlic group of pp ligand
        pKS(2)=   5.4_dp  !  >SPOHCOOH <=> >SPOCOOH- + H+ ! second phosphonate state of pp ligand
        pKSa  =   6.9_dp  !  >SPOCOOH- <=> >SPOCOO2- + H+ ! carboxylic group of pp ligand 
        pKSb  =   7.8_dp  !  >SPOHCOO- <=> >SPOCOO2- + H+ ! both 
        pKS(3)=   pKS(1)+pKSb !  >SPOHCOOH <=> >SPOCOO2- + H+ 
        

        
        do i=1,4
            KS(i)  = 10.0_dp**(-pKS(i))       ! experimental equilibruim constant surface acid
            K0S(i) = (KS(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
        enddo

        ! charges surface states
        qS(AHBH)  =  0.0_dp    !>SPOHCOOH
        qS(AHB)   = -1.0_dp    !>SPOHCOO-
        qS(ABH)   = -1.0_dp    !>SPOCOOH-
        qS(AB)    = -2.0_dp    !>SPOCOO2-
        qS(AH2BH) =  0.0_dp    ! notbound 
        qS(Su)    =  1.0_dp    !>S+ 
        qS(SuOH)  =  0.0_dp    !>SOH  
        qS(SuCl)  =  0.0_dp    !>SCl
        qS(SuNO3) =  0.0_dp    !>SNO3 
    
        ! sites density
        sigmaSurf = sigmaSurf * (4.0_dp*pi*lb)*delta ! dimensionless surface charge     


    end subroutine init_surface_pp_dynamic_cond



    subroutine init_surface_constcharge()

        use parameters,  only : delta,lb

        sigmaSurf = sigmaSurf * 4.0_dp*pi*lb*delta ! dimensionless surface charge    

    end subroutine init_surface_constcharge

    function surface_charge_quartz(psiS) result(surface_charge)
 
        use physconst
        use mathconst
        use parameters
    
        real(dp), intent(in) :: psiS
        
        real(dp) :: surface_charge

        ! .. local variables
    
        real(dp) :: xS(6)
        real(dp) :: A,avfdis
        integer :: i
    
        xS(1)= ((xbulk%Hplus/xbulk%sol)/K0S(1))*dexp(-psiS) ! SOH/SO-
        xS(2)= ((xbulk%Hplus/xbulk%sol)/K0S(2))*dexp(-psiS) ! SOH2+/SOH
        xS(3)= (((xbulk%Na/vNa)/(xbulk%sol**(vNa)))/K0S(3))*dexp(-psiS) ! SONa/SO-
        xS(4)= (((xbulk%Ca/vCa)/(xbulk%sol**(vCa)))/K0S(4))*dexp(-2.0_dp*psiS) ! SOCa+/SO-
        xS(5)= (((xbulk%Cl/vCl)/(xbulk%sol**(vCl)))/K0S(5))*dexp(psiS) ! SOH2Cl/SOH2+

        A = xS(1)*(1.0_dp+xS(2))+xS(3)+xS(4)+xS(1)*xS(2)*xS(5)
        fdisS(1)  = 1.0_dp/(1.0_dp +A) ! SO-
        fdisS(2)  = fdisS(1)*xS(1) ! SOH
        fdisS(3)  = fdisS(2)*xS(2) ! SOH2+
        fdisS(4)  = fdisS(1)*xS(3) ! SONa
        fdisS(5)  = fdisS(1)*xS(4) ! SOCa+
        fdisS(6)  = fdisS(1)*xS(5)*xS(2)*xS(1)  ! SOH2Cl  
    
        avfdis=0.0_dp
        do i=1,6       
            avfdis=avfdis +qS(i)*fdisS(i)
        enddo
    
        surface_charge=sigmaSurf*avfdis
    
    end function surface_charge_quartz


    function surface_charge_calcite(psiS) result(surface_charge)

        use physconst
        use mathconst
        use parameters

        real(dp), intent(in) :: psiS
        real(dp) :: surface_charge

        ! .. local variables                                                                                                                                          
        real(dp) :: xS(6)
        real(dp) :: A,avfdis
        integer :: i

        xS(1)= ((xbulk%Hplus/xbulk%sol)/K0S(1))*dexp(-psiS) ! CO3H/CO3- =fdisS(2)/fdisS(1)

        xS(2)= (xbulk%sol/xbulk%Hplus)*((xbulk%Ca*vsol/vCa)/(xbulk%sol**(vCa/vsol)))*K0S(2)*dexp(-psiS)   ! CO3Ca+/CO3H = fdisS(3)/fdisS(2)             
        xS(3)= ((xbulk%Hplus/xbulk%sol)/K0S(4))*dexp(-psiS) ! CaOH/CaO- =fdisS(5)/fdisS(4)   
        xS(4)= ((xbulk%Hplus/xbulk%sol)/K0S(3))*dexp(-psiS) ! CaOH2^+/CaOH = fdisS(6)/fdisS(5)

        A = xS(1)*(1.0_dp+xS(2))
        fdisS(1)  = 1.0_dp/(1.0_dp +A) ! CO3^- 
        fdisS(2)  = fdisS(1)*xS(1)   ! CO3H 
        fdisS(3)  = fdisS(2)*xS(2)   ! CO3Ca+ 
    
        A=xS(3)*(1.0_dp+xS(4))
        fdisS(4)  = 1.0_dp/(1.0_dp +A) ! CaO^- 
        fdisS(5)  = fdisS(4)*xS(3)   ! CaOH  
        fdisS(6)  = fdisS(5)*xS(4)   ! CaOH2^+                                                                                                              
        avfdis=0.0_dp
        do i=1,6
            avfdis=avfdis +qS(i)*fdisS(i)
        enddo

        surface_charge=sigmaSurf*avfdis

    end function surface_charge_calcite

  
    function surface_charge_clay(psiS) result(surface_charge)

        use physconst
        use mathconst
        use parameters

        real(dp), intent(in) :: psiS
        real(dp) :: surface_charge

        ! .. local variables                                                                                                  

        real(dp) :: xS(6)
        real(dp) :: A,avfdis
        integer :: i

        xS(1)= ((xbulk%Hplus/xbulk%sol)/K0S(1))*dexp(-psiS) ! SOH_2^0.5/SOH^-0.5                       
        xS(2)= (((xbulk%Hplus/xbulk%sol)**0.5_dp)/K0S(2))*dexp(-0.5_dp*psiS) ! SOH_1.5^0/SOH_2^0.5   
        xS(3)= (((xbulk%Na/vNa)/(xbulk%sol**(vNa)))/K0S(3))*dexp(-psiS) ! SOHNa^0.5/SOH^-0.5
        xS(4)= (((xbulk%Ca/vCa)/(xbulk%sol**(vCa)))/K0S(4))*dexp(-2.0_dp*psiS) ! SOHCa^1.5/SOH^-0.5
        xS(5)= (((xbulk%Cl/vCl)/(xbulk%sol**(vCl)))/K0S(5))*dexp(psiS) ! SOH2Cl^-0.5/SOH_2^0.5

        xS(2)=0.0 ! make zero  

        A = xS(1)*(1.0_dp+xS(2))+xS(3)+xS(4)+xS(1)*xS(5)
        fdisS(1)  = 1.0_dp/(1.0_dp +A) ! SOH^-0.5                                                                                 
        fdisS(2)  = fdisS(1)*xS(1) ! SOH_2^0.5 
        fdisS(3)  = fdisS(2)*xS(2) ! SOH_1.5^0                                                                                  
        fdisS(4)  = fdisS(1)*xS(3) ! SOHNa^0.5 
        fdisS(5)  = fdisS(1)*xS(4) ! SOHCa^1.5            
        fdisS(6)  = fdisS(1)*xS(1)*xS(5)  ! SOH_2Cl^-0.5                                                

        avfdis=0.0_dp
        do i=1,6
            avfdis=avfdis +qS(i)*fdisS(i)
        enddo

        surface_charge=sigmaSurf*avfdis

    end function surface_charge_clay


    function surface_charge_taurine(psiS) result(surface_charge)

        use physconst
        use mathconst
        use parameters

        real(dp), intent(in) :: psiS
        real(dp) :: surface_charge
    
        ! .. local variables                                                                                                  

        real(dp) :: xS(3)
        real(dp) :: A,avfdis
        integer :: i


        xS(1)= ((xbulk%Hplus/xbulk%sol)/K0Ta(1))*dexp(-psiS) ! SOH/SO^-                       
        xS(2)= (((xbulk%Na/vNa)/(xbulk%sol**(vNa)))/K0Ta(2))*dexp(-psiS) ! SONa/SO^-
        xS(3)= (((xbulk%Ca/vCa)/(xbulk%sol**(vCa)))/K0Ta(3))*dexp(-2.0_dp*psiS) ! SOCa^+/SO^-
     
        A = xS(1)+xS(2)+xS(3)
     
    
        fdisTaL(1)  = 1.0_dp/(1.0_dp + A) ! SO^-
        fdisTaL(2)  = fdisTal(1)*xS(1) ! SOH                                                                                 
        fdisTal(3)  = fdisTaL(1)*xS(2) ! SONa 
        fdisTaL(4)  = fdisTaL(1)*xS(3) ! SOCa^+                                                                                 

        avfdis=0.0_dp
        do i=1,4
            avfdis=avfdis +qTa(i)*fdisTaL(i)   
        enddo
    
        surface_charge=sigmaSurf*avfdis
    

    end function surface_charge_taurine


    function surface_charge_pp(psiS) result(surface_charge)

        use physconst
        use mathconst
        use parameters

        real(dp), intent(in) :: psiS
        real(dp) :: surface_charge
    
        ! .. local variables                                                                                                  

        real(dp) :: xS(3)
        real(dp) :: A,avfdis
        integer :: i

        ! this is a guess for the equation need to be done yet
        xS(1)= K0S(1)*(xbulk%sol/xbulk%Hplus)*dexp(psiS)        ! >SAHB-/>SAHBH                       
        xS(2)= K0S(2)*(xbulk%sol/xbulk%Hplus)*dexp(psiS)        ! >SABH-/>SAHBH 
        xS(3)= K0S(3)*((xbulk%sol/xbulk%Hplus)**2)*dexp(2.0_dp*psiS) ! >SAB2-/>SAHBH 
     
        A = xS(1)+xS(2)+xS(3)
     
        fdisS(AHBH) = 1.0_dp/(1.0_dp + A)    ! >SAHBH
        fdisS(AHB)  = fdisS(AHBH)*xS(1)      ! >SAHB-                                                                                 
        fdisS(ABH)  = fdisS(AHBH)*xS(2)      ! >SABH- 
        fdisS(AB)   = fdisS(AHBH)*xS(3)      ! >SAB2-                                                                                 
        fdisS(AH2BH) = 0.0_dp


        avfdis=0.0_dp
        do i=1,5
            avfdis=avfdis +qS(i)*fdisS(i)   
        enddo
        
        surface_charge=sigmaSurf*avfdis
        
    
    end function surface_charge_pp


    function surface_charge_pp_dynamic(psiS) result(surface_charge)

        use physconst
        use mathconst
        use parameters, only : deltaG0ads,expmu,vpp,zpp,xbulk

        real(dp), intent(in) :: psiS
        real(dp) :: surface_charge
    
        ! .. local variables                                                                                                  

        real(dp) :: xS(5),K0ads
        real(dp) :: sum_xS,avfdis
        integer :: t

        K0ads=exp(-deltaG0ads)

       ! exmpu%i := [exp(-beta(mu0_i-mu_i))v_i/v_w]exp(- beta pibulk v_i) 


        do t=1,5
            xS(t) = K0ads*exp(-qS(t)*psiS)*expmu%pp(t)/(vpp(t)*xbulk%sol**vpp(t))
        enddo    
    
        sum_xS = xS(AHBH)+xS(AHB)+xS(ABH)+xs(AB) ! do not include AB2BH assumed to not be adsorbed

        fdisS(AHBH)  = xS(AHBH)/sum_xS     ! >SAHBH
        fdisS(AHB)   = xS(AHB)/sum_xS      ! >SAHB-                                                                                 
        fdisS(ABH)   = xS(ABH)/sum_xS      ! >SABH- 
        fdisS(AB)    = xS(AB)/sum_xS       ! >SAB2-                                                                                 
        fdisS(AH2BH) = 0.0_dp

        fdisR= 1.0_dp/(1.0_dp+sum_xS)

        avfdis=0.0_dp
        do t=1,5
            avfdis=avfdis +qS(t)*fdisS(t)   
        enddo
        
        surface_charge=sigmaSurf*(1.0_dp-fdisR)*avfdis

        ! print*,"fdisS(AH2BH) =",fdisS(AH2BH)  
        ! print*,"fdisS(AHBH)  =",fdisS(AHBH)
        ! print*,"fdisS(AHB)   =",fdisS(AHB)
        ! print*,"fdisS(ABH)   =",fdisS(ABH)
        ! print*,"fdisS(AB)    =",fdisS(AB)
        ! print*,"fdisR        =",fdisR
        ! print*,"surface_charge        =",surface_charge

    end function surface_charge_pp_dynamic

    function surface_charge_pp_dynamic_cond(psiS) result(surface_charge)

        use physconst
        use mathconst
        use parameters, only : deltaG0ads,deltaG0adsSuOH,deltaG0adsSuCl,deltaG0adsSuNO3
        use parameters, only : expmu,vpp,zpp,xbulk,vCl,vNO3

        real(dp), intent(in) :: psiS
        real(dp) :: surface_charge
    
        ! .. local variables                                                                                                  

        real(dp) :: xS(9),K0ads,K0adsSuOH,K0adsSuCl,K0adsSuNO3
        real(dp) :: sum_xS,avfdis, sum_Z, sum_W, sumfdis
        integer :: t

        
        K0ads=exp(-deltaG0ads)
        if(deltaG0adsSuOH>10.0_dp) then 
            K0adsSuOH=0.0_dp
        else
            K0adsSuOH=exp(-deltaG0adsSuOH)
        endif
        if(deltaG0adsSuCl>10.0_dp) then
            K0adsSuCl=0.0_dp
        else    
            K0adsSuCl=exp(-deltaG0adsSuCl)
        endif   

        if(deltaG0adsSuNO3>10.0_dp) then
            K0adsSuNO3=0.0_dp
        else    
            K0adsSuNO3=exp(-deltaG0adsSuNO3)
        endif    
        
       ! exmpu%i := [exp(-beta(mu0_i-mu_i))v_i/v_w]exp(- beta pibulk v_i) 


        do t=1,5
            xS(t) = K0ads*exp(-(qS(t)-qS(Su))*psiS)*expmu%pp(t)/(vpp(t)*xbulk%sol**vpp(t))
        enddo    
        xS(SuOH)  = K0adsSuOH *exp(-(qS(SuOH) -qS(Su))*psiS)*expmu%OHmin/(xbulk%sol)
        xS(SuCl)  = K0adsSuCl *exp(-(qS(SuCl) -qS(Su))*psiS)*expmu%Cl/(vCl*xbulk%sol**vCl)
        xS(SuNO3) = K0adsSuNO3*exp(-(qS(SuNO3)-qS(Su))*psiS)*expmu%NO3/(vNO3*xbulk%sol**vNO3)

        sum_xS = xS(AHBH)+xS(AHB)+xS(ABH)+xs(AB) +xs(SuOH)+xs(SuCl) +xs(SuNO3) ! do not include AB2BH assumed to not be adsorbed

        ! fdisR(t) =sigma(t)/sigma0

        fdisS(Su)     = 1.0_dp/(1.0_dp+sum_xS)  ! >S+
        fdisR        = fdisS(Su)

        fdisS(AHBH)  = xS(AHBH)*fdisR     ! >SAHBH  
        fdisS(AHB)   = xS(AHB)*fdisR      ! >SAHB-                                                                                 
        fdisS(ABH)   = xS(ABH)*fdisR      ! >SABH- 
        fdisS(AB)    = xS(AB)*fdisR       ! >SAB2-                                                                                 
        fdisS(AH2BH) = 0.0_dp
        fdisS(SuOH)  = xS(SuOH)*fdisR      ! >SOH
        fdisS(SuCl)  = xS(SuCl)*fdisR      ! >SCl
        fdisS(SuNO3) = xS(SuNO3)*fdisR     ! >SNO3


        avfdis=0.0_dp
        sumfdis=0.0_dp
        do t=1,9
            avfdis=avfdis +qS(t)*fdisS(t)
            sumfdis=sumfdis+fdisS(t)   
        enddo
      !  print*,"sumfdis=", sumfdis

        surface_charge=sigmaSurf*avfdis

        ! fraction of adsorbed ligand found in state t
        sum_Z = xS(AHBH)+xS(AHB)+xS(ABH)+xs(AB) 
       

        gdisS(AHBH)  = xS(AHBH)/sum_Z     ! >SAHBH  
        gdisS(AHB)   = xS(AHB)/sum_Z      ! >SAHB-                                                                                 
        gdisS(ABH)   = xS(ABH)/sum_Z      ! >SABH- 
        gdisS(AB)    = xS(AB)/sum_Z       ! >SAB2-                                                                                 
        gdisS(AH2BH) = 0.0_dp
        
        ! fraction of ROH, RCl or R+ of non-ligand adsorbed
        sum_W = xS(SuOH)+xS(SuCl)+xS(SuNO3)+1.0_dp

        gdisS(SuOH)   = xS(SuOH)/sum_W      ! >SOH
        gdisS(SuCl)   = xS(SuCl)/sum_W      ! >SCl
        gdisS(SuNO3)  = xS(SuNO3)/sum_W      ! >SCl
        gdisS(Su)     = xS(Su)/sum_W        ! >S
        

        ! print*,"fdisS(AH2BH) =",fdisS(AH2BH)  
        ! print*,"fdisS(AHBH)  =",fdisS(AHBH)
        ! print*,"fdisS(AHB)   =",fdisS(AHB)
        ! print*,"fdisS(ABH)   =",fdisS(ABH)
        ! print*,"fdisS(AB)    =",fdisS(AB)
        ! print*,"fdisR        =",fdisR
        ! print*,"fdisS(SuOH)  =",fdisS(SuOH)
        ! print*,"fdisS(SuCl)  =",fdisS(SuCl)
        

        ! print*,"surface_charge        =",surface_charge

    end function surface_charge_pp_dynamic_cond

end module surface
