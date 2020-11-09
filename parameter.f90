module parameters

    use physconst
    use mathconst
    use random
    use volume
    use molecules
    use loopvar

    implicit none

    !  .. list of parameters

    type(moleclist) :: xbulk,expmu
    type(bornmoleclist) :: bornrad,bornbulk   
    !  .. volume 
    real(dp) :: vsol               ! volume of solvent  in nm^3       

    !  .. volume ions   
    real(dp) :: vNa                ! volume Na+ ion in units of vsol
    real(dp) :: vK                 ! volume K+ ion in units of vsol
    real(dp) :: vTB                ! volume TBA+ ion in units of vsol
    real(dp) :: vTM                ! volume TMA+ ion in units of vsol
    real(dp) :: vCl                ! volume Cl- ion in units of vsol
    real(dp) :: vNO3               ! volume NO3- ion in units of vsol      
    real(dp) :: vCa                ! volume Ca2+ positive divalent ion in units of vsol
    real(dp) :: vHplus             ! volume H+
    real(dp) :: vOHmin             ! volume OH- 
    real(dp) :: vRb                ! volume Rb+ ion in units of vsol
    real(dp) :: vIm                ! volume Im+ ion in units of vsol
    real(dp) :: vNaCl              ! volume ion pair NaCl 
    real(dp) :: vKCl               ! volume ion pair KCl
   
     ! .. volume carboxylic acid 
    real(dp) :: vAA(6)  

    ! .. volume ligand
    real(dp) :: vpp(5)             ! volume ligand 5 protonation states    
    real(dp) :: deltavpp(4)         

    !  .. radii
  
    real(dp) :: RNa
    real(dp) :: RK
    real(dp) :: RTB
    real(dp) :: RTM 
    real(dp) :: RCl
    real(dp) :: RBr 
    real(dp) :: RCa
    real(dp) :: RNO3
    real(dp) :: RRb
    real(dp) :: RIm


    ! .. segment length 
  
    ! .. valence charge 

    integer :: zNa                 ! valence charge Na+ ion 
    integer :: zK                  ! valence charge K+ ion 
    integer :: zCa                 ! valence charge Ca2+ ion 
    integer :: zCl                 ! valence charge Cl- ion 
    integer :: zTB               
    integer :: zTM              
    integer :: zNO3
    integer :: zRb
    integer :: zIm

    integer :: zpa

    integer :: zpp(5)              ! valence protonantion states        

    real(dp) :: Temp               ! temperature in K
    real(dp) :: dielectW           ! dielectric constant of water 
    real(dp) :: dielectP           ! dielectric constant of hydrocarbons/PA
    character(len=15) :: dielect_env ! selects dielectric function 

    real(dp) :: lb                 ! Bjerrum length	   
    real(dp) :: constqW            ! constant in Poisson eq dielectric constant of water 
    real(dp) :: constqE            ! electrostatic pre-factor in pdf 
   
    integer :: itmax               ! maximum number of iterations
    real(dp) :: error              ! error imposed accuaracy
    real(dp) :: fnorm              ! L2 norm of residual vector function fcn  
    integer :: infile              ! infile=1 read input files infile!=1 no input files 
    integer :: iter                ! counts number of iterations
  
    character(len=8) :: method           ! method="kinsol" or "zspow"  
    character(len=3) :: verboseflag      ! select input flag 

   
    !  .. equibrium constant
    
    real(dp) :: pKw                 ! water equilibruim constant pKw= -log[Kw] ,Kw=[H+][OH-] 
    real(dp) :: K0ionNa             ! intrinsic equilibruim constant
    real(dp) :: KionNa              ! experimemtal equilibruim constant 
    real(dp) :: pKionNa             ! experimental equilibruim constant pKion= -log[Kion]	 
    real(dp) :: K0ionK              ! intrinsic equilibruim constant
    real(dp) :: KionK               ! experimemtal equilibruim constant 
    real(dp) :: pKionK              ! experimental equilibruim constant pKion= -log[Kion]	 
    
    real(dp) :: deltaG0ads          ! adsorption energy of ligand
    real(dp) :: deltaG0adsSuOH
    real(dp) :: deltaG0adsSuCl
    real(dp) :: deltaG0adsSuNO3


    real(dp) :: K0pp(5)             ! intrinsic equilibruim constant ligand acid base equilbria  
    real(dp) :: pKpp(5)   
  
    !     .. bulk volume fractions

    real(dp), target :: cNaCl      ! concentration of NaCl in bulk in mol/liter
    real(dp) :: cKCl               ! concentration of KCl in bulk in mol/liter
    real(dp), target :: cRbCl      ! concentration of RbCl in bulk in mol/liter
    real(dp), target :: cImCl      ! concentration of ImCl in bulk in mol/liter
    real(dp) :: cCaCl2             ! concentration of CaCl2 in bulk in mol/liter
    real(dp) :: cTBCl              ! concentration of TBCl in  bulk in mol/liter
    real(dp) :: cTMNO3             ! concentration of TMNO3 in  bulk in mol/liter
    real(dp) :: cHplus             ! concentration of H+ in bulk in mol/liter
    real(dp) :: cOHmin             ! concentration of OH- in bulk in mol/liter
    real(dp) :: pHbulk             ! pH of bulk pH = -log([H+])
    real(dp) :: pOHbulk            ! p0H of bulk p0H = -log([0H-])
    real(dp), target :: cpp        ! concentration ligand in bulk in mol/liter  
    real(dp) :: rhoqbulk           ! total charge in bulk

    type (looplist), target :: pH

    logical :: isBulkHCl           ! if true adjustment of pH with HCl if false HNO3 : for ligand 
    logical :: isBulkRbOH          ! if true adjustment of pH with RbOH if false NaOH : for pa   
    
    !  .. pafiber varialbes

    real(dp) :: K0A,KA,pKA          

    real(dp) :: totalEpa
    real(dp) :: radiuspacore
    real(dp) :: totalcharge        ! equal to qres !!
    logical  :: isChargeRegularization, isCabinding
    real(dp) :: avfdispa
    real(dp) :: avfdisA(6)

    ! output varaible for charge_pa_ratio_freeRb
    integer  :: maxpalayer      ! location in layer of maximum of rhoEpa
    integer  :: maxdeltaRblayer ! maximum of layer integrated out using epsdeltaxRb tolerance 
    real(dp) :: epsdeltaxRb     ! tolerance = 0.0005_dp  
    integer  :: numlDs          ! number of Deybe length to integrate out
    real(dp) :: ratio_free_Rb,ratio_free_Rb_Debye
   

    real(dp) :: epsIm ! van der Waals interaction Im 
    real(dp) :: chiIm ! Flory-Huggins parameter between PA and Im

    !  constant for acrylic acid 
    real(dp) :: K0AA(5),pKaAA(5)
    real(dp) :: deltavA(5)
    real(dp) :: pKdRb 

contains

    ! determine total number of non linear equations

    subroutine set_size_neq()

        use globals
        use volume, only : nr

        implicit none

        integer(8) :: neq_bc

        neq_bc=0 
        if(bcflag/="cc") neq_bc=neq_bc+1

        select case (sysflag)
            case ("electnopoly") 
                neq = 2 * nr  + neq_bc  
            case ("electligand") 
                neq = 2 * nr  + neq_bc  
            case ("pafiber") 
                neq = 2 * nr  + neq_bc 
            case ("pafiberIm") 
                neq = 2 * nr  + neq_bc   
            case ("pafiberborn") 
                neq = 4 * nr  + neq_bc  
            case ("pafiberbornscf") 
                neq = 5 * nr  + neq_bc  
            case ("pafibervarelec") 
                neq = 2 * nr  + neq_bc                 
            case ("bulk water") 
                neq = 5 
            case ("bulk ligand") 
                neq = 6 
            case default
                print*,"set_size_neq: wrong value sysflag:  ",sysflag
                stop
        end select  

        neqint=int(neq)
         
    end subroutine set_size_neq

         

    ! computes the  Debye lenght for a given Bjerum lenght (lB in nm)  and 
    ! Ioinic Stenght (IS in M=mol/l 


    function DebyeLength(lB, IS) result(lD)

        use mathconst
        use physconst

        real(dp), intent(in) :: lB, IS
        real(dp) :: lD
      
        lD=1.0_dp/sqrt(8.0_dp*pi*lB*1.0e-9_dp*Na*IS*1.0e3_dp)
        lD=lD/1.0e-9_dp   
        
    end function

    ! computes the  Bjerrum lenght for a given temperature T 

    function BjerrumLength(T)result(lb)

        use mathconst
        use physconst
        
        real(dp) :: lb
        real(dp), intent(in) :: T ! temperature      

        lb=(elemcharge**2)/(4.0_dp*pi*dielectW*dielect0*kBoltzmann*T) ! bjerrum length in water=solvent in m
        lb=lb/1.0e-9_dp              ! bjerrum length in water in nm
    
    end function BjerrumLength
        
    !     purpose: initialize all constants parameter 
    !     pre: first read_inputfile has to be called   
    
    subroutine init_constants()

        use globals
        use volume
        use random
        use physconst
        
        real(dp) :: v3pp             ! local volume
        real(dp) :: lD
        
        !  .. initializations of variables
 
        pi=acos(-1.0_dp)            ! pi = arccos(-1)
        itmax=2000                  ! maximum number of iterations
        nr=nsize                    ! size of lattice in z-direction 
        
        !     .. charges
        zNa   = 1                   ! valence positive charged ion
        zK    = 1                   ! valence positive charged ion
        zCa   = 2                   ! valence divalent positive charged ion
        zCl   = -1                  ! valence negative charged ion
        zTB   = 1 
        zTM   = 1
        zNO3  = -1
        zRb   = 1
        zIm   = 1

        zpa   = -1
        
        zpp(AH2BH) = 0              ! charged states ligand
        zpp(AHBH)  = -1
        zpp(AHB)   = -2
        zpp(ABH)   = -2
        zpp(AB)    = -3

        !     .. ionic radii
        !     .. https://www.chemguide.co.uk/atoms/properties/atradius.html and http://abulafia.mt.ic.ac.uk/shannon/ptable.php
        
        RNa = 0.102_dp              ! radius of Na+ in nm
        RK  = 0.138_dp              ! radius of K+ in nm
        RBr = 0.196_dp              ! radius of Rb+ i nm 
        RCl = 0.181_dp              ! radius of Cl- in nm
            
        
        RCa = 0.106_dp              ! radius of Ca2+ in nm
        RRb = 0.152_dp              ! radius of Rb+ in nm 

        !RIm = 0.50_dp               ! radius of Imadazol ion in nm   
        RTB = 0.50_dp               ! radius of TBA+ in nm  
        RTM = 0.50_dp               ! radius of TMA+ in nm values from Wang, Nap et al in Jacs 133:2192, 2011
        RNO3= 0.30_dp               ! radius of NO3- in nm values form Kieland Jacs 59:1675, 1937
        

        !     .. volume

        if(sysflag/="neutral") then 
            vsol = 0.030_dp              ! volume water solvent molecule in (nm)^3
        elseif(sysflag=="neutral") then 
            vsol = 0.218_dp             ! volume hexane Mw=86.18 g/mol and rho=0.6548 g/ml  
        else 
            print*,"Error in call to init_constants subroutine"   
            print*,"Wrong system flag"
            stop
        endif   

        vNa  = ((4.0_dp/3.0_dp)*pi*(RNa)**3)/vsol 
        vK   = ((4.0_dp/3.0_dp)*pi*(RK)**3)/vsol 
        vCl  = ((4.0_dp/3.0_dp)*pi*(RCl)**3)/vsol 
        vCa  = ((4.0_dp/3.0_dp)*pi*(RCa)**3)/vsol 
        vRb  = ((4.0_dp/3.0_dp)*pi*(RRb)**3)/vsol
        ! vIm  = 0.09190_dp/vsol Im= C3H4N2       
        ! .. volume Im: based molecular weight  and density of v= M/(rho Na)  

        vIm  = 0.2648_dp/vsol ! Im=1-Ethyl-3-methyll imidazolium = EMIM                                                                                                       
        ! .. volume Im: based molecular weight  and density of EMIM.Cl vIM=vEMIM=vEMIMCl -vCl

        vTB  = ((4.0_dp/3.0_dp)*pi*(RTB)**3)/vsol
        vTM  = ((4.0_dp/3.0_dp)*pi*(RTM)**3)/vsol
        vNO3  = ((4.0_dp/3.0_dp)*pi*(RNO3)**3)/vsol

        vHplus = 1.0_dp
        vOHmin = 1.0_dp 

        vNaCl= (vNa+vCl)          ! contact ion pair
        vKCl = (vK+vCl)           ! contact ion pair
         
        v3pp = 0.0792_dp/vsol 

        vpp(AH2BH) = v3pp
        vpp(AHBH)  = v3pp
        vpp(AHB)   = v3pp
        vpp(ABH)   = v3pp
        vpp(AB)    = v3pp

        deltavpp(1)=vpp(AHBH)+vHplus-vpp(AH2BH) 
        deltavpp(2)=vpp(ABH)+vHplus-vpp(AHBH)    
        deltavpp(3)=vpp(AHB)+vHplus-vpp(AHBH)    
        deltavpp(4)=vpp(AB)+vHplus-vpp(AHB)     

    
        ! .. chemical equilbrium constants of ppp

        pKpp(1) =  2.26_dp        ! POH2COOH <=> POHCOOH- + H+ 
        pKpp(2) =  4.6_dp         ! POHCOOH- <=> POHCOO2- + H+ 
        pKpp(3) =  5.4_dp         ! POHCOOH- <=> POCOOH2- + H+ !
        pKpp(4) =  6.9_dp         ! POCOOH2- <=> POCOO3- + H+ 
        pKpp(5) =  7.8_dp         ! POHCOO2- <=> POCOO3- + H+ !
        pKw     = 14.0_dp         ! water equilibruim constant

        bornrad%Na  = RNa
        bornrad%Cl  = RCl
        bornrad%K   = RK             
        bornrad%Ca  = RCa
        bornrad%Hplus = radiussphere(vsol)
        bornrad%OHmin = radiussphere(vsol)
        bornrad%Rb = RRb
        bornrad%Im = radiussphere(vIm) 


        ! .. other physical variables

        Temp = 298.0_dp               ! temperature in Kelvin
        dielectW = 78.54_dp           ! dielectric constant water
        dielectP =  2.0_dp
        lb=BjerrumLength(Temp)        ! bjerrum length in water in nm
        seed  = 435672                ! seed for random number generator
        constqW = delta*delta*4.0_dp*pi*lb/vsol ! multiplicative constant Poisson Eq. 
        constqE = 1.0_dp /( 8.0_dp *constqW)      ! factor in PDF    


        call set_pa_properties()


        !  scaling of Van der Waals of Imidazolium
        if(sysflag=="pafiberIm".or.sysflag=="pafiberborn".or.sysflag=="pafibervarelec".or.&
            sysflag=="pafiberbornscf") then 
            epsIm= epsIm *((vIm*vsol)**2/vsol) 
            chiIm=chiIm*vIm
        else
            epsIm=0.0_dp
            chiIm=0.0_dp
        endif        

    end subroutine init_constants
   
    function volumesphere(radius)result(volume)

        real(dp), intent(in) :: radius
        real(dp) :: volume

        volume=(4.0_dp/3.0_dp)*pi*(radius**3)
    
    end function
      

    function radiussphere(volume)result(radius)

        real(dp), intent(in) :: volume
        real(dp) :: radius

        radius=(volume*3.0_dp/(4.0_dp*pi))**(1.0_dp/3.0_dp)
    
    end function


    !     purpose: initialize expmu needed by fcn 
    !     pre: first read_inputfile has to be called

    subroutine init_expmu_elect()
 
        use globals
        use physconst
        use dielectric_const, only : born
        
        implicit none 
        
        !     .. local variable
        
        real(dp),  dimension(:), allocatable :: x         ! volume fraction solvent iteration vector 
        real(dp),  dimension(:), allocatable :: xguess  
        real(dp) :: xNaClsalt, xKClsalt, xCaCl2salt,xRbClsalt, xImClsalt      ! volume fraction of divalent salt in bulk
        character(len=15) :: sysflag_old

        
        allocate(x(5))
        allocate(xguess(5))
        
        !     .. initializations of input dependent variables, electrostatic part 
        
        pHbulk=pH%val ! transfer pH value 

        cHplus = (10.0_dp)**(-pHbulk) ! concentration H+ in bulk
        pOHbulk = pKw -pHbulk       
        cOHmin  = (10.0_dp)**(-pOHbulk) ! concentration OH- in bulk
        
        xbulk%Hplus = (cHplus*Na/(1.0e24_dp))*(vsol) ! volume fraction H+ in bulk vH+=vsol
        xbulk%OHmin = (cOHmin*Na/(1.0e24_dp))*(vsol) ! volume fraction OH- in bulk vOH-=vsol
        

        ! NaCl in solution 
        xNaClsalt = (cNaCl*Na/(1.0d24))*((vNa+vCl)*vsol) ! volume fraction NaCl salt in mol/l
        xbulk%Na=xNaClsalt*vNa/(vNa+vCl)  
        xbulk%Cl=xNaClsalt*vCl/(vNa+vCl)
        
        ! RbCl in solution 
        xRbClsalt = (cRbCl*Na/(1.0e24_dp))*((vRb+vCl)*vsol)
        xbulk%Rb = xRbClsalt*vRb/(vRb+vCl)  
        xbulk%Cl = xbulk%Cl+xRbClsalt*vCl/(vRb+vCl)  

        print*,"pHbulk=",pHbulk
        if(pHbulk<=7) then      ! pH<= 7 
            xbulk%Cl=xbulk%Cl +(xbulk%Hplus -xbulk%OHmin)*vCl  ! NaCl+ HCl
        else     
            print*,"isbulkRbOH=",isbulkRbOH                         ! pH >7
            if(isbulkRbOH) then         
                xbulk%Rb = xbulk%Rb+ (xbulk%OHmin -xbulk%Hplus)*vRb ! RbCl +RbOH  
            else
                xbulk%Na=xbulk%Na +(xbulk%OHmin -xbulk%Hplus)*vNa ! NaCl+ NaOH    
            endif
        endif    

        ! ImCl in solution 
        xImClsalt = (cImCl*Na/(1.0e24_dp))*((vIm+vCl)*vsol)
        xbulk%Im = xImClsalt*vIm/(vIm+vCl)  
        xbulk%Cl = xbulk%Cl+xImClsalt*vCl/(vIm+vCl)  

        ! KCl in solution 
        xKClsalt = (cKCl*Na/(1.0e24_dp))*((vK+vCl)*vsol) ! volume fraction KCl salt in mol/l
        xbulk%K = xKClsalt*vK/(vK+vCl)  
        xbulk%Cl = xbulk%Cl+xKClsalt*vCl/(vK+vCl)  
        ! KCl in solution 
        xCaCl2salt = (cCaCl2*Na/(1.0e24_dp))*((vCa+2.0_dp*vCl)*vsol) ! volume fraction CaCl2 in mol/l
        xbulk%Ca=xCaCl2salt*vCa/(vCa+2.0_dp*vCl)
        xbulk%Cl=xbulk%Cl+ xCaCl2salt*2.0_dp*vCl/(vCa+2.0_dp*vCl)
        
        xbulk%NaCl=0.0_dp    ! no ion pairing
        xbulk%KCl=0.0_dp     ! no ion pairing
        
        xbulk%sol=1.0_dp -xbulk%Hplus -xbulk%OHmin -xbulk%Cl -xbulk%Na -xbulk%K-xbulk%NaCl-xbulk%KCl-xbulk%Ca &
                -xbulk%Rb-xbulk%Im  
        
        !     .. ionpairing NaCl and KCl
        !     .. only ionpairing if Kion neq 0 ion pairing !
        !     .. intrinstic equilibruim constant acid        
        !     .. Kion unit 1/M= liter per mol !

        if(sysflag=="pafiber".or.sysflag=="pafiberIm".or.sysflag=="pafiberborn".or.sysflag=="pafibervarelec"&
            .or.sysflag=="pafiberbornscf") then   ! no ion pairing
            KionNa = 0.0_dp          
            KionK  = 0.0_dp
            Ka     = 10.0_dp**(-pKa) ! experimental equilibruim constant acid 
            K0a    = (Ka*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
        endif  

        K0ionK  = KionK /(vsol*Na/1.0e24_dp) ! intrinstic equilibruim constant 
        K0ionNa = KionNa/(vsol*Na/1.0e24_dp) ! intrinstic equilibruim constant 
        
        if((KionNa.ne.0.0_dp).or.(KionK.ne.0.0_dp)) then  
            sysflag_old=sysflag 
            sysflag="bulk water"        ! set solver to fcnbulk
            call set_size_neq()         ! number of nonlinear equations
            
            x(1)=xbulk%Na
            x(2)=xbulk%Cl
            x(3)=xbulk%NaCl
            x(4)=xbulk%K
            x(5)=xbulk%KCl
            
            xguess(1)=x(1)
            xguess(2)=x(2)
            xguess(3)=x(3)
            xguess(4)=x(4)
            xguess(5)=x(5)
           
            call solver(x, xguess, error, fnorm) 
            
            !     .. return solution
            
            xbulk%Na  =x(1)
            xbulk%Cl  =x(2)
            xbulk%NaCl=x(3)
            xbulk%K   =x(4)
            xbulk%KCl =x(5)

            ! reset of flags
            iter=0
            sysflag=sysflag_old         ! switch solver back
            call set_size_neq()         ! number of non-linear  equation        
            
            xbulk%sol=1.0_dp-xbulk%Hplus-xbulk%OHmin - xbulk%Cl -xbulk%Na -xbulk%K-xbulk%NaCl-xbulk%KCl-xbulk%Ca 
            
        endif

        rhoqbulk = xbulk%Hplus -xbulk%OHmin +xbulk%Cl*zCl/vCl+xbulk%Na*zNa/vNa +xbulk%K*zK/vK+xbulk%Ca*zCa/vCa +&
            xbulk%Rb*zRb/vRb+ xbulk%Im*zIm/vIm



       
        ! pibulk = -log(xbulk%sol)  ! pressure (pi) of bulk
        ! exp(beta mu_i) = (rhobulk_i v_i) / exp(- beta pibulk v_i) 
        expmu%Na    = xbulk%Na   /(xbulk%sol**vNa) 
        expmu%K     = xbulk%K    /(xbulk%sol**vK)
        expmu%Ca    = xbulk%Ca   /(xbulk%sol**vCa) 
        expmu%Cl    = xbulk%Cl   /(xbulk%sol**vCl)
        expmu%Rb    = xbulk%Rb   /(xbulk%sol**vRb)

        expmu%Im    = xbulk%Im   /(xbulk%sol**vIm)
        expmu%NaCl  = xbulk%NaCl /(xbulk%sol**vNaCl)
        expmu%KCl   = xbulk%KCl  /(xbulk%sol**vKCl)

        expmu%Hplus = xbulk%Hplus/xbulk%sol ! vsol = vHplus 
        expmu%OHmin = xbulk%OHmin/xbulk%sol ! vsol = vOHmin 
          
        !     .. end init electrostatic part 
        
        if(sysflag=="pafiberIm".or.sysflag=="pafibervarelec") then 
           ! expmu%Im    = xbulk%Im/( exp(epsIm*(xbulk%Im/(vIm*vsol) )) * ( xbulk%sol**vIm))
            expmu%Im    = xbulk%Im/(xbulk%sol**vIm)
        endif  

        if(sysflag=="pafiberborn".or.sysflag=="pafiberbornscf") then

            bornbulk%AA   = born(lb,bornrad%AA,-1)
            bornbulk%AACa = born(lb,bornrad%AACa,1)
            
            bornbulk%Hplus = born(lb,bornrad%Hplus,1)
            bornbulk%Na    = born(lb,bornrad%Na,zNa)
            bornbulk%K     = born(lb,bornrad%K,zK)
            bornbulk%Ca    = born(lb,bornrad%Ca,zCa)
            bornbulk%Cl    = born(lb,bornrad%Cl,zCl)
            bornbulk%Rb    = born(lb,bornrad%Rb,zRb)
            bornbulk%OHmin = born(lb,bornrad%OHmin,-1)
            bornbulk%Im    = born(lb,bornrad%Im,zIm)

            expmu%Na    = (xbulk%Na   /(xbulk%sol**vNa))*exp(bornbulk%Na) 
            expmu%Cl    = (xbulk%Cl   /(xbulk%sol**vCl))*exp(bornbulk%Cl) 
            expmu%K     = (xbulk%K    /(xbulk%sol**vK) )*exp(bornbulk%K) 
            expmu%Ca    = (xbulk%Ca   /(xbulk%sol**vCa))*exp(bornbulk%Ca) 
            expmu%Rb    = (xbulk%Rb   /(xbulk%sol**vRb))*exp(bornbulk%Rb) 
            expmu%Hplus = (xbulk%Hplus/xbulk%sol) *      exp(bornbulk%Hplus)  
            expmu%OHmin = (xbulk%OHmin/xbulk%sol) *      exp(bornbulk%OHmin)  

            !expmu%Im    = xbulk%Im/( exp(epsIm*(xbulk%Im/(vIm*vsol) )) * ( xbulk%sol**vIm)) *exp(bornbulk%Im) 
          
            expmu%Im    = (xbulk%Im/(xbulk%sol**vIm)) *exp(bornbulk%Im) 
        endif  
    
        deallocate(x)
        deallocate(xguess)    
        
    end subroutine init_expmu_elect


    subroutine init_expmu_elect_qdot()
 
        use globals
        use physconst
        
        implicit none 
        
        !     .. local variable
        
        real(dp) :: xNaClsalt, xKClsalt, xCaCl2salt, xTBClsalt, xTMNO3salt            ! volume fraction of divalent salt in bulk


        !     .. initializations of input dependent variables, electrostatic part 
        
        pHbulk=pH%val ! transfer pH value 

        cHplus = (10.0_dp)**(-pHbulk) ! concentration H+ in bulk
        pOHbulk = pKw -pHbulk       
        cOHmin  = (10.0_dp)**(-pOHbulk) ! concentration OH- in bulk
        
        xbulk%Hplus = (cHplus*Na/(1.0e24_dp))*(vsol) ! volume fraction H+ in bulk vH+=vsol
        xbulk%OHmin = (cOHmin*Na/(1.0e24_dp))*(vsol) ! volume fraction OH- in bulk vOH-=vsol
        
        ! NaCl in solution 
        xNaClsalt = (cNaCl*Na/(1.0d24))*((vNa+vCl)*vsol) ! volume fraction NaCl salt in mol/l
        xbulk%Na=xNaClsalt*vNa/(vNa+vCl)
        xbulk%Cl=xNaClsalt*vCl/(vNa+vCl)

        ! adjust pH
        if(pHbulk.le.7) then      ! pH<= 7 add HCl 
            xbulk%Cl=xNaClsalt*vCl/(vNa+vCl) +(xbulk%Hplus -xbulk%OHmin)*vCl  ! NaCl+ HCl
        else                      ! pH >7
            xbulk%K= +(xbulk%OHmin -xbulk%Hplus)*vK ! NaCl+ KOH  
        endif
        
        ! KCl in solution 
        xKClsalt = (cKCl*Na/(1.0e24_dp))*((vK+vCl)*vsol) ! volume fraction KCl salt in mol/l
        xbulk%K  = xbulk%K  +  xKClsalt*vK/(vK+vCl)  
        xbulk%Cl = xbulk%Cl + xKClsalt*vCl/(vK+vCl)  
        
        ! CaCl2 in solution 
        xCaCl2salt = (cCaCl2*Na/(1.0e24_dp))*((vCa+2.0_dp*vCl)*vsol) ! volume fraction CaCl2 in mol/l
        xbulk%Ca=xCaCl2salt*vCa/(vCa+2.0_dp*vCl)
        xbulk%Cl=xbulk%Cl+ xCaCl2salt*2.0_dp*vCl/(vCa+2.0_dp*vCl)
        
        ! TBCl in solution TBCl=tetrabutyl ammonium chloride 
        xTBClsalt = (cTBCl*Na/(1.0e24_dp))*((vTB+vCl)*vsol) ! volume fraction KCl salt in mol/l
        xbulk%TB = xTBClsalt*vTB/(vTB+vCl)  
        xbulk%Cl = xbulk%Cl + xTBClsalt*vCl/(vTB+vCl)  

        ! TMNO3 in solution TMNO3=tetramethyl ammonium nitrate 
        xTMNO3salt = (cTMNO3*Na/(1.0e24_dp))*((vTM+vNO3)*vsol) ! volume fraction KCl salt in mol/l
        xbulk%TM = xTMNO3salt*vTM/(vTM+vNO3)  
        xbulk%NO3 =xTMNO3salt*vNO3/(vTM+vNO3)  



        ! no ionparing 

        xbulk%NaCl=0.0_dp  
        xbulk%KCl=0.0_dp   
        
        xbulk%sol=1.0_dp -xbulk%Hplus -xbulk%OHmin -xbulk%Cl -xbulk%Na -xbulk%K-xbulk%NaCl-xbulk%KCl-xbulk%Ca&
            -xbulk%TB -xbulk%TM - xbulk%NO3       
        
        rhoqbulk = xbulk%Hplus -xbulk%OHmin +xbulk%Cl*zCl/vCl +xbulk%Na*zNa/vNa +xbulk%K*zK/vK+&
            xbulk%TB*zTB/vTB + xbulk%Ca*zCa/vCa + xbulk%TM*zTM/vTM+xbulk%NO3*zNO3/vNO3

        !     .. if Kion neq 0 ion pairing !
        !     .. intrinstic equilibruim constant acid        
        !     .. unit 1/M= liter per mol !!!
        K0ionK  = KionK /(vsol*Na/1.0e24_dp) ! intrinstic equilibruim constant 
        K0ionNa = KionNa/(vsol*Na/1.0e24_dp) ! intrinstic equilibruim constant 
       
        if((KionNa.ne.0.0_dp).or.(KionK.ne.0.0_dp)) then  
           print*,"Input errror KionNa and KionK should be zero"
           stop
        endif

        
        ! pibulk = -log(xbulk%sol)  ! pressure (pi) of bulk
        ! exp(beta mu_i) = (rhobulk_i v_i) / exp(- beta pibulk v_i) 
        expmu%Na    = xbulk%Na   /(xbulk%sol**vNa) 
        expmu%K     = xbulk%K    /(xbulk%sol**vK)
        expmu%TB    = xbulk%TB   /(xbulk%sol**vTB)
        expmu%TM    = xbulk%TM   /(xbulk%sol**vTM)
        expmu%Ca    = xbulk%Ca   /(xbulk%sol**vCa) 
        expmu%Cl    = xbulk%Cl   /(xbulk%sol**vCl)
        expmu%NO3   = xbulk%NO3  /(xbulk%sol**vNO3)
        expmu%NaCl  = xbulk%NaCl /(xbulk%sol**vNaCl)
        expmu%KCl   = xbulk%KCl  /(xbulk%sol**vKCl)
        expmu%Hplus = xbulk%Hplus/xbulk%sol ! vsol = vHplus 
        expmu%OHmin = xbulk%OHmin/xbulk%sol ! vsol = vOHmin 

          
        !     .. end init electrostatic part 

        
    end subroutine init_expmu_elect_qdot


    ! .. computes the fraction of different charged states 
    ! .. in infinite dilution approximation ( ideal solution )
    
    function fpp_inf_dilution(xbulkHplus) result(fpp)

        use globals, only : AH2BH,AHBH,ABH,AHB,AB

        real(dp), intent(in) :: xbulkHplus
        real(dp) :: fpp(5)
    
        ! .. local variable
        real(dp) :: xA, xB, xBprime, xE, sumx, phisol

        phisol=1.0_dp

        !  .. equilibrium eq AH2BH <=> AHBH^- +H^+  A<=> B
        xA = K0pp(1)*(vpp(AHBH)/vpp(AH2BH))*(phisol)/xbulkHplus

        !  .. equilibrium eq AHBH^- <=> ABH^2- +H^+    B<=>C 
        xB = K0pp(3)*(vpp(ABH)/vpp(AHBH))*(phisol)/xbulkHplus
       
        !   .. equilibrium eq AHBH^- <=> AHB^2- +H^+    B<=>E 
        xBprime = K0pp(2)*(vpp(AHB)/vpp(AHBH))*(phisol)/xbulkHplus

        !   .. equilibrium eq AHB^2- <=> AB^3- +H^+    E<=>F
        xE = K0pp(5)*(vpp(AB)/vpp(AHB))*(phisol)/xbulkHplus

        sumx=xA + xA*xB + xA*xBprime + xA*xBprime*xE

        fpp(AH2BH) = 1.0_dp/(1.0_dp+sumx)
        fpp(AHBH)  = fpp(AH2BH) * xA
        fpp(ABH)   = fpp(AHBH)  * xB
        fpp(AHB)   = fpp(AHBH)  * xBprime 
        fpp(AB)    = fpp(AHB)   * xE

    end function


    subroutine init_expmu_elect_ligand()
 
        use globals
        use physconst
        
        implicit none 
        
        !     .. local variable
                
        real(dp),  dimension(:), allocatable :: x         ! volume fraction solvent iteration vector 
        real(dp),  dimension(:), allocatable :: xguess  
        integer :: i, t
        character(len=15) :: sysflag_old
        real(dp) :: Kpp(5), fppbulk(5)
        real(dp) :: xppbulk, rhoqppbulk, cppbulk, sumfpp
        real(dp) :: xNaClsalt, xKClsalt, xCaCl2salt, xTBClsalt ,xTMNO3salt            ! volume fraction of divalent salt in bulk

        allocate(x(6))
        allocate(xguess(6))

        !     .. initializations of input dependent variables, electrostatic part 
        
        pHbulk=pH%val ! transfer pH value 

        cHplus = (10.0_dp)**(-pHbulk) ! concentration H+ in bulk
        pOHbulk = pKw -pHbulk       
        cOHmin  = (10.0_dp)**(-pOHbulk) ! concentration OH- in bulk
        
        xbulk%Hplus = (cHplus*Na/(1.0e24_dp))*(vsol) ! volume fraction H+ in bulk vH+=vsol
        xbulk%OHmin = (cOHmin*Na/(1.0e24_dp))*(vsol) ! volume fraction OH- in bulk vOH-=vsol
        
        ! NaCl in solution 
        xNaClsalt = (cNaCl*Na/(1.0d24))*((vNa+vCl)*vsol) ! volume fraction NaCl salt in mol/l
        xbulk%Na=xNaClsalt*vNa/(vNa+vCl)
        xbulk%Cl=xNaClsalt*vCl/(vNa+vCl)
 
        ! KCl in solution 
        xKClsalt = (cKCl*Na/(1.0e24_dp))*((vK+vCl)*vsol) ! volume fraction KCl salt in mol/l
        xbulk%K = xbulk%K  +  xKClsalt*vK/(vK+vCl)  
        xbulk%Cl = xbulk%Cl + xKClsalt*vCl/(vK+vCl)  
        
        ! CaCl2 in solution 
        xCaCl2salt = (cCaCl2*Na/(1.0e24_dp))*((vCa+2.0_dp*vCl)*vsol) ! volume fraction CaCl2 in mol/l
        xbulk%Ca = xCaCl2salt*vCa/(vCa+2.0_dp*vCl)
        xbulk%Cl = xbulk%Cl+ xCaCl2salt*2.0_dp*vCl/(vCa+2.0_dp*vCl)
        
        ! TBCl in solution TBCl=tetrabutyl 
        xTBClsalt = (cTBCl*Na/(1.0e24_dp))*((vTB+vCl)*vsol) ! volume fraction KCl salt in mol/l
        xbulk%TB = xTBClsalt*vTB/(vTB+vCl)  
        xbulk%Cl = xbulk%Cl + xTBClsalt*vCl/(vTB+vCl)  


        ! TMNO3 in solution TMNO3=tetramethyl ammonium nitrate 
        xTMNO3salt = (cTMNO3*Na/(1.0e24_dp))*((vTM+vNO3)*vsol) ! volume fraction KCl salt in mol/l
        xbulk%TM  = xTMNO3salt*vTM/(vTM+vNO3)  
        xbulk%NO3 = xTMNO3salt*vNO3/(vTM+vNO3)  


        !  .. no ionpairs 
        xbulk%NaCl = 0.0_dp  
        xbulk%KCl = 0.0_dp   
    
        !  .. no ion pairing 
        KionNa = 0.0_dp
        KionK  = 0.0_dp
        K0ionK = KionK /(vsol*Na/1.0e24_dp) ! intrinstic equilibruim constant 
        K0ionNa = KionNa/(vsol*Na/1.0e24_dp) ! intrinstic equilibruim constant 
       
        !     .. intrinstic equilibruim constants      
        do i=1,5
            Kpp(i)  = 10.0_dp**(-pKpp(i)) ! experimental equilibruim constant acid 
            K0pp(i) = (Kpp(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
        enddo         
      
        ! solver non linear eq of fcnbulkligand 

        sysflag_old=sysflag 
        sysflag="bulk ligand"       ! set sysflag 
        call set_size_neq()         ! number of nonlinear equations
    
        ! .. initial guess
        if (size(fppbulk) /= 5) stop "Incorrect size of 'fppbulk'"
        fppbulk=fpp_inf_dilution(xbulk%Hplus)
        do t=1,4    
            x(t)= fppbulk(t)
        enddo

        if(isbulkHCl) then
            x(5)= xbulk%Cl
        else
            x(5)= xbulk%NO3
        endif    
        x(6)= xbulk%K
        do i=1,6
            xguess(i)=x(i)
        enddo
        ! .. end guess

        call solver(x, xguess, error, fnorm) 
        
        !     .. return solution
        sumfpp= 0.0_dp   
        do i=1,4
            fppbulk(i) = x(i)
            sumfpp=sumfpp +fppbulk(i)
        enddo
        fppbulk(5)=1.0_dp - sumfpp    
        
        cppbulk=cpp*(Na/1.0e24_dp) !  .. bulk concentration  cpp in mol/liter cppbulk in ligands/nm^3
        xppbulk=0.0_dp
        rhoqppbulk=0.0_dp
        do t=1,5
            xbulk%pp(t) = cppbulk*fppbulk(t)*vpp(t)*vsol
            xppbulk = xppbulk +xbulk%pp(t) ! total ligand volume fraction
            rhoqppbulk=rhoqppbulk +cppbulk*fppbulk(t)*zpp(t)
        enddo
        
        if(isbulkHCl) then
            xbulk%Cl =x(5)
        else
            xbulk%NO3=x(5)
        endif    
        xbulk%K  =x(6)

        xbulk%sol=1.0_dp-xbulk%Hplus-xbulk%OHmin - xbulk%Cl -xbulk%Na -xbulk%K-xbulk%TB-xbulk%Ca -xppbulk &
        -xbulk%TM-xbulk%NO3
        
        rhoqbulk = xbulk%Hplus-xbulk%OHmin +xbulk%Cl*zCl/vCl +xbulk%Na*zNa/vNa +xbulk%K*zK/vK+&
            xbulk%TB*zTB/vTB+xbulk%Ca*zCa/vCa+rhoqppbulk*vsol+xbulk%TM*zTM/vTM+xbulk%NO3*zNO3/vNO3
        

        ! reset of flags
        iter=0
        sysflag=sysflag_old         ! switch sysflag  back
        call set_size_neq()         ! number of non-linear  equation        

        ! .. make chemical potentials 

        !  pibulk = -log(xbulk%sol)  ! pressure (pi) of bulk
       
        ! expmu%i := (rhobulk_i v_i) / exp(- beta pibulk v_i) 
        ! exmpu%i := [exp(-beta(mu0_i-mu_i))v_i/v_w]exp(- beta pibulk v_i) 
       
        expmu%Na    = xbulk%Na   /(xbulk%sol**vNa) 
        expmu%K     = xbulk%K    /(xbulk%sol**vK)
        expmu%TB    = xbulk%TB   /(xbulk%sol**vTB)
        expmu%TM    = xbulk%TM   /(xbulk%sol**vTM)
        expmu%Ca    = xbulk%Ca   /(xbulk%sol**vCa) 
        expmu%Cl    = xbulk%Cl   /(xbulk%sol**vCl)
        expmu%NO3   = xbulk%NO3  /(xbulk%sol**vNO3)
        expmu%NaCl  = xbulk%NaCl /(xbulk%sol**vNaCl)
        expmu%KCl   = xbulk%KCl  /(xbulk%sol**vKCl)
        expmu%Hplus = xbulk%Hplus/xbulk%sol ! vsol = vHplus 
        expmu%OHmin = xbulk%OHmin/xbulk%sol ! vsol = vOHmin 
       
        do t=1,5
            expmu%pp(t)=xbulk%pp(t)/xbulk%sol**vpp(t)
        enddo
        

        deallocate(x)
        deallocate(xguess)
        
        
    end subroutine init_expmu_elect_ligand



    subroutine init_expmu

        use globals, only : sysflag, bcflag
      
        if(sysflag=="electnopoly") then
            if(bcflag=="pp") then 
                call init_expmu_elect_qdot()
            else
                call init_expmu_elect()
            endif    
        elseif(sysflag=="electligand") then
            call init_expmu_elect_ligand()   
        elseif(sysflag=="pafiber" ) then 
            call init_expmu_elect()
        elseif(sysflag=="pafiberIm" ) then 
            call init_expmu_elect()    
        elseif(sysflag=="pafiberborn" ) then 
            call init_expmu_elect()  
        elseif(sysflag=="pafibervarelec" ) then 
            call init_expmu_elect()    
        elseif(sysflag=="pafiberbornscf" ) then 
            call init_expmu_elect()          
        else
            print*,"Error in call to init_expmu subroutine"    
            print*,"Wrong value sysflag : ", sysflag
            stop        
        endif   

    end subroutine init_expmu


    ! dimensions pa_fiber
    ! numbers are place holder values !!!!!!!!
    ! need to be called before make_geometry 
    subroutine set_pa_properties

        integer :: i
        real(dp) :: KAA(5)
        real(dp) :: vA
 


        radiuspacore   = radius   

        ! set equilbrium constant for acrylic acid 
        pKaAA(1)=5.0_dp
        pKaAA(2)=-0.4_dp
        pKaAA(3)=1.0_dp
        pKaAA(4)=4.0_dp
        pKaAA(5)=pKdRb ! -0.6_dp

        do i=1,5
            KAA(i)=10.0_dp**(-pKaAA(i))  
            K0AA(i) = (KAA(i)*vsol)*(Na/1.0e24_dp)
        enddo
        K0AA(4) = (K0AA(4)*vsol)*(Na/1.0e24_dp)
       
        deltavA(1)=1.0_dp ! vA- + vH+ - vAH
        deltavA(2)=0.0_dp ! vA- + vNa+ - vANa+
        deltavA(3)=0.0_dp ! vA- + vCa2+ - vACa+
        deltavA(4)=0.0_dp ! 2vA- + vCa2+ -vA2Ca 
        deltavA(5)=0.0_dp ! vA- + vRb+ -vARb     
        
        vA=  0.07448_dp/vsol    ! size acrylic acid monomer !!!!
        vA=  0.038_dp/vsol      ! size COOH 
        vAA(1)= vA              ! vA-
        vAA(2)= vA              ! vAH
        vAA(3)= vA+vNa          ! vANa
        vAA(4)= vA+vCa          ! vACa
        vAA(5)= 2.0_dp*vA+vCa   ! vA2Ca     
        vAA(6)= vA+vRb          ! vARb 

        bornrad%AA   = radiussphere(vAA(1)*vsol)
        bornrad%AACa = radiussphere(vAA(4)*vsol) 
        

    end subroutine


 end module parameters
