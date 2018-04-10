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

    !  .. volume 
    real(dp) :: vsol               ! volume of solvent  in nm^3       
    !  .. volume monomers
    real(dp) :: vpolB(5)           ! volume of one polymer segment, vpol  in units of vsol
    real(dp) :: vpolA(5)           ! volume of one polymer segment, vpol  in units of vsol
    real(dp) :: vpolC              ! volume of one polymer segment hydrocarbon, vpol  in units of vsol
    real(dp) :: deltavA(4)
    real(dp) :: deltavB(4)
    !  .. volume ions   
    real(dp) :: vNa                ! volume positive ion in units of vsol
    real(dp) :: vK                 ! volume positive ion in units of vsol
    real(dp) :: vTB                ! volume positive ion in units of vsol
    real(dp) :: vCl                ! volume negative ion in units of vsol   
    real(dp) :: vCa                ! volume positive divalent ion in units of vsol
    real(dp) :: vNaCl
    real(dp) :: vKCl
    real(dp) :: vHplus
    real(dp) :: vOHmin 
    ! .. volume ligand
    real(dp) :: vpp(5)             ! volume ligand 5 protonation states    
    real(dp) :: deltavpp(4)         

    !  .. radii
  
    real(dp) :: RNa
    real(dp) :: RK
    real(dp) :: RTB 
    real(dp) :: RCl
    real(dp) :: RCa

    ! .. segment length 

    real(dp) :: lsegAB
    real(dp) :: lsegA              ! segment length of A polymer in nm
    real(dp) :: lsegB              ! segment length of B polymer in nm
    real(dp) :: lsegC              ! segment length of C polymer in nm
    real(dp) :: lsegCH2 
    real(dp) :: lsegPAA   
    real(dp) :: lsegPAMPS
    
    integer :: period              ! chain peridociy of repeat of A or B block 
  
    real(dp) :: VdWepsB            ! strenght VdW interaction in units of kT
    real(dp) :: VdWepsC            ! strenght VdW interaction in units of kT
    real(dp) :: chibulk            ! value of chibulk 
    integer :: numlayers
    integer :: VdWcutoff           ! cutoff VdW interaction in units of lseg 	
    integer :: VdWcutoffdelta      ! cutoff VdW interaction in units of delta
    integer :: layeroffset
  
    ! .. valence charge 

    integer :: zpolA(5)            ! valence charge polymer
    integer :: zpolB(5)            ! valence charge polymer
    integer :: zNa                 ! valence charge positive ion 
    integer :: zK                  ! valence charge positive ion 
    integer :: zCa                 ! valence charge divalent positive ion 
    integer :: zCl                 ! valence charge negative ion 
    integer :: zTB               
    integer :: zpp(5)              ! valence protonantion states        

    real(dp) :: Temp               ! temperature in K
    real(dp) :: dielectW           ! dielectric constant of water 
    real(dp) :: lb                 ! Bjerrum length	   
    real(dp) :: constqW            ! constant in Poisson eq dielectric constant of water 

    real(dp) :: sigmaAB            ! sigma AB polymer coated on surface
    real(dp) :: sigmaC             ! sigma C polymer coated on planar surface
  
    integer :: itmax               ! maximum number of iterations
    real(dp) :: error              ! error imposed accuaracy
    real(dp) :: fnorm              ! L2 norm of residual vector function fcn  
    integer :: infile              ! infile=1 read input files infile!=1 no input files 
    integer :: iter                ! counts number of iterations
  
    character(len=8) :: method           ! method="kinsol" or "zspow"  
    character(len=8) :: chainmethod      ! method of generating chains ="MC" or "FILE" 
    character(len=8) :: chaintype        ! type of chain: diblock,alt
    integer :: readinchains              ! nunmber of used/readin chains
    character(len=3) ::  verboseflag     ! select input flag 

    real(dp) :: heightAB           ! average height of layer
    real(dp) :: heightC            ! average height of layer 
    real(dp) :: qpolA              ! charge poly A of layer 
    real(dp) :: qpolB              ! charge poly B of layer 
    real(dp) :: qpol_tot           ! charge poly A+B of layer 
    real(dp) :: avfdisA(5)         ! average degree of dissociation 
    real(dp) :: avfdisB(5)         ! average degree of dissociation
  
    !  .. equibrium constant
  
    real(dp) :: K0A(4)              ! intrinsic equilibruim constant
    real(dp) :: KA(4)               ! experimemtal equilibruim constant 
    real(dp) :: pKA(4)              ! experimental equilibruim constant pKa= -log[Ka]
    real(dp) :: K0B(4)              ! intrinsic equilibruim constant
    real(dp) :: KB(4)               ! experimemtal equilibruim constant 
    real(dp) :: pKB(4)              ! experimental equilibruim constant pKa= -log[Ka]
    real(dp) :: pKw                 ! water equilibruim constant pKw= -log[Kw] ,Kw=[H+][OH-] 
    real(dp) :: K0ionNa             ! intrinsic equilibruim constant
    real(dp) :: KionNa              ! experimemtal equilibruim constant 
    real(dp) :: pKionNa             ! experimental equilibruim constant pKion= -log[Kion]	 
    real(dp) :: K0ionK              ! intrinsic equilibruim constant
    real(dp) :: KionK               ! experimemtal equilibruim constant 
    real(dp) :: pKionK              ! experimental equilibruim constant pKion= -log[Kion]	 
    
    real(dp) :: deltaGads           ! adsorption energy  
    real(dp) :: K0pp(5)             ! intrinsic equilibruim constant ligand acid base equilbria  
    real(dp) :: pKpp(5)   
  
    !     .. bulk volume fractions

    real(dp), target :: cNaCl      ! concentration of NaCl in bulk in mol/liter
    real(dp) :: cKCl               ! concentration of KCl in bulk in mol/liter
    real(dp) :: cCaCl2             ! concentration of CaCl2 in bulk in mol/liter
    real(dp) :: cTBCl              ! concentration of TBCl in  bulk in mol/liter 
    real(dp) :: cHplus             ! concentration of H+ in bulk in mol/liter
    real(dp) :: cOHmin             ! concentration of OH- in bulk in mol/liter
    real(dp) :: pHbulk             ! pH of bulk pH = -log([H+])
    real(dp) :: pOHbulk            ! p0H of bulk p0H = -log([0H-])
    real(dp), target :: cpp        ! concentration ligand in bulk in mol/liter  
  
    type (looplist), target :: pH
        
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
            case ("elect") 
                neq = 4 * nr + neq_bc
            case ("electdouble")  
                neq = 4 * nr
            case ("electnopoly") 
                neq = 2 * nr + neq_bc
            case ("electligand") 
                neq = 2 * nr  + neq_bc    
            case ("electHC") 
                neq = 5 * nr +neq_bc
            case ("neutral") 
                neq = 2 * nr
            case ("bulk water") 
                neq = 5 
            case ("bulk ligand") 
                neq = 6 
            case default
                print*,"Wrong value sysflag:  ",sysflag
                stop
        end select  
         
    end subroutine set_size_neq

    
    function BjerrumLenght(T)result(lb)

        use mathconst
        use physconst
        
        real(dp) :: lb
        real(dp), intent(in) :: T ! temperature      

        lb=(elemcharge**2)/(4.0_dp*pi*dielectW*dielect0*kBoltzmann*T) ! bjerrum length in water=solvent in m
        lb=lb/1.0e-9_dp              ! bjerrum length in water in nm
    
    end function BjerrumLenght
        
    !     purpose: initialize all constants parameter 
    !     pre: first read_inputfile has to be called   
    
    subroutine init_constants()

        use globals
        use volume
        use random
        use physconst
        
        implicit none      
        
        real(dp) :: vA,vB, vAA, vAMPS, v3pp
        
        !  .. initializations of variables
 
        pi=acos(-1.0_dp)          ! pi = arccos(-1)
        itmax=2000                ! maximum number of iterations
        nr=nsize                  ! size of lattice in z-direction 
        
        !     .. charges
        zNa   = 1                 ! valence positive charged ion
        zK    = 1                 ! valence positive charged ion
        zCa   = 2                 ! valence divalent positive charged ion
        zCl   =-1                 ! valence negative charged ion
        zTB   = 1 

        zpolA(1)=-1 ! A-
        zpolA(2)= 0 ! AH
        zpolA(3)= 0 ! ANa
        zpolA(4)= 1 ! ACa+
        zpolA(5)= 0 ! A2Ca
        
        zpolB(1)=-1 ! B-
        zpolB(2)= 0 ! BH
        zpolB(3)= 0 ! BNa
        zpolB(4)= 1 ! BCa+
        zpolB(5)= 0 ! B2Ca
        
        zpp(AH2BH) = 0       ! charged states ligand
        zpp(AHBH)  = -1
        zpp(AHB)   = -2
        zpp(ABH)   = -2
        zpp(AB)    = -3

        !     .. radii
        
        RNa = 0.102_dp             ! radius of Na+ in nm
        RK  = 0.138_dp             ! radius of K+ in nm
        RCl = 0.181_dp             ! radius of Cl- in nm
        RCa = 0.106_dp             ! radius of Ca2+ in nm
        RTB = 0.50_dp              ! radius of TBA+ in nm 

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
        vTB  = ((4.0_dp/3.0_dp)*pi*(RTB)**3)/vsol 
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

        !     .. volume polymer segments
        !     .. all volume scaled by vsol
        
        vAA  =  0.07448_dp/vsol ! volume based on VdW radii 
        vAMPS = 0.2134_dp/vsol
    
        vA = vAA
        vB = vAMPS

        vpolA(1)= vA              ! vA-
        vpolA(2)= vA              ! vAH
        vpolA(3)= vA+vNa          ! vANa
        vpolA(4)= vA+vCa          ! vACa
        vpolA(5)= 2.0_dp*vA+vCa   ! vA2Ca
        
        vpolB(1)= vB              ! vB-
        vpolB(2)= vB              ! vBH
        vpolB(3)= vB+vNa          ! vBNa
        vpolB(4)= vB+vCa          ! vBCa
        vpolB(5)= 2.0_dp*vB+vCa   ! vB2Ca
        
        deltavA(1)=vpolA(1)+1.0_dp-vpolA(2) ! vA-+vH+-vAH
        deltavA(2)=vpolA(1)+vNa-vpolA(3)    ! vA-+vNa+-vANa+
        deltavA(3)=vpolA(1)+vCa-vpolA(4)    ! vA- + vCa2+ -vACa+
        deltavA(4)=2.0_dp*vpolA(1)+vCa-vpolA(5) ! 2vA- + vCa2+ -vA2Ca
        
        deltavB(1)=vpolB(1)+1.0_dp-vpolB(2) ! vB-+vH+-vBH
        deltavB(2)=vpolB(1)+vNa-vpolB(3)    ! vB-+vNa+-vBNa+
        deltavB(3)=vpolB(1)+vCa-vpolB(4)    ! vB- +vCa2+ -vBCa+
        deltavB(4)=2.0_dp*vpolB(1)+vCa-vpolB(5) ! 2vB- + vCa2+ -vB2Ca+
        
        vpolC  = 0.0270_dp/vsol    ! volume CH2
       
        !  .. polymer segment lenght 
        lsegPAA   = 0.36287_dp     ! segment length in nm
        lsegPAMPS = 0.545_dp       ! segment length in nm
        lsegCH2   = 0.153_dp       ! segment length in nm od CH2 check  
        lsegAB = lsegPAMPS          
        lsegA  = lsegPAA            
        lsegB  = lsegPAMPS          
        lsegC  = lsegCH2            

        ! .. see also subroutine set_chain_properties 

        ! .. chemical equilbrium constants

        pKpp(1) =  2.26_dp        ! POH2COOH <=> POHCOOH- + H+ 
        pKpp(2) =  4.6_dp         ! POHCOOH- <=> POHCOO2- + H+ 
        pKpp(3) =  5.4_dp         ! POHCOOH- <=> POCOOH2- + H+ !
        pKpp(4) =  6.9_dp         ! POCOOH2- <=> POCOO3- + H+ 
        pKpp(5) =  7.8_dp         ! POHCOO2- <=> POCOO3- + H+ !
        pKw = 14.0_dp             ! water equilibruim constant

        ! .. other physical variables

        Temp=298.0_dp                 ! temperature in Kelvin
        dielectW=78.54_dp             ! dielectric constant water
        lb=BjerrumLenght(Temp)        ! bjerrum length in water in nm
        seed  = 435672                ! seed for random number generator
        constqW = delta*delta*4.0_dp*pi*lb/vsol ! multiplicative constant Poisson Eq. 
        
        !  .. initializations of input dependent variables 
        
        sigmaAB = sigmaAB * (1.0_dp/(delta))  ! dimensionless sigma no vpol*vsol !!!!!!!!!!!! 
        sigmaC   = sigmaC * (1.0_dp/(delta)) ! dimensionless sigma no vpol*vsol !!!!!!!!!!!!
        
        ! VdWepsC  = VdWepsC/(vpolC*vsol) ! VdW eps scaled 
        ! VdWepsB  = VdWepsB/(vpolB(3)*vsol) ! VdW eps scaled 
        
        ! .. make radius integer multiply of delta
        ! .. needed because VdW-coefficeint computed on grid 
        ! Íradius=delta*int(radius/delta)

        max_conforAB=cuantasAB
        max_conforC=cuantasC

    end subroutine init_constants
   
   
    !     purpose: initialize expmu needed by fcn 
    !     pre: first read_inputfile has to be called

    subroutine init_expmu_elect()
 
        use globals
        use physconst
        
        implicit none 
        
        !     .. local variable
        
        real(dp),  dimension(:), allocatable :: x         ! volume fraction solvent iteration vector 
        real(dp),  dimension(:), allocatable :: xguess  
        real(dp) :: xNaClsalt, xKClsalt, xCaCl2salt, xTBClsalt           ! volume fraction of divalent salt in bulk
        integer :: i
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
        
        xNaClsalt = (cNaCl*Na/(1.0d24))*((vNa+vCl)*vsol) ! volume fraction NaCl salt in mol/l
        
        if(pHbulk.le.7) then      ! pH<= 7
            xbulk%Na=xNaClsalt*vNa/(vNa+vCl)  
            xbulk%Cl=xNaClsalt*vCl/(vNa+vCl) +(xbulk%Hplus -xbulk%OHmin)*vCl  ! NaCl+ HCl
        else                      ! pH >7
            xbulk%Na=xNaClsalt*vNa/(vNa+vCl) +(xbulk%OHmin -xbulk%Hplus)*vNa ! NaCl+ NaOH  
            xbulk%Cl=xNaClsalt*vCl/(vNa+vCl)  
        endif
        
        xKClsalt = (cKCl*Na/(1.0e24_dp))*((vK+vCl)*vsol) ! volume fraction KCl salt in mol/l
        xbulk%K = xKClsalt*vK/(vK+vCl)  
        xbulk%Cl = xbulk%Cl+xKClsalt*vCl/(vK+vCl)  
        
        xCaCl2salt = (cCaCl2*Na/(1.0e24_dp))*((vCa+2.0_dp*vCl)*vsol) ! volume fraction CaCl2 in mol/l
        xbulk%Ca=xCaCl2salt*vCa/(vCa+2.0_dp*vCl)
        xbulk%Cl=xbulk%Cl+ xCaCl2salt*2.0_dp*vCl/(vCa+2.0_dp*vCl)
        
        xbulk%NaCl=0.0_dp    ! no ion pairing
        xbulk%KCl=0.0_dp     ! no ion pairing
        
        xbulk%sol=1.0_dp -xbulk%Hplus -xbulk%OHmin -xbulk%Cl -xbulk%Na -xbulk%K-xbulk%NaCl-xbulk%KCl-xbulk%Ca   
        
        !     .. if Kion neq 0 ion pairing !
        
        !     .. intrinstic equilibruim constant acid        
        !     Kion  = 0.246_dp ! unit 1/M= liter per mol !!!
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

        !     .. intrinstic equilibruim constants      
        do i=1,4
             Ka(i)  = 10.0_dp**(-pKa(i)) ! experimental equilibruim constant acid 
             Kb(i)  = 10.0_dp**(-pKb(i)) ! experimental equilibruim constant acid
             K0a(i) = (Ka(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
             K0b(i) = (Kb(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
        enddo
        !     .. rescale for i=4 2A- Ca <=> A2Ca
          
        K0a(4) = (Ka(4)*vsol)*(Na/1.0e24_dp)
        K0b(4) = (Kb(4)*vsol)*(Na/1.0e24_dp)
         

        ! pibulk = -log(xbulk%sol)  ! pressure (pi) of bulk
        ! exp(beta mu_i) = (rhobulk_i v_i) / exp(- beta pibulk v_i) 
        expmu%Na    = xbulk%Na   /(xbulk%sol**vNa) 
        expmu%K     = xbulk%K    /(xbulk%sol**vK)

        expmu%Ca    = xbulk%Ca   /(xbulk%sol**vCa) 
        expmu%Cl    = xbulk%Cl   /(xbulk%sol**vCl)
        expmu%NaCl  = xbulk%NaCl /(xbulk%sol**vNaCl)
        expmu%KCl   = xbulk%KCl  /(xbulk%sol**vKCl)
        expmu%Hplus = xbulk%Hplus/xbulk%sol ! vsol = vHplus 
        expmu%OHmin = xbulk%OHmin/xbulk%sol ! vsol = vOHmin 
          
        !     .. end init electrostatic part 
            
        VdWepsC  = VdWepsC/(vpolC*vsol) ! VdW eps scaled 
        VdWepsB  = VdWepsB/(vpolB(3)*vsol) ! VdW eps scaled 

        deallocate(x)
        deallocate(xguess)
        
    end subroutine init_expmu_elect


    subroutine init_expmu_elect_qdot()
 
        use globals
        use physconst
        
        implicit none 
        
        !     .. local variable
        
        integer :: i
        real(dp) :: xNaClsalt, xKClsalt, xCaCl2salt, xTBClsalt           ! volume fraction of divalent salt in bulk


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
        
        ! TBCl in solution TBCl=tetrabutyl 
        xTBClsalt = (cTBCl*Na/(1.0e24_dp))*((vTB+vCl)*vsol) ! volume fraction KCl salt in mol/l
        xbulk%TB = xTBClsalt*vTB/(vTB+vCl)  
        xbulk%Cl = xbulk%Cl + xTBClsalt*vCl/(vTB+vCl)  
        ! no ionparing 

        xbulk%NaCl=0.0_dp  
        xbulk%KCl=0.0_dp   
        
        xbulk%sol=1.0_dp -xbulk%Hplus -xbulk%OHmin -xbulk%Cl -xbulk%Na -xbulk%K-xbulk%NaCl-xbulk%KCl-xbulk%Ca-xbulk%TB   
        

        !     .. if Kion neq 0 ion pairing !
        !     .. intrinstic equilibruim constant acid        
        !     Kion  = 0.246_dp ! unit 1/M= liter per mol !!!
        K0ionK  = KionK /(vsol*Na/1.0e24_dp) ! intrinstic equilibruim constant 
        K0ionNa = KionNa/(vsol*Na/1.0e24_dp) ! intrinstic equilibruim constant 
       
        if((KionNa.ne.0.0_dp).or.(KionK.ne.0.0_dp)) then  
           print*,"Input errror KionNa and KionK should be zero"
           stop
        endif

        !     .. intrinstic equilibruim constants      
        do i=1,4
             Ka(i)  = 10.0_dp**(-pKa(i)) ! experimental equilibruim constant acid 
             Kb(i)  = 10.0_dp**(-pKb(i)) ! experimental equilibruim constant acid
             K0a(i) = (Ka(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
             K0b(i) = (Kb(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
        enddo
        !     .. rescale for i=4 2A- Ca <=> A2Ca
          
        K0a(4) = (Ka(4)*vsol)*(Na/1.0e24_dp)
        K0b(4) = (Kb(4)*vsol)*(Na/1.0e24_dp)
         

        ! pibulk = -log(xbulk%sol)  ! pressure (pi) of bulk
        ! exp(beta mu_i) = (rhobulk_i v_i) / exp(- beta pibulk v_i) 
        expmu%Na    = xbulk%Na   /(xbulk%sol**vNa) 
        expmu%K     = xbulk%K    /(xbulk%sol**vK)
        expmu%TB    = xbulk%TB   /(xbulk%sol**vTB)
        expmu%Ca    = xbulk%Ca   /(xbulk%sol**vCa) 
        expmu%Cl    = xbulk%Cl   /(xbulk%sol**vCl)
        expmu%NaCl  = xbulk%NaCl /(xbulk%sol**vNaCl)
        expmu%KCl   = xbulk%KCl  /(xbulk%sol**vKCl)
        expmu%Hplus = xbulk%Hplus/xbulk%sol ! vsol = vHplus 
        expmu%OHmin = xbulk%OHmin/xbulk%sol ! vsol = vOHmin 
          
        !     .. end init electrostatic part 
            
        VdWepsC  = VdWepsC/(vpolC*vsol) ! VdW eps scaled 
        VdWepsB  = VdWepsB/(vpolB(3)*vsol) ! VdW eps scaled 

        
    end subroutine init_expmu_elect_qdot


    ! .. computes the fraction of different charged states 
    ! .. in infinite dilution approximation ( ideal solution )
    
    function fpp_inf_dilution(xbulkHplus) result(fpp)

        use globals, only : AH2BH,AHBH,ABH,AHB,AB

        real(dp), intent(in) :: xbulkHplus
        real(dp) :: fpp(5)
    
        ! .. local variable

        integer :: i
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
        real(dp) :: xppbulk, rhoqppbulk, cppbulk, sumfpp , rhoqbulk
        real(dp) :: xNaClsalt, xKClsalt, xCaCl2salt, xTBClsalt           ! volume fraction of divalent salt in bulk

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

        !  .. no ionpairs 
        xbulk%NaCl = 0.0_dp  
        xbulk%KCl = 0.0_dp   
    
        !  .. no  ion pairing 
        KionNa = 0.0_dp
        KionK = 0.0_dp
        K0ionK = KionK /(vsol*Na/1.0e24_dp) ! intrinstic equilibruim constant 
        K0ionNa = KionNa/(vsol*Na/1.0e24_dp) ! intrinstic equilibruim constant 
       
        !     .. intrinstic equilibruim constants      
        do i=1,5
            Kpp(i)  = 10.0_dp**(-pKpp(i)) ! experimental equilibruim constant acid 
            K0pp(i) = (Kpp(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
        enddo         
      
        ! solver non linear eq of fcnbulkgligand 

        sysflag_old=sysflag 
        sysflag="bulk ligand"       ! set sysflag 
        call set_size_neq()         ! number of nonlinear equations
    
        ! .. initial guess
        if (size(fppbulk) /= 5) stop "Incorrect size of 'fppbulk'"
        fppbulk=fpp_inf_dilution(xbulk%Hplus)
        do t=1,4    
            x(t)= fppbulk(t)
        enddo
        x(5)= xbulk%Cl
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
        xbulk%Cl =x(5)
        xbulk%K  =x(6)

        xbulk%sol=1.0_dp-xbulk%Hplus-xbulk%OHmin - xbulk%Cl -xbulk%Na -xbulk%K-xbulk%TB-xbulk%Ca -xppbulk
        rhoqbulk= xbulk%Hplus-xbulk%OHmin +xbulk%Cl*zCl/vCl +xbulk%Na*zNa/vNa +xbulk%K*zK/vK+xbulk%TB*zTB/vTB+xbulk%Ca*zCa/vCa
        rhoqbulk=rhoqbulk+rhoqppbulk*vsol

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
        expmu%Ca    = xbulk%Ca   /(xbulk%sol**vCa) 
        expmu%Cl    = xbulk%Cl   /(xbulk%sol**vCl)
        expmu%NaCl  = xbulk%NaCl /(xbulk%sol**vNaCl)
        expmu%KCl   = xbulk%KCl  /(xbulk%sol**vKCl)
        expmu%Hplus = xbulk%Hplus/xbulk%sol ! vsol = vHplus 
        expmu%OHmin = xbulk%OHmin/xbulk%sol ! vsol = vOHmin 
       
        do t=1,5
            expmu%pp(t)=xbulk%pp(t)/xbulk%sol**vpp(t)
        enddo
          
        !     .. end init electrostatic part 
            
        VdWepsC  = VdWepsC/(vpolC*vsol) ! VdW eps scaled 
        VdWepsB  = VdWepsB/(vpolB(3)*vsol) ! VdW eps scaled 

        deallocate(x)
        deallocate(xguess)
        
        
    end subroutine init_expmu_elect_ligand



    subroutine init_expmu_neutral

        implicit none
          
        xbulk%sol=1.0_dp ! only solvent 

        VdWepsC  = VdWepsC/(vpolC*vsol) ! VdW eps scaled 
        VdWepsB  = VdWepsB/(vpolB(3)*vsol) ! VdW eps scaled

    end subroutine init_expmu_neutral

    subroutine init_expmu

        use globals, only : sysflag, bcflag
        implicit none


        if(sysflag=="elect") then 
            call init_expmu_elect()
        elseif(sysflag=="electdouble") then 
            call init_expmu_elect()
        elseif(sysflag=="electnopoly") then
            if(bcflag=="pp") then 
                call init_expmu_elect_qdot()
            else
                call init_expmu_elect()
            endif    
        elseif(sysflag=="neutral") then
            call init_expmu_neutral()
        elseif(sysflag=="electligand") then
            call init_expmu_elect_ligand()   
        else
            print*,"Error in call to init_expmu subroutine"    
            print*,"Wrong value sysflag : ", sysflag
            stop        
        endif   

    end subroutine init_expmu

 end module parameters
