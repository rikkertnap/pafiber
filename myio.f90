
module myio

    use precision_definition

    implicit none

    ! return error values

    integer, parameter ::  myio_err_sysflag   = 1
    integer, parameter ::  myio_err_runflag   = 2
    integer, parameter ::  myio_err_geometry  = 3
    integer, parameter ::  myio_err_method    = 4
    integer, parameter ::  myio_err_chaintype = 5
    integer, parameter ::  myio_err_domain    = 6
    integer, parameter ::  myio_err_inputfile = 7
    integer, parameter ::  myio_err_input     = 8
    integer, parameter ::  myio_err_bcflag    = 9 

    ! unit number 
    integer :: un_sys,un_xpolAB,un_xpolC,un_xsol,un_xNa,un_xCl,un_xK,un_xCa,un_xNaCl,un_xKCl
    integer :: un_xOHmin,un_xHplus,un_fdisA,un_fdisB,un_psi,un_charge, un_xpair, un_rhopolAB, un_xTB, un_xpp, un_cpp
   
    ! format specifiers 
    character(len=80), parameter  :: fmt = "(A8,I1,A5,ES25.16)"
    character(len=80), parameter  :: fmt2reals = "(2ES25.16E3)"   
    character(len=80), parameter  :: fmt3reals = "(3ES25.16E3)"  
    character(len=80), parameter  :: fmt4reals = "(4ES25.16E3)" 
    character(len=80), parameter  :: fmt5reals = "(5ES25.16E3)"
    character(len=80), parameter  :: fmt6reals = "(6ES25.16E3)" 
  
    integer :: num_concen      ! number of concentrations     
    real(dp), dimension(:), allocatable, target :: concen_array   ! concentrations     

    private 

    public :: read_inputfile
    public :: output
    public :: num_concen, concen_array, set_value_concen
     
contains

subroutine read_inputfile(info)

    use globals
    use parameters
    use surface 
    use myutils, only : newunit

    integer, intent(out), optional :: info

    ! .. local arguments

    integer :: info_sys, info_bc, info_run, info_geo, info_meth, info_chaintype, info_combi
    character(len=8) :: fname
    integer :: ios,un_input  ! un = unit number    
    
    !     .. reading in of variables from file
    write(fname,'(A8)')'input.in'
    open(unit=newunit(un_input),file=fname,iostat=ios,status='old')
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        if (present(info)) info = myio_err_inputfile
        return
    endif

    read(un_input,*)method
    read(un_input,*)sysflag
    read(un_input,*)bcflag
    read(un_input,*)runflag
    read(un_input,*)chainmethod
    read(un_input,*)chaintype
    read(un_input,*)sigmaAB
    read(un_input,*)sigmaC
    read(un_input,*)error             
    read(un_input,*)infile              ! guess  1==yes
    read(un_input,*)radius
    read(un_input,*)pH%val
    read(un_input,*)pH%min
    read(un_input,*)pH%max
    read(un_input,*)pH%stepsize
    read(un_input,*)pH%delta
    read(un_input,*)KionNa
    read(un_input,*)KionK
    read(un_input,*)sigmaSurf
    read(un_input,*)cNaCl
    read(un_input,*)cKCl
    read(un_input,*)cCaCl2
    if(bcflag=="pp".or.bcflag=="pd") then 
        read(un_input,*)cTBCl            !   TB=tertraButyl
        if(sysflag=="electligand") then 
            read(un_input,*)deltaG0ads
            read(un_input,*)cpp
        endif
    endif   
    read(un_input,*)pKa(1)           !   AH   <=> A- + H+ 
    read(un_input,*)pKa(2)           !   ANa  <=> A- + Na+  
    read(un_input,*)pKa(3)           !   ACa+ <=> A- + Ca2+ 
    read(un_input,*)pKa(4)           !   A2Ca <=> 2A- + Ca2+
    read(un_input,*)pKb(1)           !   BH   <=> B- + H+ 
    read(un_input,*)pKb(2)           !   BNa  <=> B- + Na+ 
    read(un_input,*)pKb(3)           !   BCa+ <=> B- + Ca2+   
    read(un_input,*)pKb(4)           !   B2Ca <=> 2B- + Ca2+   
    read(un_input,*)period
    read(un_input,*)nsize
    if(runflag=="rangenr")then
        read(un_input,*)nrmax            ! max distance
        read(un_input,*)nrmin            ! min distance
        read(un_input,*)nrstep           ! step distance  
    endif
    read(un_input,*)nsegAB
    read(un_input,*)cuantasAB
    read(un_input,*)nsegC
    read(un_input,*)cuantasC
    read(un_input,*)VdWepsC
    read(un_input,*)VdWepsB    
    read(un_input,*)VdWcutoff
    read(un_input,*)verboseflag  
    read(un_input,*)geometry
    read(un_input,*)delta   

    close(un_input)
          
    ! .. check error flag

    call check_value_sysflag(sysflag,info_sys) 
    if (info_sys == myio_err_sysflag) then
        if (present(info)) info = info_sys
        return
    endif

    call check_value_runflag(runflag,info_sys) 
    if (info_sys == myio_err_runflag) then
        if (present(info)) info = info_run
        return
    endif

    call check_value_geometry(geometry,info_geo)
    if (info_geo == myio_err_geometry) then
        print*,"Wrong geometry"   
        if (present(info)) info = info_geo
        return
    endif

    call check_value_bcflag(bcflag,info_bc) 
    if (info_bc == myio_err_bcflag) then
        if (present(info)) info = info_bc
        return
    endif

    call check_value_method(method,info_meth)
    if (info_meth == myio_err_method) then
        if (present(info)) info = info_meth
        return
    endif

    call check_value_chaintype(chaintype,info_chaintype)
    if (info_chaintype == myio_err_chaintype) then
        if (present(info)) info = info_chaintype
        return
    endif

    if (present(info)) info = 0

end subroutine read_inputfile
 

subroutine check_value_sysflag(sysflag,info)

    character(len=15), intent(in) :: sysflag
    integer, intent(out),optional :: info

    character(len=15) :: sysflagstr(6)
    integer :: i
    logical :: flag

    ! permissible values of sysflag

    sysflagstr(1)="elect"
    sysflagstr(2)="bulk water"
    sysflagstr(3)="neutral"
    sysflagstr(4)="electnopoly"
    sysflagstr(5)="electHC"
    sysflagstr(6)="electligand"

    flag=.FALSE.

    do i=1,6
        if(sysflag==sysflagstr(i)) flag=.TRUE.
    enddo

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error: value of sysflag is not permissible"
        print*,"sysflag = ",sysflag
        if (present(info)) info = myio_err_sysflag
        return
    end if

end subroutine check_value_sysflag


subroutine check_value_runflag(runflag,info)

    character(len=15), intent(in) :: runflag
    integer, intent(out),optional :: info

    character(len=15) :: runflagstr(4)
    integer :: i
    logical :: flag

    ! permissible values of runflag

    runflagstr(1)="rangepH"
    runflagstr(2)="rangepHcpp"
    runflagstr(3)="rangepHcNaCl"
    runflagstr(4)="rangenr"

    flag=.FALSE.

    do i=1,4
        if(runflag==runflagstr(i)) flag=.TRUE.
    enddo

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error: value of runflag is not permissible"
        print*,"runflag = ",runflag
        if (present(info)) info = myio_err_runflag
        return
    end if

end subroutine check_value_runflag

subroutine check_value_bcflag(bcflag,info)

    character(len=2), intent(in) :: bcflag
    integer, intent(out), optional :: info

    character(len=2) :: bcvalues(7)
    integer :: i
    logical :: flag

    ! permissible values of bcflag

    bcvalues(1)="qu"
    bcvalues(2)="cl"
    bcvalues(3)="ca"
    bcvalues(4)="ta"
    bcvalues(5)="cc"
    bcvalues(6)="pp"
    bcvalues(7)="pd"

    flag=.FALSE.
            
    do i=1,7
        if(bcflag==bcvalues(i)) flag=.TRUE.
    enddo
    if (flag.eqv. .FALSE.) then
        print*,"Error value of bcflag is not permissible"
        print*,"bcflag = ",bcflag
        if (present(info)) info = myio_err_bcflag
        return
    endif

end subroutine check_value_bcflag


subroutine check_value_geometry(geometry,info)
        

    character(len=11), intent(in) :: geometry
    integer, intent(out),optional :: info

    logical :: flag
    character(len=11) :: geometrystr(4) 
    integer :: i

    ! permissible values of geometry

    geometrystr(1)="planar"
    geometrystr(2)="spherical"
    geometrystr(3)="cylindrical"
    geometrystr(4)="invcylindrical"
    
    
    flag=.FALSE.

    do i=1,4
        if(geometry==geometrystr(i)) flag=.TRUE.
    enddo
        
    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then 
        print*,"Error: value of geometry is not permissible"
        print*,"geometry = ",geometry
        if (present(info)) info = myio_err_geometry 
        return
    endif
    
end subroutine check_value_geometry


subroutine check_value_chaintype(chaintype,info)

    character(len=8), intent(in) :: chaintype
    integer, intent(out),optional :: info

    logical :: flag
    character(len=8) :: chaintypestr(3)
    integer :: i

    ! permissible values of chaintype

    chaintypestr(1)="diblock"
    chaintypestr(2)="altA"
    chaintypestr(3)="altB"

    flag=.FALSE.

    do i=1,3
        if(chaintype==chaintypestr(i)) flag=.TRUE.
    enddo

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error: value of chaintype is not permissible"
        print*,"chaintype = ",chaintype
        if (present(info)) info = myio_err_chaintype
        return
    endif

end subroutine check_value_chaintype

subroutine check_value_method(method,info)

    character(len=8), intent(in) :: method
    integer, intent(out),optional :: info

    character(len=8) :: methodstr
    integer :: i
    logical :: flag

    ! permissible values of runflag

    methodstr="kinsol"

    flag=.FALSE.

    if (method==methodstr) flag=.TRUE.

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error: value of method is not permissible"
        print*,"method = ",method
        if (present(info)) info = myio_err_method
        return
    endif

end subroutine check_value_method



subroutine set_value_concen(runflag,info)

    use myutils, only : newunit

    character(len=15), intent(in) :: runflag
    integer, intent(out),optional :: info

    character(len=9) :: fname
    integer :: ios
    integer :: i    
    integer :: un_cs
    
    if (present(info)) info = 0

    if(runflag=="rangepHcpp".or.runflag=="rangepHcNaCl") then

       !     .. read concentrations of cpp or NaCl from file
        write(fname,'(A9)')'concen.in'
        open(unit=newunit(un_cs),file=fname,iostat=ios,status='old')
        if(ios > 0 ) then
            print*, 'Error opening file concen.in : iostat =', ios
            if (present(info)) then
                info = myio_err_inputfile
                return 
            else
                stop  
            endif        
        endif

        read(un_cs,*)num_concen ! read number of concentrations form file
        allocate(concen_array(num_concen)) 
            
        do i=1,num_concen     ! read value salt concentration
            read(un_cs,*)concen_array(i)
        enddo    
        close(un_cs)
    endif    

end subroutine  set_value_concen



subroutine output()

    use globals, only : sysflag
    implicit none

    if(sysflag=="elect") then 
        call output_elect
    elseif(sysflag=="neutral") then
        call output_neutral
    elseif(sysflag=="electnopoly") then
        call output_elect_nopoly
        call output_individualcontr_fe
    elseif(sysflag=="electligand") then
        call output_elect_nopoly
    else
        print*,"Error in output subroutine"
        print*,"Wrong value sysflag : ", sysflag
    endif     

end subroutine output


subroutine output_elect
  
    !     .. variables and constant declaractions
    use globals 
    use volume
    use parameters
    use field
    use energy
    use surface 
    use myutils, only : newunit
  
    !     .. output file names       
    
    character(len=90) :: sysfilename     
    character(len=90) :: xsolfilename 
    character(len=90) :: xpolABfilename 
    character(len=90) :: xpolCfilename 
    character(len=90) :: xpolendfilename 
    character(len=90) :: xNafilename
    character(len=90) :: xKfilename
    character(len=90) :: xTBfilename
    character(len=90) :: xppfilename
    character(len=90) :: xppfdisfilename
    character(len=90) :: xCafilename
    character(len=90) :: xNaClfilename
    character(len=90) :: xKClfilename
    character(len=90) :: xClfilename
    character(len=90) :: potentialfilename
    character(len=90) :: chargefilename
    character(len=90) :: xHplusfilename
    character(len=90) :: xOHminfilename
    character(len=90) :: densfracAfilename
    character(len=90) :: densfracBfilename
    character(len=90) :: densfracionpairfilename

    integer :: i,j,k,t     ! dummy indexes
    character(len=100) :: fnamelabel
    character(len=20) :: rstr
    logical :: isopen
    real(dp) :: xppfdis(5),cppfdis(5)
    real(dp) :: cppbulk


    ! .. executable statements 

    ! .. make label filenames 

    if((sysflag=="electnopoly".or.sysflag=="electligand").and.(bcflag=="pp".or.bcflag=="pd")) then 

         ! filelabel for qdot only                     
        
        write(rstr,'(F5.3)')sigmaSurf/(4.0_dp*pi*lb*delta)
        fnamelabel="sg"//trim(adjustl(rstr))
        write(rstr,'(F5.3)')cTBCl
        fnamelabel=trim(fnamelabel)//"cTBCl"//trim(adjustl(rstr))
        write(rstr,'(F5.3)')cNaCl
        fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))
        if(cCaCl2/=0.0_dp) then      
            write(rstr,'(F5.3)')cCaCl2
            fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
        endif    
        if(cpp/=0.0_dp) then      
            if(cpp>=0.001) then 
                write(rstr,'(F5.3)')cpp
            else
                write(rstr,'(ES8.2E2)')cpp
            endif       
            fnamelabel=trim(fnamelabel)//"cpp"//trim(adjustl(rstr))
        endif   
        write(rstr,'(F7.3)')pHbulk
        fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))
        ! nr variable in file names only in rangenr
        if(runflag/="rangenr") then 
            fnamelabel=trim(fnamelabel)//".dat"
        else 
            write(rstr,'(I4)')nr
            fnamelabel=trim(fnamelabel)//"nr"//trim(adjustl(rstr))//".dat"
        endif     

    else

        write(rstr,'(F5.3)')sigmaAB*delta 
        fnamelabel="sg"//trim(adjustl(rstr)) 
        write(rstr,'(F5.3)')cNaCl
        fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))
        write(rstr,'(F5.3)')cCaCl2
        fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
        write(rstr,'(F7.3)')pHbulk
        fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))//".dat"
    endif

    sysfilename='system.'//trim(fnamelabel)
    xpolABfilename='xpolAB.'//trim(fnamelabel)
    xpolCfilename='xpolC.'//trim(fnamelabel)
    xsolfilename='xsol.'//trim(fnamelabel)
    xNafilename='xNaions.'//trim(fnamelabel)
    xKfilename='xKions.'//trim(fnamelabel)
    xTBfilename='xTBions.'//trim(fnamelabel)
    xCafilename='xCaions.'//trim(fnamelabel)
    xNaClfilename='xNaClionpair.'//trim(fnamelabel)
    xKClfilename='xKClionpair.'//trim(fnamelabel)
    xClfilename='xClions.'//trim(fnamelabel)
    potentialfilename='potential.'//trim(fnamelabel)
    chargefilename='charge.'//trim(fnamelabel)
    xHplusfilename='xHplus.'//trim(fnamelabel)
    xOHminfilename='xOHmin.'//trim(fnamelabel)
    densfracAfilename='densityAfrac.'//trim(fnamelabel)
    densfracBfilename='densityBfrac.'//trim(fnamelabel)
    densfracionpairfilename='densityfracionpair.'//trim(fnamelabel)
    xppfilename='xppions.'//trim(fnamelabel)

    !     .. opening files        
    
    open(unit=newunit(un_sys),file=sysfilename)       
    open(unit=newunit(un_xsol),file=xsolfilename)
    open(unit=newunit(un_psi),file=potentialfilename)

    if(sysflag/="electnopoly") then          
        open(unit=newunit(un_xpolAB),file=xpolABfilename)
        open(unit=newunit(un_xpolC),file=xpolCfilename)
        open(unit=newunit(un_fdisA),file=densfracAfilename) 
        open(unit=newunit(un_fdisB),file=densfracBfilename) 
    endif   

    if(sysflag=="electligand") open(unit=newunit(un_xpp),file=xppfilename)
      
    if(verboseflag=="yes") then    
        open(unit=newunit(un_xNa),file=xNafilename)
        open(unit=newunit(un_xK),file=xKfilename)
        open(unit=newunit(un_xCa),file=xCafilename)
        open(unit=newunit(un_xTB),file=xTBfilename)
        open(unit=newunit(un_xNaCl),file=xNaClfilename)
        open(unit=newunit(un_xKCl),file=xKClfilename)
        open(unit=newunit(un_xpair),file=densfracionpairfilename)
        open(unit=newunit(un_xCl),file=xClfilename)
        open(unit=newunit(un_charge),file=chargefilename)
        open(unit=newunit(un_xHplus),file=xHplusfilename)
        open(unit=newunit(un_xOHmin),file=xOHminfilename)
    endif
    
    !   .. writting files   

    select case (geometry)
        case ("spherical")
            write(un_psi,*)radius,psiSurf
        case ("cylindrical")
            write(un_psi,*)radius,psiSurf
        case ("planar")
                write(un_psi,*)0.0,psiSurf 
        ! case invcylinder append at end file (un_psi) not begining
    end select  

    do i=1,nr
        write(un_xsol,*)rc(i),xsol(i)
        write(un_psi,*)rc(i),psi(i)
    enddo    
   
    if(geometry=="invcylindrical") write(un_psi,*)radius,psiSurf

    if(sysflag/="electnopoly") then 
        do i=1,nr
            write(un_xpolAB,fmt4reals)rc(i),xpolAB(i),rhopolA(i),rhopolB(i)
            write(un_xpolC,fmt2reals)rc(i),xpolC(i)
            write(un_fdisA,fmt6reals)rc(i),fdisA(1,i),fdisA(2,i),fdisA(3,i),fdisA(4,i),fdisA(5,i)        
            write(un_fdisB,fmt6reals)rc(i),fdisB(1,i),fdisB(2,i),fdisB(3,i),fdisB(4,i),fdisB(5,i)
        enddo
    endif   

    if(sysflag=="electligand") then 
        do i=1,nr
            write(un_xpp,fmt6reals)rc(i),xpp(i,AH2BH),xpp(i,AHBH),xpp(i,AHB),xpp(i,ABH),xpp(i,AB)
            do t=1,5    
                cppfdis(t)=(xbulk%pp(t)/(vpp(t)*vsol))/cppbulk
            enddo    
        enddo
    endif   
    
    if(verboseflag=="yes") then 
        do i=1,nr
            write(un_xNa,*)rc(i),xNa(i)
            write(un_xK,*)rc(i),xK(i)
            write(un_xCa,*)rc(i),xCa(i)
            write(un_xTB,*)rc(i),xTB(i) 
            write(un_xNaCl,*)rc(i),xNaCl(i)
            write(un_xKCl,*)rc(i),xKCl(i)
            write(un_xpair,*)rc(i),(xNaCl(i)/vNaCl)/(xNa(i)/vNa+xCl(i)/vCl+xNaCl(i)/vNaCl)
            write(un_xCl,*)rc(i),xCl(i)
            write(un_charge,*)rc(i),rhoq(i)
            write(un_xHplus,*)rc(i),xHplus(i)
            write(un_xOHmin,*)rc(i),xOHmin(i)    
        enddo    
    endif

    write(un_sys,*)'system      = planar weakpolyelectrolyte brush'
    write(un_sys,*)'version     = ',VERSION
    write(un_sys,*)'chainmethod = ',chainmethod
    write(un_sys,*)'chaintype   = ',chaintype
    if(chainmethod.eq."FILE") then
       write(un_sys,*)'readinchains = ',readinchains
    endif
    write(un_sys,*)'sysflag     = ',sysflag
    write(un_sys,*)'bcflag      = ',bcflag
    write(un_sys,*)'nsegAB      = ',nsegAB
    write(un_sys,*)'lsegAB      = ',lsegAB
    write(un_sys,*)'nsegC       = ',nsegC
    write(un_sys,*)'lsegC       = ',lsegC
    write(un_sys,*)'period      = ',period 
    write(un_sys,*)'delta       = ',delta  
    write(un_sys,*)'vsol        = ',vsol
    write(un_sys,*)'vpolA(1)    = ',vpolA(1)*vsol
    write(un_sys,*)'vpolA(2)    = ',vpolA(2)*vsol
    write(un_sys,*)'vpolA(3)    = ',vpolA(3)*vsol
    write(un_sys,*)'vpolA(4)    = ',vpolA(4)*vsol
    write(un_sys,*)'vpolA(5)    = ',vpolA(5)*vsol
    write(un_sys,*)'vpolB(1)    = ',vpolB(1)*vsol
    write(un_sys,*)'vpolB(2)    = ',vpolB(2)*vsol
    write(un_sys,*)'vpolB(3)    = ',vpolB(3)*vsol
    write(un_sys,*)'vpolB(4)    = ',vpolB(4)*vsol
    write(un_sys,*)'vpolB(5)    = ',vpolB(5)*vsol
    write(un_sys,*)'vpolC       = ',vpolC*vsol
    write(un_sys,*)'vNa         = ',vNa*vsol
    write(un_sys,*)'vCl         = ',vCl*vsol
    write(un_sys,*)'vCa         = ',vCa*vsol
    write(un_sys,*)'vK          = ',vK*vsol
    if(bcflag=="pd")then 
        write(un_sys,*)'vpp(AH2BH)  = ',vpp(AH2BH)*vsol
        write(un_sys,*)'vpp(AHBH)   = ',vpp(AHBH)*vsol
        write(un_sys,*)'vpp(AHB)    = ',vpp(AHB)*vsol
        write(un_sys,*)'vpp(ABH)    = ',vpp(ABH)*vsol
        write(un_sys,*)'vpp(AB)     = ',vpp(AB)*vsol
    endif    
    if(bcflag=="pp".or.bcflag=="pd") write(un_sys,*)'vTB         = ',vTB*vsol
    write(un_sys,*)'vNaCl       = ',vNaCl*vsol
    write(un_sys,*)'vKCl        = ',vKCl*vsol
    write(un_sys,*)'cNaCl       = ',cNaCl
    write(un_sys,*)'cKCl        = ',cKCl
    write(un_sys,*)'cCaCl2      = ',cCaCl2
    if(bcflag=="pp".or.bcflag=="pd") then 
        write(un_sys,*)'cTBCl       = ',cTBCl
        write(un_sys,*)'cpp         = ',cpp
    endif    
    write(un_sys,*)'pHbulk      = ',pHbulk
    write(un_sys,*)'pKa         = ',pKa(1)      
    write(un_sys,*)'pKaNa       = ',pKa(2)
    write(un_sys,*)'pKaACa      = ',pKa(3)
    write(un_sys,*)'pKaA2Ca     = ',pKa(4)
    write(un_sys,*)'pKb         = ',pKb(1)      
    write(un_sys,*)'pKbNa       = ',pKb(2)
    write(un_sys,*)'pKbBCa      = ',pKb(3)
    write(un_sys,*)'pKbB2Ca     = ',pKb(4)
    write(un_sys,*)'KionNa      = ',KionNa
    write(un_sys,*)'KionK       = ',KionK
    write(un_sys,*)'K0ionNa     = ',K0ionNa
    write(un_sys,*)'K0ionK      = ',K0ionK
    write(un_sys,*)'xbulk%sol   = ',xbulk%sol
    write(un_sys,*)'xbulk%Na    = ',xbulk%Na
    write(un_sys,*)'xbulk%Cl    = ',xbulk%Cl
    write(un_sys,*)'xbulk%K     = ',xbulk%K
    write(un_sys,*)'xbulk%NaCl  = ',xbulk%NaCl
    write(un_sys,*)'xbulk%KCl   = ',xbulk%KCl
    write(un_sys,*)'xbulk%Ca    = ',xbulk%Ca
    write(un_sys,*)'xbulk%Hplus = ',xbulk%Hplus
    write(un_sys,*)'xbulk%OHmin = ',xbulk%OHmin
    if(bcflag=="pp".or.bcflag=="pd") write(un_sys,*)'xbulk%TB    = ',xbulk%TB
    if(sysflag=="electligand") then
        cppbulk = (cpp*Na/(1.0e24_dp))
        write(un_sys,*)'xbulk%pp(AH2BH) = ',xbulk%pp(AH2BH)
        write(un_sys,*)'xbulk%pp(AHBH)  = ',xbulk%pp(AHBH)
        write(un_sys,*)'xbulk%pp(AHB)   = ',xbulk%pp(AHB)
        write(un_sys,*)'xbulk%pp(ABH)   = ',xbulk%pp(ABH)
        write(un_sys,*)'xbulk%pp(AB)    = ',xbulk%pp(AB)
        cppbulk = (cpp*Na/(1.0e24_dp))
        do t=1,5       
            cppfdis(t)=(xbulk%pp(t)/(vpp(t)*vsol))/cppbulk
        enddo    
        write(un_sys,*)'fdis(AH2BH) = ',cppfdis(AH2BH)
        write(un_sys,*)'fdis(AHBH)  = ',cppfdis(AHBH)
        write(un_sys,*)'fdis(AHB)   = ',cppfdis(AHB)
        write(un_sys,*)'fdis(ABH)   = ',cppfdis(ABH)
        write(un_sys,*)'fdis(AB)    = ',cppfdis(AB)
    endif    
    write(un_sys,*)'sigmaAB     = ',sigmaAB*delta
    write(un_sys,*)'sigmaC      = ',sigmaC*delta
    write(un_sys,*)'dielectW    = ',dielectW
    write(un_sys,*)'lb          = ',lb
    write(un_sys,*)'T           = ',Temp
    write(un_sys,*)'VdWepsC     = ',VdWepsC*vpolC*vsol 
    write(un_sys,*)'VdWepsB     = ',VdWepsB*vpolB(3)*vsol
    write(un_sys,*)'zpolA(1)    = ',zpolA(1)
    write(un_sys,*)'zpolA(2)    = ',zpolA(2)
    write(un_sys,*)'zpolA(3)    = ',zpolA(3)
    write(un_sys,*)'zpolA(4)    = ',zpolA(4)
    write(un_sys,*)'zpolB(1)    = ',zpolB(1)
    write(un_sys,*)'zpolB(2)    = ',zpolB(2)
    write(un_sys,*)'zpolB(3)    = ',zpolB(3)
    write(un_sys,*)'zpolB(4)    = ',zpolB(4)
    write(un_sys,*)'zpolB(5)    = ',zpolB(5)
    write(un_sys,*)'zNa         = ',zNa
    write(un_sys,*)'zCa         = ',zCa
    write(un_sys,*)'zK          = ',zK
    write(un_sys,*)'zCl         = ',zCl
    if(bcflag=="pd")then 
        write(un_sys,*)'zpp(AH2BH)  = ',zpp(AH2BH)
        write(un_sys,*)'zpp(AHBH)   = ',zpp(AHBH)
        write(un_sys,*)'zpp(AHB)    = ',zpp(AHB)
        write(un_sys,*)'zpp(ABH)    = ',zpp(ABH)
        write(un_sys,*)'zpp(AB)     = ',zpp(AB)
    endif    
    write(un_sys,*)'nr          = ',nr
    write(un_sys,*)'free energy = ',FE
    write(un_sys,*)'energy bulk = ',FEbulk 
    write(un_sys,*)'deltafenergy = ',deltaFE
    write(un_sys,*)'fnorm       = ',fnorm
    write(un_sys,*)'q residual  = ',qres
    write(un_sys,*)'error       = ',error
    write(un_sys,*)'sigmaAB     = ',sigmaAB*delta
    write(un_sys,*)'sumphiA     = ',sumphiA
    write(un_sys,*)'sumphiB     = ',sumphiB
    write(un_sys,*)'sumphiC     = ',sumphiC
    write(un_sys,*)'check phi   = ',checkphi 
    write(un_sys,*)'FEq         = ',FEq 
    write(un_sys,*)'FEpi        = ',FEpi
    write(un_sys,*)'FErho       = ',FErho
    write(un_sys,*)'FEel        = ',FEel
    write(un_sys,*)'FEelsurf    = ',FEelsurf
    write(un_sys,*)'FEbind      = ',FEbind
    write(un_sys,*)'FEVdW       = ',FEVdW 
    write(un_sys,*)'FEalt       = ',FEalt
    write(un_sys,*)'qAB         = ',qAB
    write(un_sys,*)'qC          = ',qC
    write(un_sys,*)'muAB        = ',-log(qAB)
    write(un_sys,*)'muC         = ',-log(qC)
    write(un_sys,*)'heightAB    = ',heightAB
    write(un_sys,*)'heightC     = ',heightC
    write(un_sys,*)'qpolA       = ',qpolA
    write(un_sys,*)'qpolB       = ',qpolB
    write(un_sys,*)'qpoltot     = ',qpol_tot
    write(un_sys,*)'avfdisA(1)  = ',avfdisA(1)
    write(un_sys,*)'avfdisA(2)  = ',avfdisA(2)
    write(un_sys,*)'avfdisA(3)  = ',avfdisA(3)
    write(un_sys,*)'avfdisA(4)  = ',avfdisA(4)
    write(un_sys,*)'avfdisA(5)  = ',avfdisA(5)
    write(un_sys,*)'avfdisB(1)  = ',avfdisB(1)
    write(un_sys,*)'avfdisB(2)  = ',avfdisB(2)
    write(un_sys,*)'avfdisB(3)  = ',avfdisB(3)
    write(un_sys,*)'avfdisB(4)  = ',avfdisB(4)
    write(un_sys,*)'avfdisB(5)  = ',avfdisB(5)
    write(un_sys,*)'sigmaSurf   = ',sigmaSurf/(4.0_dp*pi*lb*delta)
    write(un_sys,*)'sigmaqSurf  = ',sigmaqSurf/(4.0_dp*pi*lb*delta)
    write(un_sys,*)'psiSurf     = ',psiSurf
    if(bcflag=='ta') then
        do i=1,4   
            write(un_sys,fmt)'fdisTa(',i,')   = ',fdisTaL(i)
        enddo  
    else if(bcflag=='pp') then   
        do i=1,4   
            write(un_sys,fmt)'fdisSu(',i,')   = ',fdisS(i)
        enddo
    else if(bcflag=='pd') then   
        do i=1,4   
            write(un_sys,fmt)'fdisSu(',i,')   = ',fdisS(i)
        enddo  
        write(un_sys,*)'fdisR      = ',fdisR
        write(un_sys,*)'sigmaR     = ',fdisR*sigmaSurf/(4.0_dp*pi*lb*delta)
        write(un_sys,*)'sigmaLR    = ',(1.0_dp-fdisR)*sigmaSurf/(4.0_dp*pi*lb*delta)
    else
        do i=1,6   
            write(un_sys,fmt)' fdisSu(',i,')  = ',fdisS(i)
        enddo  
    endif
    write(un_sys,*)'nsize       = ',nsize  
    write(un_sys,*)'cuantasAB   = ',cuantasAB
    write(un_sys,*)'cuantasC    = ',cuantasC
    write(un_sys,*)'iterations  = ',iter
   
    ! .. closing files

    close(un_sys)
    close(un_xsol)
    close(un_psi)
    if(sysflag/="electnopoly") then
        close(un_xpolAB)   
        close(un_xpolC)
        close(un_fdisA)
        close(un_fdisB)
    endif
    if(sysflag=="electligand") close(un_xpp)
    if(verboseflag=="yes") then 
        close(un_xNa)   
        close(un_xK)
        close(un_xCa)
        close(un_xNaCl)
        close(un_xKCl)
        close(un_xpair)
        close(un_xCl)
        close(un_charge)
        close(un_xHplus)
        close(un_xOHmin)
        close(un_xTB)
    endif
        

end subroutine output_elect



subroutine output_neutral
  
    !     .. variables and constant declaractions
    use globals 
    use volume
    use parameters    
    use field
    use energy
    use myutils, only : newunit

    !     .. output file names         
    character(len=90) :: sysfilename     
    character(len=90) :: xsolfilename 
    character(len=90) :: xpolABfilename 
    character(len=90) :: xpolCfilename 
    character(len=90) :: xpolendfilename 
    
    character(len=80) :: fmt2reals,fmt3reals,fmt4reals,fmt5reals,fmt6reals   

    !     .. local arguments
    integer :: i
    character(len=100) :: fnamelabel
    character(len=20) :: rstr
    logical :: isopen
    !     .. executable statements 

    fmt2reals = "(2ES25.16)"  
    fmt3reals = "(3ES25.16)"  
    fmt4reals = "(4ES25.16)"  
    fmt5reals = "(5ES25.16)" 
    fmt6reals = "(6ES25.16)" 

    !     .. make label filenames 
    write(rstr,'(F5.3)')sigmaAB*delta 
    fnamelabel="sg"//trim(adjustl(rstr)) 
    write(rstr,'(F5.3)')VdWepsB
    fnamelabel=trim(fnamelabel)//"VdWepsB"//trim(adjustl(rstr))//".dat"

    !     .. make filenames 
    sysfilename='system.'//trim(fnamelabel)
    xpolABfilename='xpolAB.'//trim(fnamelabel)   
    xpolCfilename='xpolC.'//trim(fnamelabel)   
    xsolfilename='xsol.'//trim(fnamelabel)   
    xpolendfilename='xpolend.'//trim(fnamelabel)   
    
    !      .. opening files
    open(unit=newunit(un_sys),file=sysfilename)   
    open(unit=newunit(un_xpolAB),file=xpolABfilename)
    open(unit=newunit(un_xpolC),file=xpolCfilename)
    open(unit=newunit(un_xsol),file=xsolfilename)

    
    do i=1,nr    
       write(un_xpolAB,fmt4reals)rc(i),xpolAB(i),rhopolA(i),rhopolB(i)
       write(un_xpolC,fmt2reals)rc(i),xpolC(i)
       write(un_xsol,fmt2reals)rc(i),xsol(i)
    !     write(40,*)rc(i),endpol(i)
    enddo
        
    !     .. system information 

    write(un_sys,*)'system      = planar  brush' 
    write(un_sys,*)'version     = ',VERSION
    write(un_sys,*)'sysflag     = ',sysflag
    write(un_sys,*)'chainmethod = ',chainmethod
    write(un_sys,*)'chaintype   = ',chaintype
    if(chainmethod.eq."FILE") then
        write(un_sys,*)'readinchains = ',readinchains
    endif
    write(un_sys,*)'sysflag     = ',sysflag
    write(un_sys,*)'nsegAB      = ',nsegAB
    write(un_sys,*)'lsegAB      = ',lsegAB
    write(un_sys,*)'nsegC       = ',nsegC
    write(un_sys,*)'lsegC       = ',lsegC
    write(un_sys,*)'period      = ',period
    write(un_sys,*)'delta       = ',delta
    write(un_sys,*)'tol_conv    = ',error
    write(un_sys,*)'vsol        = ',vsol
    write(un_sys,*)'vpolA(1)    = ',vpolA(1)*vsol
    write(un_sys,*)'vpolA(2)    = ',vpolA(2)*vsol
    write(un_sys,*)'vpolA(3)    = ',vpolA(3)*vsol
    write(un_sys,*)'vpolA(4)    = ',vpolA(4)*vsol
    write(un_sys,*)'vpolA(5)    = ',vpolA(5)*vsol
    write(un_sys,*)'vpolB(1)    = ',vpolB(1)*vsol
    write(un_sys,*)'vpolB(2)    = ',vpolB(2)*vsol
    write(un_sys,*)'vpolB(3)    = ',vpolB(3)*vsol
    write(un_sys,*)'vpolB(4)    = ',vpolB(4)*vsol
    write(un_sys,*)'vpolB(5)    = ',vpolB(5)*vsol
    write(un_sys,*)'vpolC       = ',vpolC*vsol
    write(un_sys,*)'T           = ',Temp
    write(un_sys,*)'VdWepsC     = ',VdWepsC*vpolC*vsol
    write(un_sys,*)'VdWepsB     = ',VdWepsB*vpolB(3)*vsol
    write(un_sys,*)'cuantasAB   = ',cuantasAB
    write(un_sys,*)'cuantasC    = ',cuantasC
    
    write(un_sys,*)'distance  = ',nr*delta 
    write(un_sys,*)'nr          = ',nr
    write(un_sys,*)'free energy = ',FE  
    write(un_sys,*)'energy bulk = ',FEbulk 
    write(un_sys,*)'deltafenergy = ',deltaFE
    write(un_sys,*)'FEalt       = ',FEalt
    write(un_sys,*)'FEconfC     = ',FEconfC
    write(un_sys,*)'FEconfAB    = ',FEconfAB
    write(un_sys,*)'FEtrans%sol = ',FEtrans%sol  
    write(un_sys,*)'fnorm       = ',fnorm
    write(un_sys,*)'sumphiA     = ',sumphiA
    write(un_sys,*)'sumphiB     = ',sumphiB
    write(un_sys,*)'sumphiC     = ',sumphiC
    write(un_sys,*)'check phi   = ',checkphi 
    write(un_sys,*)'FEq         = ',FEq 
    write(un_sys,*)'FEpi        = ',FEpi
    write(un_sys,*)'FErho       = ',FErho
    write(un_sys,*)'FEVdW       = ',FEVdW
    write(un_sys,*)'qAB         = ',qAB
    write(un_sys,*)'qC          = ',qC
    write(un_sys,*)'muAB        = ',-dlog(qAB)
    write(un_sys,*)'muC         = ',-dlog(qC)
    write(un_sys,*)'heightAB    = ',heightAB
    write(un_sys,*)'heightC     = ',heightC
    write(un_sys,*)'nsize       = ',nsize  
    write(un_sys,*)'cuantasAB   = ',cuantasAB
    write(un_sys,*)'cuantasC    = ',cuantasC
    write(un_sys,*)'iterations  = ',iter
    
    ! .. closing files
 
    close(un_xsol)
    close(un_xpolAB)
    close(un_xpolC)
    close(un_sys)

  
end subroutine output_neutral



subroutine output_elect_nopoly
  
    !     .. variables and constant declaractions
    use globals 
    use volume
    use parameters
    use field
    use energy
    use surface 
    use myutils, only : newunit
  
    !     .. output file names       
    
    character(len=90) :: sysfilename  
    character(len=90) :: xsolfilename   
    character(len=90) :: xNafilename
    character(len=90) :: xKfilename
    character(len=90) :: xTBfilename
    character(len=90) :: xppfilename
    character(len=90) :: xppfdisfilename
    character(len=90) :: xCafilename
    character(len=90) :: xClfilename
    character(len=90) :: potentialfilename
    character(len=90) :: chargefilename
    character(len=90) :: xHplusfilename
    character(len=90) :: xOHminfilename

    integer :: i,j,k,t     ! dummy indexes
    character(len=100) :: fnamelabel
    character(len=20) :: rstr
    logical :: isopen
    real(dp) :: xppfdis(5),cppfdis(5)
    real(dp) :: cppbulk
    real(dp) :: eta


    ! .. executable statements 

    ! .. make label filenames 

    
    ! filelabel for qdot only                     
        
    write(rstr,'(F5.3)')sigmaSurf/(4.0_dp*pi*lb*delta)
    fnamelabel="sg"//trim(adjustl(rstr))
    write(rstr,'(F5.3)')cTBCl
    fnamelabel=trim(fnamelabel)//"cTBCl"//trim(adjustl(rstr))
    write(rstr,'(F5.3)')cNaCl
    fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))
    if(cCaCl2/=0.0_dp) then      
        write(rstr,'(F5.3)')cCaCl2
        fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
    endif    
    if(cpp/=0.0_dp) then      
        if(cpp>=0.001) then 
            write(rstr,'(F5.3)')cpp
        else
            write(rstr,'(ES8.2E2)')cpp
        endif       
        fnamelabel=trim(fnamelabel)//"cpp"//trim(adjustl(rstr))
    endif   
    write(rstr,'(F7.3)')pHbulk
    fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))
    ! nr variable in file names only in rangenr
    if(runflag/="rangenr") then 
        fnamelabel=trim(fnamelabel)//".dat"
    else 
        write(rstr,'(I4)')nr
        fnamelabel=trim(fnamelabel)//"nr"//trim(adjustl(rstr))//".dat"
    endif     




    sysfilename='system.'//trim(fnamelabel)
    xsolfilename='xsol.'//trim(fnamelabel)
    xNafilename='xNaions.'//trim(fnamelabel)
    xKfilename='xKions.'//trim(fnamelabel)
    xTBfilename='xTBions.'//trim(fnamelabel)
    xCafilename='xCaions.'//trim(fnamelabel)
    xClfilename='xClions.'//trim(fnamelabel)
    potentialfilename='potential.'//trim(fnamelabel)
    chargefilename='charge.'//trim(fnamelabel)
    xHplusfilename='xHplus.'//trim(fnamelabel)
    xOHminfilename='xOHmin.'//trim(fnamelabel)
    xppfilename='xppions.'//trim(fnamelabel)

    !     .. opening files        
    
    open(unit=newunit(un_sys),file=sysfilename)
    open(unit=newunit(un_psi),file=potentialfilename)      
    
    if(verboseflag=="yes") then  

        open(unit=newunit(un_xsol),file=xsolfilename)
        if(sysflag=="electligand") open(unit=newunit(un_xpp),file=xppfilename)
        open(unit=newunit(un_xNa),file=xNafilename)
        open(unit=newunit(un_xK),file=xKfilename)
        open(unit=newunit(un_xCa),file=xCafilename)
        open(unit=newunit(un_xTB),file=xTBfilename)
        open(unit=newunit(un_xCl),file=xClfilename)
        open(unit=newunit(un_charge),file=chargefilename)
        open(unit=newunit(un_xHplus),file=xHplusfilename)
        open(unit=newunit(un_xOHmin),file=xOHminfilename)
    endif
    
    !   .. writting files   

    select case (geometry)
        case ("spherical")
            write(un_psi,*)radius,psiSurf
        case ("cylindrical")
            write(un_psi,*)radius,psiSurf
        case ("planar")
                write(un_psi,*)0.0,psiSurf 
        ! case invcylinder append at end file (un_psi) not begining
    end select  

    do i=1,nr
        write(un_psi,*)rc(i),psi(i)
    enddo    
   
    if(geometry=="invcylindrical") write(un_psi,*)radius,psiSurf

    if(verboseflag=="yes") then 
        if(sysflag=="electligand") then 
            do i=1,nr
                write(un_xpp,fmt6reals)rc(i),xpp(i,AH2BH),xpp(i,AHBH),xpp(i,AHB),xpp(i,ABH),xpp(i,AB)
                do t=1,5    
                    cppfdis(t)=(xbulk%pp(t)/(vpp(t)*vsol))/cppbulk
                enddo    
            enddo
        endif

        do i=1,nr
            write(un_xsol,*)rc(i),xsol(i)
            write(un_xNa,*)rc(i),xNa(i)
            write(un_xK,*)rc(i),xK(i)
            write(un_xCa,*)rc(i),xCa(i)
            write(un_xTB,*)rc(i),xTB(i) 
            write(un_xCl,*)rc(i),xCl(i)
            write(un_charge,*)rc(i),rhoq(i)
            write(un_xHplus,*)rc(i),xHplus(i)
            write(un_xOHmin,*)rc(i),xOHmin(i)    
        enddo    
    endif

    write(un_sys,*)'system      = electrolyte solition qdot/NP'
    write(un_sys,*)'version     = ',VERSION
    write(un_sys,*)'sysflag     = ',sysflag
    write(un_sys,*)'bcflag      = ',bcflag
    write(un_sys,*)'delta       = ',delta   
    write(un_sys,*)'vsol        = ',vsol
    write(un_sys,*)'vNa         = ',vNa*vsol
    write(un_sys,*)'vCl         = ',vCl*vsol
    write(un_sys,*)'vCa         = ',vCa*vsol
    write(un_sys,*)'vK          = ',vK*vsol
    if(bcflag=="pd")then 
        write(un_sys,*)'vpp(AH2BH)  = ',vpp(AH2BH)*vsol
        write(un_sys,*)'vpp(AHBH)   = ',vpp(AHBH)*vsol
        write(un_sys,*)'vpp(AHB)    = ',vpp(AHB)*vsol
        write(un_sys,*)'vpp(ABH)    = ',vpp(ABH)*vsol
        write(un_sys,*)'vpp(AB)     = ',vpp(AB)*vsol
    endif    
    if(bcflag=="pp".or.bcflag=="pd") write(un_sys,*)'vTB         = ',vTB*vsol
    write(un_sys,*)'vNaCl       = ',vNaCl*vsol
    write(un_sys,*)'vKCl        = ',vKCl*vsol
    write(un_sys,*)'cNaCl       = ',cNaCl
    write(un_sys,*)'cKCl        = ',cKCl
    write(un_sys,*)'cCaCl2      = ',cCaCl2
    if(bcflag=="pp".or.bcflag=="pd") then 
        write(un_sys,*)'cTBCl       = ',cTBCl
        write(un_sys,*)'cpp         = ',cpp
        write(un_sys,*)'deltaG0ads  = ',deltaG0ads
        write(un_sys,*)'deltaGads   = ',deltaG0ads-log(Na*vsol/1.0e24_dp)
    endif    
    write(un_sys,*)'pHbulk      = ',pHbulk
    write(un_sys,*)'xbulk%sol   = ',xbulk%sol
    write(un_sys,*)'xbulk%Na    = ',xbulk%Na
    write(un_sys,*)'xbulk%Cl    = ',xbulk%Cl
    write(un_sys,*)'xbulk%K     = ',xbulk%K
    write(un_sys,*)'xbulk%NaCl  = ',xbulk%NaCl
    write(un_sys,*)'xbulk%KCl   = ',xbulk%KCl
    write(un_sys,*)'xbulk%Ca    = ',xbulk%Ca
    write(un_sys,*)'xbulk%Hplus = ',xbulk%Hplus
    write(un_sys,*)'xbulk%OHmin = ',xbulk%OHmin
    if(bcflag=="pp".or.bcflag=="pd") write(un_sys,*)'xbulk%TB    = ',xbulk%TB
    if(sysflag=="electligand") then
        cppbulk = (cpp*Na/(1.0e24_dp))
        write(un_sys,*)'xbulk%pp(AH2BH) = ',xbulk%pp(AH2BH)
        write(un_sys,*)'xbulk%pp(AHBH)  = ',xbulk%pp(AHBH)
        write(un_sys,*)'xbulk%pp(AHB)   = ',xbulk%pp(AHB)
        write(un_sys,*)'xbulk%pp(ABH)   = ',xbulk%pp(ABH)
        write(un_sys,*)'xbulk%pp(AB)    = ',xbulk%pp(AB)
        cppbulk = (cpp*Na/(1.0e24_dp))
        do t=1,5       
            cppfdis(t)=(xbulk%pp(t)/(vpp(t)*vsol))/cppbulk
        enddo    
        write(un_sys,*)'fdis(AH2BH) = ',cppfdis(AH2BH)
        write(un_sys,*)'fdis(AHBH)  = ',cppfdis(AHBH)
        write(un_sys,*)'fdis(AHB)   = ',cppfdis(AHB)
        write(un_sys,*)'fdis(ABH)   = ',cppfdis(ABH)
        write(un_sys,*)'fdis(AB)    = ',cppfdis(AB)

    endif    
    write(un_sys,*)'dielectW    = ',dielectW
    write(un_sys,*)'lb          = ',lb
    write(un_sys,*)'T           = ',Temp
    write(un_sys,*)'zNa         = ',zNa
    write(un_sys,*)'zCa         = ',zCa
    write(un_sys,*)'zK          = ',zK
    write(un_sys,*)'zCl         = ',zCl
    if(bcflag=="pd")then 
        write(un_sys,*)'zpp(AH2BH)  = ',zpp(AH2BH)
        write(un_sys,*)'zpp(AHBH)   = ',zpp(AHBH)
        write(un_sys,*)'zpp(AHB)    = ',zpp(AHB)
        write(un_sys,*)'zpp(ABH)    = ',zpp(ABH)
        write(un_sys,*)'zpp(AB)     = ',zpp(AB)
    endif    
    write(un_sys,*)'nr          = ',nr
    write(un_sys,*)'free energy = ',FE
    write(un_sys,*)'energy bulk = ',FEbulk 
    write(un_sys,*)'deltafenergy = ',deltaFE
    write(un_sys,*)'fnorm       = ',fnorm
    write(un_sys,*)'q residual  = ',qres
    write(un_sys,*)'error       = ',error
    write(un_sys,*)'FEpi        = ',FEpi
    write(un_sys,*)'FErho       = ',FErho
    write(un_sys,*)'FEel        = ',FEel
    write(un_sys,*)'FEelsurf    = ',FEelsurf
    write(un_sys,*)'FEbind      = ',FEbind
    write(un_sys,*)'FEalt       = ',FEalt
    write(un_sys,*)'sigmaSurf   = ',sigmaSurf/(4.0_dp*pi*lb*delta)
    write(un_sys,*)'sigmaqSurf  = ',sigmaqSurf/(4.0_dp*pi*lb*delta)
    write(un_sys,*)'psiSurf     = ',psiSurf
    if(bcflag=='ta') then
        do i=1,4   
            write(un_sys,fmt)'fdisTa(',i,')   = ',fdisTaL(i)
        enddo  
    else if(bcflag=='pp') then   
        do i=1,4   
            write(un_sys,fmt)'fdisSu(',i,')   = ',fdisS(i)
        enddo
    else if(bcflag=='pd') then   
        do i=1,4   
            write(un_sys,fmt)'fdisSu(',i,')   = ',fdisS(i)
        enddo  
        write(un_sys,*)'fdisR      = ',fdisR
        write(un_sys,*)'sigmaR     = ',fdisR*sigmaSurf/(4.0_dp*pi*lb*delta)
        write(un_sys,*)'sigmaLR    = ',(1.0_dp-fdisR)*sigmaSurf/(4.0_dp*pi*lb*delta)
    else
        do i=1,6   
            write(un_sys,fmt)' fdisSu(',i,')  = ',fdisS(i)
        enddo  
    endif
    write(un_sys,*)'nsize       = ',nsize  
    write(un_sys,*)'iterations  = ',iter
    if(runflag=="rangenr") then 
        eta=(radius/(nr*delta))**3
        write(un_sys,*)'eta       = ',eta
    endif  
   
    ! .. closing files

    close(un_sys)
    close(un_psi)
    
    if(verboseflag=="yes") then 
        if(sysflag=="electligand") close(un_xpp)
        close(un_xsol)
        close(un_xNa)   
        close(un_xK)
        close(un_xCa)
        close(un_xCl)
        close(un_charge)
        close(un_xHplus)
        close(un_xOHmin)
        close(un_xTB)
    endif
        

end subroutine output_elect_nopoly


subroutine output_individualcontr_fe

    use globals, only : sysflag
    use energy
    use myutils, only : newunit
    use parameters, only : sigmaAB,cNaCl,cCaCl2,pHbulk,VdWepsB
    use volume, only : delta,nr,nrmax,nrmin

    ! local arguments

    integer :: un_fe

    character(len=100) :: fenergyfilename   
    character(len=100) :: fnamelabel
    character(len=20) :: rstr

   !     .. make label filename

    if(sysflag=="elect".or.sysflag=="electdouble".or.sysflag=="electnopoly") then 
        write(rstr,'(F5.3)')sigmaAB*delta 
        fnamelabel="sg"//trim(adjustl(rstr)) 
        write(rstr,'(F5.3)')cNaCl
        fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))
        write(rstr,'(F5.3)')cCaCl2
        fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
        write(rstr,'(F7.3)')pHbulk
        fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))//".dat"
    elseif(sysflag=="neutral") then 
        write(rstr,'(F5.3)')sigmaAB*delta 
        fnamelabel="sg"//trim(adjustl(rstr)) 
        write(rstr,'(F5.3)')VdWepsB
        fnamelabel=trim(fnamelabel)//"VdWepsB"//trim(adjustl(rstr))//".dat"
    else
        print*,"Error in output_individualcontr_fe subroutine"
        print*,"Wrong value sysflag : ", sysflag
    endif    

    fenergyfilename='energy.'//trim(fnamelabel)   
        
    !     .. opening files        

    open(unit=newunit(un_fe),file=fenergyfilename) 

    write(un_fe,*)'FE              = ',FE  
    write(un_fe,*)'FEbulk          = ',FEbulk 
    write(un_fe,*)'deltaFE         = ',deltaFE
    write(un_fe,*)'FEalt           = ',FEalt  
    write(un_fe,*)'FEbulkalt       = ',FEbulkalt 
    write(un_fe,*)'deltaFEalt      = ',deltaFEalt
    
    write(un_fe,*)"FEtrans%sol     = ",FEtrans%sol   
    write(un_fe,*)"FEtrans%Na      = ",FEtrans%Na  
    write(un_fe,*)"FEtrans%Cl      = ",FEtrans%Cl  
    write(un_fe,*)"FEtrans%Ca      = ",FEtrans%Ca  
    write(un_fe,*)"FEtrans%K       = ",FEtrans%K
    write(un_fe,*)"FEtrans%KCl     = ",FEtrans%KCl
    write(un_fe,*)"FEtrans%NaCl    = ",FEtrans%NaCl  
    write(un_fe,*)"FEtrans%Hplus   = ",FEtrans%Hplus  
    write(un_fe,*)"FEtrans%OHmin   = ",FEtrans%OHmin  
    
    write(un_fe,*)"FEchempot%Na    = ",FEchempot%Na
    write(un_fe,*)"FEchempot%Cl    = ",FEchempot%Cl
    write(un_fe,*)"FEchempot%Ca    = ",FEchempot%Ca
    write(un_fe,*)"FEchempot%K     = ",FEchempot%K
    write(un_fe,*)"FEchempot%KCl   = ",FEchempot%KCl
    write(un_fe,*)"FEchempot%NaCl  = ",FEchempot%NaCl
    write(un_fe,*)"FEchempot%Hplus = ",FEchempot%Hplus
    write(un_fe,*)"FEchempot%OHmin = ",FEchempot%OHmin

    write(un_fe,*)"FEchemsurf      = ",FEchemsurf
    write(un_fe,*)"FEchemsurfalt   = ",FEchemsurfalt
    
    write(un_fe,*)"delta FEchemsurfalt= ",FEchemsurfalt-FEchemsurf-diffFEchemsurf
    
    close(un_fe)


end subroutine   output_individualcontr_fe

end module


