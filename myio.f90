
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
    integer :: un_sys,un_xpolAB,un_xpolC,un_xsol,un_xNa,un_xCl,un_xK,un_xCa,un_xNaCl,un_xKCl,un_xNO3
    integer :: un_xOHmin,un_xHplus,un_fdisA,un_fdisB,un_psi,un_charge, un_xpair, un_rhopolAB, un_xTB, un_xTM
    integer :: un_xpp, un_cpp, un_xRb, un_xIm
   
    ! format specifiers 
    character(len=80), parameter  :: fmt = "(A8,I1,A5,ES25.16)"
    character(len=80), parameter  :: fmt2reals = "(2ES25.16E3)"   
    character(len=80), parameter  :: fmt3reals = "(3ES25.16E3)"  
    character(len=80), parameter  :: fmt4reals = "(4ES25.16E3)" 
    character(len=80), parameter  :: fmt5reals = "(5ES25.16E3)"
    character(len=80), parameter  :: fmt6reals = "(6ES25.16E3)" 
  
    integer :: num_concen      ! number of concentrations     
    real(dp), dimension(:), allocatable, target :: concen_array   ! concentrations     

    real(dp),  parameter :: eps_salt = 0.00001_dp ! if salt concentration below value salt concentration not in output file
    
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

    character(len=100) :: buffer, label
    integer :: pos
    integer :: line
        
    
    !     .. reading in of variables from file
    write(fname,'(A8)')'input.in'
    open(unit=newunit(un_input),file=fname,iostat=ios,status='old')
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        if (present(info)) info = myio_err_inputfile
        return
    endif

    ios=0 
    line = 0

    ! ios<0 : if an end of record condition is encountered or if an end of file condition was detected.  
    ! ios>0 : if an error occured 
    ! ios=0 : otherwise.

    do while (ios == 0)

        read(un_input, '(A)', iostat=ios) buffer

        if (ios == 0) then
    
            line = line + 1

            !  Split buffer into label and data based on first occurence of a whitespace
            
            pos = scan(buffer, '     ')
            label = buffer(1:pos)
            buffer = buffer(pos+1:)

            select case (label) !list-directed The CHARACTER variable is treated as an 'internal file'
            case ('method')    
                read(buffer, *,iostat=ios) method
            case ('sysflag')
                read(buffer, *,iostat=ios) sysflag
            case ('runtype')
                read(buffer, *,iostat=ios) runflag
            case ('bcflag')
                read(buffer, *,iostat=ios) bcflag
            case ('error')
                read(buffer,*,iostat=ios) error 
            case ('infile')
                read(buffer,*,iostat=ios) infile              ! guess  1==yes
            case ('pH%val')    
                read(buffer,*,iostat=ios) pH%val
            case ('pH%min')
                read(buffer,*,iostat=ios) pH%min
            case ('pH%max')
                read(buffer,*,iostat=ios) pH%max
            case ('pH%stepsize')
                read(buffer,*,iostat=ios) pH%stepsize
            case ('pH%delta')
                read(buffer,*,iostat=ios) pH%delta
            case ('KionNa')
                read(buffer,*,iostat=ios) KionNa
            case ('KionK')
                read(buffer,*,iostat=ios) KionK
            case ('sigmasurf')
                read(buffer,*,iostat=ios) sigmasurf
            case ('cNaCl')
                read(buffer,*,iostat=ios) cNaCl
            case ('cKCl')
                read(buffer,*,iostat=ios) cKCl
            case ('cRbCl')
                read(buffer,*,iostat=ios) cRbCl    
            case ('cImCl')
                read(buffer,*,iostat=ios) cImCl    
            case ('cCaCl2')
                read(buffer,*,iostat=ios) cCaCl2   
            case ('cTBCl')
                read(buffer,*,iostat=ios) cTBCl  !   TB=tertraButyl 
            case ('cTMNO3')
                read(buffer,*,iostat=ios) cTMNO3
            case ('cpp')
                read(buffer,*,iostat=ios) cpp
            case ('deltaG0ads')                                        
                read(buffer,*,iostat=ios) deltaG0ads  ! only if sysflag=="electligand")
            case ('deltaG0adsSuOH ')
                read(buffer,*,iostat=ios) deltaG0adsSuOH ! only if bcflag=="pc"
            case ('deltaG0adsSuCl')
                read(buffer,*,iostat=ios) deltaG0adsSuCl ! only if bcflag=="pc"
            case ('deltaG0adsSuNO3')    
                read(buffer,*,iostat=ios) deltaG0adsSuNO3 ! only if bcflag=="pc"
            case ('nsize')
                read(buffer,*,iostat=ios) nsize  
            case ('nrmax')
                read(buffer,*,iostat=ios) nrmax ! only if runflag=="rangnr"    
            case ('nrmin')
                read(buffer,*,iostat=ios) nrmin ! only if runflag=="rangnr"    
            case ('nrstep')
                read(buffer,*,iostat=ios) nrstep ! only if runflag=="rangnr"    
            case ('verboseflag')
                read(buffer,*,iostat=ios)verboseflag 
            case ('radius')
                read(buffer,*,iostat=ios) radius  
            case ('geometry')
                read(buffer,*,iostat=ios) geometry
            case ('delta')
                read(buffer,*,iostat=ios) delta   
            case ('isbulkHCl')
                read(buffer,*,iostat=ios) isbulkHCl  
            case ('isChargeRegularization')
                read(buffer,*,iostat=ios) isChargeRegularization 
            case ('pKa')
                read(buffer,*,iostat=ios) pKa          !   AH   <=> A- + H+    
            case default
                if(pos>1) then 
                    print *, 'Invalid label at line', line  ! empty lines are skipped
                endif
            end select
        endif
    enddo

    if(ios >0 ) then
        print*, 'Error parsing file : iostat =', ios            
        if (present(info)) info = myio_err_inputfile
        return
    endif

    close(un_input)
          
     ! set geometry and bcflag for pa fiber

    if(sysflag=="pafiber") then 
        geometry="cylindrical"
        bcflag="cc" 
    endif    
    ! .. check values of certain input parametere

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


end subroutine read_inputfile
 

subroutine check_value_sysflag(sysflag,info)

    character(len=15), intent(in) :: sysflag
    integer, intent(out),optional :: info

    character(len=15) :: sysflagstr(4)
    integer :: i
    logical :: flag

    ! permissible values of sysflag

    sysflagstr(1)="bulk water"
    sysflagstr(2)="electnopoly"
    sysflagstr(3)="electligand"
    sysflagstr(4)="pafiber"


    flag=.FALSE.

    do i=1,4
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

    character(len=2) :: bcvalues(8)
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
    bcvalues(8)="pc"


    flag=.FALSE.
            
    do i=1,8
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

    if(sysflag=="electnopoly") then
        call output_ligand
        !call output_individualcontr_fe
    elseif(sysflag=="electligand") then
        call output_ligand
    else if(sysflag=="pafiber") then
        call output_pafiber
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
    !use energy
    use surface 
    use myutils, only : newunit
  
    !     .. output file names       
    
    character(len=90) :: sysfilename     
    character(len=90) :: xsolfilename 
    !character(len=90) :: xpolABfilename 
    !character(len=90) :: xpolCfilename 
    !character(len=90) :: xpolendfilename 
    character(len=90) :: xNafilename
    character(len=90) :: xKfilename
    character(len=90) :: xTBfilename
    character(len=90) :: xTMfilename
    character(len=90) :: xppfilename
    character(len=90) :: xppfdisfilename
    character(len=90) :: xCafilename
    character(len=90) :: xNaClfilename
    character(len=90) :: xKClfilename
    character(len=90) :: xClfilename
    character(len=90) :: xNO3filename
    character(len=90) :: potentialfilename
    character(len=90) :: chargefilename
    character(len=90) :: xHplusfilename
    character(len=90) :: xOHminfilename
    !character(len=90) :: densfracAfilename
    !character(len=90) :: densfracBfilename
    character(len=90) :: densfracionpairfilename

    integer :: i,j,k,t     ! dummy indexes
    character(len=100) :: fnamelabel
    character(len=20) :: rstr
    logical :: isopen
    real(dp) :: xppfdis(5),cppfdis(5)
    real(dp) :: cppbulk


    ! .. executable statements 

    ! .. make label filenames 

    if((sysflag=="electnopoly".or.sysflag=="electligand").and.&
        (bcflag=="pp".or.bcflag=="pd".or.bcflag=="pd")) then 

         ! filelabel for qdot only                     
        
        write(rstr,'(F5.3)')sigmaSurf/(4.0_dp*pi*lb*delta)
        fnamelabel="sg"//trim(adjustl(rstr))
        if(cTBCl>=0.001) then 
            write(rstr,'(F5.3)')cTBCl
        else
            write(rstr,'(ES8.2E2)')cTBCl
        endif 
        fnamelabel=trim(fnamelabel)//"cTBCl"//trim(adjustl(rstr))
        if(cTMNO3>=0.001) then 
            write(rstr,'(F5.3)')cTMNO3
        else
            write(rstr,'(ES8.2E2)')cTMNO3
        endif 
        fnamelabel=trim(fnamelabel)//"cTMNO3"//trim(adjustl(rstr))
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
 
        write(rstr,'(F5.3)')cNaCl
        fnamelabel="cNaCl"//trim(adjustl(rstr))
        write(rstr,'(F5.3)')cCaCl2
        fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
        write(rstr,'(F7.3)')pHbulk
        fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))//".dat"
    endif

    sysfilename='system.'//trim(fnamelabel)
    !xpolABfilename='xpolAB.'//trim(fnamelabel)
    !xpolCfilename='xpolC.'//trim(fnamelabel)
    xsolfilename='xsol.'//trim(fnamelabel)
    xNafilename='xNaions.'//trim(fnamelabel)
    xKfilename='xKions.'//trim(fnamelabel)
    xTBfilename='xTBions.'//trim(fnamelabel)
    xTMfilename='xTMions.'//trim(fnamelabel)
    xCafilename='xCaions.'//trim(fnamelabel)
    xNaClfilename='xNaClionpair.'//trim(fnamelabel)
    xKClfilename='xKClionpair.'//trim(fnamelabel)
    xClfilename='xClions.'//trim(fnamelabel)
    xNO3filename='xNO3ions.'//trim(fnamelabel)
    potentialfilename='potential.'//trim(fnamelabel)
    chargefilename='charge.'//trim(fnamelabel)
    xHplusfilename='xHplus.'//trim(fnamelabel)
    xOHminfilename='xOHmin.'//trim(fnamelabel)
    !densfracAfilename='densityAfrac.'//trim(fnamelabel)
    !densfracBfilename='densityBfrac.'//trim(fnamelabel)
    densfracionpairfilename='densityfracionpair.'//trim(fnamelabel)
    xppfilename='xppions.'//trim(fnamelabel)

    !     .. opening files        
    
    open(unit=newunit(un_sys),file=sysfilename)       
    open(unit=newunit(un_xsol),file=xsolfilename)
    open(unit=newunit(un_psi),file=potentialfilename)

    !if(sysflag/="electnopoly") then          
    !    open(unit=newunit(un_xpolAB),file=xpolABfilename)
    !    open(unit=newunit(un_xpolC),file=xpolCfilename)
    !   open(unit=newunit(un_fdisA),file=densfracAfilename) 
    !    open(unit=newunit(un_fdisB),file=densfracBfilename) 
    ! endif   

    if(sysflag=="electligand") open(unit=newunit(un_xpp),file=xppfilename)
      
    if(verboseflag=="yes") then    
        open(unit=newunit(un_xNa),file=xNafilename)
        open(unit=newunit(un_xK),file=xKfilename)
        open(unit=newunit(un_xCa),file=xCafilename)
        open(unit=newunit(un_xTB),file=xTBfilename)
        open(unit=newunit(un_xTM),file=xTMfilename)
        open(unit=newunit(un_xNaCl),file=xNaClfilename)
        open(unit=newunit(un_xKCl),file=xKClfilename)
        open(unit=newunit(un_xpair),file=densfracionpairfilename)
        open(unit=newunit(un_xCl),file=xClfilename)
        open(unit=newunit(un_xNO3),file=xNO3filename)
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

  !  if(sysflag/="electnopoly") then 
  !      do i=1,nr
  !          write(un_xpolAB,fmt4reals)rc(i),xpolAB(i),rhopolA(i),rhopolB(i)
  !          write(un_xpolC,fmt2reals)rc(i),xpolC(i)
  !          write(un_fdisA,fmt6reals)rc(i),fdisA(1,i),fdisA(2,i),fdisA(3,i),fdisA(4,i),fdisA(5,i)        
  !          write(un_fdisB,fmt6reals)rc(i),fdisB(1,i),fdisB(2,i),fdisB(3,i),fdisB(4,i),fdisB(5,i)
  !      enddo
  !  endif   

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
    if(bcflag=="pp".or.bcflag=="pd".or.bcflag=="pc") write(un_sys,*)'vTB         = ',vTB*vsol
    write(un_sys,*)'vNaCl       = ',vNaCl*vsol
    write(un_sys,*)'vKCl        = ',vKCl*vsol
    write(un_sys,*)'cNaCl       = ',cNaCl
    write(un_sys,*)'cKCl        = ',cKCl
    write(un_sys,*)'cCaCl2      = ',cCaCl2
    if(bcflag=="pp".or.bcflag=="pd".or.bcflag=="pc") then 
        write(un_sys,*)'cTBCl       = ',cTBCl
        write(un_sys,*)'cpp         = ',cpp 
        write(un_sys,*)'deltaG0ads  = ',deltaG0ads
        write(un_sys,*)'deltaGads   = ',deltaG0ads-log(Na*vsol/1.0e24_dp)
    endif  
    if(bcflag=="pc")then
        write(un_sys,*)'deltaG0adsSuOH  = ',deltaG0adsSuOH
        write(un_sys,*)'deltaGadsSuOH   = ',deltaG0adsSuOH-log(Na*vsol/1.0e24_dp)
        write(un_sys,*)'deltaG0adsSuCl  = ',deltaG0adsSuCl
        write(un_sys,*)'deltaGadsSuCl   = ',deltaG0adsSuCl-log(Na*vsol/1.0e24_dp)
    endif      
    write(un_sys,*)'pHbulk      = ',pHbulk
!    write(un_sys,*)'pKa         = ',pKa(1)      
!    write(un_sys,*)'pKaNa       = ',pKa(2)
!    write(un_sys,*)'pKaACa      = ',pKa(3)
!    write(un_sys,*)'pKaA2Ca     = ',pKa(4)
!    write(un_sys,*)'pKb         = ',pKb(1)      
!    write(un_sys,*)'pKbNa       = ',pKb(2)
!    write(un_sys,*)'pKbBCa      = ',pKb(3)
!    write(un_sys,*)'pKbB2Ca     = ',pKb(4)
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
    if(bcflag=="pp".or.bcflag=="pd".or.bcflag=="pc") write(un_sys,*)'xbulk%TB    = ',xbulk%TB
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
    !write(un_sys,*)'free energy = ',FE
    !write(un_sys,*)'energy bulk = ',FEbulk 
    !write(un_sys,*)'deltafenergy = ',deltaFE
    write(un_sys,*)'fnorm       = ',fnorm
    !write(un_sys,*)'q residual  = ',qres
    write(un_sys,*)'error       = ',error
    ! write(un_sys,*)'check phi   = ',checkphi 
    !write(un_sys,*)'FEq         = ',FEq 
    !write(un_sys,*)'FEpi        = ',FEpi
    !write(un_sys,*)'FErho       = ',FErho
    !write(un_sys,*)'FEel        = ',FEel
    !write(un_sys,*)'FEelsurf    = ',FEelsurf
    !write(un_sys,*)'FEbind      = ',FEbind
    !write(un_sys,*)'FEVdW       = ',FEVdW 
    !write(un_sys,*)'FEalt       = ',FEalt
    
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
    else if(bcflag=='pc') then   
        do i=1,8   
            write(un_sys,fmt)'fdisSu(',i,')   = ',fdisS(i)
        enddo
        do i=1,8   
            write(un_sys,fmt)'gdisSu(',i,')   = ',gdisS(i)
        enddo   
        write(un_sys,*)'fdisR      = ',fdisR
        write(un_sys,*)'sigmaR     = ',fdisR*sigmaSurf/(4.0_dp*pi*lb*delta)
        write(un_sys,*)'sigmaLR    = ',(1.0_dp-fdisS(Su)-fdisS(SuOH)-fdisS(SuCl))*sigmaSurf/(4.0_dp*pi*lb*delta)
    else
        do i=1,6   
            write(un_sys,fmt)' fdisSu(',i,')  = ',fdisS(i)
        enddo  
    endif
    write(un_sys,*)'nsize       = ',nsize  
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



subroutine output_ligand
  
    !     .. variables and constant declaractions
    use globals 
    use volume
    use parameters
    use field
   ! use energy
    use surface 
    use myutils, only : newunit
  
    !     .. output file names       
    
    character(len=90) :: sysfilename  
    character(len=90) :: xsolfilename   
    character(len=90) :: xNafilename
    character(len=90) :: xKfilename
    character(len=90) :: xTBfilename
    character(len=90) :: xTMfilename
    character(len=90) :: xppfilename
    character(len=90) :: xppfdisfilename
    character(len=90) :: xCafilename
    character(len=90) :: xClfilename
    character(len=90) :: xNO3filename
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
    real(dp) :: eta, sigmaSurf0


    ! .. executable statements 

    ! .. make label filenames 

    sigmaSurf0=sigmaSurf/(4.0_dp*pi*lb*delta)  
    ! filelabel for qdot only                     
    if(sigmaSurf0<10.0_dp) then 
        write(rstr,'(F5.3)')sigmaSurf0    
    else    
        write(rstr,'(F6.3)')sigmaSurf0
    endif
        
    fnamelabel="sg"//trim(adjustl(rstr))
    
    if(cTBCl>=0.001) then 
        write(rstr,'(F5.3)')cTBCl
    else
        write(rstr,'(ES8.2E2)')cTBCl
    endif 
    fnamelabel=trim(fnamelabel)//"cTBCl"//trim(adjustl(rstr))
    if(cTMNO3>=0.001) then 
        write(rstr,'(F5.3)')cTMNO3
    else
        write(rstr,'(ES8.2E2)')cTMNO3
    endif 
    fnamelabel=trim(fnamelabel)//"cTMNO3"//trim(adjustl(rstr))
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
    xTMfilename='xTMions.'//trim(fnamelabel)
    xCafilename='xCaions.'//trim(fnamelabel)
    xClfilename='xClions.'//trim(fnamelabel)
    xNO3filename='xNO3ions.'//trim(fnamelabel)
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
        open(unit=newunit(un_xTM),file=xTMfilename)
        open(unit=newunit(un_xCl),file=xClfilename)
        open(unit=newunit(un_xNO3),file=xNO3filename)
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
            write(un_xTM,*)rc(i),xTM(i) 
            write(un_xCl,*)rc(i),xCl(i)
            write(un_xNO3,*)rc(i),xNO3(i)
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
    if(bcflag=="pd".or.bcflag=="pc")then 
        write(un_sys,*)'vpp(AH2BH)  = ',vpp(AH2BH)*vsol
        write(un_sys,*)'vpp(AHBH)   = ',vpp(AHBH)*vsol
        write(un_sys,*)'vpp(AHB)    = ',vpp(AHB)*vsol
        write(un_sys,*)'vpp(ABH)    = ',vpp(ABH)*vsol
        write(un_sys,*)'vpp(AB)     = ',vpp(AB)*vsol
    endif    
    if(bcflag=="pp".or.bcflag=="pd".or.bcflag=="pc") write(un_sys,*)'vTB         = ',vTB*vsol
    write(un_sys,*)'vNaCl       = ',vNaCl*vsol
    write(un_sys,*)'vKCl        = ',vKCl*vsol
    write(un_sys,*)'cNaCl       = ',cNaCl
    write(un_sys,*)'cKCl        = ',cKCl
    write(un_sys,*)'cCaCl2      = ',cCaCl2
    if(bcflag=="pp".or.bcflag=="pd".or.bcflag=="pc") then 
        write(un_sys,*)'cTBCl       = ',cTBCl
        write(un_sys,*)'cTMNO3      = ',cTMNO3
        write(un_sys,*)'cpp         = ',cpp
        write(un_sys,*)'deltaG0ads  = ',deltaG0ads
        write(un_sys,*)'deltaGads   = ',deltaG0ads-log(Na*vsol/1.0e24_dp)
    endif    
    if(bcflag=="pc")then
        write(un_sys,*)'deltaG0adsSuOH  = ',deltaG0adsSuOH
        write(un_sys,*)'deltaGadsSuOH   = ',deltaG0adsSuOH-log(Na*vsol/1.0e24_dp)
        write(un_sys,*)'deltaG0adsSuCl  = ',deltaG0adsSuCl
        write(un_sys,*)'deltaGadsSuCl   = ',deltaG0adsSuCl-log(Na*vsol/1.0e24_dp)
    endif    
    write(un_sys,*)'pHbulk      = ',pHbulk
    write(un_sys,*)'rhoqbulk    = ',rhoqbulk
    write(un_sys,*)'xbulk%sol   = ',xbulk%sol
    write(un_sys,*)'xbulk%Na    = ',xbulk%Na
    write(un_sys,*)'xbulk%Cl    = ',xbulk%Cl
    write(un_sys,*)'xbulk%K     = ',xbulk%K
    write(un_sys,*)'xbulk%NaCl  = ',xbulk%NaCl
    write(un_sys,*)'xbulk%KCl   = ',xbulk%KCl
    write(un_sys,*)'xbulk%Ca    = ',xbulk%Ca
    write(un_sys,*)'xbulk%Hplus = ',xbulk%Hplus
    write(un_sys,*)'xbulk%OHmin = ',xbulk%OHmin
    if(bcflag=="pp".or.bcflag=="pd".or.bcflag=="pc") then 
        write(un_sys,*)'xbulk%TB    = ',xbulk%TB
        write(un_sys,*)'xbulk%TM    = ',xbulk%TM
        write(un_sys,*)'xbulk%NO3   = ',xbulk%NO3
    endif    
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
    write(un_sys,*)'zTB         = ',zTB
    write(un_sys,*)'zTM         = ',zTM
    
    
    if(bcflag=="pd".or.bcflag=="pc") then 
        write(un_sys,*)'zpp(AH2BH)  = ',zpp(AH2BH)
        write(un_sys,*)'zpp(AHBH)   = ',zpp(AHBH)
        write(un_sys,*)'zpp(AHB)    = ',zpp(AHB)
        write(un_sys,*)'zpp(ABH)    = ',zpp(ABH)
        write(un_sys,*)'zpp(AB)     = ',zpp(AB)
    endif    
    write(un_sys,*)'nr          = ',nr
    !write(un_sys,*)'free energy = ',FE
    !write(un_sys,*)'energy bulk = ',FEbulk 
    !write(un_sys,*)'deltafenergy = ',deltaFE
    write(un_sys,*)'fnorm       = ',fnorm
    !write(un_sys,*)'q residual  = ',qres
    write(un_sys,*)'error       = ',error
    ! write(un_sys,*)'FEpi        = ',FEpi
    ! write(un_sys,*)'FErho       = ',FErho
    ! write(un_sys,*)'FEel        = ',FEel
    ! write(un_sys,*)'FEelsurf    = ',FEelsurf
    ! write(un_sys,*)'FEbind      = ',FEbind
    ! write(un_sys,*)'FEalt       = ',FEalt
    ! write(un_sys,*)'sigmaSurf   = ',sigmaSurf/(4.0_dp*pi*lb*delta)
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
    else if(bcflag=='pc') then   
        do i=1,9   
            write(un_sys,fmt)'fdisSu(',i,')   = ',fdisS(i)
        enddo
        do i=1,8   
            write(un_sys,fmt)'gdisSu(',i,')   = ',gdisS(i)
        enddo   
        write(un_sys,*)'fdisR      = ',fdisR
        write(un_sys,*)'sigmaR     = ',fdisR*sigmaSurf/(4.0_dp*pi*lb*delta)
        write(un_sys,*)'sigmaLR    = ',(1.0_dp-fdisS(Su)-fdisS(SuOH)-fdisS(SuCl)-fdisS(SuNO3))&
            *sigmaSurf/(4.0_dp*pi*lb*delta)
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
        close(un_xTM)
        close(un_xNO3)
    endif
        

end subroutine output_ligand

subroutine output_pafiber
  
    !     .. variables and constant declaractions
    use globals 
    use volume
    use parameters
    use field
    use surface 
    use myutils, only : newunit
  
    !     .. output file names       
    
    character(len=90) :: sysfilename     
    character(len=90) :: xsolfilename 
    character(len=90) :: xNafilename
    character(len=90) :: xRbfilename
    character(len=90) :: xImfilename
    character(len=90) :: xKfilename
    character(len=90) :: xCafilename
    character(len=90) :: xNaClfilename
    character(len=90) :: xKClfilename
    character(len=90) :: xClfilename
    character(len=90) :: xHplusfilename
    character(len=90) :: xOHminfilename
    character(len=90) :: potentialfilename
    character(len=90) :: chargefilename
    character(len=90) :: densfracionpairfilename
    character(len=90) :: fdispafilename

    integer :: i,j,k,t     ! dummy indexes
    character(len=100) :: fnamelabel
    character(len=20) :: rstr
    logical :: isopen


    ! .. executable statements 

    ! .. make label filenames 
    
    fnamelabel=""

    if(cNaCl>eps_salt) then 
        write(rstr,'(F5.3)')cNaCl
        fnamelabel="cNaCl"//trim(adjustl(rstr))
    endif

    if(cKCl>eps_salt) then 
        write(rstr,'(F5.3)')cKCl
        fnamelabel=trim(fnamelabel)//"cKCl"//trim(adjustl(rstr))
    endif 
    
    if(cCaCl2>eps_salt) then 
        write(rstr,'(F5.3)')cCaCl2
        fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
    endif 

    if(cRbCl>eps_salt) then 
        write(rstr,'(F5.3)')cRbCl
        fnamelabel=trim(fnamelabel)//"cRbCl"//trim(adjustl(rstr))
    endif

    if(cImCl>eps_salt) then 
        write(rstr,'(F5.3)')cImCl
        fnamelabel=trim(fnamelabel)//"cImCl"//trim(adjustl(rstr))
    endif 

    write(rstr,'(F7.3)')pHbulk
    fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))//".dat"


    sysfilename='system.'//trim(fnamelabel)
    xsolfilename='xsol.'//trim(fnamelabel)
    xNafilename='xNaions.'//trim(fnamelabel)
    xRbfilename='xRbions.'//trim(fnamelabel)
    xImfilename='xImions.'//trim(fnamelabel)
    xKfilename='xKions.'//trim(fnamelabel)
    xCafilename='xCaions.'//trim(fnamelabel)
    xNaClfilename='xNaClionpair.'//trim(fnamelabel)
    xKClfilename='xKClionpair.'//trim(fnamelabel)
    xClfilename='xClions.'//trim(fnamelabel)
    potentialfilename='potential.'//trim(fnamelabel)
    chargefilename='charge.'//trim(fnamelabel)
    xHplusfilename='xHplus.'//trim(fnamelabel)
    xOHminfilename='xOHmin.'//trim(fnamelabel)
    densfracionpairfilename='densityfracionpair.'//trim(fnamelabel)
    fdispafilename='fdispa.'//trim(fnamelabel)
    
    !     .. opening files        
    
    open(unit=newunit(un_sys),file=sysfilename)       
    open(unit=newunit(un_xsol),file=xsolfilename)
    open(unit=newunit(un_psi),file=potentialfilename)
    

    if(verboseflag=="yes") then   
        open(unit=newunit(un_xNa),file=xNafilename)
        open(unit=newunit(un_xRb),file=xRbfilename)
        open(unit=newunit(un_xIm),file=xImfilename)
        open(unit=newunit(un_xK),file=xKfilename)
        open(unit=newunit(un_xCa),file=xCafilename)
        open(unit=newunit(un_xNaCl),file=xNaClfilename)
        open(unit=newunit(un_xKCl),file=xKClfilename)
        open(unit=newunit(un_xpair),file=densfracionpairfilename)
        open(unit=newunit(un_xCl),file=xClfilename)
        open(unit=newunit(un_charge),file=chargefilename)
        open(unit=newunit(un_xHplus),file=xHplusfilename)
        open(unit=newunit(un_xOHmin),file=xOHminfilename)
        open(unit=newunit(un_fdisA),file=fdispafilename)  
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
   
    
    if(verboseflag=="yes") then 
        do i=1,nr
            write(un_xNa,*)rc(i),xNa(i)
            write(un_xRb,*)rc(i),xRb(i)
            write(un_xIm,*)rc(i),xIm(i)
            write(un_xK,*)rc(i),xK(i)
            write(un_xCa,*)rc(i),xCa(i)
            write(un_xNaCl,*)rc(i),xNaCl(i)
            write(un_xKCl,*)rc(i),xKCl(i)
            write(un_xpair,*)rc(i),(xNaCl(i)/vNaCl)/(xNa(i)/vNa+xCl(i)/vCl+xNaCl(i)/vNaCl)
            write(un_xCl,*)rc(i),xCl(i)
            write(un_charge,*)rc(i),rhoq(i)
            write(un_xHplus,*)rc(i),xHplus(i)
            write(un_xOHmin,*)rc(i),xOHmin(i)    
            write(un_fdisA,*)rc(i),fdispa(i)    
        enddo    
    endif

    write(un_sys,*)'system      = pafiber'
    write(un_sys,*)'version     = ',VERSION
    write(un_sys,*)'sysflag     = ',sysflag
    write(un_sys,*)'bcflag      = ',bcflag
    write(un_sys,*)'delta       = ',delta  
    write(un_sys,*)'vsol        = ',vsol
    write(un_sys,*)'vNa         = ',vNa*vsol
    write(un_sys,*)'vCl         = ',vCl*vsol
    write(un_sys,*)'vCa         = ',vCa*vsol
    write(un_sys,*)'vK          = ',vK*vsol
    write(un_sys,*)'vRb         = ',vRb*vsol
    write(un_sys,*)'vIm         = ',vIm*vsol
    write(un_sys,*)'vNaCl       = ',vNaCl*vsol
    write(un_sys,*)'vKCl        = ',vKCl*vsol
    write(un_sys,*)'cNaCl       = ',cNaCl
    write(un_sys,*)'cRbCl       = ',cRbCl
    write(un_sys,*)'cImCl       = ',cImCl
    write(un_sys,*)'cKCl        = ',cKCl
    write(un_sys,*)'cRbCl       = ',cRbCl
    write(un_sys,*)'cImCl       = ',cImCl
    write(un_sys,*)'cCaCl2      = ',cCaCl2
    write(un_sys,*)'pHbulk      = ',pHbulk
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
    write(un_sys,*)'dielectW    = ',dielectW
    write(un_sys,*)'lb          = ',lb
    write(un_sys,*)'T           = ',Temp
   
    write(un_sys,*)'zNa         = ',zNa
    write(un_sys,*)'zCa         = ',zCa
    write(un_sys,*)'zK          = ',zK
    write(un_sys,*)'zRb         = ',zRb
    write(un_sys,*)'zIm         = ',zIm
    write(un_sys,*)'zpa         = ',zpa
    write(un_sys,*)'nr          = ',nr

    !write(un_sys,*)'free energy = ',FE
    !write(un_sys,*)'energy bulk = ',FEbulk 
    !write(un_sys,*)'deltafenergy = ',deltaFE
    
    write(un_sys,*)'fnorm       = ',fnorm
    write(un_sys,*)'totalcharge = ',totalcharge
    write(un_sys,*)'totalEpa    = ',totalEpa
    write(un_sys,*)'avfdispa    = ',avfdispa
     
   
    !write(un_sys,*)'q residual  = ',qres
    
    write(un_sys,*)'error       = ',error
    ! write(un_sys,*)'check phi   = ',checkphi 
    !write(un_sys,*)'FEq         = ',FEq 
    !write(un_sys,*)'FEpi        = ',FEpi
    !write(un_sys,*)'FErho       = ',FErho
    !write(un_sys,*)'FEel        = ',FEel
    !write(un_sys,*)'FEelsurf    = ',FEelsurf
    !write(un_sys,*)'FEbind      = ',FEbind
    !write(un_sys,*)'FEVdW       = ',FEVdW 
    !write(un_sys,*)'FEalt       = ',FEalt
    
    write(un_sys,*)'sigmaSurf   = ',sigmaSurf/(4.0_dp*pi*lb*delta)
    write(un_sys,*)'sigmaqSurf  = ',sigmaqSurf/(4.0_dp*pi*lb*delta)
    write(un_sys,*)'psiSurf     = ',psiSurf
    if(bcflag=='ta') then
        do i=1,4   
            write(un_sys,fmt)'fdisTa(',i,')   = ',fdisTaL(i)
        enddo  
    else if(bcflag=='pp'.or.bcflag=='pd'.or.bcflag=='pc') then
        print*,"output_pafiber: wrong bcflag : ",bcflag   
    else
        do i=1,6   
            write(un_sys,fmt)' fdisSu(',i,')  = ',fdisS(i)
        enddo  
    endif
    write(un_sys,*)'nsize       = ',nsize  
    write(un_sys,*)'geometry    = ',geometry  
    write(un_sys,*)'iterations  = ',iter
    write(un_sys,*)'radius      = ',radius 
    
   
    ! .. closing files

    close(un_sys)
    close(un_xsol)
    close(un_psi)

    
    if(verboseflag=="yes") then 
        close(un_xNa)
        close(un_xRb)
        close(un_xIm)
        close(un_xK)
        close(un_xCa)
        close(un_xNaCl)
        close(un_xKCl)
        close(un_xpair)
        close(un_xCl)
        close(un_charge)
        close(un_xHplus)
        close(un_xOHmin)
        close(un_fdisA)
    endif
        

end subroutine output_pafiber


! subroutine output_individualcontr_fe

!     use globals, only : sysflag
!     use energy
!     use myutils, only : newunit
!     use parameters, only : sigmaAB,cNaCl,cCaCl2,pHbulk,VdWepsB
!     use volume, only : delta,nr,nrmax,nrmin

!     ! local arguments

!     integer :: un_fe

!     character(len=100) :: fenergyfilename   
!     character(len=100) :: fnamelabel
!     character(len=20) :: rstr

!    !     .. make label filename

!     if(sysflag=="elect".or.sysflag=="electdouble".or.sysflag=="electnopoly") then 
!         write(rstr,'(F5.3)')sigmaAB*delta 
!         fnamelabel="sg"//trim(adjustl(rstr)) 
!         write(rstr,'(F5.3)')cNaCl
!         fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))
!         write(rstr,'(F5.3)')cCaCl2
!         fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
!         write(rstr,'(F7.3)')pHbulk
!         fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))//".dat"
!     elseif(sysflag=="neutral") then 
!         write(rstr,'(F5.3)')sigmaAB*delta 
!         fnamelabel="sg"//trim(adjustl(rstr)) 
!         write(rstr,'(F5.3)')VdWepsB
!         fnamelabel=trim(fnamelabel)//"VdWepsB"//trim(adjustl(rstr))//".dat"
!     else
!         print*,"Error in output_individualcontr_fe subroutine"
!         print*,"Wrong value sysflag : ", sysflag
!     endif    

!     fenergyfilename='energy.'//trim(fnamelabel)   
        
!     !     .. opening files        

!     open(unit=newunit(un_fe),file=fenergyfilename) 

!     write(un_fe,*)'FE              = ',FE  
!     write(un_fe,*)'FEbulk          = ',FEbulk 
!     write(un_fe,*)'deltaFE         = ',deltaFE
!     write(un_fe,*)'FEalt           = ',FEalt  
!     write(un_fe,*)'FEbulkalt       = ',FEbulkalt 
!     write(un_fe,*)'deltaFEalt      = ',deltaFEalt
    
!     write(un_fe,*)"FEtrans%sol     = ",FEtrans%sol   
!     write(un_fe,*)"FEtrans%Na      = ",FEtrans%Na  
!     write(un_fe,*)"FEtrans%Cl      = ",FEtrans%Cl  
!     write(un_fe,*)"FEtrans%Ca      = ",FEtrans%Ca  
!     write(un_fe,*)"FEtrans%K       = ",FEtrans%K
!     write(un_fe,*)"FEtrans%KCl     = ",FEtrans%KCl
!     write(un_fe,*)"FEtrans%NaCl    = ",FEtrans%NaCl  
!     write(un_fe,*)"FEtrans%Hplus   = ",FEtrans%Hplus  
!     write(un_fe,*)"FEtrans%OHmin   = ",FEtrans%OHmin  
    
!     write(un_fe,*)"FEchempot%Na    = ",FEchempot%Na
!     write(un_fe,*)"FEchempot%Cl    = ",FEchempot%Cl
!     write(un_fe,*)"FEchempot%Ca    = ",FEchempot%Ca
!     write(un_fe,*)"FEchempot%K     = ",FEchempot%K
!     write(un_fe,*)"FEchempot%KCl   = ",FEchempot%KCl
!     write(un_fe,*)"FEchempot%NaCl  = ",FEchempot%NaCl
!     write(un_fe,*)"FEchempot%Hplus = ",FEchempot%Hplus
!     write(un_fe,*)"FEchempot%OHmin = ",FEchempot%OHmin

!     write(un_fe,*)"FEchemsurf      = ",FEchemsurf
!     write(un_fe,*)"FEchemsurfalt   = ",FEchemsurfalt
    
!     write(un_fe,*)"delta FEchemsurfalt= ",FEchemsurfalt-FEchemsurf-diffFEchemsurf
    
!     close(un_fe)


! end subroutine   output_individualcontr_fe

end module


