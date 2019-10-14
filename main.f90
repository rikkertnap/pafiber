! ---------------------------------------------------------------|
! Solves the SCMFT eqs in planar/spheric/cylindrical geometry    |
! surface NP/colloid/fiber/ has a surface charge                 |
! Different surface charge (bcflags)                             |
! 	 1) "qu" = quartz                                            |
!    2) "cl" = clay                                              |
!    3) "ca" = calcite                                           |
!    4) "ta" = taurine                                           | 
!    5) "cc" = constant charge                                   | 
!    6) "pp" = phoshonatepropionate                              |  
!    7) "pd" ==pp ligand in solution dynamics adsorption         |
!    input/output: see myio.f90                                  | 
! ---------------------------------------------------------------|
      
program main
    
    !     .. variable and constant declaractions 
    use globals       ! parameters definitions 
    use physconst
    use mathconst
    use volume
    use random
    use field
    use parameters
    !use energy
    use listfcn
    use initxvector
    use surface
    use myio
    use myutils
    use fcnpointer
    
    implicit none  
    
    real(dp),  dimension(:), allocatable :: x         ! iteration vector 
    real(dp),  dimension(:), allocatable :: xguess    ! guess iteration vector
    real(dp),  dimension(:), allocatable :: xstored   ! stored iteration vector
    real(dp),  dimension(:), allocatable :: fvec      ! stored iteration vector


    integer :: i, c, info             ! dummy indices       
    logical :: use_xstored       
    logical :: isfirstguess 
    logical :: issolution  

    character(len=lenText) :: text
    character(len=20) :: rstr
    
    real(dp), dimension(:),  pointer :: list
    real(dp), pointer :: list_val
    real(dp), target :: list_first


    ! .. executable statements 
    ! .. init 

    LogName='status.log'
    call open_logfile(LogUnit,LogName) 
    text='program begins'
    call print_to_log(LogUnit,text)


    call read_inputfile()
    call init_constants()
    call allocate_geometry(nsize)
    call make_geometry()            ! generate volume elements lattice 
    call allocate_field(nsize) 
    call set_size_neq()             ! number of non-linear equation neq    
    call init_expmu()
    call init_surface(bcflag)
    
    call set_fcn

    if(sysflag=="pafiber".or.sysflag=="pafiberIm".or.sysflag=="pafiberborn".or.sysflag=="pafibervarelec") then 
        call init_xpa_volume_dist
        call init_rhoEpa_dist
    endif   

    !  .. computation starts
                     
    ! .. select variable with which list_array associated
    
    if (runflag=="rangepHcpp"  .or. runflag=="rangepHcNaCl" .or. runflag=="rangepHcRbCl") then
        call set_value_concen(runflag,info)
        if(info/=0) then
            print*,"Error in input file: info = ",info," : end program." 
            stop
        endif
        list=>concen_array 
        if(runflag=="rangepHcpp")   list_val => cpp    
        if(runflag=="rangepHcNaCl") list_val => cNaCl
        if(runflag=="rangepHcRbCl") list_val => cRbCl
                
    else
        if(associated(list)) nullify(list) 
        if(associated(list_val)) nullify(list_val) 
        num_concen=1
        allocate(list(1))
        list_val => list_first ! need to point to a valid variable
    endif    

    if(runflag/="rangenr") then  

        allocate(xstored(neq))
        allocate(x(neq))
        allocate(xguess(neq))   
        allocate(fvec(neq))   
    
        isfirstguess = .true.    
        use_xstored = .false.             
        iter = 0
        list_first = list(1)

        do c=1,num_concen        ! loop 

            iter = 0                  ! iteration counter 
            list_val = list(c)
            isfirstguess = .true.
            pH%val = pH%min    

            do while (pH%min<=pH%val.and.pH%val<=pH%max.and.(abs(pH%stepsize)>=pH%delta)) 
                          
                call init_expmu()
                call make_guess(x, xguess, isfirstguess)     
                call solver(x, xguess, error, fnorm, issolution)
                call fcnptr(x,fvec,neq)
                
                if(isSolution) then
                    ! call fcnenergy()        
                    totalcharge=total_charge(rhoq,sigmaqSurf)
                    avfdispa=average_charge_pa() 
                    avfdisA=average_charge_pa_Ca()
                    call output()           ! writing of output
                    write(rstr,'(F7.3)')pH%val
                    text="solution pH="//trim(adjustl(rstr))
                    pH%val=pH%val+pH%stepsize
                else
                     call output() 
                    text="no solution: backstep"
                    call print_to_log(LogUnit,text)
                    pH%stepsize=pH%stepsize/2.0_dp  ! smaller 
                    pH%val=pH%val-pH%stepsize       ! sttep back

                    do i=1,neqint
                        x(i)=xguess(i)
                    enddo       
                endif 

                isfirstguess= .false.
                iter  = 0              ! reset of iteration counter 

            enddo 

        enddo

        deallocate(x)
        deallocate(xguess)   
        deallocate(fvec)   

    else  ! runflag==rangenr

        allocate(xstored(neq))    
      
        nr=nrmax
        isfirstguess = .true.    
        use_xstored = .false.
        iter = 0 

        do while (nr>=nrmin)        ! loop distances
             
          
            call set_size_neq()  
            
            allocate(x(neq))
            allocate(xguess(neq))

            call make_guess(x, xguess, isfirstguess, use_xstored, xstored)
            call solver(x, xguess, error, fnorm)
            call output()           ! writing of output

            isfirstguess =.false.    
            use_xstored = .true.
            iter = 0                ! reset of iteration counter 
            nr = nr-nrstep          ! reduce distance 
            do i=1,neqint
                xstored(i)=x(i)
            enddo

            deallocate(x)   
            deallocate(xguess)
        enddo  

    endif    
 
    deallocate(xstored)

    call deallocate_field()

    text="program end"
    call print_to_log(LogUnit,text)
    call close_logfile(LogUnit)

end program main
