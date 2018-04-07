! ---------------------------------------------------------------|
! Solves the SCMFT eqs planar/spherical/(inv)cylidrical surface, |
! coated with weak poyelectrolytes and or have a surface charge  |
! Different surface charge (bcflags)                             |
! 	 1) "qu" = quartz                                            |
!    2) "cl" = clay                                              |
!    3) "ca" = calcite                                           |
!    4) "ta" = taurine                                           | 
!    5) "cc" = constant charge                                   | 
!    6) "pp" = phoshonatepropionate                              |  
!    input/output: see myio.f90                                  | 
! ---------------------------------------------------------------|
      
program brushweakpolyelectrolyte 
    
    !     .. variable and constant declaractions 
    use globals       ! parameters definitions 
    use physconst
    use mathconst
    use volume
    use random
    use field
    use parameters
    use matrices
    use energy
    use chains
    use VdW
    use listfcn
    use initxvector
    use surface
    use myio
    use myutils
    use chaingenerator
    use fcnpointer
    
    implicit none  
    
    real(dp),  dimension(:), allocatable :: x         ! iteration vector 
    real(dp),  dimension(:), allocatable :: xguess    ! guess iteration vector
    real(dp),  dimension(:), allocatable :: xstored   ! stored iteration vector
    real(dp),  dimension(:), allocatable :: fvec      ! stored iteration vector


    integer :: i             ! dummy indices       
    logical :: use_xstored       
    logical :: isfirstguess   
    character(len=lenText) :: text
    character(len=20) :: rstr
    
    ! .. executable statements 
    ! .. init 

    LogName='status.log'
    call open_logfile(LogUnit,LogName) 
    text='program begins'
    call print_to_log(LogUnit,text)


    call read_inputfile()
    call init_constants()
    call init_matrices()            ! init matrices for chain generation
    call allocate_chains(cuantasAB,nsegAB,cuantasC,nsegC)  
    call make_sequence_chain(period,chaintype)
    call set_properties_chain(period,chaintype)  
    call make_chains(chainmethod)   ! generate polymer configurations 
   

    call allocate_geometry(nsize)
    call make_geometry()            ! generate volume elements lattice 
    call allocate_field(nsize) 

    !    call read_VdWCoeff() ! THIS NEED TO BE CHANGED !!!
    
    call set_size_neq()             ! number of non-linear equation neq    
    call init_expmu()

    call init_surface(bcflag)

    !  .. computation starts
    
    ! nr=nrmax                        
    allocate(xstored(neq))
    allocate(x(neq))
    allocate(xguess(neq))   
    allocate(fvec(neq))   
    

    isfirstguess = .true.    
    use_xstored = .false.             
    iter = 0
    
    pH%val=pH%min


    ! print*,"-main->sysflag=",sysflag

    if(runflag=="rangepH") then 

        !  .. first increase pH value

        do while (pH%min<=pH%val.and.pH%val<=pH%max.and.(abs(pH%stepsize)>=pH%delta)) 
           
            call init_expmu()
            call make_guess(x, xguess, isfirstguess) 
            call solver(x, xguess, error, fnorm) 
            
            if(isNaN(fnorm)) then  
                text="no solution: backstep"
                call print_to_log(LogUnit,text)
                pH%stepsize=pH%stepsize/2.0_dp ! smaller 
                pH%val=pH%val-pH%stepsize ! step back
                do i=1,neq
                    x(i)=xguess(i)
                enddo       
            else 
                ! call fcnenergy()        
                ! call average_height()      
                ! call charge_polymer()
                ! call average_charge_polymer()
                call output()           ! writing of output

                write(rstr,'(F7.3)')pH%val
                text="solution pH="//trim(adjustl(rstr))
                pH%val=pH%val+pH%stepsize
            endif 
            isfirstguess= .false.
            iter  = 0              ! reset of iteration counter 

        enddo 

    endif



    deallocate(xstored)
    call deallocate_field()

    text="program end"
    call print_to_log(LogUnit,text)
    call close_logfile(LogUnit)

end program brushweakpolyelectrolyte
