module initxvector

    use globals
    use parameters
    use field
    use volume 
    use surface

    implicit none
      
contains



subroutine init_guess_electnopoly(x, xguess)
      
    implicit none
  
    real(dp), intent(inout) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp), intent(out) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: i
    character(len=8) :: fname(2)
    integer :: ios,nfile(4)
    integer :: neq_bc 
  
  
    ! .. init guess all xbulk     

    do i=1,nr
        x(i)=xbulk%sol
        x(i+nr)=0.000_dp
    enddo

    if(sysflag=="pafiber") then 
        do i=1,18
            x(i)=abs(1.0_dp-xpa(i))
            x(i+nr)=0.000_dp
        enddo
    endif    

    neq_bc=0
    if(bcflag/="cc") then
        neq_bc=1 
        x(2*nr+neq_bc)=0.00_dp
    endif        

    if (infile.eq.1) then   ! infile is read in from file/stdio  
    
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'psi.in'
     
        nfile(1)=100
        nfile(2)=200
     
        do i=1,2 ! loop files
            open(unit=nfile(i),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then    
                print*, 'file num ber =',nfile(i),' file name =',fname(i)
                print*, 'Error opening file : iostat =', ios
                stop
            endif
        enddo
     
        if(bcflag/="cc") read(200,*)psisurf     ! surface potential 
        do i=1,nr
            read(100,*)xsol(i)    ! solvent
            read(200,*)psi(i)     ! potential 
            x(i)      = xsol(i)    ! placing xsol  in vector x
            x(i+nr)   = psi(i)     ! placing xsol  in vector x
        enddo

        do i=1,2
            close(nfile(i))
        enddo

    endif
    !     .. end init from file 
  
    do i=1,neqint
        xguess(i)=x(i)
    enddo

end subroutine init_guess_electnopoly
!     purpose: initalize x and xguess


subroutine init_guess_pafiberIm(x, xguess)
      
    use parameters, only : xbulk
      
    implicit none
  
    real(dp), intent(inout) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp), intent(out) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: i
    character(len=8) :: fname(2)
    integer :: ios,nfile(4)
    integer :: neq_bc 
    real(dp) :: rhoIm_bulk
  
    ! .. init guess all xbulk     


    rhoIm_bulk=xbulk%Im/(vIm*vsol)

    do i=1,nr
        x(i)=xbulk%sol
        x(i+nr)=0.000_dp
        x(i+2*nr)=rhoIm_bulk
    enddo

    if(sysflag=="pafiber") then 
        do i=1,18
            x(i)=abs(1.0_dp-xpa(i))
            x(i+nr)=0.000_dp
        enddo
    endif    

    neq_bc=0
    if(bcflag/="cc") then
        neq_bc=1 
        x(2*nr+neq_bc)=0.00_dp
    endif        

    if (infile.eq.1) then   ! infile is read in from file/stdio  
    
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'psi.in'
     
        nfile(1)=100
        nfile(2)=200
     
        do i=1,2 ! loop files
            open(unit=nfile(i),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then    
                print*, 'file num ber =',nfile(i),' file name =',fname(i)
                print*, 'Error opening file : iostat =', ios
                stop
            endif
        enddo
     
        if(bcflag/="cc") read(200,*)psisurf     ! surface potential 
        do i=1,nr
            read(100,*)xsol(i)    ! solvent
            read(200,*)psi(i)     ! potential 
            x(i)      = xsol(i)    ! placing xsol  in vector x
            x(i+nr)   = psi(i)     ! placing xsol  in vector x
        enddo

        do i=1,2
            close(nfile(i))
        enddo

    endif
    !     .. end init from file 
  
    do i=1,neqint
        xguess(i)=x(i)
    enddo

end subroutine init_guess_pafiberIm



subroutine init_guess_pafiberborn(x, xguess)
      
    use parameters, only : xbulk
      
    implicit none
  
    real(dp), intent(inout) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp), intent(out) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: i
    character(len=8) :: fname(5)
    integer :: ios,nfile(5)
    integer :: neq_bc 
    real(dp) :: xAA,xAACa
  
    ! .. init guess all xbulk     

    do i=1,nr
        x(i)=xbulk%sol
        x(i+nr)=0.0_dp
        x(i+2*nr)=0.0_dp
        x(i+3*nr)=0.0_dp
    !    x(i+4*nr)=xbulk%Im/vIm
    enddo

    
    do i=1,18
        x(i)=abs(1.0_dp-xpa(i)) 
    enddo
        

    neq_bc=0
    if(bcflag/="cc") then
        neq_bc=1 
        x(5*nr+neq_bc)=0.00_dp
    endif        

    if (infile.eq.1) then   ! infile is read in from file/stdio  
    
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'psi.in'
        write(fname(3),'(A6)')'xAA.in'
        write(fname(4),'(A8)')'xAACa.in' 
        !write(fname(5),'(A8)')'rhoIm.in'
        
        nfile(1)=100
        nfile(2)=200
        nfile(3)=300
        nfile(4)=400
        !nfile(5)=500
     
        do i=1,4 ! loop files
            open(unit=nfile(i),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then    
                print*, 'file num ber =',nfile(i),' file name =',fname(i)
                print*, 'Error opening file : iostat =', ios
                stop
            endif
        enddo
     
        if(bcflag/="cc") read(200,*)psisurf     ! surface potential 
        do i=1,nr
            read(100,*)xsol(i)    ! solvent
            read(200,*)psi(i)     ! potential 
            read(300,*)xAA        
            read(400,*)xAACa 
            x(i)      = xsol(i)   
            x(i+nr)   = psi(i)   
            x(i+2*nr) = xAA   
            x(i+3*nr) = xAACa
        enddo
            
        do i=1,4
            close(nfile(i))
        enddo

    endif
    !     .. end init from file 
  
    do i=1,neqint
        xguess(i)=x(i)
    enddo

end subroutine init_guess_pafiberborn


subroutine init_guess_pafiberbornscf(x, xguess)
      
    use parameters, only : xbulk
    use myutils, only :  newunit
    implicit none
  
    real(dp), intent(inout) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp), intent(out) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: i
    character(len=8) :: fname(5)
    integer :: ios,unfile(5)
    integer :: neq_bc 
    real(dp) :: xAA,xAACa,phi
  
    ! .. init guess all xbulk     

    do i=1,nr
        x(i)=xbulk%sol
        x(i+nr)=0.0_dp
        x(i+2*nr)=0.0_dp
        x(i+3*nr)=0.0_dp
        x(i+4*nr)=xpa(i)
    enddo

    do i=1,18
        x(i)=abs(1.0_dp-xpa(i)) 
    enddo
        
    neq_bc=0
    if(bcflag/="cc") then
        neq_bc=1 
        x(5*nr+neq_bc)=0.00_dp
    endif        

    if (infile.eq.1) then   ! infile is read in from file/stdio  
    
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'psi.in'
        write(fname(3),'(A6)')'xAA.in'
        write(fname(4),'(A8)')'xAACa.in' 
        write(fname(5),'(A6)')'phi.in' 
      
        do i=1,5 ! loop files
            open(unit=newunit(unfile(i)),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then    
                print*, 'file num ber =',unfile(i),' file name =',fname(i)
                print*, 'Error opening file : iostat =', ios
                stop
            endif
        enddo
     
        if(bcflag/="cc") read(200,*)psisurf     ! surface potential 
        do i=1,nr
            read(unfile(1),*)xsol(i)    ! solvent
            read(unfile(2),*)psi(i)     ! potential 
            read(unfile(3),*)xAA        
            read(unfile(4),*)xAACa 
            read(unfile(5),*)phi
            x(i)      = xsol(i)   
            x(i+nr)   = psi(i)   
            x(i+2*nr) = xAA   
            x(i+3*nr) = xAACa
            x(i+4*nr) = phi
        enddo
            
        do i=1,5
            close(unfile(i))
        enddo

    endif
    !     .. end init from file 
  
    do i=1,neqint
        xguess(i)=x(i)
    enddo

end subroutine init_guess_pafiberbornscf

!     purpose: initalize x and xguess
!     .. copy solution of previous solution ( distance ) to create new guess
!     .. data x=(pi,psi) and pi and psi order and split into  blocks
!     .. of size (nptso,nptsi,nptsb,nptss) =( outside, inside, boundary, on sphere )


subroutine make_guess_from_xstored(xguess,xstored)

    use globals
    use parameters
    use volume

    implicit none

    real(dp), intent(out) :: xguess(:)    ! guess volume fraction solvent and potentia 
    real(dp), intent(in) :: xstored(:)  ! assumed-shape array
 

    !   .. local variables
    integer :: i,neq_bc

    neq_bc=0    
    if(bcflag/="cc") neq_bc=neq_bc+1
    
    if (sysflag=="electnopoly".or.sysflag=="electligand") then 
        do i=1,nr/2
            xguess(i)=xstored(i)                    ! volume fraction solvent 
            xguess(i+nr)=xstored(i+nr+nrstep)       ! potential
        enddo
        do i=nr/2+1,nr  ! shift by nrstep
            xguess(i)=xstored(i+nrstep)    
            xguess(i+nr)=xstored(i+nr+2*nrstep)      
        enddo       

        do i=1,neq_bc
            xguess(2*nr+i)=xstored(2*(nr+nrstep)+i) 
        enddo   
    else
        print*,"Error : make_guess_from_xstored wrong sysflag"
        print*,"sysflag",sysflag
        stop
    endif    

end subroutine make_guess_from_xstored

subroutine make_guess(x, xguess,isfirstguess,flagstored,xstored)
  
    use globals
    use parameters
    use volume

    implicit none

    real(dp), intent(inout) :: x(neq)          ! iteration vector 
    real(dp), intent(out)   :: xguess(neq)     ! guess volume fraction solvent and potential 
    logical,  intent(in)    :: isfirstguess    ! first guess   
    logical,  optional, intent(in) :: flagstored
    real(dp), optional, intent(in) :: xstored(:) ! assumed-shaped array

    !     ..local variables 
    integer :: i,neq_bc

!    print*,"value isfirstguess=",isfirstguess    
  
    neq_bc=0    
    if(bcflag/="cc") neq_bc=neq_bc+1 


    if(present(flagstored)) then
        if(present(xstored)) then
            if(flagstored) then  
                ! print*,"flagstored==true"   
                ! copy solution xstored in xguess
                call make_guess_from_xstored(xguess,xstored)

            else if(isfirstguess) then       ! first guess

                if(sysflag=="electnopoly".or.sysflag=="electligand".or.sysflag=="pafiber") then 
                    call init_guess_electnopoly(x,xguess)
                else if(sysflag=="pafiberIm") then 
                    call init_guess_pafiberIm(x,xguess)
                else if(sysflag=="pafibervarelec") then 
                    call init_guess_electnopoly(x,xguess)
                else if(sysflag=="pafiberborn") then 
                    call init_guess_pafiberborn(x,xguess)           
                else if(sysflag=="pafiberbornscf") then 
                    call init_guess_pafiberbornscf(x,xguess)     
                else     
                    print*,"make_guess: wrong value sysflag : ", sysflag
                endif

            else  
                do i=1,neqint
                    xguess(i)=x(i)      ! volume fraction solvent 
                enddo
            endif
        else
            print*,"Error: argument xstored not present, while flagstored present"
            stop 
        endif 
    else if(isfirstguess) then       ! first guess

        if(sysflag=="electnopoly".or.sysflag=="electligand".or.sysflag=="pafiber") then 
            call init_guess_electnopoly(x,xguess)   
        else if(sysflag=="pafiberIm") then 
            call init_guess_pafiberIm(x,xguess)  
        else if(sysflag=="pafibervarelec") then 
            call init_guess_electnopoly(x,xguess) 
        else if(sysflag=="pafiberborn") then 
            call init_guess_pafiberborn(x,xguess)   
        else if(sysflag=="pafiberbornscf") then 
            call init_guess_pafiberbornscf(x,xguess)     
        else
            print*,"make_guess: wrong value sysflag : ", sysflag
        endif
    else      
        do i=1,neqint
            xguess(i)=x(i)      ! volume fraction solvent 
        enddo
    endif

end subroutine make_guess

end module initxvector
