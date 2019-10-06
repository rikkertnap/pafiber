! module dielectfcn.f90:                                                 | 
! computes the effectice dielectic function epsfcn and          |
! the derivate of the dielectric function  with respect to      |
! the polymer volume fraction                                   |
! both are in units of the dielectric constant of ionic solvent |
! pre : phi = 1-xsol -xpos -xneg xHplus- xOHmin                 |
! post: espfcn and Despfcn                                      |

module dielectric_const

    use precision_definition  
    implicit none
    
    private                    ! default all routines in this module private 
    public  ::  dielectfcn,born 

contains

subroutine dielectfcn(phi,epsfcn,Depsfcn,dielectP,dielectW,n) 
    
    integer, intent(in) :: n
    real(dp), intent(in) ::  dielectP, dielectW 
    real(dp), intent(inout) :: epsfcn(:),Depsfcn(:)
    real(dp), intent(in)  :: phi(:)

    call dielectfcnAV(phi,epsfcn,Depsfcn,dielectP,dielectW,n) 
 !   call dielectfcnConst(epsfcn,Depsfcn,n) 

end subroutine

subroutine dielectfcnConst(epsfcn,Depsfcn,n) 

    integer, intent(in)  :: n 
    real(dp), intent(inout) :: epsfcn(:),Depsfcn(:)
 
    integer :: i  
    do i=1,n  
        epsfcn(i)=  1.0d0 ! dielectric function
        Depsfcn(i)= 0.0d0 ! derivative dielectric function    
    enddo
                                
end subroutine

! volume fraction weighted average of dielectric constant

subroutine dielectfcnAV(phi,epsfcn,Depsfcn,dielectP,dielectW ,n) 

    integer, intent(in)  :: n
    real(dp), intent(in) ::  dielectP, dielectW 
    real(dp), intent(inout) :: epsfcn(:),Depsfcn(:)
    real(dp), intent(in) :: phi(:)
    

    integer :: i
    real(dp)  :: ratioeps

    ratioeps = dielectP/dielectW 
   
    do i=1,n  
        epsfcn(i)= 1.0_dp-phi(i) + ratioeps * phi(i) ! dieletric function
        Depsfcn(i)= -1.0_dp+ratioeps ! derivative dieletric function    
    enddo
                                
end subroutine

! Maxwell-Garnett mixing rule for dielectic constant
! valid only for low phi<10-5 : 
subroutine dielectfcnMG(phi,epsfcn,Depsfcn,dielectP, dielectW, n)  
    use mathconst

    integer, intent(in) :: n
    real(dp), intent(in) ::  dielectP, dielectW 
    real(dp) , intent(inout):: epsfcn(:),Depsfcn(:)
    real(dp), intent(in) :: phi(:)

    !     .. local variables

    real(dp) :: Kalpha, Kmossotti, gamma
    real(dp) :: epspol,epssol
    integer ::i
  
    epssol =dielectW          ! permittivty ionic solution 
    epspol =dielectP          ! permittivty polymer  
   
    Kmossotti = (epssol-epspol)/(2.0_dp*epssol+epspol)  ! Mossotti factor
    gamma=1.5_dp         ! ratio between electric radius and matter radius
    Kalpha= Kmossotti * 4.0_dp* pi /( 3.0_dp* gamma**3)
      
    do i=1,n  
        epsfcn(i)= 1.0_dp-3.0_dp*Kalpha*phi(i)/(1.0_dp+Kalpha*phi(i)) ! dieletric function
        Depsfcn(i)= -(3.0_dp*Kalpha/((1.0_dp+Kalpha*phi(i))**2))   ! derivative dieletric function
    enddo
                                
end subroutine

                                                
! computes the born energy for a given  Bjerrum lenght and size 
! U_B = (beta z^2e^2)/(8 pi eps eps0 a )= lB z /(2 *a)                       


function born(lB, radius, z) result(born_energy)

    real(dp), intent(in) :: lB, radius
    integer, intent(in) :: z
    real(dp) :: born_energy

    born_energy=(lB*z*z/(2.0_dp*radius))

end function
      


end module

