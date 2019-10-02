module molecules

    use precision_definition
    implicit none

    type moleclist
        real(dp) :: sol
        real(dp) :: Na
        real(dp) :: Cl
        real(dp) :: NO3 
        real(dp) :: K
        real(dp) :: Ca
        real(dp) :: NaCl
        real(dp) :: KCl
        real(dp) :: Hplus
        real(dp) :: OHmin
        real(dp) :: TB
        real(dp) :: TM
        real(dp) :: Rb 
        real(dp) :: Im  ! Imidazolium
        real(dp), dimension(5) :: pp
    end type moleclist
  

    type bornmoleclist
        real(dp) :: AA
        real(dp) :: AACa
        real(dp) :: Na
        real(dp) :: Cl
        real(dp) :: K
        real(dp) :: Ca
        real(dp) :: Hplus
        real(dp) :: OHmin
        real(dp) :: Rb 
        real(dp) :: Im  
    end type bornmoleclist 

end module molecules
