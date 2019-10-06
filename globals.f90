                                                                  
    !     .. module file of global variables 

module globals

    use precision_definition
    use mathconst
  
    implicit none
  
    !     .. variables

    integer  :: nsize         ! size lattice, numer of layers
    integer(8)  :: neq        ! number of non-linear equations
    integer(8)  :: neqmax     ! maximum number of non-linear equations
    integer  :: neqint        ! number of non-linear equations, integer instead of integer(8)

    character(len=15) :: sysflag   ! sysflag selects fcn    
    character(len=15) :: runflag   ! runflag
    character(len=2)  :: bcflag    ! bcflag selects bc surface 

    integer, parameter :: AH2BH = 5
    integer, parameter :: AHBH  = 1
    integer, parameter :: AHB   = 2
    integer, parameter :: ABH   = 3
    integer, parameter :: AB    = 4
    integer, parameter :: SuOH  = 6
    integer, parameter :: SuCl  = 7  
    integer, parameter :: Su    = 8
    integer, parameter :: SuNO3 = 9
    


end module globals

