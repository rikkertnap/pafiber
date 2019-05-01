subroutine solver(x, xguess, accuracy, residual, issolution)
    use globals
    use parameters
    use listfcn
    use fcnpointer

    implicit none
  
    !     .. arguments
    real(dp) :: x(neq)
    real(dp) :: xguess(neq)  
    real(dp) :: accuracy
    real(dp) :: residual
    logical  :: issolution
    
    call set_size_neq  
    call set_fcn
    
    if(method.eq."kinsol") then
     
        call kinsol_gmres_solver(x, xguess, error, fnorm, issolution)
 
    else  
        print*,"Solver method incorrect"
        stop
    endif
  
end subroutine solver
