subroutine solver(x, xguess, accuracy, residual, issolution)
    
    use precision_definition
    use globals, only : neq    
    use parameters, only : method, set_size_neq
    use listfcn, only : set_fcn
  
  
    !     .. arguments
    real(dp) :: x(neq)
    real(dp) :: xguess(neq)  
    real(dp) :: accuracy
    real(dp) :: residual
    logical  :: issolution
    
    call set_size_neq  
    call set_fcn
    
    if(method.eq."kinsol") then
     
        call kinsol_gmres_solver(x, xguess, accuracy, residual, issolution)
 
    else  
        print*,"Solver method incorrect"
        stop
    endif
  
end subroutine solver
