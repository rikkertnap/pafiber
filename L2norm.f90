
module vectornorm

! computes L2norm
    use precision_definition
    implicit none

contains

    function l2norm (f,n)result(norm)
      
        implicit none

        integer, intent(in)  :: n 
        real(dp), intent(in) :: f(n)
        real(dp)             :: norm ! output

        integer :: i ! dummy index

        norm=0.0_dp
        do i=1,n
            norm = norm + f(i)*f(i)
        enddo
        norm=sqrt(norm)
        
    end function l2norm

    function l2norm_part(f,n,m,k)result(norm)
      
        implicit none

        integer, intent(in)  :: n,m,k 
        real(dp), intent(in) :: f(n)
        real(dp)             :: norm ! output

        integer :: i ! dummy index

        norm=0.0_dp
        do i=m,k
            norm = norm + f(i)*f(i)
        enddo
        norm=sqrt(norm)
        
    end function l2norm_part

end module vectornorm

