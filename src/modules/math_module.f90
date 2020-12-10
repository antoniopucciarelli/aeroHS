module math_module
    implicit none

    real(kind=8), parameter :: pi = 4 * atan(1.0)  

    contains 

    function cross(x, y)
    ! this function compute the cross product between 2 vectors
        implicit none

        real(kind=8)              :: cross
        real(kind=8),dimension(2) :: x
        real(kind=8),dimension(2) :: y

        cross = x(1)*y(2) - x(2)*y(1)

    end function cross

    function norm(x)
    ! this function computes the norm of a vector

        implicit none 

        real(kind=8)                         :: norm
        real(kind=8),dimension(2),intent(in) :: x

        norm = sqrt(x(1)**2 + x(2)**2)

    end function norm

end module math_module