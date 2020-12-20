program beta_angle

    use math_module

    implicit none 

    ! initializing variables
    real(kind=8)                :: beta, acosbeta
    real(kind=8)                :: theta
    real(kind=8)                :: r1mod, r2mod
    real(kind=8),dimension(2)   :: N1, N2, T, r1, r2, midpoint
    real(kind=8),dimension(2,2) :: ROT

    ! getting points coords
    read(*,*) N1(1), N1(2)
    read(*,*) N2(1), N2(2)
    read(*,*) T(1), T(2)

    midpoint = (N1 + N2)/2.0

    ! compute panel inclination angle 
    theta = atan2((N2(2)-N1(2)),(N2(1)-N1(1)))

    ! declaring rotation matrix
    ROT(1,1) =  cos(theta)
    ROT(1,2) =  sin(theta)
    ROT(2,1) = -sin(theta)
    ROT(2,2) =  cos(theta)

    ! compute distance vectors between target point and N1, N2
    r1 = T - N1
    r2 = T - N2
    ! compute r1, r2 module
    r1mod = norm(r1)
    r2mod = norm(r2)
    ! compute beta angle
    acosbeta = dot_product(r1,r2)/(r1mod * r2mod)
    beta = acos(acosbeta)

    ! determine beta angle sign
    ! 1 step - rotation of points by the panel inclination angle
    T  = matmul(ROT,  T)
    midpoint = matmul(ROT, midpoint)

    ! 2 step - points translation in order to have the midpoint of the panel at point (0, 0)
    T  = T  - midpoint 

    ! 3 step - determine if the angle is > 0 < by the position of T with respect to the panel
    if(T(2)>0)then
        beta = abs(beta)
    else
        beta = -abs(beta)
    end if

    print*, 'theta = ', theta
    print*, 'beta  = ', beta/pi*180

end program beta_angle
