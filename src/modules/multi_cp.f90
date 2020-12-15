module multi_cp 

    contains 

    subroutine generate_airfoils(selection,selection_type,alpha1,alpha2,PANELsize1,PANELsize2,PANEL_array1, & 
                                 PANEL_array2,MEANLINEarray1,MEANLINEarray2)
        use PANEL_object 
        use MEANline_object    
        use AIRFOIL_object
        use airfoilgenerator
             
        implicit none 

        type(panel),dimension(:),allocatable    :: PANEL_array1 
        type(panel),dimension(:),allocatable    :: PANEL_array2 
        type(MEANline),dimension(:),allocatable :: MEANLINEarray1 
        type(MEANline),dimension(:),allocatable :: MEANLINEarray2    
        integer(kind=4)                         :: PANELsize1
        integer(kind=4)                         :: PANELsize2 
        type(NACA_airfoil)                      :: airfoil1
        type(NACA_airfoil)                      :: airfoil2
        real(kind=8)                            :: alpha1
        real(kind=8)                            :: alpha2 
        real(kind=8)                            :: V              = 0.0
        real(kind=8)                            :: P0             = 0.0
        real(kind=8)                            :: rho            = 0.0
        integer(kind=4)                         :: start_angle    = 0
        integer(kind=4)                         :: end_angle      = 0
        integer(kind=4)                         :: dim            = 0
        integer(kind=4)                         :: selection      = 0
        integer(kind=4)                         :: selection_type = 0

        ! setting properties 
        call setting_properties(P0,V,rho,alpha1,alpha2,start_angle,end_angle,dim,selection,selection_type) 

        ! setting first airfoil 
        call make_airfoil(PANELsize1,MEANLINEarray1,PANEL_array1,airfoil1,alpha1)
        
        ! setting second airfoil 
        call make_airfoil(PANELsize2,MEANLINEarray2,PANEL_array2,airfoil2,alpha2) 

    end subroutine generate_airfoils 
    
    subroutine compute_multi_matrix() 
        use PANEL_object
        use cp 
        
        implicit none 
        
        integer(kind=4),intent(in)                   :: PANELsize1 
        integer(kind=4),intent(in)                   :: PANELsize2 
        type(panel),dimension(PANELsize1),intent(in) :: PANEL_array1 
        type(panel),dimension(PANELsize2),intent(in) :: PANEL_array2
        real(kind=8),dimension(:,:),allocatable      :: matrix 
        real(kind=8),dimension(2)                    :: normal_ith 
        real(kind=8)                                 :: vortex_value
        integer(kind=4)                              :: i, j
        
        ! matrix allocation 
        allocate(matrix(PANELsize1+PANELsize2+2,PANELsize1+PANELsize2+2))
        
        ! computing 1st airfoil self induction 
        do i=1,PANELsize1
            
            normal_ith   = PANEL_array1(i)%normal 
            vortex_value = 0.0

            do j=1,PANELsize1  
                matrix(i,j) = dot_product(normal_ith,integral(PANEL_array1(i),PANEL_array1(j),'source')) 
                vortex_value = vortex_value + dot_product(normal_ith,interal(PANEL_array1(i),PANEL_array1,'vortex'))
            end do 
            
            matrix(i,PANELsize1+1) = vortex_value 

        end do 

        ! computing 1st airfoil induced by the 2nd airfoil 
        do i=1,PANELsize1 
            
            normal_ith   = PANEL_array1(i)%normal 
            vortex_value = 0.0
            
            do j=1,PANELsize2 
                matrix(i,PANELsize1+1) = dot_product(normal_ith,integral(PANEL_array1(i),PANEL_array2(j),'source'))
                vortex_value           = vortex_value + dot_product(normal_ith,integral(PANEL_array1(i),PANEL_array2(j),'vortex')) 
            end do 
            
            matrix(i,PANELsize1+PANELsize+2) = vortex_value

        end do 
        
        ! computing 2nd airfoil induced by the 1st airfoil
        do i=1,PANELsize2 
            
            normal_ith   = PANEL_array2(i)%normal 
            vortex_value = 0.0 
            
            do j=1,PANELsize1
                matrix(PANELsize1+1+i,j) = dot_product(normal_ith,integral(PANEL_array2(i),PANEL_array1(j),'source'))
                vortex_value             = vortex_value + dot_product(normal_ith,integral(PANEL_array2(i),PANEL_array1(j),'vortex'))
            end do 
            
            matrix(PANELsize1+1+j,PANELsize1+1)

        end do 

        ! computing 2nd airfoil self induction 
        do i=1,PANELsize2 
            
            normal_ith   = PANEL_array2(i)%normal 
            vortex_value = 0.0 
            
            do j=1,PANELsize2
                matrix(PANELsize1+1+i,PANELsize1+1+j) = dot_product(normal_ith,integral(PANEL_array2(i),PANEL_array2(j),'source'))
                vortex_value             = vortex_value + dot_product(normal_ith,integral(PANEL_array2(i),PANEL_array2(j),'vortex'))
            end do 
            
            matrix(PANELsize1+1+j,PANELsize1+1+j)

        end do 
    
    end subroutine compute_multi_matrix

    subroutine compute_multi_vector(vector,PANELsize1,PANELsize2,PANEL_array1,PANEL_array2,alpha,V) 
        use PANEL_object 
        implicit none 

        integer(kind=4),intent(in)                   :: PANELsize1
        integer(kind=4),intent(in)                   :: PANELsize2 
        type(panel),dimension(PANELsize1),intent(in) :: PANEL_array1 
        type(panel),dimension(PANELsize2),intent(in) :: PANEL_array2  
        integer(kind=4)                              :: i
        real(kind=8),intent(in)                      :: V
        real(kind=8),intent(in)                      :: alpha 
        real(kind=8),dimension(2)                    :: vel 
        real(kind=8),dimension(:),allocatable        :: vector 

        ! determining vec --> more general description 
        vel(1) = V*cos(alpha)
        vel(2) = V*sin(alpha)
            
        ! vector variable allocation 
        allocate(vector(PANELsize1+PANELsize2+2))
        
        ! velocity on 1st airfoil 
        do i=1,PANELsize  
            vector(i) = dot_product(PANEL_array1(i)%normal,vel)
        end do 
        
        ! KUTTA condition on 1st airfoil 
        vector(PANELsize+1) = dot_product(PANEL_array1(1)%tangent,vel) + dot_product(PANEL_array1(PANELsize1)%tangent,vel)

        ! velocity on 2nd airfoil
        do i=1,PANELsize2
            vector(PANELsize1+1+j) = dot_product(PANEL_array2(i)%normal,vel) 
        end do 

        ! KUTTA condition on 2nd airfoil 
        vector(PANELsize1+PANELsize2+2) =dot_product(PANEL_array2(1)%tangent,vel) + & 
                                          dot_product(PANEL_array2(PANELsize2)%tangent,vel)

    end subroutine compute_multi_vector 

    subroutine solve_multi() 

    end subroutine solve_multi 

    subroutine compute_multi_field()

    end subroutine compute_multi_field

end module multi_cp
