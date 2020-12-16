module multi_cp 

    contains 

    subroutine generate_airfoils(selection,selection_type,alpha1,alpha2,PANELsize1,PANELsize2,PANEL_array1, & 
                                 PANEL_array2,MEANLINEarray1,MEANLINEarray2)
        use PANEL_object 
        use MEANline_object    
        use AIRFOIL_object
        use airfoilgenerator
        use ask_module 
        use discretization_module 

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
        integer(kind=4)                         :: selection     
        integer(kind=4)                         :: selection_type
        character(len=30)                       :: filename 
        character(len=30)                       :: filename1
        character(len=30)                       :: filename2

        ! setting properties 
        call setting_properties(P0,V,rho,alpha1,alpha2,start_angle,end_angle,dim,selection,selection_type) 

        ! setting first airfoil 
        filename  = 'GNUplot_coord_data1.dat'
        filename1 = 'GNUplot_mean_data1.dat'
        filename2 = 'GNUplot_tg_norm1.dat'
        call make_airfoil(PANELsize1,MEANLINEarray1,PANEL_array1,airfoil1,alpha1,filename,filename1,filename2)
            
        ! saving 1st airfoil data
        call GNUplot_saving(PANEL_array1,MEANLINEarray1,airfoil1%get_npoints(),filename,filename1,filename2)

        ! setting second airfoil 
        filename  = 'GNUplot_coord_data2.dat'
        filename1 = 'GNUplot_mean_data2.dat'
        filename2 = 'GNUplot_tg_norm2.dat'
        call make_airfoil(PANELsize2,MEANLINEarray2,PANEL_array2,airfoil2,alpha2,filename,filename1,filename2) 
            
        ! savin 2nd airfoil data
        call GNUplot_saving(PANEL_array2,MEANLINEarray2,airfoil2%get_npoints(),filename,filename1,filename2)

    end subroutine generate_airfoils 
    
    subroutine compute_multi_matrix(PANELsize1,PANELsize2,PANEL_array1,PANEL_array2,matrix) 
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
                vortex_value = vortex_value + dot_product(normal_ith,integral(PANEL_array1(i),PANEL_array1(j),'vortex'))
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
            
            matrix(i,PANELsize1+PANELsize2+2) = vortex_value

        end do 
        
        ! computing 2nd airfoil induced by the 1st airfoil
        do i=1,PANELsize2 
            
            normal_ith   = PANEL_array2(i)%normal 
            vortex_value = 0.0 
            
            do j=1,PANELsize1
                matrix(PANELsize1+1+i,j) = dot_product(normal_ith,integral(PANEL_array2(i),PANEL_array1(j),'source'))
                vortex_value             = vortex_value + dot_product(normal_ith,integral(PANEL_array2(i),PANEL_array1(j),'vortex'))
            end do 
            
            matrix(PANELsize1+1+i,PANELsize1+1) = vortex_value

        end do 

        ! computing 2nd airfoil self induction 
        do i=1,PANELsize2 
            
            normal_ith   = PANEL_array2(i)%normal 
            vortex_value = 0.0 
            
            do j=1,PANELsize2
                matrix(PANELsize1+1+i,PANELsize1+1+j) = dot_product(normal_ith,integral(PANEL_array2(i),PANEL_array2(j),'source'))
                vortex_value             = vortex_value + dot_product(normal_ith,integral(PANEL_array2(i),PANEL_array2(j),'vortex'))
            end do 
            
            matrix(PANELsize1+1+i,PANELsize1+1+j) = vortex_value 

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
        
        ! allocating vector variable 
        allocate(vector(PANELsize1+PANELsize2+2))

        ! determining vec --> more general description 
        vel(1) = V*cos(alpha)
        vel(2) = V*sin(alpha)
            
        ! velocity on 1st airfoil --> no penetration conditions 
        do i=1,PANELsize1 
            vector(i) = dot_product(PANEL_array1(i)%normal,vel)
        end do 
        
        ! KUTTA condition on 1st airfoil 
        vector(PANELsize1+1) = dot_product(PANEL_array1(1)%tangent,vel) + dot_product(PANEL_array1(PANELsize1)%tangent,vel)

        ! velocity on 2nd airfoil --> no penetration conditions
        do i=1,PANELsize2
            vector(PANELsize1+1+i) = dot_product(PANEL_array2(i)%normal,vel) 
        end do 

        ! KUTTA condition on 2nd airfoil 
        vector(PANELsize1+PANELsize2+2) =dot_product(PANEL_array2(1)%tangent,vel) + & 
                                          dot_product(PANEL_array2(PANELsize2)%tangent,vel)

    end subroutine compute_multi_vector 

    subroutine solvemulti(PANELsize1,PANELsize2,matrix,vector,solution) 
        use cp
        implicit none 
        
        integer(kind=4),intent(in)                                                         :: PANELsize1
        integer(kind=4),intent(in)                                                         :: PANELsize2
        real(kind=8),dimension(PANELsize1+PANELsize2+2,PANELsize1+PANELsize2+2),intent(in) :: matrix 
        real(kind=8),dimension(PANELsize1+PANELsize2+2,PANELsize1+PANELsize2+2)            :: M 
        real(kind=8),dimension(PANELsize1+PANELsize2+2),intent(in)                         :: vector 
        real(kind=8),dimension(:),allocatable                                              :: solution 
        integer(kind=4),dimension(PANELsize1+PANELsize2+2)                                 :: pivoting_vec 
        integer(kind=4),dimension(PANELsize1+PANELsize1+2)                                 :: info
        
        ! allocating solution 
        allocate(solution(PANELsize1+PANELsize2+2))

        ! allocating matrix and vector --> during DGESV they are changed 
        M        = matrix 
        solution = vector 
        
        ! compute solution 
        call DGESV(PANELsize1+PANELsize2+2, 1, M, PANELsize1+PANELsize2+2, pivoting_vec, solution, PANELsize1+PANELsize2+2, info)
        
        ! test solution 
        call test_matrix(matrix,solution,vector,PANELsize1+PANELsize2+1)

    end subroutine solvemulti 

    subroutine compute_multi_field(solution,PANEL_array1,PANEL_array2,PANELsize1,PANELsize2,P0,alpha,V,rho,filename)
        use PANEL_object 
        use math_module
        use FOUL
        use cp 
        
        implicit none 

        ! generating a grid of measure points 
        ! this grid can be expressed in a monodimensional array that contains 2 dimensional vectors 
        !   the procedure can be possible if it's used a type variable that describes a vector
        type :: array_type  
            real(kind=8),pointer :: coords(:)
        end type array_type 
        
        integer(kind=4),intent(in)                                 :: PANELsize1
        integer(kind=4),intent(in)                                 :: PANELsize2 
        type(panel),dimension(PANELsize1),intent(in)               :: PANEL_array1
        type(panel),dimension(PANELsize2),intent(in)               :: PANEL_array2
        type(panel)                                                :: dummy_panel
        real(kind=8),dimension(PANELsize1+PANELsize2+2),intent(in) :: solution 
        real(kind=8),dimension(2)                                  :: velocity
        real(kind=8),dimension(2)                                  :: vel
        integer(kind=4)                                            :: i, j, k       
        real(kind=8)                                               :: deltax, deltay
        real(kind=8)                                               :: x_start, x_end
        real(kind=8)                                               :: y_start, y_end
        real(kind=8),intent(in)                                    :: V
        real(kind=8),intent(in)                                    :: alpha
        real(kind=8)                                               :: pressure
        real(kind=8)                                               :: norm_vel
        real(kind=8),intent(in)                                    :: P0
        real(kind=8),intent(in)                                    :: rho
        type(array_type),dimension(:,:),allocatable                :: grid
        integer(kind=4)                                            :: nrows, ncols
        character(len=30),intent(in)                               :: filename
        integer(kind=4)                                            :: x
        
        ! declaring grid dimensions
        nrows = 100
        ncols = 300

        ! grid allocation process in memory
        allocate(grid(nrows,ncols))
        do i=1,nrows
            do j=1,ncols
                allocate(grid(i,j)%coords(2))
            end do
        end do
           
        ! allocating grid dimensions
        call write_formatted('SETTING GRID','yellow')
        print*, 'type x coords:'
        print*, '   start point'
        read*,  x_start
        print*, '   end point'
        read*,  x_end
        print*, 'type y coords:'
        print*, '   start point'
        read*,  y_start 
        print*, '   end point'
        read*,  y_end

        ! grid dimension allocation process
        deltax = (x_end - x_start)/ncols
        deltay = (y_end - y_start)/nrows        
        
        ! grid data allocation process
        ! there may be some errors with the computation of the velocity values in the system 
        ! this depends on the position of the grid points and it can generate some errors in the field
        ! -- this depends on the position of the study points in the grid 
        ! -- a wrong placed study point can generate some troubles in the field study  
        do i=1,nrows 
            do j=1,ncols
            ! adding deltax/2 and deltay/2 to avoid computation errors such as NaN
                grid(i,j)%coords(1) = x_start + (j-1)*deltax + deltax/2
                grid(i,j)%coords(2) = y_start + (i-1)*deltay + deltay/2
            end do
        end do

        call write_formatted('COMPUTING FIELD','yellow')

        ! open file to save data 
        open(unit=1, file=filename, status='replace')
        
        ! allocating velocity vector --> in this case the alpha angle should be 0
        ! --> description of the vel variable is made to render it more flexible for future changes 
        vel(1) = V*cos(alpha)
        vel(2) = V*sin(alpha) 

        do i=1,nrows 
                
             do j=1,ncols   
                
                ! velocity vector starts from null vector 
                velocity = (/0.0, 0.0/)
                ! this panel object allows to compute the velocity at each point of the field 
                ! described by the dummy_panel midpoint 
                dummy_panel%midpoint = grid(i,j)%coords
                
!                if(grid(i,j)%coords(1) >= PANEL_array(PANELsize/2)%coords1(1) .and. & 
!                   grid(i,j)%coords(1) <= PANEL_array(1)%coords1(1))then 
!                    
!                    x = 0
!                    k = 1
!
!                    do while(x == 0 .and. k<=PANELsize/2)
!                             
!                        if(grid(i,j)%coords(2) >= PANEL_array(PANELsize/2-k+1)%midpoint(2) .and. &
!                           grid(i,j)%coords(2) <= PANEL_array(PANELsize/2+k)%midpoint(2))then
!                            
!                           if(grid(i,j)%coords(1) <= PANEL_array(PANELsize/2+k)%coords2(1) .and. & 
!                              grid(i,j)%coords(1) >= PANEL_array(PANELsize/2+k-1)%coords1(1))then
!                                                                  
!                                x = 1                             
!                                                                  
!                           end if                                 
!                                                                  
!                        end if
!                        
!                        k = k + 1
!
!                    end do
!                end if
                
                ! computing induction by the 1st airfoil panels 
                ! source distribution induction 
                ! vortex distribution induction 
                do k=1,PANELsize1 
                    velocity = velocity + integral(dummy_panel,PANEL_array1(k),'source')*solution(k)
                    velocity = velocity + integral(dummy_panel,PANEL_array1(k),'vortex')*solution(PANELsize1+1)
                end do

                ! computing induction by the 2nd airfoil panels
                ! source distribution induction
                ! vortex distribution induction 
                do k=1,PANELsize2
                    velocity = velocity + integral(dummy_panel,PANEL_array2(k),'source')*solution(PANELsize1+1+k)
                    velocity = velocity + integral(dummy_panel,PANEL_array2(k),'vortex')*solution(PANELsize1+PANELsize2+2)
                end do
                
                ! adding the flow stream velocity 
                velocity = velocity + vel
                
                ! computing the norm of velocity 
                norm_vel = norm(velocity)
                
                ! computing the pressure with BERNOULLI's theorem 
                pressure = P0 + 0.5*rho*(V**2 - norm_vel**2)  
                
!                if(x == 1)then 
!                    velocity = (/0.0, 0.0/)
!                    norm_vel = 0
!                end if

                ! writing data in file 
                write(1,*) dummy_panel%midpoint(1), dummy_panel%midpoint(2), &
                           velocity(1), velocity(2), norm_vel, & 
                           pressure 
            end do 
        end do
        
        call write_formatted('[','normal','OK','green','] -- velocity and pressure field computed                ','normal', &
                              '--> ','normal',filename,'normal') 
        
        ! plot data 
        call system('gnuplot -p PRESSUREfieldMULTI.plt')

        close(1)
        deallocate(grid)
    end subroutine compute_multi_field

end module multi_cp
