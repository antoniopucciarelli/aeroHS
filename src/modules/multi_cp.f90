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
        use FOUL 

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
        character(len=1)                        :: answer

        ! setting properties 
        call setting_properties(P0,V,rho,alpha1,alpha2,start_angle,end_angle,dim,selection,selection_type) 
        
        ! setting first airfoil 
        filename  = 'GNUplot_coord_data1.dat'
        filename1 = 'GNUplot_mean_data1.dat'
        filename2 = 'GNUplot_tg_norm1.dat'
        call write_formatted('FIRST AIRFOIL','yellow')
        call make_airfoil(PANELsize1,MEANLINEarray1,PANEL_array1,airfoil1,alpha1,filename,filename1,filename2)
            
        ! saving 1st airfoil data
        call GNUplot_saving(PANEL_array1,MEANLINEarray1,airfoil1%get_npoints(),filename,filename1,filename2)

        ! setting second airfoil 
        filename  = 'GNUplot_coord_data2.dat'
        filename1 = 'GNUplot_mean_data2.dat'
        filename2 = 'GNUplot_tg_norm2.dat'
        call write_formatted('SECOND AIRFOIL','yellow')
        call make_airfoil(PANELsize2,MEANLINEarray2,PANEL_array2,airfoil2,alpha2,filename,filename1,filename2) 
        
        ! asking to invert airfoil --> useful for ground effect study
        print*, 'do you want to invert the 2nd airfoil coords? [Y\n]'
        read*, answer

        if(answer == 'Y' .or. answer == 'y')then 
            call invert(PANELsize2,PANEL_array2)
        end if 

        ! savin 2nd airfoil data
        call GNUplot_saving(PANEL_array2,MEANLINEarray2,airfoil2%get_npoints(),filename,filename1,filename2)

    end subroutine generate_airfoils 
    
    subroutine compute_multi_matrix(PANELsize1,PANELsize2,PANEL_array1,PANEL_array2,matrix) 
        use PANEL_object
        use cp 
        
        implicit none 
        
        integer(kind=4),intent(in)              :: PANELsize1 
        integer(kind=4),intent(in)              :: PANELsize2 
        type(panel),dimension(PANELsize1)       :: PANEL_array1 
        type(panel),dimension(PANELsize2)       :: PANEL_array2
        real(kind=8),dimension(:,:),allocatable :: matrix 
        real(kind=8),dimension(2)               :: normal_ith 
        real(kind=8),dimension(2)               :: tangent_first
        real(kind=8),dimension(2)               :: tangent_last
        real(kind=8)                            :: vortex_value
        integer(kind=4)                         :: i, j
        
        ! matrix allocation 
        allocate(matrix(PANELsize1+PANELsize2+2,PANELsize1+PANELsize2+2))
        
        ! computing 1st airfoil self induction 
        do i=1,PANELsize1
            
            normal_ith   = PANEL_array1(i)%normal 
            vortex_value = 0.0

            do j=1,PANELsize1  
                matrix(i,j)  = dot_product(normal_ith,integral(PANEL_array1(i),PANEL_array1(j),'source')) 
                vortex_value = vortex_value + dot_product(normal_ith,integral(PANEL_array1(i),PANEL_array1(j),'vortex'))
            end do 
                       
            matrix(i,PANELsize1+1) = vortex_value 

        end do 

        ! computing 1st airfoil induced by the 2nd airfoil 
        do i=1,PANELsize1 
            
            normal_ith   = PANEL_array1(i)%normal 
            vortex_value = 0.0
            
            do j=1,PANELsize2 
                matrix(i,PANELsize1+1+j) = dot_product(normal_ith,integral(PANEL_array1(i),PANEL_array2(j),'source'))
                vortex_value             = vortex_value + dot_product(normal_ith,integral(PANEL_array1(i),PANEL_array2(j),'vortex')) 
            end do 
            
            matrix(i,PANELsize1+PANELsize2+2) = vortex_value

        end do 

        ! KUTTA condition on 1st airfoil 
        ! tangent declaration for 1st airfoil
        tangent_first = PANEL_array1(1)%tangent
        tangent_last  = PANEL_array1(PANELsize1)%tangent
        
        ! 1st airfoil on itself
        vortex_value = 0.0

        do j=1,PANELsize1
            matrix(PANELsize1+1,j) = dot_product(tangent_first,integral(PANEL_array1(1),PANEL_array1(j),'source')) + & 
                                     dot_product(tangent_last,integral(PANEL_array1(PANELsize1),PANEL_array1(j),'source'))
            vortex_value = vortex_value + dot_product(tangent_first,integral(PANEL_array1(1),PANEL_array1(j),'vortex')) + &
                           dot_product(tangent_last,integral(PANEL_array1(PANELsize1),PANEL_array1(j),'vortex'))                
        end do

        matrix(PANELsize1+1,PANELsize1+1) = vortex_value

        ! 2nd airfoil on 1st airfoil
        vortex_value = 0.0

        do j=1,PANELsize2
            matrix(PANELsize1+1,PANELsize1+1+j) = & 
                                     dot_product(tangent_first,integral(PANEL_array1(1),PANEL_array2(j),'source')) + & 
                                     dot_product(tangent_last,integral(PANEL_array1(PANELsize1),PANEL_array2(j),'source'))
            vortex_value = vortex_value + dot_product(tangent_first,integral(PANEL_array1(1),PANEL_array2(j),'vortex')) + &
                           dot_product(tangent_last,integral(PANEL_array1(PANELsize1),PANEL_array2(j),'vortex'))                
        end do

        matrix(PANELsize1+1,PANELsize1+PANELsize2+2) = vortex_value

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

        ! KUTTA condition on 2nd airfoil 
        ! tangent declaration for 2nd airfoil
        tangent_first = PANEL_array2(1)%tangent
        tangent_last  = PANEL_array2(PANELsize2)%tangent
        
        ! 1st airfoil on 2nd airfoil 
        vortex_value = 0.0

        do j=1,PANELsize1
            matrix(PANELsize1+PANELsize2+2,j) = &  
                                       dot_product(tangent_first,integral(PANEL_array2(1),PANEL_array1(j),'source')) + & 
                                       dot_product(tangent_last,integral(PANEL_array2(PANELsize2),PANEL_array1(j),'source'))
            vortex_value = vortex_value + dot_product(tangent_first,integral(PANEL_array2(1),PANEL_array1(j),'vortex')) + &
                           dot_product(tangent_last,integral(PANEL_array2(PANELsize2),PANEL_array1(j),'vortex'))                
        end do

        matrix(PANELsize1+PANELsize2+2,PANELsize1+1) = vortex_value

        ! 2nd airfoil itself                
        vortex_value = 0.0

        do j=1,PANELsize2
            matrix(PANELsize1+PANELsize2+2,PANELsize1+1+j) = & 
                                     dot_product(tangent_first,integral(PANEL_array2(1),PANEL_array2(j),'source')) + & 
                                     dot_product(tangent_last,integral(PANEL_array2(PANELsize2),PANEL_array2(j),'source'))
            vortex_value = vortex_value + dot_product(tangent_first,integral(PANEL_array2(1),PANEL_array2(j),'vortex')) + &
                           dot_product(tangent_last,integral(PANEL_array2(PANELsize2),PANEL_array2(j),'vortex'))                
        end do

        matrix(PANELsize1+PANELsize2+2,PANELsize1+PANELsize2+2) = vortex_value
    
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
        vector(PANELsize1+PANELsize2+2) = dot_product(PANEL_array2(1)%tangent,vel) + & 
                                          dot_product(PANEL_array2(PANELsize2)%tangent,vel)
        
        vector = - vector

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
        integer(kind=4)                                            :: x1, x2
        
        ! declaring grid dimensions
        nrows = 150
        ncols = 350

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
                
                ! checking internal airfoil area --> 1st arifoil                 
                if(grid(i,j)%coords(1) >= PANEL_array1(PANELsize1/2)%coords1(1)-0.001 .and. & 
                   grid(i,j)%coords(1) <= PANEL_array1(1)%coords1(1)+0.001)then 
                    
                    x1 = 0
                    k = 1

                    do while(x1 == 0 .and. k<=PANELsize1/2)
                             
                        if(grid(i,j)%coords(2) >= PANEL_array1(PANELsize1/2-k+1)%midpoint(2)-0.001 .and. &
                           grid(i,j)%coords(2) <= PANEL_array1(PANELsize1/2+k)%midpoint(2)+0.001)then
                            
                           if(grid(i,j)%coords(1) <= PANEL_array1(PANELsize1/2+k)%coords2(1)+0.001 .and. & 
                              grid(i,j)%coords(1) >= PANEL_array1(PANELsize1/2+k-1)%coords1(1)-0.001)then
                                                                  
                                x1 = 1                             
                                                                  
                           end if                                 
                                                                  
                        end if
                        
                        k = k + 1

                    end do
                end if
                
                ! checking internal airfoil area --> 2nd airfoil
                if(grid(i,j)%coords(1) >= PANEL_array2(PANELsize2/2)%coords1(1) .and. & 
                   grid(i,j)%coords(1) <= PANEL_array2(1)%coords1(1))then 
                    
                    x2 = 0
                    k = 1

                    do while(x2 == 0 .and. k<=PANELsize2/2)
                             
                        if(grid(i,j)%coords(2) >= PANEL_array2(PANELsize2/2-k+1)%midpoint(2) .and. &
                           grid(i,j)%coords(2) <= PANEL_array2(PANELsize2/2+k)%midpoint(2))then
                            
                           if(grid(i,j)%coords(1) <= PANEL_array2(PANELsize2/2+k)%coords2(1) .and. & 
                              grid(i,j)%coords(1) >= PANEL_array2(PANELsize2/2+k-1)%coords1(1))then
                                                                  
                                x2 = 1                             
                                                                  
                           end if                                 
                                                                  
                        end if
                        
                        k = k + 1

                    end do
                end if
               
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
                
                if(grid(i,j)%coords(1) >= PANEL_array1(PANELsize1/2)%coords1(1)-0.001 .and. & 
                   grid(i,j)%coords(1) <= PANEL_array1(1)%coords1(1)+0.001)then 
                    if(x1 == 1)then 
                        velocity = (/0.0, 0.0/)
                        norm_vel = 0
                    end if
                end if

                if(grid(i,j)%coords(1) >= PANEL_array2(PANELsize2/2)%coords1(1) .and. & 
                   grid(i,j)%coords(1) <= PANEL_array2(1)%coords1(1))then 
                    if(x2 == 1)then 
                        velocity = (/0.0, 0.0/)
                        norm_vel = 0
                    end if
                end if

                ! computing the pressure with BERNOULLI's theorem 
                pressure = P0 + 0.5*rho*(V**2 - norm_vel**2)  
                
                ! writing data in file 
                write(1,*) dummy_panel%midpoint(1), dummy_panel%midpoint(2), &
                           velocity(1), velocity(2), norm_vel, & 
                           pressure 
            end do 
        end do
        
        call write_formatted('[','normal','OK','green','] -- velocity and pressure field computed                ','normal', &
                              '--> ','normal',filename,'normal') 
        
        close(1)

        ! plot data 
        call system('gnuplot -p PRESSUREfieldMULTI.plt')
        
        ! grid deallocation
        do i=1,nrows
            do j=1,ncols
                deallocate(grid(i,j)%coords)
            end do
        end do 
        deallocate(grid)
    end subroutine compute_multi_field

    subroutine compute_MULTIairfoilFIELD(solution,PANEL_array1,PANEL_array2,PANELsize1,PANELsize2,cp_vec1,cp_vec2,P0,alpha,V,rho)  
        use PANEL_object 
        use math_module
        use FOUL
        use cp 
        
        implicit none 
        
        integer(kind=4),intent(in)                                  :: PANELsize1
        integer(kind=4),intent(in)                                  :: PANELsize2 
        type(panel),dimension(PANELsize1),intent(in)                :: PANEL_array1
        type(panel),dimension(PANELsize2),intent(in)                :: PANEL_array2
        real(kind=8),dimension(PANELsize1+PANELsize2+2),intent(in)  :: solution
        real(kind=8),dimension(PANELsize1),intent(out)              :: cp_vec1
        real(kind=8),dimension(PANELsize2),intent(out)              :: cp_vec2       
        real(kind=8),intent(in)                                     :: V
        real(kind=8),intent(in)                                     :: alpha
        real(kind=8),intent(in)                                     :: P0
        real(kind=8),intent(in)                                     :: rho
        real(kind=8),dimension(2,2)                                 :: ROT  
        real(kind=8),dimension(2)                                   :: velocity
        real(kind=8),dimension(2)                                   :: vel
        real(kind=8)                                                :: pressure
        real(kind=8)                                                :: norm_vel
        real(kind=8),dimension(PANELsize1)                          :: cp_coeff1
        real(kind=8),dimension(PANELsize2)                          :: cp_coeff2
        real(kind=8)                                                :: normal_velocity
        real(kind=8)                                                :: circulation
        integer(kind=4)                                             :: j, k       

        ! initializing circulation
        circulation = 0.0

        ! open file to save data 
        open(unit=1, file='FLOWfield_airfoilMULTI1.dat', status='replace')
        open(unit=2, file='cp_dataMULTI1.dat', status='replace')
        
        ! exernal velocity 
        vel(1) = V*cos(alpha)
        vel(2) = V*sin(alpha) 
        
        ! rotation matrix --> alpha angle > 0 if it's counterclockwise
        ! -- the sine and cosine terms respect the direct rotation matrix description 
        ROT(1,1) =  cos(alpha)
        ROT(2,1) =  sin(alpha) 
        ROT(1,2) = -sin(alpha)
        ROT(2,2) =  cos(alpha)      

        do j=1,PANELsize1 

            ! updating circulation 
            circulation = circulation + solution(PANELsize1+1)*PANEL_array1(j)%get_length()
                
            velocity = (/0.0, 0.0/)
            
            ! computing airfoil's panels induction 
            do k=1,PANELsize1 
                velocity = velocity + integral(PANEL_array1(j),PANEL_array1(k),'source')*solution(k)
                velocity = velocity + integral(PANEL_array1(j),PANEL_array1(k),'vortex')*solution(PANELsize1+1)
            end do
            
            ! computing ground panels induction 
            do k=1,PANELsize2
                velocity = velocity + integral(PANEL_array1(j),PANEL_array2(k),'source')*solution(PANELsize1+1+k)
                velocity = velocity + integral(PANEL_array1(j),PANEL_array2(k),'vortex')*solution(PANELsize1+PANELsize2+2)
            end do

            velocity = velocity + vel
            
            norm_vel = norm(velocity)

            pressure = P0 + 0.5*rho*(V**2 - norm_vel**2)

            cp_coeff1(j) = 1 - (norm_vel/V)**2  
                
            normal_velocity = dot_product(PANEL_array1(j)%normal,velocity)

            write(1,*) PANEL_array1(j)%midpoint, matmul(ROT,PANEL_array1(j)%midpoint), &
                       velocity(1), velocity(2), norm_vel, & 
                       pressure   , cp_coeff1(j)   , normal_velocity       
            
        end do
        
        call write_formatted('[','normal','OK','green','] -- velocity and pressure field around airfoil computed ','normal', & 
                             '--> FLOWfield_airfoil.dat','normal') 
        
        print*, 'total circulation (GAMMA) 1st airfoil = ', circulation 
        
        ! writing Cp data in a separate file in order to separate the UPPER and LOWER surface in plotting procedure 
        do j=1,PANELsize1/2
            write(2,*) PANEL_array1(j)%midpoint(1),              - cp_coeff1(j), & 
                       PANEL_array1(j+PANELsize1/2)%midpoint(1), - cp_coeff1(PANELsize1/2+j)
        end do
                
        close(1)
        close(2)

        cp_vec1 = cp_coeff1
        
        ! initializing circulation
        circulation = 0.0

        ! open file to save data 
        open(unit=1, file='FLOWfield_airfoilMULTI2.dat', status='replace')
        open(unit=2, file='cp_dataMULTI2.dat', status='replace')
        
        ! exernal velocity 
        vel(1) = V*cos(alpha)
        vel(2) = V*sin(alpha) 
        
        ! rotation matrix --> alpha angle > 0 if it's counterclockwise
        ! -- the sine and cosine terms respect the direct rotation matrix description 
        ROT(1,1) =  cos(alpha)
        ROT(2,1) =  sin(alpha) 
        ROT(1,2) = -sin(alpha)
        ROT(2,2) =  cos(alpha)      

        do j=1,PANELsize2

            ! updating circulation 
            circulation = circulation + solution(PANELsize1+PANELsize2+2)*PANEL_array2(j)%get_length()
                
            velocity = (/0.0, 0.0/)
            
            ! computing airfoil's panels induction 
            do k=1,PANELsize2
                velocity = velocity + integral(PANEL_array2(j),PANEL_array2(k),'source')*solution(PANELsize1+1+k)
                velocity = velocity + integral(PANEL_array2(j),PANEL_array2(k),'vortex')*solution(PANELsize1+PANELsize2+2)
            end do
            
            ! computing ground panels induction 
            do k=1,PANELsize1
                velocity = velocity + integral(PANEL_array2(j),PANEL_array1(k),'source')*solution(k)
                velocity = velocity + integral(PANEL_array2(j),PANEL_array1(k),'vortex')*solution(PANELsize1+1)
            end do

            velocity = velocity + vel
            
            norm_vel = norm(velocity)

            pressure = P0 + 0.5*rho*(V**2 - norm_vel**2)

            cp_coeff2(j) = 1 - (norm_vel/V)**2  
                
            normal_velocity = dot_product(PANEL_array2(j)%normal,velocity)

            write(1,*) PANEL_array2(j)%midpoint, matmul(ROT,PANEL_array2(j)%midpoint), &
                       velocity(1), velocity(2), norm_vel, & 
                       pressure   , cp_coeff2(j)   , normal_velocity       
            
        end do
        
        call write_formatted('[','normal','OK','green','] -- velocity and pressure field around airfoil computed ','normal', & 
                             '--> FLOWfield_airfoil.dat','normal') 
        
        print*, 'total circulation (GAMMA) 2nd airfoil = ', circulation 
        
        ! writing Cp data in a separate file in order to separate the UPPER and LOWER surface in plotting procedure 
        do j=1,PANELsize2/2
            write(2,*) PANEL_array2(j)%midpoint(1),              - cp_coeff2(j), & 
                       PANEL_array2(j+PANELsize2/2)%midpoint(1), - cp_coeff2(PANELsize2/2+j)
        end do
                
        close(1)
        close(2)

        cp_vec2 = cp_coeff2

        ! plot data 
        call system('gnuplot -p CP_plot_MULTI.plt')

    end subroutine compute_MULTIairfoilFIELD
        
    subroutine invert(PANELsize,PANEL_array)
        use PANEL_object 
    
        implicit none 
        
        integer(kind=4),intent(in)                     :: PANELsize
        integer(kind=4)                                :: i 
        type(panel),dimension(PANELsize),intent(inout) :: PANEL_array
        type(panel),dimension(PANELsize)               :: PANEL_array_temp 
        real(kind=8),dimension(2)                      :: temp 

        do i=1,PANELsize
            
            PANEL_array(i)%coords1(2)  = - PANEL_array(i)%coords1(2)
            PANEL_array(i)%coords2(2)  = - PANEL_array(i)%coords2(2)
            PANEL_array(i)%midpoint(2) = - PANEL_array(i)%midpoint(2) 
            
            ! swapping panel ID
            PANEL_array(i)%id = PANELsize - (i-1) 

            ! swapping coords to match the requirements
            temp                   = PANEL_array(i)%coords1 
            PANEL_array(i)%coords1 = PANEL_array(i)%coords2
            PANEL_array(i)%coords2 = temp

            call PANEL_array(i)%set_angle() 

            if(PANEL_array(i)%POS == 'UP')then 
                PANEL_array(i)%POS = 'DW'
            else 
                PANEL_array(i)%POS = 'UP'
            end if
            
            call PANEL_array(i)%compute_tangent_and_normal() 

            call PANEL_array(i)%compute_ROT()

            PANEL_array_temp(PANELsize-(i-1)) = PANEL_array(i) 

        end do  

        PANEL_array = PANEL_array_temp

    end subroutine invert  

end module multi_cp
