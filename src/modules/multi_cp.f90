module ground_cp

    contains 
    
    subroutine generate_ground(GROUNDpanel,GROUNDsize) 
        use PANEL_object
        use FOUL 

        implicit none 
        
        integer(kind=4),intent(inout)                      :: GROUNDsize
        type(panel),dimension(:),allocatable,intent(inout) :: GROUNDpanel
        real(kind=8)                                       :: GROUNDstart
        real(kind=8)                                       :: GROUNDend 
        integer(kind=4)                                    :: j
        real(kind=8)                                       :: theta 
        real(kind=8)                                       :: deltax    
        real(kind=8)                                       :: Y
        
        theta = 0.0

        ! asking # of panels to discretize the ground
        print*, 'type the number of panels for discretize the ground effet'
        read*, GROUNDsize
        
        ! allocate GROUNpanel
        allocate(GROUNDpanel(GROUNDsize))
        
        ! GROUNpanel extremes
        GROUNDend   = -0.2
        GROUNDstart =  1.5    

        ! panel length
        deltax = abs(GROUNDstart - GROUNDend)/GROUNDsize
        
        ! panel Y position 
        print*, 'type ground position --> it should be a negative number (depending on the posiotion of the airfoil)'
        read*, Y
        

        do j=1,GROUNDsize
            call GROUNDpanel(j)%set_coords(GROUNDstart + (j-1)*deltax, Y, GROUNDstart + j*deltax, Y)
            GROUNDpanel(j)%angle   = theta 
            GROUNDpanel(j)%normal  = (/0.0, 1.0/)
            GROUNDpanel(j)%tangent = (/1.0, 0.0/)
            GROUNDpanel(j)%length  = deltax
            call GROUNDpanel(j)%compute_ROT()
            call GROUNDpanel(j)%compute_midpoint('noprt')
            call GROUNDpanel(j)%set_id(j)
         end do  
        
         ! saving data in file
         open(unit=1, file='GROUNDdata.dat', status='replace')
            
         write(1,*) 'GROUND panel objects'

         do j=1,GROUNDsize 
            call GROUNDpanel(j)%saving(1)            
         end do 

         call write_formatted('[','normal','OK','green','] -- ground panels generated','normal') 
        
         close(1)
    end subroutine generate_ground 
        
    subroutine computeGROUNDmatrix(matrix,PANEL_array,GROUNDpanel,PANELsize,GROUNDsize)
        use PANEL_object
        use cp 
        use FOUL
        
        implicit none 
        
        integer(kind=4),intent(in)                   :: PANELsize
        integer(kind=4),intent(in)                   :: GROUNDsize
        integer(kind=4)                              :: i, j
        real(kind=8)                                 :: vortex_value 
        real(kind=8),dimension(2)                    :: normal_ith
        real(kind=8),dimension(2)                    :: tangent_first, tangent_last          
        real(kind=8),dimension(:,:),allocatable      :: matrix
        type(panel),dimension(GROUNDsize),intent(in) :: GROUNDpanel
        type(panel),dimension(PANELsize),intent(in)  :: PANEL_array
        
        ! allocation matrix 
        allocate(matrix(PANELsize+GROUNDsize+1,PANELsize+GROUNDsize+1))

        ! declaring tangent vectors for the KUTTA-JOUKOWSKY condition
        tangent_first = PANEL_array(1)%tangent 
        tangent_last  = PANEL_array(PANELsize)%tangent 
        
        ! computing matrix process
        do i=1,PANELsize
            
            normal_ith   = PANEL_array(i)%normal
            vortex_value = 0.0
                
            ! no penetration condition on the airfoil panels induced by itself
            do j=1,PANELsize
                matrix(i,j)  = dot_product(normal_ith,integral(PANEL_array(i),PANEL_array(j),'source'))
                vortex_value = vortex_value + dot_product(normal_ith,integral(PANEL_array(i),PANEL_array(j),'vortex'))
            end do
           
            matrix(i,PANELsize+1) = vortex_value
            
            ! no penetration condition on the airfoil panels induced by the ground panels
            do j=PANELsize+2,PANELsize+GROUNDsize+1
                matrix(i,j) = dot_product(normal_ith,integral(PANEL_array(i),GROUNDpanel(j-PANELsize-1),'source')) 
            end do 

        end do 

        ! KUTTA-JOUKOWSKY condition
        vortex_value = 0.0

        do j=1,PANELsize
            matrix(PANELsize+1,j)  = dot_product(tangent_first,integral(PANEL_array(1),PANEL_array(j),'source')) + &
                                     dot_product(tangent_last,integral(PANEL_array(PANELsize),PANEL_array(j),'source'))
            vortex_value = vortex_value + dot_product(tangent_first,integral(PANEL_array(1),PANEL_array(j),'vortex')) + &
                           dot_product(tangent_last,integral(PANEL_array(PANELsize),PANEL_array(j),'vortex'))     
        end do     

        matrix(PANELsize+1,PANELsize+1) = vortex_value 

        do j=1,GROUNDsize 
            matrix(PANELsize+1,j+PANELsize+1) = dot_product(tangent_first,integral(PANEL_array(1),GROUNDpanel(j),'source')) + &
                                            dot_product(tangent_last,integral(PANEL_array(PANELsize),GROUNDpanel(j),'source'))
        end do
        
        ! no penetration conditions on ground panels
        do i=PANELsize+2,PANELsize+GROUNDsize+1
            
            normal_ith   = GROUNDpanel(i)%normal
            vortex_value = 0.0  
           
            do j=1,PANELsize
                matrix(i,j)  = dot_product(normal_ith,integral(GROUNDpanel(i-PANELsize-1),PANEL_array(j),'source'))
                vortex_value = vortex_value + dot_product(normal_ith,integral(GROUNDpanel(i-PANELsize-1),PANEL_array(j),'vortex'))
            end do
            
            matrix(i,PANELsize+1) = vortex_value 
            
            do j=PANELsize+2,PANELsize+GROUNDsize+1
                matrix(i,j) = dot_product(normal_ith,integral(GROUNDpanel(i-PANELsize-1), GROUNDpanel(j-PANELsize-1),'source'))
            end do 

        end do 
        
        call write_formatted('[','normal','OK','green','] -- matrix generated','normal') 

    end subroutine computeGROUNDmatrix

    subroutine computeGROUNDvector(vector,PANEL_array,GROUNDpanel,PANELsize,GROUNDsize,alpha,V)
        use PANEL_object        
        use FOUL 
        use cp
        implicit none             
        
        integer(kind=4),intent(in)                   :: PANELsize
        integer(kind=4),intent(in)                   :: GROUNDsize 
        type(panel),dimension(PANELsize),intent(in)  :: PANEL_array
        type(panel),dimension(GROUNDsize),intent(in) :: GROUNDpanel
        real(kind=8),intent(in)                      :: alpha
        real(kind=8),intent(in)                      :: V
        real(kind=8),dimension(2)                    :: velocity    
        real(kind=8),dimension(:),allocatable        :: vector
        integer(kind=4)                              :: i, j
        
        ! allocation of vector
        allocate(vector(PANELsize+GROUNDsize+1))

        velocity(1) = V*cos(alpha)
        velocity(2) = V*sin(alpha)
        
        ! known conditions on the airfoil panels
        do j=1,PANELsize 
            vector(j) = dot_product(PANEL_array(j)%normal,velocity)        
        end do 
        
        ! known condition for the KUTTA_JOUKOWSKY theorem
        vector(PANELsize+1) = dot_product(PANEL_array(1)%tangent,velocity) + dot_product(PANEL_array(PANELsize)%tangent,velocity)     

        ! known conditions for the ground panels
        do j=PANELsize+2,PANELsize+GROUNDsize+1
            vector(j) = dot_product(GROUNDpanel(j-PANELsize)%normal,velocity) 
        end do  

        call write_formatted('[','normal','OK','green','] -- vector generated','normal') 
        
        vector = - vector

    end subroutine computeGROUNDvector
    
    subroutine solveGROUND(solution,GROUNDmatrix,GROUNDvector,PANELsize,GROUNDsize)
        use cp
        use PANEL_object 
        implicit none 
        
        integer(kind=4),intent(in)                                            :: PANELsize
        integer(kind=4),intent(in)                                            :: GROUNDsize 
        real(kind=8),dimension(PANELsize+GROUNDsize+1,PANELsize+GROUNDsize+1) :: GROUNDmatrix
        real(kind=8),dimension(PANELsize+GROUNDsize+1)                        :: GROUNDvector
        real(kind=8),dimension(PANELsize+GROUNDsize+1,PANELsize+GROUNDsize+1) :: M    
        integer(kind=4),dimension(PANELsize+GROUNDsize+1)                     :: info
        integer(kind=4),dimension(PANELsize+GROUNDsize+1)                     :: pivoting_vec
        real(kind=8),dimension(:),allocatable                                 :: solution  
        
        ! allocation of solution 
        allocate(solution(PANELsize+GROUNDsize+1))

        ! allocating new variables that stores initial data
        solution = GROUNDvector  
        M        = GROUNDmatrix

        ! solving problem with lapack libraries 
        call DGESV(PANELsize+GROUNDsize+1, 1, M, PANELsize+GROUNDsize+1, pivoting_vec, solution, PANELsize+GROUNDsize+1, info)
        
        ! testing results
        call test_matrix(GROUNDmatrix,solution,GROUNDvector,PANELsize+GROUNDsize+1)

    end subroutine solveGROUND

    subroutine computeGROUNDfield(solution,PANEL_array,GROUNDpanel,PANELsize,GROUNDsize,P0,alpha,V,rho)
        use PANEL_object 
        use cp 
        use math_module
        use FOUL
        implicit none 

        ! generating a grid of measure points 
        ! this grid can be expressed in a monodimensional array that contains 2 dimensional vectors 
        !   the procedure can be possible if it's used a type variable that describes a vector
        type :: array_type  
            real(kind=8),pointer :: coords(:)
        end type array_type 
        
        integer(kind=4),intent(in)                                :: PANELsize
        integer(kind=4),intent(in)                                :: GROUNDsize 
        type(panel),dimension(PANELsize),intent(in)               :: PANEL_array
        type(panel),dimension(GROUNDsize),intent(in)              :: GROUNDpanel
        type(panel)                                               :: dummy_panel
        real(kind=8),dimension(PANELsize+GROUNDsize+1),intent(in) :: solution 
        real(kind=8),dimension(2)                                 :: velocity
        real(kind=8),dimension(2)                                 :: vel
        integer(kind=4)                                           :: i, j, k       
        real(kind=8)                                              :: deltax, deltay
        real(kind=8)                                              :: x_start, x_end
        real(kind=8)                                              :: y_start, y_end
        real(kind=8),intent(in)                                   :: V
        real(kind=8),intent(in)                                   :: alpha
        real(kind=8)                                              :: pressure
        real(kind=8)                                              :: norm_vel
        real(kind=8),intent(in)                                   :: P0
        real(kind=8),intent(in)                                   :: rho
        type(array_type),dimension(:,:),allocatable               :: grid
        integer(kind=4)                                           :: nrows, ncols
        
        ! declaring grid dimensions
        nrows = 1e+2
        ncols = 2e+2

        ! grid allocation process in memory
        allocate(grid(nrows,ncols))
        do i=1,nrows
            do j=1,ncols
                allocate(grid(i,j)%coords(2))
            end do
        end do
            
        ! allocating grid dimensions
        x_start = -0.2
        x_end   =  1.5
        y_start = -0.2 
        y_end   =  0.2

        ! grid dimension allocation process
        deltax = (x_end - x_start)/ncols
        deltay = (y_end - y_start)/nrows        
        ! grid data allocation process
        do i=1,nrows 
            do j=1,ncols
                grid(i,j)%coords(1) = x_start + (j-1)*deltax
                grid(i,j)%coords(2) = y_start + (i-1)*deltay
            end do
        end do

        ! open file to save data 
        open(unit=1, file='FLOWfield.dat', status='replace')
        
        vel(1) = V*cos(alpha)
        vel(2) = V*sin(alpha) 

        do i=1,nrows 
                
             do j=1,ncols   
           
                 velocity = (/0.0, 0.0/)
                 dummy_panel%midpoint = grid(i,j)%coords
                
                do k=1,PANELsize 
                    velocity = velocity + integral(dummy_panel,PANEL_array(k),'source')*solution(k)
                    velocity = velocity + integral(dummy_panel,PANEL_array(k),'vortex')*solution(PANELsize+1)
                end do
                
                do k=1,GROUNDsize 
                    velocity = velocity + integral(dummy_panel,GROUNDpanel(k),'source')*solution(PANELsize+1+k)
                end do

                velocity = velocity + vel
                
                norm_vel = norm(velocity)

                pressure = P0 + 0.5*rho*(V**2 - norm_vel**2)  
            
                write(1,*) dummy_panel%midpoint(1), dummy_panel%midpoint(2), &
                           velocity(1), velocity(2), norm_vel, & 
                           pressure 
            end do 
        end do
        
        call write_formatted('[','normal','OK','green','] -- velocity and pressure field computed','normal') 

        close(1)
        deallocate(grid)

    end subroutine computeGROUNDfield

end module ground_cp
