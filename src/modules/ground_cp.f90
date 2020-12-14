module ground_cp

    contains 
    
    subroutine generate_ground(GROUNDpanel,GROUNDsize) 
    ! this subroutine compute the ground panelization 
        use PANEL_object
        use FOUL 

        implicit none 
        
        integer(kind=4),intent(inout)                      :: GROUNDsize  ! # of panel to discretize the ground 
        type(panel),dimension(:),allocatable,intent(inout) :: GROUNDpanel ! ground panels --> objects that discritize the ground 
        real(kind=8)                                       :: GROUNDstart ! start of the ground 
        real(kind=8)                                       :: GROUNDend   ! end of the ground 
        real(kind=8)                                       :: theta       ! inclination of each panel wrt the x axis 
        real(kind=8)                                       :: deltax      ! length of each panel
        real(kind=8)                                       :: Y           ! y position wrt the origin 
        integer(kind=4)                                    :: j
        
        ! allocating angle 
        theta = 0.0
        
        call write_formatted('GENERATE GROUND PANELIZATION','yellow')

        ! asking # of panels to discretize the ground
        print*, 'type the number of panels for discretize the ground effect'
        read*,   GROUNDsize
        
        ! allocate GROUNpanel
        allocate(GROUNDpanel(GROUNDsize))
        
        ! GROUNDpanel extremes
        print*, 'type the x coordinate for the start point of the ground panelization'
        read*,   GROUNDstart
        print*, 'type the x coordinate for the end point of the ground panelization'
        read*,   GROUNDend     

        ! panel length
        deltax = (GROUNDend - GROUNDstart)/GROUNDsize
        
        ! ground panels position 
        print*, 'type ground position'
        read*,   Y
        
        ! allocation process 
        do j=1,GROUNDsize
            call GROUNDpanel(j)%set_id(j)
            call GROUNDpanel(j)%set_coords(GROUNDstart + (j-1)*deltax, Y, GROUNDstart + j*deltax, Y)
            call GROUNDpanel(j)%compute_midpoint('noprt')
            GROUNDpanel(j)%normal   = (/0.0, 1.0/)
            GROUNDpanel(j)%tangent  = (/1.0, 0.0/)
            GROUNDpanel(j)%length   = deltax
            GROUNDpanel(j)%angle    = theta 
            GROUNDpanel(j)%ROT(1,1) = 1     
            GROUNDpanel(j)%ROT(1,2) = 0   
            GROUNDpanel(j)%ROT(2,1) = 0   
            GROUNDpanel(j)%ROT(2,2) = 1    
         end do  
        
         ! saving data in file
         open(unit=1, file='GROUNDdata.dat',   status='replace')
         open(unit=2, file='GROUNDpanels.dat', status='replace')        

         write(1,*) 'GROUND panel objects'
         
         ! saving process 
         do j=1,GROUNDsize 
            call GROUNDpanel(j)%saving(1)
            write(2,*) GROUNDpanel(j)%coords1            
         end do 

         write(2,*) GROUNDpanel(GROUNDsize)%coords2

         call write_formatted('[','normal','OK','green','] -- ground panels generated','normal') 
        
         close(1)
         close(2)
    end subroutine generate_ground 
        
    subroutine computeGROUNDmatrix(matrix,PANEL_array,GROUNDpanel,PANELsize,GROUNDsize,panel_type)
    ! this subroutine compute the system matrix for the airfoil + ground panelization 
        use PANEL_object
        use FOUL
        use cp 
        
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
        character(len=6),intent(in)                  :: panel_type
        
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
            
            ! induced velocity by vortex distribution          
            matrix(i,PANELsize+1) = vortex_value
            
            ! no penetration condition on the airfoil panels induced by the ground panels
            do j=PANELsize+2,PANELsize+GROUNDsize+1
                matrix(i,j) = dot_product(normal_ith,integral(PANEL_array(i),GROUNDpanel(j-(PANELsize+1)),panel_type)) 
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
            matrix(PANELsize+1,j+PANELsize+1) = dot_product(tangent_first,integral(PANEL_array(1),GROUNDpanel(j),panel_type)) + &
                                            dot_product(tangent_last,integral(PANEL_array(PANELsize),GROUNDpanel(j),panel_type))
        end do
        
        ! no penetration conditions on ground panels
        do i=PANELsize+2,PANELsize+GROUNDsize+1
            
            normal_ith   = GROUNDpanel(i-(PANELsize+1))%normal
            vortex_value = 0.0  
            
            ! induction by airfoil panels on the ground panels 
            do j=1,PANELsize
                matrix(i,j)  = dot_product(normal_ith,integral(GROUNDpanel(i-(PANELsize+1)),PANEL_array(j),'source'))
                vortex_value = vortex_value + dot_product(normal_ith,integral(GROUNDpanel(i-(PANELsize+1)),PANEL_array(j),'vortex'))
            end do
            
            ! induction by airfoil distribution of vorticity on the ground panels 
            matrix(i,PANELsize+1) = vortex_value 
            
            ! induction by gruond panels on themselves
            do j=PANELsize+2,PANELsize+GROUNDsize+1
                matrix(i,j) = dot_product(normal_ith,integral(GROUNDpanel(i-(PANELsize+1)),GROUNDpanel(j-(PANELsize+1)),panel_type))
            end do 

        end do 
        
        call write_formatted('[','normal','OK','green','] -- matrix generated','normal') 

    end subroutine computeGROUNDmatrix

    subroutine computeGROUNDvector(vector,PANEL_array,GROUNDpanel,PANELsize,GROUNDsize,alpha,V)
    ! this subroutine computes the vector velocity on the panels induced by the streamflow 
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
        
        ! the following variable is computed with alpha in order to be more general and to allow future changes in the code 
        velocity(1) = V*cos(alpha)
        velocity(2) = V*sin(alpha)
        
        ! known conditions on the airfoil panels
        do j=1,PANELsize 
            vector(j) = dot_product(PANEL_array(j)%normal,velocity)        
        end do 
        
        ! known condition for the KUTTA-JOUKOWSKY theorem
        vector(PANELsize+1) = dot_product(PANEL_array(1)%tangent,velocity) + dot_product(PANEL_array(PANELsize)%tangent,velocity)     

        ! known conditions for the ground panels
        do j=PANELsize+2,PANELsize+GROUNDsize+1
            vector(j) = dot_product(GROUNDpanel(j-PANELsize-1)%normal,velocity) 
        end do  
        
        ! the vector variable is the known vector --> in order to solve the system it must be multiplied for -1
        vector = - vector

        call write_formatted('[','normal','OK','green','] -- vector generated','normal') 
        
    end subroutine computeGROUNDvector
    
    subroutine solveGROUND(solution,GROUNDmatrix,GROUNDvector,PANELsize,GROUNDsize)
    ! this subroutine computes the solution for the system and tests the results 
        use PANEL_object 
        use cp
        implicit none 
        
        integer(kind=4),intent(in)                                            :: PANELsize
        integer(kind=4),intent(in)                                            :: GROUNDsize 
        real(kind=8),dimension(PANELsize+GROUNDsize+1,PANELsize+GROUNDsize+1) :: GROUNDmatrix
        real(kind=8),dimension(PANELsize+GROUNDsize+1,PANELsize+GROUNDsize+1) :: M    
        real(kind=8),dimension(PANELsize+GROUNDsize+1)                        :: GROUNDvector
        integer(kind=4),dimension(PANELsize+GROUNDsize+1)                     :: pivoting_vec
        integer(kind=4),dimension(PANELsize+GROUNDsize+1)                     :: info
        real(kind=8),dimension(:),allocatable                                 :: solution  
        
        ! allocation of solution 
        allocate(solution(PANELsize+GROUNDsize+1))

        ! allocating new variables that stores initial data
        solution = GROUNDvector  
        M        = GROUNDmatrix

        ! solving problem with LAPACK libraries 
        call DGESV(PANELsize+GROUNDsize+1, 1, M, PANELsize+GROUNDsize+1, pivoting_vec, solution, PANELsize+GROUNDsize+1, info)
        
        ! testing results allows to know the quality of results --> it's physically the normal velocity from the panels 
        call test_matrix(GROUNDmatrix,solution,GROUNDvector,PANELsize+GROUNDsize)

    end subroutine solveGROUND

    subroutine computeGROUNDfield(solution,PANEL_array,GROUNDpanel,PANELsize,GROUNDsize,P0,alpha,V,rho,panel_type)
    ! this suboruitne computes the velocity and pressure field in the system 
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
        character(len=6),intent(in)                               :: panel_type 

        ! declaring grid dimensions
        ncols = 2.5e+2
        nrows = 1.5e+2

        ! grid allocation process in memory
        allocate(grid(nrows,ncols))
        do i=1,nrows
            do j=1,ncols
                allocate(grid(i,j)%coords(2))
            end do
        end do
           
        ! allocating grid dimensions
        x_start = GROUNDpanel(1)%coords1(1) 
        x_end   = GROUNDpanel(GROUNDsize)%coords2(1)
        y_start = GROUNDpanel(1)%coords1(2) 
        y_end   = PANEL_array(PANELsize/2)%coords1(2) + 0.2

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
        open(unit=1, file='FLOWfieldGE.dat', status='replace')
        
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
                
                ! computing induction by the airfoil panels 
                ! source distribution induction 
                ! vortex distribution induction 
                do k=1,PANELsize 
                    velocity = velocity + integral(dummy_panel,PANEL_array(k),'source')*solution(k)
                    velocity = velocity + integral(dummy_panel,PANEL_array(k),'vortex')*solution(PANELsize+1)
                end do
                
                ! computing induction by the ground panels 
                ! source distribuiton induction 
                do k=1,GROUNDsize 
                    velocity = velocity + integral(dummy_panel,GROUNDpanel(k),panel_type)*solution(PANELsize+1+k)
                end do
                
                ! adding the flow stream velocity 
                velocity = velocity + vel
                
                ! computing the norm of velocity 
                norm_vel = norm(velocity)
                
                ! computing the pressure with BERNOULLI's theorem 
                pressure = P0 + 0.5*rho*(V**2 - norm_vel**2)  
                
              !  ! plotting condition --> avoid to display non physical values for pressure 
              !  if(pressure < 0)then 
              !      pressure = - 0.5
              !  end if

                ! writing data in file 
                write(1,*) dummy_panel%midpoint(1), dummy_panel%midpoint(2), &
                           velocity(1), velocity(2), norm_vel, & 
                           pressure 
            end do 
        end do
        
        call write_formatted('[','normal','OK','green','] -- velocity and pressure field computed                ','normal', &
                              '--> FLOWfieldGE.dat','normal') 
        
        ! plot data 
        call system('gnuplot -p PRESSUREfieldGE.plt')

        close(1)
        deallocate(grid)
    end subroutine computeGROUNDfield
    
    subroutine compute_airfoilFIELD(solution,PANEL_array,GROUNDpanel,PANELsize,GROUNDsize,cp_vec,P0,alpha,V,rho,panel_type)  
        use PANEL_object 
        use math_module
        use FOUL
        use cp 
        
        implicit none 
        
        integer(kind=4),intent(in)                                :: PANELsize
        integer(kind=4),intent(in)                                :: GROUNDsize 
        type(panel),dimension(PANELsize),intent(in)               :: PANEL_array
        type(panel),dimension(GROUNDsize),intent(in)              :: GROUNDpanel
        real(kind=8),dimension(PANELsize+GROUNDsize+1),intent(in) :: solution
        real(kind=8),dimension(PANELsize),intent(out)             :: cp_vec            
        real(kind=8),intent(in)                                   :: V
        real(kind=8),intent(in)                                   :: alpha
        real(kind=8),intent(in)                                   :: P0
        real(kind=8),intent(in)                                   :: rho
        real(kind=8),dimension(2,2)                               :: ROT  
        real(kind=8),dimension(2)                                 :: velocity
        real(kind=8),dimension(2)                                 :: vel
        real(kind=8)                                              :: pressure
        real(kind=8)                                              :: norm_vel
        real(kind=8),dimension(PANELsize)                         :: cp_coeff
        real(kind=8)                                              :: normal_velocity
        real(kind=8)                                              :: circulation
        integer(kind=4)                                           :: j, k       
        character(len=6),intent(in)                               :: panel_type 
        
        ! initializing circulation
        circulation = 0.0

        ! open file to save data 
        open(unit=1, file='FLOWfield_airfoilGE.dat', status='replace')
        open(unit=2, file='cp_dataGE.dat', status='replace')
        
        ! exernal velocity 
        vel(1) = V*cos(alpha)
        vel(2) = V*sin(alpha) 
        
        ! rotation matrix --> alpha angle > 0 if it's counterclockwise
        ! -- the sine and cosine terms respect the direct rotation matrix description 
        ROT(1,1) =  cos(alpha)
        ROT(2,1) =  sin(alpha) 
        ROT(1,2) = -sin(alpha)
        ROT(2,2) =  cos(alpha)      

        do j=1,PANELsize 

            ! updating circulation 
            circulation = circulation + solution(PANELsize+1)*PANEL_array(j)%get_length()
                
            velocity = (/0.0, 0.0/)
            
            ! computing airfoil's panels induction 
            do k=1,PANELsize 
                velocity = velocity + integral(PANEL_array(j),PANEL_array(k),'source')*solution(k)
                velocity = velocity + integral(PANEL_array(j),PANEL_array(k),'vortex')*solution(PANELsize+1)
            end do
            
            ! computing ground panels induction 
            do k=1,GROUNDsize 
                velocity = velocity + integral(PANEL_array(j),GROUNDpanel(k),panel_type)*solution(PANELsize+1+k)
            end do

            velocity = velocity + vel
            
            norm_vel = norm(velocity)

            pressure = P0 + 0.5*rho*(V**2 - norm_vel**2)

            cp_coeff(j) = 1 - (norm_vel/V)**2  
                
            normal_velocity = dot_product(PANEL_array(j)%normal,velocity)

            write(1,*) PANEL_array(j)%midpoint, matmul(ROT,PANEL_array(j)%midpoint), &
                       velocity(1), velocity(2), norm_vel, & 
                       pressure   , cp_coeff   , normal_velocity       
            
        end do
        
        call write_formatted('[','normal','OK','green','] -- velocity and pressure field around airfoil computed ','normal', & 
                             '--> FLOWfield_airfoil.dat','normal') 
        
        print*, 'total circulation (GAMMA) = ', circulation 
        
        ! writing Cp data in a separate file in order to separate the UPPER and LOWER surface in plotting procedure 
        do j=1,PANELsize/2
            write(2,*) PANEL_array(j)%midpoint(1),             - cp_coeff(j), & 
                       PANEL_array(j+PANELsize/2)%midpoint(1), - cp_coeff(PANELsize/2+j)
        end do
                
        close(1)
        close(2)
        
        cp_vec = cp_coeff

        ! plot data 
        call system('gnuplot -p CP_plotGE.plt')

    end subroutine compute_airfoilFIELD

end module ground_cp
