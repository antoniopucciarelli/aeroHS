module cp 

    contains 

    subroutine manage_data(airfoil,MEAN_array,PANEL_array,MEANsize,PANELsize,filename)
    ! this subroutine reads data from .dat file created with airfoilgenerator.f90 and stores into new 
    ! meanline and panel object arrays
        use MEANline_object
        use PANEL_object
        use AIRFOIL_object
        use FOUL
        use discretization_module
        implicit none

        type(NACA_airfoil),intent(inout)                    :: airfoil
        type(panel),dimension(:),allocatable,intent(out)    :: PANEL_array
        type(MEANline),dimension(:),allocatable,intent(out) :: MEAN_array
        character(len=30)                                   :: filename
        character(len=100)                                  :: str 
        character(len=1)                                    :: delimiter
        character(len=15)                                   :: str_num
        integer(kind=4)                                     :: num, i 
        integer(kind=4),intent(out)                         :: MEANsize, PANELsize
        real(kind=8)                                        :: panel_tangent_module
        character(len=30)                                   :: filename1
        character(len=30)                                   :: filename2
        character(len=30)                                   :: filename3

        print*, 'type airfoil name to analize -- NACA-4digits'
        read*, filename
        print*, new_line('(A)')
        
        filename(9:30) = '_airfoil_MEAN.dat' 

        ! setting airfoil name
        airfoil%airfoilname = filename(1:8)
        
        open(unit=1, file=filename, status='old', action='read')
        
        ! getting started with number of mean-line objects
        ! number of discretization points process
        read(1,'(A)') str
        delimiter = ':'
        read(str(index(str,delimiter)+1:),*) num 
        
        allocate(MEAN_array(1:num))
        MEANsize = num

        airfoil%n_points = num

        read(1,'(A)') str
        read(1,'(A)') str

        do i=1,num 
            
            ! ID reading process
            read(1,'(A)') str
            read(1,'(T30, I4)') MEAN_array(i)%id
            
            ! POINT COORDINATES reading process
            read(1,'(T28, 2F8.4)') MEAN_array(i)%coords(1), MEAN_array(i)%coords(2)
            
            ! THICKNESS reading process
            read(1,'(T28, F8.4)') MEAN_array(i)%thickness
            
            ! THETA reading process
            read(1,'(T28, F8.4)') MEAN_array(i)%theta
            
            read(1,'(A)') str
            read(1,'(A)') str
            
        end do

        write(str_num,*) num
        call write_formatted('[','normal','OK','green','] -- mean-line object loaded => # of objects','normal',str_num,'yellow')
        
        close(1)
        
        filename(9:) = '_airfoil_PANEL.dat'
        
        open(unit=1, file=filename, status='old', action='read')
        
        ! getting started with number of panel objects
        ! number of discretization points process
        read(1,'(A)') str
        delimiter = ':'
        read(str(index(str,delimiter)+1:),*) num 
        
        allocate(PANEL_array(1:num))
        PANELsize = num 

        ! setting airfoil number of division points
        airfoil%n_points = (PANELsize + 2)/2 

        read(1,'(A)') str
        read(1,'(A)') str

        do i=1,num 

            ! ID reading process
            read(1,'(A)') str
            read(1,'(T30, I4)') PANEL_array(i)%id
            
            ! LENGTH reading process
            read(1,'(T28, F12.8)') PANEL_array(i)%length
            
            ! MID-POINT COORDINATES reading process
            read(1,'(T28, 2F12.8)') PANEL_array(i)%midpoint(1), PANEL_array(i)%midpoint(2)
            
            ! STARTING-POINT COORDINATES reading process
            read(1,'(T28, 2F12.8)') PANEL_array(i)%coords1(1), PANEL_array(i)%coords1(2)
            
            ! ENDING-POINT COORDINATES reading process
            read(1,'(T28,2F12.8)') PANEL_array(i)%coords2(1), PANEL_array(i)%coords2(2)
            
            ! TANGENT VECTOR COORDINATES reading process
            read(1,'(T28, 2F12.8)') PANEL_array(i)%tangent(1), PANEL_array(i)%tangent(2)
            
            ! NORMAL VECTOR COORDINATES reading process
            read(1,'(T28, 2F12.8)') PANEL_array(i)%normal(1), PANEL_array(i)%normal(2)

            ! PANEL INCLINATION reading process
            read(1,'(T28,  F12.8)') PANEL_array(i)%angle   
            
            ! PANEL ROTATION MATRIX reading process 
            read(1,'(T28, 4F12.8)') PANEL_array(i)%ROT(1,1), PANEL_array(i)%ROT(2,1), & 
                                    PANEL_array(i)%ROT(1,2), PANEL_array(i)%ROT(2,2)

            ! PANEL POSITION reading process
            read(1, '(T28,    A2)') PANEL_array(i)%POS

            ! run check on the tangent vector
            panel_tangent_module = sqrt((PANEL_array(i)%tangent(1))**2 + (PANEL_array(i)%tangent(2))**2)
            if (panel_tangent_module < (1-1e-3) .or. panel_tangent_module > (1+1e-3)) then
                write(str_num,*) i
                call write_formatted('[','normal','WARNING','red','] -- panel problem => # ','normal',str_num,'yellow')
                print*, 'tangentx = ',PANEL_array(i)%tangent(1), 'tangenty = ', PANEL_array(i)%tangent(2)
                print*, 'panel # ', i ,' tangent vector module =', panel_tangent_module
            end if
            
            read(1,'(A)') str
            read(1,'(A)') str
            
        end do
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!! GRAPHICS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! already exists filename variable
        filename1  = 'GNUplot_coord_data.dat'
        filename2 = 'GNUplot_mean_data.dat'
        filename3 = 'GNUplot_tg_norm.dat'
        call GNUplot_print(airfoil,PANEL_array,MEAN_array,filename1,filename2,filename3)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!! GRAPHICS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        write(str_num,*) num
        call write_formatted('[','normal','OK','green','] -- initial condition and geometry','normal')
        call write_formatted('[','normal','OK','green','] -- panel object loaded     => # of objects','normal',str_num,'yellow')
        
        close(1)
        
    end subroutine manage_data 

    recursive subroutine compute_COEFF(beta,r1mod,r2mod,PANEL_ith,PANEL_jth)
    ! this suboutine computes the variables needed to compute the integral for source and vortex panels on the airfoil

        use math_module
        use PANEL_object

        implicit none 
    
        ! initializing variables 
        ! in the previous version of this subroutine the [N1, N2, T, midpoint] vectors were
        ! computed from 0 --> this time the suborutine uses the values of PANEL_jth and PANEL_ith
        ! in orfer to describe the [N1, N2, T, midpoint] vectors 
        type(panel),intent(in)      :: PANEL_ith, PANEL_jth
        real(kind=8),dimension(2)   :: N1, N2, T
        real(kind=8),dimension(2)   :: r1, r2, midpoint
        real(kind=8),dimension(2,2) :: ROT                             
        real(kind=8)                :: beta
        real(kind=8)                :: r1mod, r2mod
        real(kind=8)                :: acosbeta
        real(kind=8)                :: theta

        ! compute panel ith midpoint and declaring N1 and N2 vectors
        N1 = PANEL_jth%coords1
        N2 = PANEL_jth%coords2
        T  = PANEL_ith%midpoint

        ! compute panel jth midpoint
        midpoint = PANEL_jth%midpoint
        
        ! compute panel inclination angle
        theta = PANEL_jth%get_angle()

        ! declaring rotation matrix --> have a look at compute_ROT in PANEL_object.f90
        ROT = PANEL_jth%ROT

        ! compute distance vectors between target point and N1, N2
        r1 = T - N1
        r2 = T - N2
        ! compute r1, r2 module
        r1mod = norm(r1)
        r2mod = norm(r2)
        ! compute beta angle
        acosbeta = dot_product(r1,r2)/(r1mod * r2mod)
        if(acosbeta >= 1)then
            acosbeta =  1
        else if(acosbeta <= -1)then
            acosbeta = -1
        end if

        beta = acos(acosbeta)
    
        ! determine beta angle sign
        ! 1 step - rotation of points by the panel inclination angle
        T        = matmul(ROT, T)
        midpoint = matmul(ROT, midpoint)
    
        ! 2 step - points translation in order to have the midpoint of the panel at point (0, 0)
        T  = T  - midpoint 
    
        ! 3 step - determine if the angle is > 0 < by the position of T with respect to the panel
        if(T(2)>=0)then
            beta =  abs(beta)
        else
            beta = -abs(beta)
        end if

    end subroutine compute_COEFF

    function integral(PANEL_ith,PANEL_jth,input_type)
    ! this funciton computes the integral for the source and vortex panel induced on PANEL_ith by PANEL_jth
        use PANEL_object
        use math_module
        
        implicit none 

        character(len=6),intent(in) :: input_type
        type(panel),intent(in)      :: PANEL_ith
        type(panel),intent(in)      :: PANEL_jth
        real(kind=8)                :: beta
        real(kind=8)                :: r1mod, r2mod
        real(kind=8),dimension(2)   :: integral
        real(kind=8),dimension(2)   :: tangent, normal
        
        call compute_COEFF(beta,r1mod,r2mod,PANEL_ith,PANEL_jth)
        
        ! declaration of jth PANEL normal and tangent vectors 
        tangent = PANEL_jth%tangent
        normal  = PANEL_jth%normal

        if(input_type == 'source')then
            integral = - 1/(2*pi)*log(r2mod/r1mod)*tangent + 1/(2*pi)*beta*normal
        else if(input_type == 'vortex')then
            integral = - 1/(2*pi)*beta*tangent - 1/(2*pi)*log(r2mod/r1mod)*normal
        end if

    end function integral

    subroutine compute_matrix(matrix,PANEL_array,PANELsize)
    ! this subroutine generates the problem matrix that its used for compute the vortices strength 
    ! every matrix element can be seen as a weight of the panel vorticity into the airfoil system
    ! its is important to keep in mind that we use the midpoint of every panel to compute the weigths/elements of the matrix    
        use PANEL_object
        use math_module
        use FOUL
        implicit none 

        integer(kind=4),intent(in)                          :: PANELsize
        integer(kind=4)                                     :: i, j 
        type(panel),dimension(PANELsize),intent(in)         :: PANEL_array
        real(kind=8),dimension(:,:),allocatable,intent(out) :: matrix
        real(kind=8)                                        :: vortex_value
        real(kind=8),dimension(2)                           :: normal_ith
        real(kind=8),dimension(2)                           :: tangent_first
        real(kind=8),dimension(2)                           :: tangent_last

        allocate(matrix(PANELsize+1,PANELsize+1))

        do i=1,PANELsize
            
            ! declaring normal vector
            normal_ith  = PANEL_array(i)%normal

            vortex_value = 0.0

            do j=1,PANELsize
                matrix(i,j)  = dot_product(normal_ith,integral(PANEL_array(i),PANEL_array(j),'source'))
                vortex_value = vortex_value + dot_product(normal_ith,integral(PANEL_array(i),PANEL_array(j),'vortex'))
            end do

            matrix(i,PANELsize+1) = vortex_value

        end do 
        
        ! declaring tangent and normal vectors
        ! first panel -- AIRFOIL UPPER PART 
        tangent_first = PANEL_array(1)%tangent
        ! last panel  -- AIRFOIL LOWER PART 
        tangent_last  = PANEL_array(PANELsize)%tangent 
        
        vortex_value = 0.0

        ! KUTTA-JOUKOSKY condition        
        do j=1,PANELsize

            matrix(PANELsize+1,j) = dot_product(tangent_first,integral(PANEL_array(1),PANEL_array(j),'source')) + &
                                    dot_product(tangent_last,integral(PANEL_array(PANELsize),PANEL_array(j),'source'))
            
            vortex_value = vortex_value + dot_product(tangent_first,integral(PANEL_array(1),PANEL_array(j),'vortex')) + &
                           dot_product(tangent_last,integral(PANEL_array(PANELsize),PANEL_array(j),'vortex'))

        end do 

        matrix(PANELsize+1,PANELsize+1) = vortex_value

        call write_formatted('[','normal','OK','green','] -- matrix created','normal')    
        
    end subroutine compute_matrix

    subroutine compute_vector(V,alpha,vector,PANEL_array,PANELsize)
    ! this subroutine computes the known vector from the free stram velocity and the panel geometry of the airfoil
        use PANEL_object
        use math_module
        use FOUL

        implicit none

        integer(kind=4),intent(in)                    :: PANELsize
        type(panel),dimension(PANELsize),intent(in)   :: PANEL_array
        real(kind=8),dimension(PANELsize),intent(out) :: vector
        real(kind=8),dimension(2)                     :: V_vec
        real(kind=8),intent(in)                       :: V
        real(kind=8)                                  :: alpha
        integer(kind=4)                               :: i
        real(kind=8),dimension(2)                     :: tangent_ith, normal_ith
        real(kind=8),dimension(2)                     :: tangent_first, normal_first
        real(kind=8),dimension(2)                     :: tangent_last, normal_last

        ! alpha is already expressed in radiants -- setting_properties subroutine
        V_vec(1) = V*cos(alpha)
        V_vec(2) = V*sin(alpha)

        do i=1,PANELsize

            ! declaring normal vector
            normal_ith = PANEL_array(i)%normal

            vector(i)  = dot_product(V_vec,normal_ith)
        
        end do

        ! declaring tangent vectors
        tangent_first = PANEL_array(1)%tangent
        tangent_last  = PANEL_array(PANELsize)%tangent
       
        vector(PANELsize+1) = dot_product(V_vec,tangent_first) + dot_product(V_vec,tangent_last)

        call write_formatted('[','normal','OK','green','] -- velocity vector created','normal')

        vector = - vector 

    end subroutine compute_vector

    function solve_matrix(matrix,vector,PANELsize)
    ! this subroutine computes the solution => vorticity and source strenght of each panel
        implicit none 
        
        integer(kind=4),intent(in)                                 :: PANELsize
        integer(kind=4),dimension(PANELsize+1)                     :: pivoting_vec
        integer(kind=4),dimension(PANELsize+1)                     :: info
        real(kind=8),dimension(PANELsize+1)                        :: solve_matrix
        real(kind=8),dimension(PANELsize+1),intent(in)             :: vector 
        real(kind=8),dimension(PANELsize+1,PANELsize+1),intent(in) :: matrix
        real(kind=8),dimension(PANELsize+1,PANELsize+1)            :: M

        solve_matrix = vector
        M            = matrix

        call DGESV(PANELsize+1, 1, M, PANELsize+1, pivoting_vec, solve_matrix, PANELsize+1, info)

        call test_matrix(matrix,solve_matrix,vector,PANELsize)

    end function solve_matrix

    subroutine test_matrix(matrix,solution,vector,PANELsize)
    ! this subroutine tests the solution of the vortex-source problem
        use FOUL
        implicit none

        integer(kind=4)                                 :: PANELsize 
        real(kind=8),dimension(PANELsize+1,PANELsize+1) :: matrix
        real(kind=8),dimension(PANELsize+1)             :: solution
        real(kind=8),dimension(PANELsize+1)             :: vector
        real(kind=8),dimension(PANELsize+1)             :: test
        integer(kind=4)                                 :: i
        real(kind=8)                                    :: maxerr
        character(len=30)                               :: char_error

        test   = matmul(matrix,solution) - vector

        maxerr = maxval(abs(test), PANELsize+1)

        write(char_error,'(EN15.5)') maxerr

        if(maxerr < 1e-10)then
            call write_formatted('[','normal','OK','green','] -- matrix error  =','normal',char_error,'normal')
        else
            call write_formatted('[','normal','WARNING','red','] -- matrix error  =','normal',char_error,'normal')
        end if

    end subroutine test_matrix   

    function compute_vel(solution, PANELsize, PANEL_array, V, alpha)
    ! this subroutine computes the velocity tangent to the airfoil-body
        use PANEL_object
        use math_module
        implicit none 
    
        integer(kind=4)                                :: PANELsize
        type(panel),dimension(PANELsize),intent(in)    :: PANEL_array
        real(kind=8),dimension(PANELsize)              :: compute_vel
        real(kind=8),dimension(PANELsize+1),intent(in) :: solution
        real(kind=8),dimension(PANELsize,PANELsize+1)  :: matrix
        real(kind=8)                                   :: vortex_value = 0.0
        real(kind=8)                                   :: alpha
        real(kind=8),intent(in)                        :: V
        real(kind=8),dimension(2)                      :: V_vec
        real(kind=8),dimension(PANELsize)              :: vector
        integer(kind=4)                                :: i, j
        real(kind=8),dimension(2)                      :: tangent_ith, normal_ith
        real(kind=8)                                   :: circulation 
        
        ! setting circualtion null value 
        circulation = 0.0  

        do i=1,PANELsize

            ! declaring tangent vector
            tangent_ith = PANEL_array(i)%tangent

            do j=1,PANELsize
        
                matrix(i,j)  = dot_product(tangent_ith,integral(PANEL_array(i),PANEL_array(j),'source')) 
        
                vortex_value = vortex_value + & 
                               dot_product(tangent_ith,integral(PANEL_array(i),PANEL_array(j),'vortex'))
        
            end do

            matrix(i,PANELsize+1) = vortex_value
            vortex_value          = 0.0
        
        end do

        ! alpha is already expressed in radiants -- setting_properties subroutine
        V_vec(1) = V*cos(alpha)
        V_vec(2) = V*sin(alpha)

        do i=1,PANELsize

            ! declaring tangent vector for the panel ith
            tangent_ith = PANEL_array(i)%tangent

            vector(i)   = dot_product(V_vec,tangent_ith)
        
        end do

        do i=1,PANELsize 
            circulation = circulation + solution(PANELsize+1)*PANEL_array(i)%get_length() 
        end do

        compute_vel = matmul(matrix, solution) + vector 
        
        print*, 'total circulation (GAMMA) = ',   circulation 
        print*, 'Cl                        = ', - circulation * 2 

    end function compute_vel

    function compute_pressure(P0,V,rho,Vvec,PANELsize) 
    ! this subroutine computes the pressure distribution at the midpoins of every panel in the airfoil frame
        implicit none   

        integer(kind=4),intent(in)        :: PANELsize
        real(kind=8),intent(in)           :: P0
        real(kind=8),intent(in)           :: V
        real(kind=8),intent(in)           :: rho
        real(kind=8),dimension(PANELsize) :: Vvec
        real(kind=8),dimension(PANELsize) :: compute_pressure
        
        compute_pressure = P0 + 0.5*rho*(V**2 - Vvec**2)

    end function compute_pressure

    function compute_cp(Vvec,V,PANELsize)
    ! this subroutine computes the Cp distribution at the midpoins of every panel in the airfoil frame
        implicit none 
        
        integer(kind=4),intent(in)        :: PANELsize
        real(kind=8),intent(in)           :: V
        real(kind=8),dimension(PANELsize) :: Vvec
        real(kind=8),dimension(PANELsize) :: compute_cp

        compute_cp = 1 - (Vvec / V)**2

    end function compute_cp

    function compute_cl(cp_vec,PANEL_array,PANELsize)
    ! this function compute the cl coefficient at a defined angle of attack 
        use PANEL_object
        use math_module
        implicit none 

        integer(kind=4),intent(in)                   :: PANELsize 
        type(panel),dimension(PANELsize),intent(in)  :: PANEL_array
        real(kind=8),dimension(PANELsize),intent(in) :: cp_vec
        real(kind=8)                                 :: compute_cl 
        real(kind=8)                                 :: panel_length
        real(kind=8)                                 :: S 
        real(kind=8),dimension(2)                    :: N1 
        real(kind=8),dimension(2)                    :: N2
        integer(kind=4)                              :: i

        compute_cl = 0.0
        
        do i=1,PANELsize 

            panel_length = PANEL_array(i)%get_length()

            if(PANEL_array(i)%get_position() == 'UP')then
                panel_length = - panel_length
            end if

            compute_cl = compute_cl + panel_length * cp_vec(i)

        end do 

        S = norm(PANEL_array(PANELsize/2+1)%coords1 - PANEL_array(1)%coords1)

        compute_cl = compute_cl / S

    end function compute_cl

    subroutine compute_angle(alpha_angle,start_angle,end_angle,angle_num)
    ! this subroutine computes the angle vector for the computation of Cl    
        use math_module 

        implicit none 

        integer(kind=4),intent(in)                        :: start_angle
        integer(kind=4),intent(in)                        :: end_angle
        real(kind=8),dimension(:),allocatable,intent(out) :: alpha_angle
        integer(kind=4),intent(in)                        :: angle_num
        integer(kind=4)                                   :: i
        real(kind=8)                                      :: k
        
        allocate(alpha_angle(angle_num))

        k = real((end_angle - start_angle),8) /real(angle_num,8)

        do i=1,angle_num
            
            alpha_angle(i) = (real(start_angle,8) + k*real(i,8))/180*pi  

        end do 

        alpha_angle(1)         = real(start_angle,8)/180*pi
        alpha_angle(angle_num) = real(end_angle,8)/180*pi
    
    end subroutine compute_angle
    
    subroutine CLalpha(start_angle,end_angle,angle_num)
    ! this subroutine computes the cl alpha diagram for the airfoil 
        use AIRFOIL_object
        use FOUL
        use PANEL_object
        use MEANline_object
        use plot
        use math_module
        use ask_module

        implicit none 

        integer(kind=4)                                 :: start_angle
        integer(kind=4)                                 :: end_angle
        integer(kind=4)                                 :: angle_num                             
        type(NACA_airfoil)                              :: airfoil
        real(kind=8),dimension(:),allocatable           :: alpha_angle 
        integer(kind=4)                                 :: PANELsize
        real(kind=8),dimension(:,:),allocatable         :: matrix
        type(panel),dimension(:),allocatable            :: PANEL_array
        real(kind=8),dimension(:),allocatable           :: vector
        real(kind=8),dimension(:,:),allocatable         :: cl_alpha 
        real(kind=8),dimension(:),allocatable           :: solution
        real(kind=8),dimension(:),allocatable           :: Vvec
        real(kind=8),dimension(:),allocatable           :: pressure
        real(kind=8),dimension(:),allocatable           :: cp_vec
        integer(kind=4)                                 :: i, j
        real(kind=8)                                    :: k
        
        ! computing alpha angle set        
        call compute_angle(alpha_angle,start_angle,end_angle,angle_num)

        ! generating airfoil object 
        call airfoil%set_AIRFOILname()
        call airfoil%set_npoints()
        airfoil%scaling = 1.0
        
        PANELsize = airfoil%n_points*2 - 2        
        allocate(PANEL_array(PANELsize))

        do i=1,angle_num
            
            print*, 'alpha angle = ', alpha_angle(i)

            ! making geometry --> it uses the airfoil object and then varies the alpha angle 
            call airfoil_CL(alpha_angle(i),airfoil,PANELsize,PANEL_array)
            
            if(i == 1)then
                ! allocating variables
                allocate(vector(PANELsize+1))
                allocate(solution(PANELsize+1))
                allocate(cp_vec(PANELsize))
                allocate(cl_alpha(angle_num,2))
            end if 

            ! computing matrix process
            call compute_matrix(matrix,PANEL_array,PANELsize)

            ! computing known vector 
            call compute_vector(real(1.0,8),real(0.0,8),vector,PANEL_array,PANELsize)
            
            ! solving system
            solution = solve_matrix(matrix,vector,PANELsize)
            
           ! ! computing velocity 
           ! Vvec     = compute_vel(solution,PANELsize,PANEL_array,real(1.0,8),real(0.0,8))
           ! 
           ! ! computing Cp
           ! cp_vec   = compute_cp(Vvec,real(1.0,8),PANELsize)

            ! computing Cl wrt alpha
            cl_alpha(i,1) = alpha_angle(i)*180/pi                    ! conversion to deg angle
           ! cl_alpha(i,2) = compute_cl(cp_vec,PANEL_array,PANELsize)
                
            cl_alpha(i,2) = 0.0

            do j=1,PANELsize 

                cl_alpha(i,2) = cl_alpha(i,2) - 2 * PANEL_array(j)%length * solution(PANELsize+1)

            end do   

            ! deallocation of matrix --> this because in the compute_matrix subroutine it allocates the matrix in the memory            
            deallocate(matrix)

        end do 

        call plot_cl(cl_alpha,angle_num)
        
        ! deallocation process            
        deallocate(vector)
        deallocate(solution)
        deallocate(cp_vec)
        deallocate(cl_alpha)

    end subroutine CLalpha

    subroutine airfoil_CL(alpha,airfoil,PANELsize,PANELarray)

        ! module declaration
        use AIRFOIL_object
        use PANEL_object
        use MEANline_object
        use discretization_module
        use FOUL
        use math_module  

        implicit none

        ! variable declaration
        type(NACA_airfoil)                                  :: airfoil            ! airfoil object
        type(panel),dimension(:),intent(inout)              :: PANELarray         ! panel object array
        type(MEANline),allocatable,dimension(:)             :: MEANLINEarray      ! mean line point object array
        real(kind=8),allocatable,dimension(:)               :: coordyUP           ! x-coords array -> describes the UPPER airfoil geometry
        real(kind=8),allocatable,dimension(:)               :: coordxUP           ! y-coords array -> describes the UPPER airfoil geometry
        real(kind=8),allocatable,dimension(:)               :: coordyDW           ! x-coords array -> describes the LOWER airfoil geometry
        real(kind=8),allocatable,dimension(:)               :: coordxDW           ! y-coords array -> describes the LOWER airfoil geometry
        integer(kind=4)                                     :: counter = 1        ! # of airfoil studied counter
        integer(kind=4)                                     :: x = 1              ! checking variable in while loop
        integer(kind=4)                                     :: dim                ! auxiliary variable -> # of discretisation points for the airfoil
        integer(kind=4)                                     :: k                  ! auxiliary variable
        integer(kind=4)                                     :: j                  ! auxiliary variable
        integer(kind=4)                                     :: i                  ! auxiliary variable
        real(kind=8)                                        :: airfoil_data1      ! easy-access variable -> 1st number of airfoil%data
        real(kind=8)                                        :: airfoil_data2      ! easy-access variable -> 2nd number of arifoil%data
        real(kind=8),intent(in)                             :: alpha              ! AOA 
        integer(kind=4),intent(inout)                       :: PANELsize          ! number of discretization panels
        real(kind=8),dimension(2)                           :: transl             ! translation vector        
        !character(len=30),intent(in)                        :: GNUplot_coord_data ! filename of the airfoil coords data container 
        !character(len=30),intent(in)                        :: GNUplot_mean_data  ! filename of the airfoil mean data container
        !character(len=30),intent(in)                        :: GNUplot_tg_norm    ! filename of the airfoil coords, tangent and normal data container
        
        ! setting properties 
        ! PAY ATTENTION the properties are already set in the setting_properties() subroutine
        airfoil%AOA     = alpha

        dim = airfoil%get_npoints()

        ! data allocation procedure
        allocate(coordxUP(1:dim)) 
        allocate(coordyUP(1:dim))
        allocate(coordxDW(1:dim)) 
        allocate(coordyDW(1:dim))
        allocate(MEANLINEarray(1:dim))


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MEAN LINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        ! calculate airfoil mean line point discretisation along x-axis
        do k=1,dim
            call MEANLINEarray(k)%set_id(k)
            call MEANLINEarray(k)%set_coordx(dim,k)
        end do
    
        ! extracting information from NACA digits 
        airfoil_data1 = real(airfoil%data(1),8)/100.0
        airfoil_data2 = real(airfoil%data(2),8)/10.0
    
        ! calculate airfoil mean line point discretisation along y-axis
        !     y coord
        !     gradient [dy] => theta
        k = 1 ! initialising data for the loop
    
        do while(MEANLINEarray(k)%coords(1)<airfoil_data2)
            call MEANLINEarray(k)%set_coordy_leading(airfoil_data1,airfoil_data2)
            call MEANLINEarray(k)%set_gradient_leading(airfoil_data1,airfoil_data2)
            k = k + 1
        end do
    
        do i=k,dim ! loop starts counting from k
            call MEANLINEarray(i)%set_coordy_trailing(airfoil_data1,airfoil_data2)
            call MEANLINEarray(i)%set_gradient_trailing(airfoil_data1,airfoil_data2)
        end do
    
        MEANLINEarray(dim)%coords = (/ 1.0, 0.0 /) ! mean-line end point coords
    
        ! thickness
        do k=1,dim
            call MEANLINEarray(k)%set_thickness(airfoil%data(3)) 
        end do
        ! coordx|coordy
        do k=1,dim
            call MEANLINEarray(k)%compute_UPcoords(coordxUP(k),coordyUP(k))
        end do
        do k=1,dim
            call MEANLINEarray(k)%compute_DOWNcoords(coordxDW(k),coordyDW(k))
        end do
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MEAN LINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! closure coords
        coordxUP(dim)                = 1
        coordyUP(dim)                = 0
        coordxDW(dim)                = 1
        coordyDW(dim)                = 0
        MEANLINEarray(dim)%coords(1) = 1
        MEANLINEarray(dim)%coords(2) = 0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PANEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!!!!!!!!!!!!!!!!!!!!!!!! PANEL DATA ALLOCATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! panel array allocation data -> LOWER AIRFOIL PART
                do k=1,dim-1
                    call PANELarray(k)%set_id(k)
                    call PANELarray(k)%set_coords(coordxDW(dim-k+1),coordyDW(dim-k+1),&
                                                  coordxDW(dim-k),coordyDW(dim-k))
                    call PANELarray(k)%set_position('DW')
                end do
                ! panel array allocation data -> UPPER AIRFOIL PART
                do k=1,dim-1
                    call PANELarray(dim+k-1)%set_id(dim+k-1)
                    call PANELarray(dim+k-1)%set_coords(coordxUP(k),coordyUP(k),&
                                                  coordxUP(k+1),coordyUP(k+1))
                    call PANELarray(dim+k-1)%set_position('UP')
                end do
            !!!!!!!!!!!!!!!!!!!!!!!!! PANEL DATA ALLOCATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!! PANEL PROPERTIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! compute segments length -> length doesn't vary with rotation and translation 
                do k=1,2*dim-2
                    call PANELarray(k)%compute_length()
                end do
                ! compute panel angle between x-axis and the tangent vector direction
                do k=1,2*dim-2
                    call PANELarray(k)%set_angle()
                end do
                ! compute segments tangent and normal versors
                do j=1,2*dim-2
                    call PANELarray(j)%compute_tangent_and_normal()  
                end do
                ! check on leading edge panels
                call check_LE_panels(PANELarray,dim)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!! PANEL PROPERTIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PANEL ROTATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if(alpha == 0.0)then 
                    ! compute inverse rotation matrix for each panel
                    do j=1,2*dim-2
                        call PANELarray(j)%compute_ROT()
                    end do 
                else  
                ! this procedure computes the rotation of the airfoil in space 
                    ! this procedure rotates the tangent vector of the airfoil
                    do j=1,2*dim-2                   
                        call rot(PANELarray(j)%tangent(1),PANELarray(j)%tangent(2),alpha)
                    end do
                    ! this procedure rotates the normal vector of the airfoil 
                    do j=1,2*dim-2                   
                        call rot(PANELarray(j)%normal(1),PANELarray(j)%normal(2),alpha)
                    end do
                    ! this procedure rotates the airfoil points of the airfoil
                    do j=1,2*dim-2
                        call rot(PANELarray(j)%coords1(1),PANELarray(j)%coords1(2),alpha)
                        call rot(PANELarray(j)%coords2(1),PANELarray(j)%coords2(2),alpha)
                    end do
                    ! this procedure rotates the airfoil meanline 
                    do j=1,dim
                        call rot(MEANLINEarray(j)%coords(1),MEANLINEarray(j)%coords(2),alpha)
                    end do
                    ! this procedure computes the variation of the panel angle
                    ! !!! PAY ATTENTION !!! the angle to vary is in the opposite direction of the convention
                    ! ---> an angle alpha > 0 will result on the system as the different of such angle
                    do j=1,2*dim-2
                        PANELarray(j)%angle = PANELarray(j)%angle - alpha
                    end do
                    ! this procedure computes the inverse rotation matrix used in the integral function
                    do j=1,2*dim-2
                        call PANELarray(j)%compute_ROT()
                    end do
                end if
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PANEL ROTATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!! PANEL TRANSLATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                airfoil%transl = transl
                
                ! translation process                        
                do j=1,2*dim-2
                    PANELarray(j)%coords1 = PANELarray(j)%coords1 + transl
                    PANELarray(j)%coords2 = PANELarray(j)%coords2 + transl
                end do

                do j=1,dim
                    MEANLINEarray(j)%coords = MEANLINEarray(j)%coords + transl
                end do
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!! PANEL TRANSLATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PANEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!! MIDPOINT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! compute segments middle points 
            !  -> middle-points depend on actual wing-points coordinates
            do j=1,2*dim-2
                call PANELarray(j)%compute_midpoint('noprt')
            end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!! MIDPOINT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !!!!!!!!!!!!!!!!!!!!!!! SAVING & GRAPHICS !!!!!!!!!!!!!!!!!!!!!!!!!
            ! storing data
            ! call ask_and_save(airfoil,PANELarray,MEANLINEarray)

            ! printing data
            ! call GNUplot_print(airfoil,PANELarray,MEANLINEarray,GNUplot_coord_data,GNUplot_mean_data,GNUplot_tg_norm) 
        !!!!!!!!!!!!!!!!!!!!!!! SAVING & GRAPHICS !!!!!!!!!!!!!!!!!!!!!!!!!
        
        ! data deallocation procedure
        deallocate(coordxUP) 
        deallocate(coordyUP)
        deallocate(coordxDW)
        deallocate(coordyDW)
        deallocate(MEANLINEarray)

        PANELsize = 2*dim-2

        call write_formatted('[','normal','OK','green','] -- airfoil geometry generated','normal')

    end subroutine airfoil_CL

    subroutine compute_field(PANEL_array,PANELsize,solution,V,P0,rho,alpha,filename)
    ! this subroutine computes the velocity field of the system after have computed the values of every singularity [gamma; sigma(i)]     
        use PANEL_object
        use math_module
        use FOUL
        implicit none 
        
        ! generating a grid of measure points 
        ! this grid can be expressed in a monodimensional array that contains 2 dimensional vectors 
        !   the procedure can be possible if it's used a type variable that describes a vector
        type :: array_type  
            real(kind=8),pointer :: coords(:)
        end type array_type 
        
        integer(kind=4),intent(in)                     :: PANELsize
        integer(kind=4)                                :: ncols 
        integer(kind=4)                                :: nrows 
        integer(kind=4)                                :: i, j, k
        real(kind=8)                                   :: deltax, deltay
        real(kind=8)                                   :: x_start, x_end
        real(kind=8)                                   :: y_start, y_end
        real(kind=8),dimension(2)                      :: velocity
        real(kind=8),intent(in)                        :: V
        real(kind=8),intent(in)                        :: alpha
        real(kind=8),dimension(PANELsize+1),intent(in) :: solution
        real(kind=8)                                   :: pressure
        real(kind=8)                                   :: norm_vel
        real(kind=8),intent(in)                        :: P0
        real(kind=8),intent(in)                        :: rho
        type(panel)                                    :: dummy_panel
        type(panel),dimension(PANELsize),intent(in)    :: PANEL_array
        type(array_type),dimension(:,:),allocatable    :: grid
        character(len=30),intent(in)                   :: filename
        integer(kind=4)                                :: x

        ! # of rows and columns 
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

        call write_formatted('COMPUTING FIELD','yellow')

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

        ! computing velocity 
        ! -- saving results in FLOWfield.dat
        ! -- using integral, computeCOEFF to compute velocity
        
        ! open file to save dat
        open(unit=1, file= filename, status='replace')
        
        ! compute velocity process
        do i=1,nrows
            do j=1,ncols
                ! allocating midpoint in dummy_panel
                ! it's the only point necessary to compute the velocity
                dummy_panel%midpoint(1) = grid(i,j)%coords(1)
                dummy_panel%midpoint(2) = grid(i,j)%coords(2)
                
                ! initializing velocity
                velocity(1) = 0.0
                velocity(2) = 0.0
                
                if(grid(i,j)%coords(1) >= PANEL_array(PANELsize/2+1)%coords1(1) .and. & 
                   grid(i,j)%coords(1) <= PANEL_array(1)%coords1(1))then 
                    
                    x = 0
                    k = 1

                    do while(x == 0 .and. k<=PANELsize/2)
                             
                        if(grid(i,j)%coords(2) >= PANEL_array(PANELsize/2-k+1)%midpoint(2) .and. &
                           grid(i,j)%coords(2) <= PANEL_array(PANELsize/2+k)%midpoint(2))then
                            
                           if(grid(i,j)%coords(1) <= PANEL_array(PANELsize/2+k)%coords2(1) .and. & 
                              grid(i,j)%coords(1) >= PANEL_array(PANELsize/2+k-1)%coords1(1))then
                                                                  
                                x = 1                             
                                                                  
                           end if                                 
                                                                  
                        end if
                        
                        k = k + 1

                    end do
                end if

                ! computing velocity through integral funcition
                do k=1,PANELsize
                    velocity = velocity + integral(dummy_panel,PANEL_array(k),'source') * solution(k)                    
                    velocity = velocity + integral(dummy_panel,PANEL_array(k),'vortex') * solution(PANELsize+1)
                end do

                ! imposing external velocity
                velocity(1) = velocity(1) + V*cos(alpha)
                velocity(2) = velocity(2) + V*sin(alpha)
                    
                ! computing norm of velocity
                norm_vel = norm(velocity)
                
                if(grid(i,j)%coords(1) >= PANEL_array(PANELsize/2)%coords2(1) .and. & 
                   grid(i,j)%coords(1) <= PANEL_array(1)%coords1(1))then                 
                    if(x == 1)then 
                        velocity = (/0.0, 0.0/)
                        norm_vel = 0
                    end if
                end if 
                ! computing pressure
                pressure = P0 + 0.5*rho*(V**2 - norm_vel**2)  
                
                if(pressure < -1e-1)then
                    pressure = P0 + 0.5*rho*V**2
                end if

                ! writing data in FLOWfield.dat
                ! x_coord, y_coords, x_vel, y_vel, velocity_norm, pressure    
                write(1,*) dummy_panel%midpoint(1), dummy_panel%midpoint(2), &
                           velocity(1)            , velocity(2),             & 
                           norm_vel               , pressure

            end do
        end do

        close(1)
        ! deallocation process
        do i=1,nrows
            do j=1,ncols
                deallocate(grid(i,j)%coords)
            end do 
        end do 
        deallocate(grid)
    end subroutine compute_field

end module cp   
