module cp 

    contains 

    subroutine setting_properties(P0,V,rho,alpha,start_angle,end_angle,dim,selection)
    ! this subroutine sets the flow's external conditions (at infinity) 
        use math_module
        implicit none
         
        real(kind=8),intent(inout)    :: P0
        real(kind=8),intent(inout)    :: V
        real(kind=8),intent(inout)    :: rho
        real(kind=8),intent(inout)    :: alpha       ! [deg]
        integer(kind=4),intent(inout) :: start_angle ! [deg]
        integer(kind=4),intent(inout) :: end_angle   ! [deg]
        integer(kind=4),intent(inout) :: dim 
        integer(kind=4),intent(inout) :: selection
        integer(kind=4)               :: x 

        x = 1

        do while(x==1)
            print*, 'to compute Cp       vs X --> type(1)'
            print*, 'to compute pressure vs X --' 
            print*, '           velocity vs X --> type(2)'
            print*, '           Cp       vs X --'
            print*, 'to compute Cl       vs X --> type(3)'
            read*, selection

            select case (selection)
                case(1)
                    print*, 'type airfoil angle of attack -- AOA [deg]'
                    read*, alpha
                    V = 1
                    x = 0

                case(2)
                    print*, 'type ambient pressure P0 [Pa] '
                    read*, P0
                    print*, 'type air velocity at infinity [m/s]'
                    read*, V
                    print*, 'type air density [kg/m**3]   -- method hp: constant along the airfoil'
                    read*, rho
                    print*, 'type airfoil angle of attack -- AOA [deg]'
                    read*, alpha

                    x = 0
                
                case(3)
                    call ask_angle(start_angle,end_angle,dim)
                    V = 1
                    x = 0

                case default 
                    print*, 'you have selected an invalid action', new_line('(A)'), 'type again'

            end select 

        end do

        alpha = alpha/180.0*pi

    end subroutine setting_properties

    subroutine ask_angle(start_angle,end_angle,dim)
        implicit none 

        integer(kind=4),intent(inout) :: start_angle
        integer(kind=4),intent(inout) :: end_angle
        integer(kind=4),intent(inout) :: dim

        print*, 'type the initial AOA of the sequence'
        read*, start_angle
        print*, 'type the end AOA of the sequence'
        read*, end_angle
        print*, 'type # of discretization points'
        read*, dim

    end subroutine ask_angle

    subroutine ask_to_continue_cp(i)
        implicit none
        character(len=1)              :: resp
        integer(kind=4),intent(inout) :: i 
        
        print*, 'do you want to create a new study? [Y\n]'
        read*, resp
        
        if(resp=='Y' .or. resp=='y')then
            i = 1
        else 
            i = 0
        end if
    end subroutine ask_to_continue_cp

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
        call GNUplot_print(airfoil,PANEL_array,MEAN_array)
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
        real(kind=8)                                        :: vortex_value = 0.0
        real(kind=8),dimension(2)                           :: tangent_ith, normal_ith
        real(kind=8),dimension(2)                           :: tangent_first, normal_first
        real(kind=8),dimension(2)                           :: tangent_last, normal_last

        allocate(matrix(PANELsize+1,PANELsize+1))

        do i=1,PANELsize
            
            ! declaring normal vector
            normal_ith  = PANEL_array(i)%normal

            do j=1,PANELsize
                matrix(i,j)  = dot_product(normal_ith,integral(PANEL_array(i),PANEL_array(j),'source'))
                vortex_value = vortex_value + dot_product(normal_ith,integral(PANEL_array(i),PANEL_array(j),'vortex'))
            end do

            matrix(i,PANELsize+1) = vortex_value
            vortex_value          = 0.0

        end do 
        
        ! declaring tangent and normal vectors
        ! first panel -- AIRFOIL UPPER PART 
        normal_first  = PANEL_array(1)%normal
        tangent_first = PANEL_array(1)%tangent
        ! last panel  -- AIRFOIL LOWER PART 
        normal_last   = PANEL_array(PANELsize)%normal
        tangent_last  = PANEL_array(PANELsize)%tangent 
        
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
            normal_ith =   PANEL_array(i)%normal

            vector(i)  = - dot_product(V_vec,normal_ith)
        
        end do

        ! declaring tangent vectors
        tangent_first = PANEL_array(1)%tangent
        tangent_last  = PANEL_array(PANELsize)%tangent
       
        vector(PANELsize+1) = - (dot_product(V_vec,tangent_first) + &
                                 dot_product(V_vec,tangent_last))

        call write_formatted('[','normal','OK','green','] -- velocity vector created','normal')

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

        test = matmul(matrix,solution) - vector

        maxerr = maxval(abs(test), PANELsize)

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

        compute_vel = matmul(matrix, solution) + vector 
        
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
        S          = 0.0
        
        do i=1,PANELsize 

            panel_length = PANEL_array(i)%get_length()

            S = S + panel_length

            if(PANEL_array(i)%get_position() == 'UP')then
                panel_length = - panel_length
            end if

            compute_cl = compute_cl + panel_length * cp_vec(i)

        end do 

        S = abs(norm(PANEL_array(PANELsize/2+1)%coords1 - PANEL_array(1)%coords1))

        compute_cl = compute_cl / S

    end function compute_cl

    subroutine CLalpha(cl_alpha,start_angle,end_angle,dim,matrix,PANEL_array,PANELsize)
    ! this suborutine computes the cl alpha diagram for the airfoil 
        use FOUL 
        use PANEL_object
        use plot
        use math_module
        implicit none 

        integer(kind=4)                                 :: start_angle
        integer(kind=4)                                 :: end_angle
        integer(kind=4)                                 :: dim                             
        real(kind=8),dimension(dim)                     :: alpha_angle 
        integer(kind=4),intent(in)                      :: PANELsize
        real(kind=8),dimension(PANELsize+1,PANELsize+1) :: matrix
        type(panel),dimension(PANELsize)                :: PANEL_array
        real(kind=8),dimension(PANELsize+1)             :: vector
        real(kind=8),dimension(dim,2),intent(inout)     :: cl_alpha 
        real(kind=8),dimension(PANELsize)               :: solution
        real(kind=8),dimension(PANELsize)               :: Vvec
        real(kind=8),dimension(PANELsize)               :: pressure
        real(kind=8),dimension(PANELsize)               :: cp_vec
        integer(kind=4)                                 :: i
        real(kind=8)                                    :: k

        k = real((end_angle - start_angle),8) /real(dim,8)

        do i=1,dim
            
            alpha_angle(i) = (real(start_angle,8) + k*real(i,8))/180*pi  

        end do 

        alpha_angle(1)   = real(start_angle,8)/180*pi
        alpha_angle(dim) = real(end_angle,8)/180*pi

        do i=1,dim

            call compute_vector(real(1.0,8),alpha_angle(i),vector,PANEL_array,PANELsize)
            
            solution      = solve_matrix(matrix,vector,PANELsize)

            Vvec          = compute_vel(solution,PANELsize,PANEL_array,real(1.0,8),alpha_angle(i))
            
            cp_vec        = compute_cp(Vvec,real(1.0,8),PANELsize)

            cl_alpha(i,1) = alpha_angle(i)*180/pi
            cl_alpha(i,2) = compute_cl(cp_vec,PANEL_array,PANELsize)

        end do 

        call plot_cl(cl_alpha,dim)

    end subroutine CLalpha
        
    subroutine compute_vel_field(PANEL_array,PANELsize,solution,V,alpha)
    ! this subroutine computes the velocity field of the system after have computed the values of every singularity [gamma; sigma(i)]     
        use PANEL_object
        use math_module
        implicit none 
        
        ! generating a grid of measure points 
        ! this grid can be expressed in a monodimensional array that contains 2 dimensional vectors 
        !   the procedure can be possible if it's used a type variable that describes a vector
        type :: array_type  
            real(kind=8),pointer :: coords(:)
        end type array_type 
        
        integer(kind=4),intent(in)                     :: PANELsize
        integer(kind=4)                                :: ncols = 5e+2
        integer(kind=4)                                :: nrows = 2e+2
        integer(kind=4)                                :: i, j, k
        real(kind=8)                                   :: deltax, deltay
        real(kind=8)                                   :: x_start, x_end
        real(kind=8)                                   :: y_start, y_end
        real(kind=8),dimension(2)                      :: velocity
        real(kind=8),intent(in)                        :: V
        real(kind=8),intent(in)                        :: alpha
        real(kind=8),dimension(PANELsize+1),intent(in) :: solution
        type(panel)                                    :: dummy_panel
        type(panel),dimension(PANELsize),intent(in)    :: PANEL_array
        type(array_type),dimension(:,:),allocatable    :: grid

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

        ! computing velocity 
        ! -- saving results in VELfield.dat
        ! -- using integral, computeCOEFF to compute velocity
        
        ! open file to save dat
        open(unit=1, file='VELfield.dat', status='replace')
        write(1,*) 'x_coord, y_coords, x_vel, y_vel, velocity_norm'    
        
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
                
                ! computing velocity through integral funcition
                do k=1,PANELsize
                    velocity = velocity + integral(dummy_panel,PANEL_array(k),'source') * solution(k)                    
                    velocity = velocity + integral(dummy_panel,PANEL_array(k),'vortex') * solution(PANELsize+1)
                end do

                ! imposing external velocity
                velocity(1) = velocity(1) + V*cos(alpha)
                velocity(2) = velocity(2) + V*sin(alpha)
                
                ! writing data in VELfield.dat
                write(1,*) dummy_panel%midpoint(1), dummy_panel%midpoint(2), &
                           velocity(1)            , velocity(2), norm(velocity)

            end do
        end do

        close(1)
   
        deallocate(grid)
    end subroutine compute_vel_field

end module cp   
