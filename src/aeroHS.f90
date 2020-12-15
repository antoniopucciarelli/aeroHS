!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! author           : antonio pucciarelli                                  !
! date of creation : 11/20/2020                                           !
! written with     : vim                                                  ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program aeroHS
    
    use airfoilgenerator
    use MEANline_object
    use AIRFOIL_object
    use PANEL_object
    use print_save
    use ask_module
    use ground_cp
    use discretization_module  
    use cp
    use plot

    implicit none
    
    type(NACA_airfoil)                      :: airfoil
    type(MEANline),dimension(:),allocatable :: MEAN_array    
    type(panel),dimension(:),allocatable    :: PANEL_array
    type(panel),dimension(:),allocatable    :: GROUNDpanel
    integer(kind=4)                         :: MEANsize
    real(kind=8),dimension(:,:),allocatable :: matrix
    real(kind=8),dimension(:),allocatable   :: vector     
    real(kind=8),dimension(:),allocatable   :: solution
    real(kind=8),dimension(:),allocatable   :: Vvec
    real(kind=8),dimension(:),allocatable   :: pressure
    real(kind=8),dimension(:),allocatable   :: cp_vec 
    real(kind=8),dimension(:,:),allocatable :: cl_alpha
    integer(kind=4)                         :: PANELsize      = 0
    integer(kind=4)                         :: maxsize        = 0
    real(kind=8)                            :: alpha          = 0.0
    real(kind=8)                            :: V              = 0.0
    real(kind=8)                            :: P0             = 0.0
    real(kind=8)                            :: rho            = 0.0
    integer(kind=4)                         :: start_angle    = 0
    integer(kind=4)                         :: end_angle      = 0
    integer(kind=4)                         :: dim            = 0
    integer(kind=4)                         :: selection      = 0
    integer(kind=4)                         :: selection_type = 0
    integer(kind=4)                         :: i              = 1
    integer(kind=4)                         :: GROUNDsize
    real(kind=8)                            :: CL
    character(len=30)                       :: filename
    character(len=6)                        :: panel_type     

    do while(i==1)
        
        print*, 'there are 3 different options:'
        print*, '   option 1 --> analize 1 airfoil'
        print*, '                   1 --> compute  Cp, Cl'
        print*, '                   2 --> compute [Cp; Cl; velocity; pressure]'
        print*, '   option 2 --> analize 2 airfoils'
        print*, '                   1 --> [L; M]'
        print*, '   option 3 --> analize ground effect'
        print*, '                   1 --> compute  Cp, Cl'
        print*, 'type an option'
        read*,  selection_type 

        if(selection_type == 1)then
            
            call setting_properties(P0,V,rho,alpha,start_angle,end_angle,dim,selection,selection_type)

            if(selection == 1)then
            
                call ask_geometry(PANELsize,PANEL_array,MEAN_array,airfoil,alpha,selection)
                
                ! autosaving geometry for plots
                call GNUplot_saving(PANEL_array,MEAN_array,airfoil%get_npoints())
                
                ! computing matrix process
                call compute_matrix(matrix,PANEL_array,PANELsize)
                
                ! compute known vector properties from geometry and external flow 
                allocate(vector(1:PANELsize))
                ! because the alpha angle is already expressed in the airfoil geometry
                !    there is no need to put alpha into the next function 
                ! alpha = real(0.0,8)
                call compute_vector(V,real(0.0,8),vector,PANEL_array,PANELsize)

                ! allocation of solution array's dimensions 
                allocate(solution(1:PANELsize+1))
                ! computing solution
                solution = solve_matrix(matrix,vector,PANELsize)

                allocate(Vvec(1:PANELsize))
                ! alpha = real(0.0,8)
                Vvec = compute_vel(solution,PANELsize,PANEL_array,V,real(0.0,8))

                allocate(cp_vec(1:PANELsize))
                cp_vec = compute_cp(Vvec,V,PANELsize)
                
                ! compute CL value 
                CL = compute_cl(cp_vec,PANEL_array,PANELsize)
                print*, 'CL value = ', CL

                ! asking to save matrices, known vector and solution
                ! call ask_to_save_matrix_vector(PANELsize,matrix,vector,solution)

                !!!!!!!!!!!!!!!!!!! COMPUTING VELOCITY FIELD !!!!!!!!!!!!!!!!!!!
                ! high demanding process 
                ! -- it depends on the dimension of the system and its discrtization
                ! alpha = real(0.0,8)
                call compute_field(PANEL_array,PANELsize,solution,V,P0,rho,real(0.0,8))
                !!!!!!!!!!!!!!!!!!! COMPUTING VELOCITY FIELD !!!!!!!!!!!!!!!!!!!

                !!!!!!!!!!!!!!!!!!!!!! PLOTTING RESULTS !!!!!!!!!!!!!!!!!!!!!!!!
                call plot_cp(cp_vec,PANEL_array,PANELsize)
                call plot_vel_field()
                !!!!!!!!!!!!!!!!!!!!!! PLOTTING RESULTS !!!!!!!!!!!!!!!!!!!!!!!!

                ! deallocation process
                deallocate(Vvec)
                deallocate(cp_vec)
                deallocate(vector)
                deallocate(solution)

            else if(selection == 2)then

                call ask_geometry(PANELsize,PANEL_array,MEAN_array,airfoil,alpha,selection)

                ! autosaving geometry for plots
                call GNUplot_saving(PANEL_array,MEAN_array,airfoil%get_npoints())
                
                ! computing matrix process
                call compute_matrix(matrix,PANEL_array,PANELsize)
            
                ! compute known vector properties from geometry and external flow 
                allocate(vector(1:PANELsize))
                ! alpha = real(0.0,8)
                call compute_vector(V,real(0.0,8),vector,PANEL_array,PANELsize)

                ! allocation of solution array's dimensions 
                allocate(solution(1:PANELsize+1))
                ! computing solution
                solution = solve_matrix(matrix,vector,PANELsize)
            
                allocate(Vvec(1:PANELsize))
                ! because the alpga angle is already expressed in the airfoil geometry
                !    there is no need to put alpha into the next function 
                ! alpha = real(0.0,8)
                Vvec = compute_vel(solution,PANELsize,PANEL_array,V,real(0.0,8))

                allocate(pressure(1:PANELsize))
                pressure = compute_pressure(P0,V,rho,Vvec,PANELsize)

                allocate(cp_vec(1:PANELsize))
                cp_vec = compute_cp(Vvec,V,PANELsize)

                ! compute CL value 
                CL = compute_cl(cp_vec,PANEL_array,PANELsize)
                print*, 'CL value = ', CL

                ! asking to save matrices, known vector and solution
                ! call ask_to_save_matrix_vector(PANELsize,matrix,vector,solution)
                
                !!!!!!!!!!!!!!!!!!! COMPUTING VELOCITY FIELD !!!!!!!!!!!!!!!!!!!
                ! high demanding process 
                ! -- it depends on the dimension of the system and its discrtization
                call compute_field(PANEL_array,PANELsize,solution,V,P0,rho,real(0.0,8))
                !!!!!!!!!!!!!!!!!! COMPUTING VELOCITY FIELD !!!!!!!!!!!!!!!!!!!!

                !!!!!!!!!!!!!!!!!!!!!! PLOTTING RESULTS !!!!!!!!!!!!!!!!!!!!!!!!
                call plot_cp(cp_vec,PANEL_array,PANELsize)
                call plot_pressure(pressure,PANEL_array,PANELsize)
                call plot_vel(Vvec,PANEL_array,PANELsize)
                !!!!!!!!!!!!!!!!!!!!!! PLOTTING RESULTS !!!!!!!!!!!!!!!!!!!!!!!!

                ! deallocation process
                deallocate(Vvec)
                deallocate(pressure)
                deallocate(cp_vec)
                deallocate(vector)
                deallocate(solution)

            else if(selection == 3)then

                call ask_geometry(PANELsize,PANEL_array,MEAN_array,airfoil,real(0.0,8),selection)

                ! computing matrix process
                call compute_matrix(matrix,PANEL_array,PANELsize)
            
                allocate(cl_alpha(dim,2))
                ! modify alpha angle = real(0.0,8) because the airfoil geometry inclination
                call CLalpha(cl_alpha,start_angle,end_angle,dim,matrix,PANEL_array,PANELsize)

                ! deallocation process
                deallocate(cl_alpha)
            
            end if

            call ask_to_continue_cp(i)

            ! deallocation process
            deallocate(PANEL_array)
            deallocate(MEAN_array)
            deallocate(matrix)
            
        else if(selection_type == 2)then 
                    
            ! call 2 times airfoil generator
            ! modify extension at every suborutine and function regarding # of panels
            ! use multi_cp to made this changements -- you have cp that works fine
            ! enable difference dimension between airfoils
            ! modify matrix generation
            ! compute moment and lift
            
        else if(selection_type == 3)then 

            ! asking method to adopt for the computation
            call ask_method(panel_type)            

            ! generate flow properties
            call setting_properties(P0,V,rho,alpha,start_angle,end_angle,dim,selection,selection_type)
            
            ! generate airfoil geometry 
            call ask_geometry(PANELsize,PANEL_array,MEAN_array,airfoil,alpha,selection)

            ! autosaving geometry for plots
            call GNUplot_saving(PANEL_array,MEAN_array,airfoil%get_npoints())
            
            ! ground panels generation             
            call generate_ground(GROUNDpanel,GROUNDsize)
                 
            ! compute system matrix 
            call computeGROUNDmatrix(matrix,PANEL_array,GROUNDpanel,PANELsize,GROUNDsize,panel_type)

            ! compute known vector 
            call computeGROUNDvector(vector,PANEL_array,GROUNDpanel,PANELsize,GROUNDsize,real(0.0,8),V)

            ! computing solution
            call solveGROUND(solution,matrix,vector,PANELsize,GROUNDsize)
            
            ! asking to save matrices, known vector and solution
            call ask_to_save_matrix_vector(PANELsize+GROUNDsize,matrix,vector,solution)
            
            !!!!!!!!!!!!!!!!!! COMPUTING VELOCITY FIELD !!!!!!!!!!!!!!!!!!!!
            ! high demanding process 
            ! -- it depends on the number of point for the field discretization 
            call computeGROUNDfield(solution,PANEL_array,GROUNDpanel,PANELsize,GROUNDsize,P0,real(0.0,8),V,rho,panel_type)
            !!!!!!!!!!!!!!!!!! COMPUTING VELOCITY FIELD !!!!!!!!!!!!!!!!!!!!

            ! computing cp, pressure and velocity around the airfoil 
            allocate(cp_vec(PANELsize))
            call compute_airfoilFIELD(solution,PANEL_array,GROUNDpanel,PANELsize,GROUNDsize,cp_vec,P0,real(0.0,8),V,rho,panel_type)
                        
            ! compute CL value 
            CL = compute_cl(cp_vec,PANEL_array,PANELsize)
            print*, 'CL value = ', CL    

            ! asking user to continue
            call ask_to_continue_cp(i)

            ! deallocation process
            deallocate(PANEL_array)
            deallocate(GROUNDpanel)
            deallocate(MEAN_array)
            deallocate(solution)
            deallocate(matrix)
            deallocate(vector)
            deallocate(cp_vec)        

        end if 
    end do 

end program aeroHS
