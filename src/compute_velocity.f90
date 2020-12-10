!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program compute_velocity
    
    use airfoilgenerator
    use MEANline_object
    use AIRFOIL_object
    use PANEL_object
    use print_save
    use plot
    use cp

    implicit none
    
    type(NACA_airfoil)                      :: airfoil
    type(MEANline),dimension(:),allocatable :: MEAN_array    
    type(panel),dimension(:),allocatable    :: PANEL_array
    integer(kind=4)                         :: MEANsize
    real(kind=8),dimension(:,:),allocatable :: matrix
    real(kind=8),dimension(:),allocatable   :: vector     
    real(kind=8),dimension(:),allocatable   :: solution
    real(kind=8),dimension(:),allocatable   :: Vvec
    real(kind=8),dimension(:),allocatable   :: pressure
    real(kind=8),dimension(:),allocatable   :: cp_vec 
    real(kind=8),dimension(:,:),allocatable :: cl_alpha
    integer(kind=4)                         :: PANELsize   = 0
    integer(kind=4)                         :: maxsize     = 0
    real(kind=8)                            :: alpha       = 0.0
    real(kind=8)                            :: V           = 0.0
    real(kind=8)                            :: P0          = 0.0
    real(kind=8)                            :: rho         = 0.0
    integer(kind=4)                         :: start_angle = 0
    integer(kind=4)                         :: end_angle   = 0
    integer(kind=4)                         :: dim         = 0
    integer(kind=4)                         :: selection   = 0
    character(len=30)                       :: filename
    integer(kind=4)                         :: i           = 1
    
    do while(i==1)

        ! setting up velocity, ambient pressure and air density
        call setting_properties(P0,V,rho,alpha,start_angle,end_angle,dim,selection)

        call ask_geometry(PANELsize,PANEL_array,MEAN_array,airfoil)

        ! computing matrix process
        call compute_matrix(matrix,PANEL_array,PANELsize)

        if(selection == 1)then
        
            ! compute known vector properties from geometry and external flow 
            allocate(vector(1:PANELsize))
            call compute_vector(real(1.0,8),alpha,vector,PANEL_array,PANELsize)

            ! allocation of solution array's dimensions 
            allocate(solution(1:PANELsize+1))
            ! computing solution
            solution = solve_matrix(matrix,vector,PANELsize)

            allocate(Vvec(1:PANELsize))
            Vvec = compute_vel(solution,PANELsize,PANEL_array,real(1.0,8),alpha)

            allocate(cp_vec(1:PANELsize))
            cp_vec = compute_cp(Vvec,real(1.0,8),PANELsize)

            call ask_to_save_matrix_vector(PANELsize,matrix,vector,solution)

            !!!!!!!!!!!!!!!!!!!!!! PLOTTING RESULTS !!!!!!!!!!!!!!!!!!!!!!!!
            call plot_cp(cp_vec,PANEL_array,PANELsize)
            !!!!!!!!!!!!!!!!!!!!!! PLOTTING RESULTS !!!!!!!!!!!!!!!!!!!!!!!!

            ! deallocation process
            deallocate(Vvec)
            deallocate(cp_vec)
            deallocate(vector)
            deallocate(solution)

        else if(selection == 2)then

            ! compute known vector properties from geometry and external flow 
            allocate(vector(1:PANELsize))
            call compute_vector(V,alpha,vector,PANEL_array,PANELsize)

            ! allocation of solution array's dimensions 
            allocate(solution(1:PANELsize+1))
            ! computing solution
            solution = solve_matrix(matrix,vector,PANELsize)
        
            allocate(Vvec(1:PANELsize))
            Vvec = compute_vel(solution,PANELsize,PANEL_array,V,alpha)

            allocate(pressure(1:PANELsize))
            pressure = compute_pressure(P0,V,rho,Vvec,PANELsize)

            allocate(cp_vec(1:PANELsize))
            cp_vec = compute_cp(Vvec,V,PANELsize)

            call ask_to_save_matrix_vector(PANELsize,matrix,vector,solution)

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

            allocate(cl_alpha(dim,2))
            call CLalpha(cl_alpha,start_angle,end_angle,dim,matrix,PANEL_array,PANELsize)

            ! dealloctaion process
            deallocate(cl_alpha)
        
        end if

        call ask_to_continue_cp(i)

        ! deallocation process
        deallocate(PANEL_array)
        deallocate(MEAN_array)
        deallocate(matrix)

    end do 

end program compute_velocity