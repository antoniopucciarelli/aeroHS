!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! author           : antonio pucciarelli                                  !
! date of creation : 11/20/2020                                           !
! written with     : vim                                                  ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program aeroHS
    
    use discretization_module  
    use airfoilgenerator
    use MEANline_object
    use AIRFOIL_object
    use PANEL_object
    use ask_module
    use print_save
    use ground_cp
    use multi_cp
    use plot
    use cp

    implicit none
    
    type(NACA_airfoil)                      :: airfoil
    type(NACA_airfoil)                      :: airfoil1
    type(NACA_airfoil)                      :: airfoil2
    type(MEANline),dimension(:),allocatable :: MEAN_array   
    type(MEANline),dimension(:),allocatable :: MEAN_array1
    type(MEANline),dimension(:),allocatable :: MEAN_array2 
    type(panel),dimension(:),allocatable    :: PANEL_array
    type(panel),dimension(:),allocatable    :: PANEL_array1
    type(panel),dimension(:),allocatable    :: PANEL_array2 
    type(panel),dimension(:),allocatable    :: GROUNDpanel
    real(kind=8),dimension(:,:),allocatable :: matrix
    real(kind=8),dimension(:),allocatable   :: vector     
    real(kind=8),dimension(:),allocatable   :: solution
    real(kind=8),dimension(:),allocatable   :: Vvec
    real(kind=8),dimension(:),allocatable   :: pressure
    real(kind=8),dimension(:),allocatable   :: cp_vec
    real(kind=8),dimension(:),allocatable   :: cp_vec1
    real(kind=8),dimension(:),allocatable   :: cp_vec2 
    real(kind=8),dimension(:,:),allocatable :: cl_alpha
    real(kind=8)                            :: alpha          = 0.0
    real(kind=8)                            :: alpha1         = 0.0
    real(kind=8)                            :: alpha2         = 0.0
    real(kind=8)                            :: V              = 0.0
    real(kind=8)                            :: P0             = 0.0
    real(kind=8)                            :: rho            = 0.0
    real(kind=8)                            :: CL             = 0.0
    integer(kind=4)                         :: MEANsize       = 0
    integer(kind=4)                         :: PANELsize      = 0
    integer(kind=4)                         :: PANELsize1     = 0
    integer(kind=4)                         :: PANELsize2     = 0
    integer(kind=4)                         :: maxsize        = 0
    integer(kind=4)                         :: start_angle    = 0
    integer(kind=4)                         :: end_angle      = 0
    integer(kind=4)                         :: dim            = 0
    integer(kind=4)                         :: selection      = 0
    integer(kind=4)                         :: selection_type = 0
    integer(kind=4)                         :: i              = 1
    integer(kind=4)                         :: GROUNDsize
    character(len=30)                       :: filename
    character(len=30)                       :: filename1
    character(len=30)                       :: filename2
    character(len=6)                        :: panel_type     
        
    print*, 'there are 3 different options:'
    print*, '   option 1 --> analize 1 airfoil'
    print*, '                   1 --> compute  Cp, Cl'
    print*, '                   2 --> compute  Cl/alpha'
    print*, '   option 2 --> analize 2 airfoils'
    print*, '                   1 --> compute  Cp, Cl'
    print*, '   option 3 --> analize ground effect [with ground panels] '
    print*, '                   1 --> compute  Cp, Cl'
    print*, 'type an option'
    read*,  selection_type 

    if(selection_type == 1)then
        
        call setting_properties(P0,V,rho,alpha1,alpha2,start_angle,end_angle,dim,selection,selection_type)

        if(selection == 1)then
            
            ! describing file name
            filename  = 'GNUplot_coord_data.dat'
            filename1 = 'GNUplot_mean_data.dat'
            filename2 = 'GNUplot_tg_norm.dat'
            call ask_geometry(PANELsize,PANEL_array,MEAN_array,airfoil,alpha1,selection, & 
                              filename,filename1,filename2)
            
            ! autosaving geometry for plots
            call GNUplot_saving(PANEL_array,MEAN_array,airfoil%get_npoints(),filename,filename1,filename2)
            
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
            
            ! allocation of velocity vector
            allocate(Vvec(1:PANELsize))
            ! alpha = real(0.0,8)
            Vvec = compute_vel(solution,PANELsize,PANEL_array,V,real(0.0,8))
            
            ! allocation of cp
            allocate(cp_vec(1:PANELsize))
            cp_vec = compute_cp(Vvec,V,PANELsize)
            
            ! compute CL value 
            CL = compute_cl(solution(PANELsize+1),PANEL_array,PANELsize)
            print*, 'CL value                  = ', CL

            !!!!!!!!!!!!!!!!!!! COMPUTING VELOCITY FIELD !!!!!!!!!!!!!!!!!!!
            ! high demanding process 
            ! -- it depends on the dimension of the system and its discrtization
            ! alpha = real(0.0,8)
            filename = 'FLOWfield.dat'
            call compute_field(PANEL_array,PANELsize,solution,V,P0,rho,real(0.0,8),filename)
            !!!!!!!!!!!!!!!!!!! COMPUTING VELOCITY FIELD !!!!!!!!!!!!!!!!!!!

            !!!!!!!!!!!!!!!!!!!!!! PLOTTING RESULTS !!!!!!!!!!!!!!!!!!!!!!!!
            call plot_cp(cp_vec,PANEL_array,PANELsize,alpha1)
            call plot_vel_field()
            !!!!!!!!!!!!!!!!!!!!!! PLOTTING RESULTS !!!!!!!!!!!!!!!!!!!!!!!!

            ! deallocation process
            deallocate(Vvec)
            deallocate(cp_vec)
            deallocate(vector)
            deallocate(matrix)
            deallocate(solution)
            deallocate(MEAN_array)
            deallocate(PANEL_array)

        else if(selection == 2)then

            call CLalpha(start_angle,end_angle,dim)
        
        end if
        
    else if(selection_type == 2)then 
        
        ! known variables 
        V     = 1.0 
        alpha = 0.0
        rho   = 1.0
        P0    = 1.0
        
        ! generating airfoils 
        call generate_airfoils(selection,selection_type,alpha1,alpha2,PANELsize1,PANELsize2,PANEL_array1,PANEL_array2, &
                               MEAN_array1,MEAN_array2)
        
        ! plotting generated airfoils             
        call plot_airfoils()

        ! computing system matrix 
        call compute_multi_matrix(PANELsize1,PANELsize2,PANEL_array1,PANEL_array2,matrix)
         
        ! computing known vector  
        call compute_multi_vector(vector,PANELsize1,PANELsize2,PANEL_array1,PANEL_array2,alpha,V)
        
        ! saving matrix, vectors
        filename = 'matrix_data.dat' 
        call save_matrix(matrix,PANELsize1+PANELsize2+2,filename) 
        filename = 'vector_data.dat'
        call save_vector(vector,PANELsize1+PANELsize2+2,filename)

        ! computing and testing solution
        call solvemulti(PANELsize1,PANELsize2,matrix,vector,solution)
        
        ! computing Cp vs X
        allocate(cp_vec1(PANELsize1))
        allocate(cp_vec2(PANELsize2))
        call compute_MULTIairfoilFIELD(solution,PANEL_array1,PANEL_array2,PANELsize1,PANELsize2,cp_vec1,cp_vec2,P0,alpha,V,rho)
        
        ! computing Cl
        ! 1st airfoil 
        CL = compute_cl(solution(PANELsize1+1),PANEL_array1,PANELsize1)
        print*, '1st airfoil Cl                = ', CL
        ! 2nd airfoil
        CL = compute_cl(solution(PANELsize1+PANELsize2+2),PANEL_array2,PANELsize2)
        print*, '2nd airfoil Cl                = ', CL
        
        ! computing velocity at y = 0
        call  compute_midflow(PANELsize1,PANELsize2,PANEL_array1,PANEL_array2,solution,V,real(0.0,8))
        
        ! computing velocity and pressure fields
        ! plotting pressure field 
        filename = 'FLOWfieldMULTI.dat'
        call compute_multi_field(solution,PANEL_array1,PANEL_array2,PANELsize1,PANELsize2,P0,alpha,V,rho,filename) 

        ! deallocation process
        deallocate(matrix)
        deallocate(vector)
        deallocate(cp_vec1)
        deallocate(cp_vec2)
        deallocate(solution) 
        deallocate(MEAN_array1)
        deallocate(MEAN_array2)
        deallocate(PANEL_array1)
        deallocate(PANEL_array2)

    else if(selection_type == 3)then 

        ! asking method to adopt for the computation
        call ask_method(panel_type)            

        ! generate flow properties
        call setting_properties(P0,V,rho,alpha1,alpha2,start_angle,end_angle,dim,selection,selection_type)
        
        ! generate airfoil geometry 
        ! describing file name
        filename  = 'GNUplot_coord_dataGE.dat'
        filename1 = 'GNUplot_mean_dataGE.dat'
        filename2 = 'GNUplot_tg_normGE.dat'
        call ask_geometry(PANELsize,PANEL_array,MEAN_array,airfoil,alpha1,selection, & 
                          filename,filename1,filename2)

        ! autosaving geometry for plots
        call GNUplot_saving(PANEL_array,MEAN_array,airfoil%get_npoints(),filename,filename1,filename2)
        
        ! ground panels generation             
        call generate_ground(GROUNDpanel,GROUNDsize)
             
        ! compute system matrix 
        call computeGROUNDmatrix(matrix,PANEL_array,GROUNDpanel,PANELsize,GROUNDsize,panel_type)

        ! compute known vector 
        call computeGROUNDvector(vector,PANEL_array,GROUNDpanel,PANELsize,GROUNDsize,real(0.0,8),V)

        ! computing solution
        call solveGROUND(solution,matrix,vector,PANELsize,GROUNDsize)

        ! computing ground velocity 
        call compute_GROUNDpanelVEL(PANELsize,GROUNDsize,PANEL_array,GROUNDpanel,solution,real(0.0,8),V,panel_type) 

        !!!!!!!!!!!!!!!!!! COMPUTING VELOCITY FIELD !!!!!!!!!!!!!!!!!!!!
        ! high demanding process 
        ! -- it depends on the number of point for the field discretization 
        filename = 'FLOWfieldGE.dat'
        call computeGROUNDfield(solution,PANEL_array,GROUNDpanel,PANELsize,GROUNDsize,P0,real(0.0,8),V,rho,panel_type,filename)
        !!!!!!!!!!!!!!!!!! COMPUTING VELOCITY FIELD !!!!!!!!!!!!!!!!!!!!

        ! computing cp, pressure and velocity around the airfoil 
        allocate(cp_vec(PANELsize))
        call compute_airfoilFIELD(solution,PANEL_array,GROUNDpanel,PANELsize,GROUNDsize,cp_vec,P0,real(0.0,8),V,rho,panel_type)
                    
        ! compute CL value 
        CL = compute_cl(solution(PANELsize+1),PANEL_array,PANELsize)
        print*, 'CL value = ', CL    

        ! deallocation process
        deallocate(matrix)
        deallocate(vector)
        deallocate(cp_vec)        
        deallocate(solution)
        deallocate(MEAN_array)
        deallocate(GROUNDpanel)
        deallocate(PANEL_array)

    end if 

end program aeroHS
