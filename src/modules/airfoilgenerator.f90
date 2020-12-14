!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the aim of the module is to calculate the geometry of a NACA**** airfoil  !
!   -> the airfoil must be NACA 4-digits                                    !
!       -> there are different options for saving and displaying data       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module airfoilgenerator

    contains

        subroutine make_airfoil(PANELsize,MEANLINEarray,PANELarray,airfoil,alpha)

        ! module declaration
        use AIRFOIL_object
        use PANEL_object
        use MEANline_object
        use discretization_module
        use FOUL
        use math_module  

        implicit none

        ! variable declaration
        type(NACA_airfoil),intent(out)                      :: airfoil       ! airfoil object
        type(panel),allocatable,dimension(:),intent(out)    :: PANELarray    ! panel object array
        type(MEANline),allocatable,dimension(:),intent(out) :: MEANLINEarray ! mean line point object array
        real(kind=8),allocatable,dimension(:)               :: coordyUP      ! x-coords array -> describes the UPPER airfoil geometry
        real(kind=8),allocatable,dimension(:)               :: coordxUP      ! y-coords array -> describes the UPPER airfoil geometry
        real(kind=8),allocatable,dimension(:)               :: coordyDW      ! x-coords array -> describes the LOWER airfoil geometry
        real(kind=8),allocatable,dimension(:)               :: coordxDW      ! y-coords array -> describes the LOWER airfoil geometry
        integer(kind=4)                                     :: counter = 1   ! # of airfoil studied counter
        integer(kind=4)                                     :: x = 1         ! checking variable in while loop
        integer(kind=4)                                     :: dim           ! auxiliary variable -> # of discretisation points for the airfoil
        integer(kind=4)                                     :: k             ! auxiliary variable
        integer(kind=4)                                     :: j             ! auxiliary variable
        integer(kind=4)                                     :: i             ! auxiliary variable
        real(kind=8)                                        :: airfoil_data1 ! easy-access variable -> 1st number of airfoil%data
        real(kind=8)                                        :: airfoil_data2 ! easy-access variable -> 2nd number of arifoil%data
        real(kind=8),intent(in)                             :: alpha         ! AOA 
        integer(kind=4),intent(out)                         :: PANELsize     ! number of discretization panels
        real(kind=8),dimension(2)                           :: transl        ! translation vector        

        ! this program allows you to create multiple NACA**** profile each run 

        call airfoil%set_AIRFOILname() ! setting airfoil name

        call airfoil%set_npoints()     ! setting airfoil # of discretisation points
        
        ! setting properties 
        ! PAY ATTENTION the properties are already set in the setting_properties() subroutine
        airfoil%AOA     = alpha*180/pi
        airfoil%scaling = 1.0

        dim = airfoil%get_npoints()

        ! data allocation procedure
        allocate(coordxUP(1:dim)) 
        allocate(coordyUP(1:dim))
        allocate(coordxDW(1:dim)) 
        allocate(coordyDW(1:dim))
        allocate(PANELarray(1:2*dim-2)) 
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
                    ! asking LE position
                    print*, 'type leading edge position'
                    read*, transl(1), transl(2)

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

            ! print airfoil data 
            call airfoil%print_data()
            
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
                ! call GNUplot_print(airfoil,PANELarray,MEANLINEarray)
            !!!!!!!!!!!!!!!!!!!!!!! SAVING & GRAPHICS !!!!!!!!!!!!!!!!!!!!!!!!!
            
            ! data deallocation procedure
            deallocate(coordxUP) 
            deallocate(coordyUP)
            deallocate(coordxDW)
            deallocate(coordyDW)

            PANELsize = 2*dim-2

            call write_formatted('[','normal','OK','green','] -- airfoil geometry generated','normal')

    end subroutine make_airfoil

    subroutine ask_geometry(PANELsize,PANEL_array,MEAN_array,airfoil,alpha,selection)
    ! this subroutine allows to ask the user the action to take
        use cp 
        use AIRFOIL_object
        use PANEL_object
        use MEANline_object
        use FOUL

        implicit none 
        
        type(NACA_airfoil),intent(out)                      :: airfoil
        type(panel),dimension(:),allocatable,intent(out)    :: PANEL_array
        type(MEANline),dimension(:),allocatable,intent(out) :: MEAN_array
        integer(kind=4),intent(inout)                       :: PANELsize
        integer(kind=4)                                     :: MEANsize
        character(len=30)                                   :: filename
        integer(kind=4)                                     :: x
        integer(kind=4),intent(in)                          :: selection
        integer(kind=4)                                     :: selection_type 
        real(kind=8),intent(in)                             :: alpha 
        
        x = 1

        do while(x==1)
            
            call write_formatted('SETTING GEOMETRY','yellow')

            if(selection == 3)then 
                print*, 'creaing a new 4 digits airfoil'

                ! making new geometry from scratch
                call make_airfoil(PANELsize,MEAN_array,PANEL_array,airfoil,alpha)
                
                x = 0

            else 
                print*, 'load data from an existing .dat file created by previous runs  --'
                print*, '--- PAY ATTENTION --- in case of loading file, the airfoil AOA --'
                print*, '    is referred to the saved airfoil NOT the AOA just typed    --> type(1)'
                print*, 'create a new 4 digits airfoil from scratch                     --> type(2)'
                read*, selection_type

                select case(selection_type)
                    case(1)
                        ! extracting data from .dat files from data created by airfoilgenerator.f90
                        call manage_data(airfoil,MEAN_array,PANEL_array,MEANsize,PANELsize,filename)
                        
                        x = 0

                    case(2)
                        ! making new geometry form scratch 
                        call make_airfoil(PANELsize,MEAN_array,PANEL_array,airfoil,alpha)

                        x = 0 

                    case default 
                        print*, 'you have selected an invalid action', new_line('(A)'), 'type again'

                end select 
            end if 
        end do 
        
        print*, new_line('(A)')

    end subroutine ask_geometry 

end module airfoilgenerator
