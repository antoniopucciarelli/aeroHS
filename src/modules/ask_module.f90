module ask_module 
    
    contains 


    subroutine setting_properties(P0,V,rho,alpha,start_angle,end_angle,dim,selection,selection_type)
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
        integer(kind=4),intent(in)    :: selection_type
        integer(kind=4)               :: x 

        x = 1
        
        if(selection_type == 1)then  
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
                        V   = 1
                        rho = 1 
                        P0  = 1
                        x   = 0

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
        else 
            do while(x==1)
                print*, 'to compute Cp       vs X --> type(1)'
                print*, 'to compute pressure vs X --' 
                print*, '           velocity vs X --> type(2)'
                print*, '           Cp       vs X --'
                read*, selection

                select case (selection)
                    case(1)
                        print*, 'type airfoil angle of attack -- AOA [deg]'
                        read*, alpha
                        V   = 1
                        rho = 1 
                        P0  = 1
                        x   = 0

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
                    
                    case default 
                        print*, 'you have selected an invalid action', new_line('(A)'), 'type again'

                end select 

            end do
        end if 

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

end module ask_module 
