module PANEL_object

    type panel
        !
        ! panel object
        ! id                                = panel identification number
        ! midpoint = [midpointx, midpointy] = panel middle-point position vector 
        ! tangent  = [tangentx, tangenty]   = panel tangent vector coordinates
        ! normal   = [normalx, normaly]     = panel normal vector coordinates
        ! length                            = panel length
        ! coords1  = [coodrs1_x, coords1_y] = panel starting point
        ! coords2  = [coodrs2_x, coords2_y] = panel ending point
        !
        integer(kind=4)           :: id       = 0
        real(kind=8),dimension(2) :: midpoint = (/0.0, 0.0/) ! [x, y]
        real(kind=8),dimension(2) :: tangent  = (/0.0, 0.0/) ! [x, y]
        real(kind=8),dimension(2) :: normal   = (/0.0, 0.0/) ! [x, y]
        real(kind=8)              :: length   = 0.0
        real(kind=8),dimension(2) :: coords1  = (/0.0, 0.0/) ! [x, y]
        real(kind=8),dimension(2) :: coords2  = (/0.0, 0.0/) ! [x, y]
        real(kind=8)              :: angle    = 0.0          ! rads
        character(len=2)          :: POS      = 'ND'         ! [UP => upper surface; DW => lower surface; ND => not defined]

        contains

        !!!!!!!!!!!!!!! GET FUNCTION - PASS PROCEDURE !!!!!!!!!!!!!!!
        procedure, pass(this) :: get_id
        procedure, pass(this) :: get_midpointx
        procedure, pass(this) :: get_midpointy
        procedure, pass(this) :: get_tangentx
        procedure, pass(this) :: get_tangenty
        procedure, pass(this) :: get_normalx
        procedure, pass(this) :: get_normaly
        procedure, pass(this) :: get_length
        procedure, pass(this) :: get_coords1
        procedure, pass(this) :: get_coords2
        procedure, pass(this) :: get_angle
        procedure, pass(this) :: get_position
        !!!!!!!!!!!!!!! GET FUNCTION - PASS PROCEDURE !!!!!!!!!!!!!!!    
        procedure, pass(this) :: set_coords
        procedure, pass(this) :: set_id
        procedure, pass(this) :: set_angle
        procedure, pass(this) :: set_position
        procedure, pass(this) :: SCALINGfunc
        procedure, pass(this) :: compute_length
        procedure, pass(this) :: compute_tangent_and_normal
        procedure, pass(this) :: compute_transl
        procedure, pass(this) :: compute_midpoint
        procedure, pass(this) :: saving

    end type panel

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! GET functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer(kind=4) function get_id(this)
            implicit none
            class(panel), intent(in) :: this
            get_id = this%id
        end function get_id
    
        real(kind=8) function get_midpointx(this)
            implicit none
            class(panel), intent(in) :: this

            get_midpointx = this%midpoint(1)
        end function get_midpointx
    
        real(kind=8) function get_midpointy(this)
            implicit none
            class(panel), intent(in) :: this

            get_midpointy = this%midpoint(2)
        end function get_midpointy
    
        real(kind=8) function get_tangentx(this)
            implicit none 
            class(panel), intent(in) :: this

            get_tangentx = this%tangent(1)
        end function get_tangentx
    
        real(kind=8) function get_tangenty(this)
            implicit none
            class(panel), intent(in) :: this

            get_tangenty = this%tangent(2)
        end function get_tangenty
    
        real(kind=8) function get_normalx(this)
            implicit none
            class(panel), intent(in) :: this

            get_normalx = this%normal(1)
        end function get_normalx
    
        real(kind=8) function get_normaly(this)
            implicit none 
            class(panel), intent(in) :: this

            get_normaly = this%normal(2)
        end function get_normaly
    
        real(kind=8) function get_length(this)
            implicit none
            class(panel), intent(in) :: this

            get_length = this%length
        end function get_length
    
        function get_coords1(this)
            implicit none
            real(kind=8),dimension(2) :: get_coords1
            class(panel), intent(in) :: this
        
            get_coords1 = this%coords1
        end function get_coords1
    
        function get_coords2(this)
            implicit none
            real(kind=8),dimension(2) :: get_coords2
            class(panel), intent(in) :: this
        
            get_coords2 = this%coords2
        end function get_coords2

        real(kind=8) function get_angle(this)
            implicit none
            class(panel), intent(in) :: this
            get_angle = this%angle
        end function get_angle

        character(len=2) function get_position(this)
            implicit none
            class(panel), intent(in) :: this
            get_position = this%POS
        end function get_position
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! GET functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine set_coords(this,coordx1,coordy1,coordx2,coordy2)
        implicit none

        class(panel),intent(inout) :: this
        real(kind=8),intent(in)    :: coordx1
        real(kind=8),intent(in)    :: coordx2
        real(kind=8),intent(in)    :: coordy1 
        real(kind=8),intent(in)    :: coordy2

        this%coords1(1) = coordx1
        this%coords1(2) = coordy1
        this%coords2(1) = coordx2
        this%coords2(2) = coordy2

    end subroutine set_coords

    subroutine set_id(this,ID)
        implicit none

        class(panel),intent(inout) :: this
        integer(kind=4),intent(in) :: ID  

        this%id = ID

    end subroutine set_id

    subroutine set_angle(this)
        implicit none

        class(panel),intent(inout) :: this
        real(kind=8)               :: dx
        real(kind=8)               :: dy

        dx = this%coords2(1) - this%coords1(1)
        dy = this%coords2(2) - this%coords1(2)

        this%angle = atan2(dy,dx)  
        
    end subroutine set_angle

    subroutine set_position(this,flag)
        implicit none 

        class(panel),intent(inout)  :: this
        character(len=2),intent(in) :: flag

        this%POS = flag
        
    end subroutine set_position

    subroutine SCALINGfunc(this,scale) 
        implicit none

        class(panel),intent(inout) :: this
        real(kind=4),intent(in)    :: scale

        ! scaling 
        this%coords1(:) = this%coords1(:) * scale
        this%coords2(:) = this%coords2(:) * scale

    end subroutine SCALINGfunc

    subroutine compute_length(this)
        use math_module
        implicit none

        class(panel),intent(inout) :: this
        real(kind=8),dimension(2)  :: coord
        
        coord = this%get_coords1() - this%get_coords2() 
        
        this%length = norm(coord) 
    
    end subroutine compute_length

    subroutine compute_tangent_and_normal(this)
    ! this subroutine compute the normal and tangent vectors of a panel 
    ! it changes automatically the direction of each vector with the help
    !   of the panel object type POS attribute
    !   
    ! OUTPUT VECTORS:  
    !   UP position:      
    !       normal vector points upward
    !       tangent vector points rearward
    !   DW position:    
    !       normal vector points downward
    !       tangent vector points rearward
    !
    ! !!! KEEP IN MIND !!!
    ! there are some vectors that do not coincide with the panel source-vorticity convenction regarding their direction
    ! 
        use math_module
        use FOUL
        implicit none

        class(panel),intent(inout)  :: this
        real(kind=8)                :: coordx1
        real(kind=8)                :: coordx2
        real(kind=8)                :: coordy1
        real(kind=8)                :: coordy2
        real(kind=8)                :: dx
        real(kind=8)                :: dy
        real(kind=8)                :: theta
        real(kind=8)                :: cross_value
        real(kind=8),dimension(2)   :: temp
        character(len=2)            :: flag           

        ! starting point coords
        coordx1 = this%coords1(1)
        coordy1 = this%coords1(2)
        ! ending point coords
        coordx2 = this%coords2(1)
        coordy2 = this%coords2(2)

        flag = this%get_position()

        if (this%length /= 0.0) then 
            ! description of panel vector
            dx = coordx2 - coordx1
            dy = coordy2 - coordy1

            theta = atan2(dy,dx)
            
            this%tangent(1) =   cos(theta)
            this%tangent(2) =   sin(theta)
            this%normal(1)  = - sin(theta)
            this%normal(2)  =   cos(theta)
            
            if (flag == 'UP') then
                if (this%normal(2) < 0) then
                    this%normal  = - this%normal
                end if
            
                cross_value = cross(this%tangent,this%normal)
            
                if (cross_value < 0) then
                    this%tangent = - this%tangent
                end if

            else if (flag == 'DW') then
                if (this%normal(2) > 0 ) then
                    this%normal = - this%normal
                end if

                cross_value = cross(this%tangent,this%normal)
                
                if (cross_value > 0) then 
                    this%tangent = - this%tangent
                end if

            end if
        else 
            call write_formatted('[','normal','WARNING','red','] -- panel object with 0 length','normal')
        end if

    end subroutine compute_tangent_and_normal

    subroutine compute_midpoint(this,flag)
        use FOUL
        implicit none

        class(panel),intent(inout)  :: this
        real(kind=8)                :: norm_1
        real(kind=8)                :: norm_2
        character(len=5),intent(in) :: flag

        this%midpoint = (this%coords1 + this%coords2)/2

        if(flag == 'print')then 
            norm_1 = sqrt((this%midpoint(1)- this%coords1(1))**2 + (this%midpoint(2)- this%coords1(2))**2)
            norm_2 = sqrt((this%midpoint(1)- this%coords2(1))**2 + (this%midpoint(2)- this%coords2(2))**2)
            if(norm_1 == norm_2)then
                call write_formatted('[','normal','OK','green','] -- norm1 = norm2','normal')
            else
                call write_formatted('[','normal','WARNING','red','] -- norm1 = norm2','normal')
                print*, norm_1 - norm_2
                print*, log(norm_1/norm_2)
            end if
        end if
    
    end subroutine compute_midpoint

    subroutine compute_transl(this,airfoil)
        use AIRFOIL_object
        implicit none

        class(panel),intent(inout)     :: this
        class(NACA_airfoil),intent(in) :: airfoil

        this%coords1  = this%coords1  + airfoil%transl
        this%coords2  = this%coords2  + airfoil%transl

    end subroutine compute_transl

    subroutine saving(this,writing_file)
        implicit none 

        class(panel),intent(in) :: this
        integer,intent(in)      :: writing_file

        write(writing_file,'(A12)')         'PANEL OBJECT'
        write(writing_file,'(A27,T29, I4)') '    id                   : ', this%get_id()
        write(writing_file,'(A27,  F12.8)') '    length               : ', this%get_length()
        write(writing_file,'(A27, 2F12.8)') '    midpoint       [x,y] : ', this%get_midpointx(), this%get_midpointy() 
        write(writing_file,'(A27, 2F12.8)') '    starting point [x,y] : ', this%get_coords1()
        write(writing_file,'(A27, 2F12.8)') '    ending point   [x,y] : ', this%get_coords2()
        write(writing_file,'(A27, 2F12.8)') '    tangent vector [x,y] : ', this%get_tangentx(), this%get_tangenty()
        write(writing_file,'(A27, 2F12.8)') '    normal vector  [x,y] : ', this%get_normalx(), this%get_normaly()
        write(writing_file,'(A27,  F12.8)') '    angle          [rad] : ', this%get_angle()
        write(writing_file,'(A27,     A2)') '    position             : ', this%get_position()
        write(writing_file,*) new_line('A')   

    end subroutine saving

end module PANEL_object 