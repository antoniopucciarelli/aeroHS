module MEANline_object
    
    type MEANline
        !
        ! mean line object 
        ! id        = mean line identification number
        ! coords    = [coords_x, coords_y] = position vector
        ! thickness = thickness of the airfoil relative the MEANline object
        ! dy        = mean line curvature gradient
        ! theta     = mean line 
        !
        integer(kind=4)           :: id        = 0
        real(kind=8),dimension(2) :: coords    = (/0.0, 0.0/)
        real(kind=8)              :: thickness = 0.0
        real(kind=8)              :: dy        = 0.0
        real(kind=8)              :: theta     = 0.0

        contains
        
        !!!!!!!!!!!!!!! GET FUNCTION - PASS PROCEDURE !!!!!!!!!!!!!!!
            procedure, pass(this) :: get_id
            procedure, pass(this) :: get_coords
            procedure, pass(this) :: get_dy
            procedure, pass(this) :: get_theta
            procedure, pass(this) :: get_thickness
        !!!!!!!!!!!!!!! GET FUNCTION - PASS PROCEDURE !!!!!!!!!!!!!!!
        procedure, pass(this) :: set_id
        procedure, pass(this) :: set_coordx
        procedure, pass(this) :: set_coordy_leading
        procedure, pass(this) :: set_coordy_trailing
        procedure, pass(this) :: set_gradient_leading
        procedure, pass(this) :: set_gradient_trailing
        procedure, pass(this) :: set_thickness
        procedure, pass(this) :: compute_UPcoords
        procedure, pass(this) :: compute_DOWNcoords
        procedure, pass(this) :: compute_transl
        procedure, pass(this) :: saving

    end type MEANline

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! GET functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer(kind=4) function get_id(this)
            implicit none 
            class(MEANline),intent(in) :: this

            get_id = this%id
        end function get_id

        real(kind=8) function get_dy(this)
            implicit none
            class(MEANline),intent(in) :: this
            
            get_dy = this%dy
        end function get_dy

        function get_coords(this)
            implicit none 
            real(kind=8),dimension(2)  :: get_coords
            class(MEANline),intent(in) :: this
            
            get_coords = this%coords
        end function get_coords

        real(kind=8) function get_theta(this)
            implicit none 
            class(MEANline),intent(in) :: this
            
            get_theta = this%theta
        end function get_theta

        real(kind=8) function get_thickness(this)
            implicit none 
            class(MEANline),intent(in) :: this
            
            get_thickness = this%thickness
        end function get_thickness
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! GET functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine set_id(this,x)
        implicit none 
        integer(kind=4),intent(in)    :: x
        class(MEANline),intent(inout) :: this

        this%id = x
    end subroutine set_id

    subroutine set_coordx(this,dim,i)
        use math_module
        implicit none 
        class(MEANline),intent(inout) :: this
        integer(kind=4),intent(in)    :: dim 
        integer(kind=4),intent(in)    :: i

        this%coords(1) = 0.5 * (1.0 - cos(pi*(real(i,8)-1.0)/real(dim,8)))
    end subroutine set_coordx

    subroutine set_coordy_leading(this,airfoildata1,airfoildata2)
        implicit none
        class(MEANline),intent(inout) :: this
        real(kind=8),intent(in)       :: airfoildata1, airfoildata2
        real(kind=8)                  :: x

        x = this%coords(1)
        
        this%coords(2) = airfoildata1*(2.0*airfoildata2*x - x**2)/(airfoildata2**2)
    end subroutine set_coordy_leading

    subroutine set_coordy_trailing(this,airfoildata1,airfoildata2)
        implicit none
        class(MEANline),intent(inout) :: this
        real(kind=8),intent(in)       :: airfoildata1, airfoildata2
        real(kind=8)                  :: x
        
        x = this%coords(1)

        this%coords(2) = airfoildata1*(1.0 - 2.0*airfoildata2 + 2.0*airfoildata2*x - x**2)/((1.0-airfoildata2)**2)
    end subroutine set_coordy_trailing

    subroutine set_gradient_leading(this,airfoildata1,airfoildata2)
        implicit none 
        class(MEANline),intent(inout) :: this
        real(kind=8),intent(in)       :: airfoildata1, airfoildata2
        real(kind=8)                  :: x

        x = this%coords(1)
    
        this%dy = 2.0*airfoildata1*(airfoildata2 - x)/(airfoildata2**2)
        this%theta = atan(this%dy)
    end subroutine set_gradient_leading

    subroutine set_gradient_trailing(this,airfoildata1,airfoildata2)
        implicit none
        class(MEANline),intent(inout) :: this
        real(kind=8),intent(in)       :: airfoildata1, airfoildata2
        real(kind=8)                  :: x

        x = this%coords(1)
        
        this%dy = 2.0*airfoildata1*(airfoildata2 - x)/((1.0 - airfoildata2)**2)
        this%theta = atan(this%dy)
    end subroutine set_gradient_trailing

    subroutine set_thickness(this,th)
        implicit none
        class(MEANline),intent(inout) :: this
        integer(kind=4),intent(in)    :: th  
    
        this%thickness = real(th,8)/20 * (0.2969*this%coords(1)**0.5 - 0.126*this%coords(1) &
                        -0.3516*this%coords(1)**2 + 0.28430*this%coords(1)**3 - 0.1036*this%coords(1)**4)
    end subroutine set_thickness

    subroutine compute_UPcoords(this,coordx,coordy)
        implicit none 
        class(MEANline),intent(in) :: this
        real(kind=8),intent(out)   :: coordx, coordy
        
        coordx = this%coords(1) - this%thickness * sin(this%theta) 
        coordy = this%coords(2) + this%thickness * cos(this%theta)  
    end subroutine compute_UPcoords

    subroutine compute_DOWNcoords(this,coordx,coordy)
        implicit none 
        class(MEANline),intent(in) :: this
        real(kind=8),intent(out)   :: coordx, coordy
        
        coordx = this%coords(1) + this%thickness * sin(this%theta) 
        coordy = this%coords(2) - this%thickness * cos(this%theta)  
    end subroutine compute_DOWNcoords

    subroutine compute_transl(this,airfoil)
        use AIRFOIL_object
        implicit none
        class(NACA_airfoil),intent(in) :: airfoil
        class(MEANline),intent(inout)  :: this

        this%coords = this%coords + airfoil%transl
    end subroutine compute_transl

    subroutine saving(this,writing_file)
        implicit none 
        class(MEANline),intent(in) :: this
        integer,intent(in)         :: writing_file

        write(writing_file,'(A16)')        'MEAN-LINE OBJECT'
        write(writing_file,'(A27,T30,I4)') '    id                 : ', this%get_id()
        write(writing_file,'(A27, 2F8.4)') '    point coords [x,y] : ', this%get_coords()
        write(writing_file,'(A27,  F8.4)') '    thickness          : ', this%get_thickness()
        write(writing_file,'(A27,  F8.4)') '    theta              : ', this%get_theta()
        write(writing_file,*) new_line('A')

    end subroutine saving
    
end module MEANline_object  