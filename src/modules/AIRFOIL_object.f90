module AIRFOIL_object    

    type NACA_airfoil
        !
        ! NACA airfoil object
        ! id       = airfoil identification number 
        ! n_points = # of points for the airfoil discretisation 
        ! scaling  = airfoil scale model 
        ! [traslx, trasly] = translation vector in xy plane
        ! AOA      = angle of attack
        ! data     = array with all the information about the airfoil's dimensions
        !
        integer(kind=4)              :: id          = 0
        integer(kind=4)              :: n_points    = 0
        real(kind=8)                 :: scaling     = 0.0
        real(kind=8),dimension(2)    :: transl      = (/0.0, 0.0/) 
        real(kind=8)                 :: AOA         = 0.0
        integer(kind=4),dimension(3) :: data        = 0
        character(len=8)             :: airfoilname = 'NACA0000'

        contains

        procedure, pass(this) :: get_id
        procedure, pass(this) :: get_scaling
        procedure, pass(this) :: get_transl
        procedure, pass(this) :: get_translx
        procedure, pass(this) :: get_transly
        procedure, pass(this) :: get_AOA
        procedure, pass(this) :: get_data
        procedure, pass(this) :: get_npoints
        procedure, pass(this) :: get_airfoilname
        procedure, pass(this) :: set_AIRFOILname
        procedure, pass(this) :: set_npoints
        procedure, pass(this) :: set_scaling
        procedure, pass(this) :: set_AOA
        procedure, pass(this) :: set_transl
        procedure, pass(this) :: saving
        procedure, pass(this) :: print_data

    end type NACA_airfoil

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! GET functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer(kind=4) function get_id(this)
            implicit none
            class(NACA_airfoil),intent(in) :: this
            get_id = this%id
        end function get_id

        real(kind=8) function get_AOA(this)
            implicit none
            class(NACA_airfoil),intent(in) :: this
            get_AOA = this%AOA
        end function get_AOA

        function get_transl(this)
            implicit none
            real(kind=8),dimension(2)       :: get_transl
            class(NACA_airfoil),intent(in) :: this
            get_transl = this%transl
        end function get_transl

        real(kind=8) function get_translx(this)
            implicit none
            class(NACA_airfoil),intent(in) :: this
            get_translx = this%transl(1)
        end function get_translx

        real(kind=8) function get_transly(this)
            implicit none
            class(NACA_airfoil),intent(in) :: this
            get_transly = this%transl(2)
        end function get_transly

        real(kind=8) function get_scaling(this)
            implicit none
            class(NACA_airfoil),intent(in) :: this
            get_scaling = this%scaling
        end function get_scaling

        integer(kind=4) function get_npoints(this)
            implicit none
            class(NACA_airfoil),intent(in) :: this
            get_npoints = this%n_points
        end function get_npoints

        function get_data(this)
            implicit none
            integer(kind=4),dimension(3)   :: get_data 
            class(NACA_airfoil),intent(in) :: this
            get_data = this%data
        end function 

        character(len=8) function get_airfoilname(this)
            implicit none 
            class(NACA_airfoil),intent(in) :: this 

            get_airfoilname = this%airfoilname
        end function get_airfoilname
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! GET functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine set_AIRFOILname(this)
    ! this subroutine gets the airfoil name and converts its data into a vector
        implicit none
        class(NACA_airfoil),intent(inout) :: this
        character(len=8)                  :: airfoil_name
        integer(kind=4)                   :: i, j, k
        integer(kind=4),dimension(3)      :: airfoil_data
        
        print*, 'type 4 digits NACA airfoil'
        read*, airfoil_name 
        read(airfoil_name(5:8),*) i

        this%airfoilname = airfoil_name

        ! input check -> the airfoil must be NACA 4-digits
        ! if the numerical digits are more than 4, the array runs out
        if(airfoil_name(1:4) .eq. 'naca' .or.  airfoil_name(1:4) .eq. 'NACA') then
        else 
          print*, 'error with airfoil input'
          return
        end if
      
        ! airfoil data declaration  
        do j=1,2
            k = 10**(4-j)
            airfoil_data(j) = i/k
            i = i-airfoil_data(j)*k
        end do
        airfoil_data(3) = i
      
        ! 1st integer: max value of the ordinate in percent of the cord
        ! 2nd integer: absissa of respective max ordinate point 
        ! 3rd+4th integer: thickness in percentage of the cord
        
        this%data(:) = airfoil_data(:)
    end subroutine set_AIRFOILname

    subroutine set_npoints(this)
    ! this subroutine gets the number of point for the geometry discretisation
        implicit none
        class(NACA_airfoil),intent(inout) :: this
        integer(kind=4)                   :: dim
        ! number of reference points
        print*, 'type number of coords for the airfoil'
        read*, dim
        if(dim<3) then
            print*, 'invalid input'
            return
        end if
        this%n_points = dim
    end subroutine set_npoints

    subroutine set_scaling(this)
    ! this subroutine gets the scaling factor for the airfoil
        implicit none
        class(NACA_airfoil),intent(inout) :: this
        real(kind=4)                      :: scale

        ! wing scaling
        print*, 'type the scaling factor'
        read*, scale
        ! scaling condition
        if(scale<=0) then
            print*, 'error with scaling factor'
            return
        end if

        this%scaling = scale
    end subroutine set_scaling

    subroutine set_AOA(this)
        implicit none 
        class(NACA_airfoil),intent(inout) :: this
            
        ! wing rotation
        print*, 'type AOA, clockwise positive rotation'
        read*, this%AOA
    end subroutine set_AOA

    subroutine set_transl(this)
        implicit none
        class(NACA_airfoil),intent(inout) :: this

        ! wing translation
        print*, 'type (x,y) translation coordinates'
        read*, this%transl(1), this%transl(2)        
    end subroutine set_transl

    subroutine saving(this)
        implicit none
        class(NACA_airfoil),intent(in) :: this
        character(len=25)              :: filename

        read(this%airfoilname,*) filename
        filename(9:25) = '_airfoil.dat'

        open(unit=60, file=filename, status='replace')

        write(60,*) 'AIRFOIL OBJECT'
        write(60,*) '    name                  :', this%airfoilname
        write(60,*) '    discretisation points :', this%get_npoints()
        write(60,*) '    scaling factor        :', this%get_scaling()
        write(60,*) '    AOA                   :', this%get_AOA() 
        write(60,*) '    translation [x,y]     :', this%get_transl()
        write(60,*) new_line('A')

        close(60)
    end subroutine saving

    subroutine print_data(this)
        implicit none
        class(NACA_airfoil),intent(in) :: this

        print*, 'AIRFOIL OBJECT'
        print*, '    name                  : ', this%airfoilname
        print*, '    discretisation points : ', this%get_npoints()
        print*, '    scaling factor        : ', this%get_scaling()
        print*, '    AOA                   : ', this%get_AOA() 
        print*, '    translation [x,y]     : ', this%get_transl()
        print*, new_line('A')
    end subroutine print_data

end module AIRFOIL_object
