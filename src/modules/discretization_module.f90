module discretization_module

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!! PANEL & MEANline RELATED !!!!!!!!!!!!!!!!!!!!!!!!!
        recursive subroutine SCALINGfunction(x,y,scale)
        ! this subroutine modifies the panel points position  
            implicit none
            real(kind=8),intent(in)    :: scale
            real(kind=8),intent(inout) :: x
            real(kind=8),intent(inout) :: y

            ! scaling 
            x = x * scale
            y = y * scale
        end subroutine SCALINGfunction

        recursive subroutine rot(coordx,coordy,alpha)
        ! this subroutine computes the rotation of the vector of coords [coordx, coordy]
            use math_module
            implicit none
            real(kind=8),intent(inout) :: coordx, coordy ! input coordinates
            real(kind=8)               :: angle          ! angle between x axis and [coordx, coordy]
            real(kind=8)               :: radius         ! lenght of [coordx, coordy] vector 
            real(kind=8),intent(in)    :: alpha          ! AOA => describes the angle of rotation of the vector

            angle  = atan2(coordy,coordx)

            angle  = angle - alpha
            radius = abs(sqrt(coordx**2 + coordy**2))

            coordx = radius * cos(angle)
            coordy = radius * sin(angle)

        end subroutine rot

        subroutine check_LE_panels(PANELarray,dim)
        ! this subroutine checks the normal and tangent vectors in the proximity of the leading edge
            use PANEL_object
            implicit none

            integer(kind=4),intent(in)                   :: dim
            integer(kind=4)                              :: interval
            integer(kind=4)                              :: i
            type(panel),dimension(2*dim-2),intent(inout) :: PANELarray

            ! interval variable description
            ! this variable describes the length of the check interval related to the total panels used 
            ! to discretize the airfoil geometry
            interval = dim/10

            if(interval == 0)then
                interval = 1 
            end if

            ! this loop checks all the panel at the leading edge subourb
            do i=dim-interval,dim+interval 
                if(PANELarray(i)%get_normalx() > 0)then 
                    PANELarray(i)%normal  = - PANELarray(i)%normal
                    PANELarray(i)%tangent = - PANELarray(i)%tangent
                end if
            end do

        end subroutine check_LE_panels
    !!!!!!!!!!!!!!!!!!!!!!!!! PANEL & MEANline RELATED !!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!! GNUplot PROCESS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine ask_to_continue(i)
            implicit none
            character(len=1)              :: resp
            integer(kind=4),intent(inout) :: i 
            
            print*, 'do you want to create a new geometry? [Y\n]'
            read*, resp
            
            if(resp=='Y' .or. resp=='y')then
                i = 1
            else 
                i = 0
            end if
        end subroutine ask_to_continue

        subroutine ask_and_save(airfoil,PANELarray,MEANLINEarray)
            use AIRFOIL_object
            use PANEL_object
            use MEANline_object

            implicit none
            class(NACA_airfoil),intent(in)          :: airfoil    
            class(panel),intent(in),dimension(:)    :: PANELarray
            class(MEANline),intent(in),dimension(:) :: MEANLINEarray
            character(len=1)                        :: resp
            character(len=30)                       :: filename
            integer(kind=4)                         :: dim, k, writing_file

            print*, 'do you want to store data? [Y\n]'
            read*, resp

            if(resp=='Y' .or. resp=='y')then
                call airfoil%saving()

                dim = airfoil%get_npoints()
                read(airfoil%airfoilname,*) filename
                filename(9:30) = '_airfoil_MEAN.dat'

                writing_file = 1

                open(unit=writing_file, file=filename, status='replace')

                write(1,*) '# of mean line objects : ', airfoil%get_npoints()
                write(1,*) new_line('A')

                do k=1,dim
                    call MEANLINEarray(k)%saving(writing_file)
                end do

                close(writing_file)

                filename(9:30) = '_airfoil_PANEL.dat'
                writing_file = 1
                open(unit=writing_file, file=filename, status='replace')

                write(1,*) '# of panel objects : ', 2*airfoil%get_npoints() - 2
                write(1,*) new_line('A')

                do k=1,2*dim-2 
                    call PANELarray(k)%saving(writing_file)
                end do

                close(writing_file)
            end if

        end subroutine ask_and_save

        subroutine GNUplot_print(airfoil,PANELarray,MEANLINEarray,GNUplot_coord_data,GNUplot_mean_data,GNUplot_tg_norm)
            use AIRFOIL_object
            use MEANline_object
            use PANEL_object
            implicit none   
            class(NACA_airfoil),intent(in)          :: airfoil
            class(MEANline),intent(in),dimension(:) :: MEANLINEarray
            class(panel),intent(in),dimension(:)    :: PANELarray
            character(len=1)                        :: resp
            character(len=30)                       :: GNUplot_coord_data
            character(len=30)                       :: GNUplot_mean_data
            character(len=30)                       :: GNUplot_tg_norm
            
            ! printing option via gnuplot
            print*, 'do you want to print ',airfoil%get_airfoilname(),' ? [Y\n]'
            read*, resp
            if(resp=='Y' .or. resp=='y')then
                call GNUplot_saving(PANELarray,MEANLINEarray,airfoil%get_npoints(), & 
                                    GNUplot_coord_data,GNUplot_mean_data,GNUplot_tg_norm)
                call system('gnuplot -p AIRFOILgnuplot.plt')
                call system('gnuplot -p AIRFOILelement_plot.plt')
            end if
            
        end subroutine GNUplot_print

        subroutine GNUplot_saving(PANELarray,MEANLINEarray,dim,GNUplot_coord_data,GNUplot_mean_data,GNUplot_tg_norm)
            use MEANline_object
            use PANEL_object
            implicit none
            class(MEANline),intent(in),dimension(:) :: MEANLINEarray
            class(panel),intent(in),dimension(:)    :: PANELarray
            integer(kind=4),intent(in)              :: dim
            integer(kind=4)                         :: k
            character(len=30)                       :: GNUplot_coord_data
            character(len=30)                       :: GNUplot_mean_data
            character(len=30)                       :: GNUplot_tg_norm
        
            open(unit=1, file=GNUplot_coord_data, status='replace')
            open(unit=2, file=GNUplot_mean_data,  status='replace')
            open(unit=3, file=GNUplot_tg_norm,    status='replace')

            do k=1,dim
                write(1,'(2F8.4)') MEANLINEarray(k)%get_coords()
            end do

            do k=1,2*dim-2
                write(2,'(2F8.4)') PANELarray(k)%get_coords1()
            end do
            write(2,'(2F8.4)') PANELarray(1)%get_coords1()

            write(2,'(2F8.4)') PANELarray(1)%get_coords2()

            do k=1,2*dim-2
                write(3,'(6F8.4)') PANELarray(k)%get_midpointx(), &
                                PANELarray(k)%get_midpointy(), &
                                PANELarray(k)%get_tangentx() , &
                                PANELarray(k)%get_tangenty() , & 
                                PANELarray(k)%get_normalx()  , &
                                PANELarray(k)%get_normaly()
            end do

            close(1)
            close(2)
            close(3)
        end subroutine GNUplot_saving
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! GNUplot PROCESS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!! GNUplot -excess- PROCESS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine savevector(MEANarray,dim)
            use MEANline_object
            implicit none

            integer(kind=4)               :: dim, j
            type(MEANline),dimension(dim) :: MEANarray

            open(unit=1, file='mean.dat', status='replace')

            do j=1,dim
                write(1,'(2F8.4)') MEANarray(j)%coords(1), MEANarray(j)%coords(2)
            end do

            close(1)

        end subroutine savevector

        subroutine savevectors(vec1,vec2,vec3,vec4,dim)
            implicit none

            integer(kind=4)             :: dim, j
            real(kind=8),dimension(dim) :: vec1, vec2, vec3, vec4

            open(unit=1, file='UP.dat', status='replace')
            open(unit=2, file='DOWN.dat', status='replace')

            do j=1,dim
                write(1,'(2F8.4)') vec1(j), vec2(j)
            end do
            do j=1,dim
                write(1,'(2F8.4)') vec3(j), vec4(j)
            end do

            close(1)
            close(2)

        end subroutine savevectors
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! GNUplot -excess- PROCESS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
end module discretization_module
