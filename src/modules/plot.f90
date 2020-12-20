module plot 

    contains

    subroutine plot_cp(cp_vec,PANEL_array,PANELsize,alpha)
    ! this subroutine plots the cp numbers over the aifoil geometry
        use FOUL
        use PANEL_object
        use discretization_module
        implicit none

        integer(kind=4),intent(in)                   :: PANELsize
        real(kind=8),dimension(PANELsize),intent(in) :: cp_vec
        type(panel),dimension(PANELsize),intent(in)  :: PANEL_array
        integer(kind=4)                              :: i
        real(kind=8),dimension(2)                    :: position_vector1
        real(kind=8),dimension(2)                    :: position_vector2
        real(kind=8),intent(in)                      :: alpha

        open(unit=1,file='cp_gnuplot.dat',status='replace')

        ! saving vectors
        do i=1,PANELsize/2
            
            position_vector1 = PANEL_array(i)%midpoint
            position_vector2 = PANEL_array(i+PANELsize/2)%midpoint

            call rot(position_vector1(1),position_vector1(2),-alpha)  
            call rot(position_vector2(1),position_vector2(2),-alpha)  

            write(1,*) position_vector1(1), - cp_vec(i), &
                       position_vector2(1), - cp_vec(PANELsize/2+i)
        end do

        call system('gnuplot -p CP_plot.plt')

        call write_formatted('[','normal','OK','green','] -- Cp plotted ','normal')

        close(1)

    end subroutine plot_cp
    
    subroutine plot_vel(Vvec,PANEL_array,PANELsize)
    ! this subroutine plots the velocity numbers over the aifoil geometry
        use FOUL
        use PANEL_object
        implicit none
        
        integer(kind=4),intent(in)                   :: PANELsize
        real(kind=8),dimension(PANELsize),intent(in) :: Vvec
        type(panel),dimension(PANELsize),intent(in)  :: PANEL_array
        integer(kind=4)                              :: i

        open(unit=1,file='VEL_gnuplot.dat',status='replace')

        ! saving vectors
        do i=1,PANELsize/2
            write(1,*) PANEL_array(i)%get_midpointx(), Vvec(i), &
                       PANEL_array(PANELsize/2+i)%get_midpointx(), Vvec(PANELsize/2+i)
        end do

        call system('gnuplot -p vel_plot.plt')

        call write_formatted('[','normal','OK','green','] -- velocity plotted ','normal')

        close(1)

    end subroutine plot_vel

    subroutine plot_pressure(pressure,PANEL_array,PANELsize)
    ! this subroutine plots the pressure numbers over the aifoil geometry
        use FOUL
        use PANEL_object
        implicit none

        integer(kind=4),intent(in)                   :: PANELsize
        real(kind=8),dimension(PANELsize),intent(in) :: pressure
        type(panel),dimension(PANELsize),intent(in)  :: PANEL_array
        integer(kind=4)                              :: i

        open(unit=1,file='PRESSURE_gnuplot.dat',status='replace')

        ! saving vectors
        do i=1,PANELsize/2
            write(1,*) PANEL_array(i)%get_midpointx(), pressure(i), &
                       PANEL_array(PANELsize/2+i)%get_midpointx(), pressure(PANELsize/2+i)
        end do

        call system('gnuplot -p pressure_plot.plt')

        call write_formatted('[','normal','OK','green','] -- pressure plotted ','normal')

        close(1)

    end subroutine plot_pressure

    subroutine plot_cl(cl_alpha,array_size)
    ! this subroutine plots the cl coefficients versus the angle of attack
        use FOUL
        use PANEL_object
        implicit none 

        real(kind=8),dimension(array_size,2),intent(in) :: cl_alpha  
        integer(kind=4)                                 :: array_size
        integer(kind=4)                                 :: i

        open(unit=1,file='cl_gnuplot.dat',status='replace')

        ! saving vectors
        do i=1,array_size
            write(1,*) cl_alpha(i,:)
        end do

        close(1)

        call system('gnuplot -p CL_plot.plt')

        call write_formatted('[','normal','OK','green','] -- Cl plotted ','normal')

    end subroutine plot_cl
    
    subroutine plot_vel_field()
        use FOUL
        implicit none 
 
        call system('gnuplot -p PRESSUREfield.plt')

        call write_formatted('[','normal','OK','green','] -- velocity field plotted ','normal')
    end subroutine plot_vel_field 
        
    subroutine plot_airfoils()
        use FOUL 
        implicit none 

        call system('gnuplot -p PLOTairfoils.plt')

        call write_formatted('[','normal','OK','green','] -- velocity field plotted ','normal')
    
    end subroutine plot_airfoils

end module plot 
