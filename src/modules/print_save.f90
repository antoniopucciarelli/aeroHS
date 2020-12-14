module print_save

    contains

    subroutine print_mean_panel(MEAN_array,PANEL_array,MEANsize,PANELsize)
    ! this subroutine prints the airfoil data on the screen
            use PANEL_object
            use MEANline_object
            implicit none
    
            type(panel),dimension(:),intent(in)    :: PANEL_array
            type(MEANline),dimension(:),intent(in) :: MEAN_array
            integer(kind=4)                        :: i
            integer(kind=4),intent(in)             :: MEANsize, PANELsize
    
            do i=1,MEANsize
                print('(A18,    I4)'), 'MEANline object # ', i
                print('(A27,    I4)'), '   id                    = ', MEAN_array(i)%get_id()
                print('(A27, 2F8.4)'), '   coordinates           = ', MEAN_array(i)%get_coords()
                print('(A27, 2F8.4)'), '   thickness             = ', MEAN_array(i)%get_thickness()
                print('(A27,  F8.4)'), '   theta                 = ', MEAN_array(i)%get_theta()
                print*, new_line('A')
            end do
    
            do i=1,PANELsize
                print('(A15,    I4)'), 'PANEL object # ', i
                print('(A27,    I4)'), '   id                    = ', PANEL_array(i)%get_id()
                print('(A27,  F8.4)'), '   length                = ', PANEL_array(i)%get_length()
                print('(A27, 2F8.4)'), '   mid-point coordinates = ', PANEL_array(i)%get_midpointx(), PANEL_array(i)%get_midpointy()
                print('(A27, 2F8.4)'), '   starting-coordinates  = ', PANEL_array(i)%get_coords1()
                print('(A27, 2F8.4)'), '   ending-coordinates    = ', PANEL_array(i)%get_coords2()
                print('(A27, 2F8.4)'), '   tangent               = ', PANEL_array(i)%get_tangentx(),  PANEL_array(i)%get_tangenty()
                print('(A27, 2F8.4)'), '   normal                = ', PANEL_array(i)%get_normalx(),   PANEL_array(i)%get_normaly()
                print*, new_line('A')
            end do 
      
    end subroutine print_mean_panel
        
    subroutine print_matrix(matrix,matrix_size,maxsize)
    ! this subroutine prints the matrix computed by subroutine => compute_matrix
    ! maxsize => it allows printng a squared matrix [1:maxsize , 1:maxsize]
        use FOUL
        implicit none 

        real(kind=8),dimension(:,:),intent(in) :: matrix
        integer(kind=4)                        :: i, j, maxsize
        integer(kind=4),intent(in)             :: matrix_size

        j = 0

        if (maxsize>matrix_size)then
            maxsize = matrix_size
            j = 1
        end if 

        do i=1,maxsize
                print*, matrix(i,1:maxsize)
        end do 

        if (j.ne.0)then
            call write_formatted('YOU HAVE INSERTED AS INPUT A NUMBER','red', &
            ' > THAN THE MATRIX DIMENSION','red')
        end if

    end subroutine print_matrix
    
    subroutine print_vector(vector,vector_size,maxsize)
    ! this subroutine prints the vector element
    ! maxsize => it allows printing the first maxsize elements => [1:maxsize] 
        use FOUL
        implicit none 

        real(kind=8),dimension(:),intent(in) :: vector 
        integer(kind=4)                      :: i, j, vector_size, maxsize

        j = 0 

        if (maxsize>vector_size)then
            maxsize = vector_size
            j = 1
        end if 

        do i=1,maxsize
            print*, vector(i)
        end do

        if (j .ne. 0)then
            call write_formatted('YOU HAVE INSERTED AS INPUT A NUMBER','red', &
            ' > THAN THE VECTOR DIMENSION','red')
        end if 

    end subroutine print_vector
    
    subroutine save_matrix(matrix,matrix_size,filename)
    ! this subroutine saves the matrix data computed by compute_matrix in matrix.dat file
        use FOUL
        implicit none 
        
        real(kind=8),dimension(:,:),intent(in) :: matrix 
        integer(kind=4),intent(in)             :: matrix_size
        integer(kind=4)                        :: i
        character(len=30)                      :: filename

        open(unit=1, file=filename, status='replace')

        do i=1,matrix_size
            write(1,*) matrix(i,:)
        end do

        call write_formatted('[','normal','OK','green','] -- matrix saved --> ','normal',filename,'normal')
        
        close(1)

    end subroutine save_matrix 
    
    subroutine save_vector(vector,vector_size,filename)
    ! this subroutine saves the vector data computed by compute_vector in vector.dat file
        use FOUL
        implicit none 
        
        real(kind=8),dimension(:),intent(in) :: vector
        integer(kind=4),intent(in)           :: vector_size
        integer(kind=4)                      :: i
        character(len=30)                    :: filename
    
        open(unit=1, file=filename, status='replace')

        do i=1,vector_size
            write(1,*) vector(i)
        end do

        call write_formatted('[','normal','OK','green','] -- vector saved --> ','normal',filename,'normal')

        close(1)

    end subroutine save_vector 

    subroutine print_vortex_data(PANEL_ith,PANEL_jth,ROT,vortex1,vortex2)
        use PANEL_object
        use FOUL
        implicit none
        type(panel),intent(in)      :: PANEL_ith, PANEL_jth
        real(kind=8),dimension(2,2) :: ROT
        real(kind=8),dimension(2)   :: vortex1, vortex2
    
        print*, 'PANEL_ith = ', PANEL_ith%get_id()
        print*, 'PANEL_jth = ', PANEL_jth%get_id()
        print*, new_line('(A)')
        print*, 'TANGENT VECTOR ith-panel => tg(ith) = ', PANEL_ith%get_tangentx(), PANEL_ith%get_tangenty()
        print*, 'NORMAL  VECTOR ith-panel => nr(ith) = ', PANEL_ith%get_normalx(),  PANEL_ith%get_normaly()
        print*, 'TANGENT VECTOR jth-panel => tg(ith) = ', PANEL_jth%get_tangentx(), PANEL_jth%get_tangenty()
        print*, 'NORMAL  VECTOR jth-panel => nr(ith) = ', PANEL_jth%get_normalx(),  PANEL_jth%get_normaly()
        print*, new_line('(A)'), 'ROTATION MATRIX'
        print('(2(F8.4), / , 2(F8.4))'), ROT 
        print*, 'VORTEX IJ'
        print('(2(F8.4))'), vortex1
        print*, 'VORTEX ij'
        print('(2(F8.4))'),  vortex2
        print*, new_line('(A)')
    
    end subroutine print_vortex_data

    subroutine ask_to_save_matrix_vector(PANELsize,matrix,vector,solution)
        implicit none 
        character(len=30)                                          :: filename
        character(len=1)                                           :: resp
        integer(kind=4),intent(in)                                 :: PANELsize
        real(kind=8),dimension(PANELsize+1,PANELsize+1),intent(in) :: matrix 
        real(kind=8),dimension(PANELsize+1),intent(in)             :: vector 
        real(kind=8),dimension(PANELsize+1),intent(in)             :: solution
        real(kind=8),dimension(PANELsize+1)                        :: error 

        print*, 'do you want to save [matrix, known_vector, solution_vector]? [Y\n]'
        read*, resp
        
        if(resp=='Y' .or. resp=='y')then
            !!!!!!!!!!!!!!!!!!!!!!! SAVING VECTORS & MATRIX !!!!!!!!!!!!!!!!!!!!!!!!!
            filename = 'matrix.dat'
            call save_matrix(matrix,PANELsize+1,filename)
            filename = 'vector.dat'
            call save_vector(vector,PANELsize+1,filename)
            filename = 'solution.dat'
            call save_vector(solution,PANELsize+1,filename)
            filename = 'error.dat'
            error = matmul(matrix,solution) - vector 
            call save_vector(error,PANELsize+1,filename)    
            !!!!!!!!!!!!!!!!!!!!!!! SAVING VECTORS & MATRIX !!!!!!!!!!!!!!!!!!!!!!!!!
        end if

    end subroutine ask_to_save_matrix_vector

end module print_save
