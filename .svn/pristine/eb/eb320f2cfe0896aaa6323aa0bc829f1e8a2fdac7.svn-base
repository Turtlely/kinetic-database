module kdbinsert
    implicit none

    contains


    subroutine atoms_dimensions(pos, ibox, n, dimensions)
        !! Returns side lengths (in direct coordinates) of a
        !! parallelepiped that can enclose the atoms provided
        !!
        !! Intended to measure size of an already-clumped atoms object,
        !! and thus does NOT apply PBC
        !!
        !! Parameters:
        !!      pos:  positions (n x 3 array) of n chosen atoms
        !!      ibox: inverse of the cell from which the atoms originate

        ! inputs
        integer, intent(in)          :: n
        double precision, intent(in) :: pos(n, 3)
        double precision, intent(in) :: ibox(3, 3)
        !output
        double precision :: dimensions(3)

!f2py   intent(in)   :: pos, ibox
!f2py   intent(hide) :: n=len(pos)
!f2py   intent(out)  :: dimensions

        ! other variables
        double precision :: vdir(3)
        integer          :: i, j

        dimensions = [0, 0, 0]  ! intialize to smallest values possible
        do i = 1, n-1
            do j = i+1, n
                vdir = matmul(pos(j, :) - pos(i, :), ibox)
                dimensions = max(dimensions, vdir)  ! get the larger length
            end do
        end do
    end subroutine atoms_dimensions


end module kdbinsert
