module kdb
    use fhelper
    implicit none

    contains
 

    subroutine atom_atom_distance(v, dist)
        double precision, intent(in) :: v(3)
        double precision             :: dist
!f2py   intent(in)  :: v
!f2py   intent(out) :: dist
        dist = norm2(v)
    end subroutine atom_atom_distance


   subroutine atom_atom_perm_distance(v, box, ibox, d)
       !! Returns array of distances of vector v permuted in the 8 most-adjacent PBC cells
       !!
       !! Parameters:
       !!      v:     input vector (rank 1 array of length 3)
       !!      box:   cell that defines the boundary conditions (3x3 array)
       !!      ibox:  inverse of box

       ! inputs
       double precision, intent(in) :: v(3), box(3, 3), ibox(3, 3)
       ! output
       double precision, dimension(0:7) :: d

!f2py   intent(in)  :: v, box, ibox, ortho
!f2py   intent(out) :: d

       ! other variables
       double precision :: vpbc(3), vdir(3)
       integer          :: i, j, k, sgn(3)

       ! pbc
       vdir = matmul(v, ibox)
       do i = 1, 3
           vdir(i) = modulo(vdir(i), 1.d0) + 1.5d0
           vdir(i) = modulo(vdir(i), 1.d0) - 0.5d0
       end do
       vpbc = matmul(vdir, box)

       sgn = -1 * sign([1.d0, 1.d0, 1.d0], vdir)

       ! permuting vdir in 8 neighboring cells
       do i = 0, 1
           do j = 0, 1
               do k = 0, 1
                   d(4*i + 2*j + k) = norm2(vpbc + i*sgn(1)*box(1,:) + j*sgn(2)*box(2,:) + k*sgn(3)*box(3,:))
               end do
           end do
       end do
    
   end subroutine atom_atom_perm_distance


    subroutine atom_atom_pseudo_mic_distance(v, box, ibox, dmic)
        !! Returns the distance^2 according to the Minimum Image Convention (almost)
        !!
        !! For unit cells without any non-orthogonal angles, result is identical
        !! to PBC distance^2
        !!
        !! Parameters:
        !!      v:     input vector (rank 1 array of length 3)
        !!      box:   cell that defines the boundary conditions (3x3 array)
        !!      ibox:  inverse of box

        ! inputs
        double precision, intent(in) :: v(3), box(3, 3), ibox(3, 3)
        ! output
        double precision :: dmic

!f2py   intent(in)  :: v, box, ibox, ortho
!f2py   intent(out) :: dmic

        ! other variables
        double precision :: vdir(3), vnew(3), dnew
        integer          :: i, j, k

        ! pbc
        vdir = matmul(v, ibox)
        do i = 1, 3
            vdir(i) = modulo(vdir(i), 1.d0) + 1.5d0
            vdir(i) = modulo(vdir(i), 1.d0) - 0.5d0
        end do
        dmic = 1.d10
        ! permuting vdir in 8 neighboring cells (only permute if ortho(i) is 0)
        do i = 0, 1
            vdir(1) = sign(1.d0, vdir(1)) * (abs(vdir(1)) - 1.d0)
            do j = 0, 1
                vdir(2) = sign(1.d0, vdir(2)) * (abs(vdir(2)) - 1.d0)
                do k = 0, 1
                    vdir(3) = sign(1.d0, vdir(3)) * (abs(vdir(3)) - 1.d0)
                    vnew = matmul(vdir, box)
                    dnew = sum(vnew * vnew)
                    if (dnew < dmic) then
                        dmic = dnew
                    end if
                end do
            end do
        end do
        dmic = sqrt(dmic)
       
    end subroutine atom_atom_pseudo_mic_distance


    function pbc(r, box, ibox) result(v)
        !! Applies periodic boundary conditions.
        !! Parameters:
        !!      r:    the vector the boundary conditions are applied to (rank 1 array of length 3)
        !!      box:  the box that defines the boundary conditions (3x3 array)
        !!      ibox: the inverse of the box (required because f2py's optional argument feature is messed up)

        ! inputs
        double precision, intent(in) :: r(3), box(3, 3), ibox(3, 3)
        ! output
        double precision :: v(3)
        ! other variables
        integer :: i

        v = matmul(r, ibox)
        do i = 1, 3
            v(i) = modulo(v(i), 1.d0) + 1.5d0
            v(i) = modulo(v(i), 1.d0) - 0.5d0
        end do
        v = matmul(v, box)

    end function pbc


    function pbcs(r, box, ibox, n) result(v)
        !! Applies periodic boundary conditions.
        !! Parameters:
        !!      r:    the row vector(s) the boundary conditions are applied to (nx3 array)
        !!      box:  the box that defines the boundary conditions (3x3 array)
        !!      ibox: the inverse of the box (required because f2py's optional argument feature is messed up)

        ! inputs
        integer, intent(in)          :: n
        double precision, intent(in) :: r(n, 3), box(3, 3), ibox(3, 3)
        ! output
        double precision :: v(n, 3)
        ! other variables
        integer :: i, j

        v = matmul(r, ibox)
        do j = 1, 3
            do i = 1, n
                v(i, j) = modulo(v(i, j), 1.d0) + 1.5d0
                v(i, j) = modulo(v(i, j), 1.d0) - 0.5d0
            end do
        end do
        v = matmul(v, box)

    end function pbcs

    
    !!========================================================================================!!
    !!                                      Clump Methods                                     !!
    !!========================================================================================!!

    subroutine find_clump_order(pos, n, clump_order)
        !! Generates clump order for all atoms in pos, based on closest pairs (non-PBC).
        !! Assumes that atoms already had PBC removed using kdb.blind_clump(). Intended to be
        !! used with kdb.ordered_clump()
        !! Parameters:
        !!      pos: array (n x 3) of atoms positions
        !! Returns:
        !!      clump_order: integer list (length 2*n) of the following format:
        !!                   [source_1, target_1, source_2, target_2, ...]

        ! inputs
        integer, intent(in)                                 :: n
        double precision, intent(in), dimension(0:n-1, 0:2) :: pos  ! indices start at 0
        ! output
        double precision :: clump_order(2*n - 2)

!f2py   intent(in)  :: pos
!f2py   intent(out) :: clump_order

        ! other variables
        integer, dimension(0:n-1) :: done, atoms
        integer                   :: length_done = 1, i, j, min_i = 0, min_j = 0, temp, index
        double precision          :: mindist2, dist2  ! the "2" means "squared"
        double precision          :: v(3)
        integer, parameter        :: null_val = -1  ! value to set when atom is "removed" from array
        
        atoms = [(i, i=0,n-1, 1)]  ! [0, 1, 2, ..., n-2, n-1]
        done(0) = atoms(0)
        atoms(0) = null_val  ! nullified
        
        do length_done = 1, n-1  ! number of atoms done at beginning of loop
            mindist2 = 1.d100
            do i = 1, n-1  ! remaining elements in atoms
                if (atoms(i) > null_val) then
                    do j = 0, length_done-1  ! goes through atoms whose order has been determined
                        v = pos( atoms(i) , : ) - pos( done(j) , : )  
                        dist2 = sum(v*v)
                        if (dist2 < mindist2) then
                            mindist2 = dist2
                            min_i = i
                            min_j = j
                        end if
                    end do
                end if
            end do

            index = 2*length_done
            temp = atoms(min_i)
            clump_order(index - 1) = done(min_j)
            clump_order(index) = temp
            done(length_done) = temp
            atoms(min_i) = null_val  ! nullified

        end do 

    end subroutine find_clump_order


    subroutine ordered_clump(pos, clump_order, box, ibox, n_pos, n_order, clumped_pos)
        !! Clump performed with an order determined beforehand ( from kdb.find_clump_order() )
        !! To be used on reactant during query (assumes PBC still applies)
        !! Parameters:
        !!      pos        : positions (n_pos x 3 array) of atoms, whose subset will be clumped
        !!      clump_order: order to clump (length n_order). Format should be identical to
        !!                   return list of kdb.find_clump_order()
        !!      box        : cell from which the atoms came from
        !!      ibox       : inverse of box (required because f2py's optional argument feature
        !!                   is messed up)
        !! Returns:
        !!      clumped_pos: clumped positions of pos (n_pos x 3 array)

        ! inputs
        integer, intent(in)                                     :: n_pos, n_order
        double precision, intent(in), dimension(0:n_pos-1, 0:2) :: pos  ! indices start at 0
        integer, intent(in)                                     :: clump_order(n_order)
        double precision, intent(in)                            :: box(3, 3), ibox(3, 3)
        ! output
        double precision, dimension(0:n_pos-1, 0:2) :: clumped_pos

!f2py   intent(in)  :: pos, clump_order, box, ibox
!f2py   intent(out) :: clumped_pos

        ! other variables
        integer          :: i, source, successor
        double precision :: v(3)

        clumped_pos = pos

        do i = 1, n_order, 2
            source = clump_order(i)
            successor = clump_order(i + 1)
            v = clumped_pos(successor, :) - clumped_pos(source, :)
            v = pbc(v, box, ibox)
            clumped_pos(successor, :) = clumped_pos(source, :) + v
        end do

    end subroutine ordered_clump


    subroutine blind_clump(pos, atoms, box, ibox, n, p, clumped_pos)
        !! Clump performed without an order determined beforehand.
        !! To be used to remove PBC from unstripped (or stripped) atoms before kdbinsert.
        !! Parameters:
        !!      pos  : *positions* (n x 3 array) of ASE's Atoms object (NOT the object itself)
        !!      atoms: list of atoms to clump with length p (elements are all non-negative
        !!             integers)
        !!      box  : Atoms object's cell (3 x 3 array)
        !!      ibox : inverse of box (required because f2py's optional argument feature is
        !!             messed up)

        ! inputs
        integer, intent(in)                                 :: n, p
        double precision, intent(in), dimension(0:n-1, 0:2) :: pos  ! indices start at 0
        double precision, intent(in)                        :: box(3, 3), ibox(3, 3)
        integer                                             :: atoms(p)
        ! output
        double precision, dimension(0:n-1, 0:2) :: clumped_pos

!f2py   intent(in)  :: pos, box, atoms, ibox
!f2py   intent(out) :: clumped_pos

        ! other variables
        integer            :: done(p)
        integer            :: length_done = 1, i, j, min_i = 0, min_j = 0
        double precision   :: mindist2, dist2  ! the "2" means "squared"
        double precision   :: v(3), vmin(3)
        integer, parameter :: null_val = -1  ! value to set when atom is "removed" from array
        
        clumped_pos = pos
        done(1) = atoms(1)
        atoms(1) = null_val  ! nullified
        
        do length_done = 1, p-1  ! number of atoms done at beginning of loop
            mindist2 = 1.d100
            do i = 2, p  ! remaining elements in atoms
                if (atoms(i) > null_val) then
                    do j = 1, length_done
                        ! SJ: better to do column-wise, but found no good way to interface it with python
                        v = clumped_pos( atoms(i) , : ) - clumped_pos( done(j) , : )  
                        v = pbc(v , box, ibox)
                        dist2 = sum(v*v)
                        if (dist2 < mindist2) then
                            vmin = v
                            mindist2 = dist2
                            min_i = i
                            min_j = j
                        end if
                    end do
                end if
            end do

            clumped_pos( atoms(min_i) , : ) = clumped_pos( done(min_j) , : ) + vmin
            done(length_done + 1) = atoms(min_i)
            atoms(min_i) = null_val  ! nullified
        end do

    end subroutine blind_clump

end module kdb
