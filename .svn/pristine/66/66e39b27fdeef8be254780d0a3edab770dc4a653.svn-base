! helper functions for KDB's fortran code

module fhelper
    implicit none

    contains


    pure function matinv3(A) result(B)
        !! https://fortranwiki.org/fortran/show/Matrix+inversion
        !! Performs a direct calculation of the inverse of a 3×3 matrix.
        double precision, intent(in) :: A(3,3)   !! Matrix
        double precision             :: B(3,3)   !! Inverse matrix
        double precision             :: detinv
    
        ! Calculate the inverse determinant of the matrix
        detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                  - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                  + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
    
        ! Calculate the inverse of the matrix
        B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
        B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
        B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
        B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
        B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
        B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    end function matinv3


end module fhelper
