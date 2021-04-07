subroutine frobenius(A,B,M,eps,equal)
        ! Computes the Frobenius norm of the difference 
        ! of two matrices A and B of size M*M 
        ! Returns the boolean variable 'equal' accounting for
        ! matrices equality according to the given threshold 'eps'

        implicit none
        integer :: M, i, j
        double precision :: eps, norm
        double precision, dimension(M,M) :: A, B
        logical :: equal

        !#### COMPUTING THE FROBENIUS NORM OF A-B ####
        norm = 0.d0
        do i = 1, M
                do j = 1, M
                        norm = norm + abs(A(i,j) - B(i,j))**2
                end do
        end do
        norm = sqrt(norm)

        !#### COMPARING THE FROBENIUS NORM OF A-B TO THE THRESHOLD ####
        equal = norm .lt. eps

end subroutine frobenius
