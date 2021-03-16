! DEADLINE = 4 weeks from 10/03/21

program helium
        !### general settings ###
        implicit none
        !########################
        !### variables ###
        double precision, dimension(:,:), allocatable :: &
                S, T, Sbar, Sm12, H, F, C, Ds, R, G
        double precision, dimension(:), allocatable :: alpha, d, e, c_1
        double precision :: pqrs, tmp, e1, Et, eps = 10e-8
        integer :: M, i, j, Z = 2, p, q, k, it = 1
        logical :: condition
        !#################

! 1. Read of the data and construction of the integral matrices

        !### reading input parameters ###
        ! M is the dimension of the basis set
        ! alpha is the matrix of Slater coefficients
        open(unit=1, file="input_params.txt", status="old")
        read(1,*) M
        allocate(alpha(M), S(M,M), T(M,M), Sbar(M,M), Sm12(M,M), d(M), e(M), &
                H(M,M), c_1(M), F(M,M), C(M,M), Ds(M,M), R(M,M), G(M,M))
        read(1,*) alpha
        close(1)
        !################################

        !### building overlap matrix S ###
        do i = 1, M
                do j = 1, M
                        if (i .eq. j) then
                                S(i,i) = 1
                        else
                                S(i,j) = (2 * sqrt(alpha(i) * alpha(j))/ (alpha(i) + alpha(j))) ** 3
                        end if
                end do
        end do
        !#################################
        !### printing S matrix ###
        !write(*,*) ''
        !!write(*,*) 'The overlap matrix S is :'
        !write(*,*) ''
        !do i = 1, M
        !        write(*,'(10F10.4)') (S(i,j), j=1,M)
        !end do
        !write(*,*) ''
        !#########################
        !### diagonalization of S ###
        T = S
        call tred2(T,M,M,d,e)
        call tqli(d,e,M,M,T)
        
        !############################
        !### printing the matrix of eigenvectors T ###
        !write(*,*) 'The matrix of eigenvectors T is :'
        !write(*,*) ''
        !do i = 1, M
        !        write(*,'(10F10.4)') (T(i,j), j=1,M)
        !end do
        !write(*,*) ''
        !#############################################
        !### printing the diagonalized matrix Sbar ###
        do i = 1, M
                Sbar(i,i) = d(i)
        end do
        !write(*,*) 'The diagonalized matrix Sbar is :'
        !write(*,*) ''
        !do i = 1, M
        !        write(*,'(10F10.4)') (Sbar(i,j), j=1,M)
        !end do
        !write(*,*) ''
        !#############################################
        !### computing Sbar**-1/2 ###
        do i = 1, M
                Sbar(i,i) = 1 / sqrt(Sbar(i,i))
        end do
        !############################
        !### computing and printing S**-1/2 ###
        Sm12 = matmul(matmul(T,Sbar),transpose(T))
        !write(*,*) 'The S**-1/2 matrix is :'
        !write(*,*) ''
        !do i = 1, M
        !        write(*,'(10F10.4)') (Sm12(i,j), j=1,M)
        !end do
        !write(*,*) ''
        !#########################
        !### computing the hpq elements of the H matrix ###
        do i = 1, M
                do j = 1, M
                        H(i,j) = 4 * (alpha(i) * alpha(j) - Z * (alpha(i) + alpha(j))) * &
                        (sqrt(alpha(i) * alpha(j)) / (alpha(i) + alpha(j))) ** 3
                end do
        end do
        !write(*,*) 'The H matrix is :'
        !write(*,*) ''
        !do i = 1, M
        !        write(*,'(10F10.6)') (H(i,j), j=1,M)
        !end do
        !write(*,*) ''
        !##################################################
        !### computing the pqrs integrals ###
        ! call pqrsint(pqrs,p,q,r,s,M,alpha)
        ! it returns the analytical value of the pqrs integral in the pqrs variable
        !####################################

! Iterative SCF process
        
        !### initialization of the LCAO coeffs of phi1 ###
        ! use a matrix of size Mx(number of AO) (?) to allow for expansion over more AO ? (rather than a vectorial notation which
        ! allows for expansion over a single AO
        do i = 1, M
                if (i .eq. 1) then
                        c_1(i) = 1.d0
                else 
                        c_1(i) = 0.d0
                end if
        end do
        !### printing of the LCAO coeffs of phi1 ###
        !write(*,*) 'The LCAO coeffs of phi1 are'
        !write(*,*) ''
        !write(*,*) c_1
        !write(*,*) ''
        !###########################################
        !#################################################
        !### computation of the Fock matrix F ###
        G = 0.d0
        do i = 1, M
                do j = 1, M
                        tmp = 0
                        do p = 1, M
                                do q = 1, M
                                        if (c_1(p) * c_1(q) .ne. 0) then
                                                call pqrsint(pqrs,i,j,p,q,M,alpha)
                                                tmp = tmp + c_1(p) * c_1(q) * pqrs
                                        end if
                                end do
                        end do
                        G(i,j) = tmp
                end do
        end do
        F = H + G
        !########################################
        !### printing of the G matrix ###
        !write(*,*) 'The G matrix is'
        !write(*,*) ''
        !do i = 1, M
        !        write(*,'(10F10.6)') (G(i,j), j=1,M)
        !end do
        !write(*,*) ''
        !################################
        !### printing of the Fock matrix F ###
        !write(*,*) 'The Fock matrix F is :'
        !write(*,*) ''
        !do i = 1, M
        !        write(*,'(10F10.6)') (F(i,j), j=1,M)
        !end do
        !write(*,*) ''
        !#####################################
        !### computation of the total energy Et ###
        Et = 0.d0
        Et = 2.d0 * dot_product(matmul(c_1,H),c_1) + dot_product(matmul(c_1,G),c_1)
        !##########################################
        !### printing of the total energy Et ###
        write(*,*) 'The total energy Et in a.u. is'
        write(*,*) ''
        write(*,*) Et
        write(*,*) ''
        !#######################################
        !### proceding to the Lowdin orthogonalization for F ###
        F = matmul(matmul(Sm12,F),Sm12) ! it's actually F'
        d = 0.d0
        e = 0.d0
        call tred2(F,M,M,d,e)
        call tqli(d,e,M,M,F)
        C = matmul(Sm12,F)
        !write(*,*) 'The matrix C of the Fock matrix eigenvectors is :'
        !write(*,*) ''
        !do i = 1, M
        !       write(*,'(10F10.6)') (C(i,j), j=1,M)
        !end do
        !write(*,*) ''
        F = 0.d0
        do i = 1, M
                F(i,i) = d(i)
        end do
        !write(*,*) 'The matrix of the Fock matrix eigenvalues is :'
        !write(*,*) ''
        !do i = 1, M
        !        write(*,'(10F10.6)') (F(i,j), j=1,M)
        !end do
        !write(*,*) ''
        !#######################################################
        !### Storing lowest energy eigenvalue as e1 ###
        e1 = minval(d)
        !##############################################
        !### Storing the e1 associated eigenvectors as the new c_1 ###
        c_1 = C(:,minloc(d,DIM=1))
        !write(*,*) 'The new LCAO coeffs for phi1 are'
        !write(*,*) ''
        !write(*,*) c_1
        !write(*,*)
        !#############################################################
        !### Building the density matrix R ###
        do i = 1, M
                do j = 1, M
                        R(i,j) = c_1(i) * c_1(j)
                end do
        end do
        !#####################################
        !### Verifying equation 8 ###
        tmp = 0.d0
        do i = 1, M
                do j = 1, M
                        tmp = tmp + 2 * R(i,j)*S(i,j)
                end do
        end do
        !write(*,*) 'Equation 8 is :', nint(tmp) .eq. Z
        !############################






        deallocate(alpha, S, T, Sbar, Sm12, d, e, H, c_1, F, C, Ds, R, G)
end program helium
