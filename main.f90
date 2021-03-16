! DEADLINE = 4 weeks from 10/03/21

program helium
        !### general settings ###
        implicit none
        !########################
        !### variables ###
        double precision, dimension(:,:), allocatable :: &
                S, T, Sbar, Sm12, H, F, C, Ds, R, G, Fprime
        double precision, dimension(:), allocatable :: alpha, d, e, c_1
        double precision :: pqrs, tmp, e1, eps = 10e-8
        double precision, dimension(2) :: Et
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
                H(M,M), c_1(M), F(M,M), C(M,M), Ds(M,M), R(M,M), G(M,M), Fprime(M,M))
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
        !### diagonalization of S ###
        T = S
        call tred2(T,M,M,d,e)
        call tqli(d,e,M,M,T)
        
        !############################
        !### computing Sbar**-1/2 ###
        do i = 1, M
                Sbar(i,i) = d(i)
        end do
        do i = 1, M
                Sbar(i,i) = 1 / sqrt(Sbar(i,i))
        end do
        !############################
        !### computing S**-1/2 ###
        Sm12 = matmul(matmul(T,Sbar),transpose(T))
        !#########################
        !### building the H matrix ###
        do i = 1, M
                do j = 1, M
                        H(i,j) = 4 * (alpha(i) * alpha(j) - Z * (alpha(i) + alpha(j))) * &
                        (sqrt(alpha(i) * alpha(j)) / (alpha(i) + alpha(j))) ** 3
                end do
        end do
        !##################################################

! Iterative SCF process
        
        !### initialization of the LCAO coeffs of phi1 ###
        do i = 1, M
                if (i .eq. 1) then
                        c_1(i) = 1.d0
                else 
                        c_1(i) = 0.d0
                end if
        end do
        !#################################################
        

        write(*,'(1X,11A)') 'Iteration',char(9),'c1',char(9),&
                'c2',char(9),'F11',char(9),'F22',char(9),'F22',char(9),'e',char(9),'E'

        !### iterating SCF process ###
        condition = .true.
        do while (condition .eqv. .true.)
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
                !### computation of the total energy Et ###
                Et(2) = 2.d0 * dot_product(matmul(c_1,H),c_1) + dot_product(matmul(c_1,G),c_1)
                !##########################################
                !### proceding to the Lowdin orthogonalization for F ###
                Fprime = matmul(matmul(Sm12,F),Sm12) ! it's actually F'
                d = 0.d0
                e = 0.d0
                call tred2(Fprime,M,M,d,e)
                call tqli(d,e,M,M,Fprime)
                C = matmul(Sm12,Fprime)
                Fprime = 0.d0
                do i = 1, M
                        Fprime(i,i) = d(i)
                end do
                !#######################################################
                !### Storing lowest energy eigenvalue and associated eigenvector ###
                e1 = minval(d)
                c_1 = C(:,minloc(d,DIM=1))
                !###################################################################
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
                !############################




                it = it + 1
                write(*,'(I5.1,7F10.6)') it, c_1(1), c_1(2), F(1,1), F(2,2), F(1,2), e1, Et(2)
                condition = (abs(Et(1) - Et(2)) .gt. eps) .and. (nint(tmp) .eq. Z)
                Et(1) = Et(2)
        end do
        !#############################


        deallocate(alpha, S, T, Sbar, Sm12, d, e, H, c_1, F, C, Ds, R, G, Fprime)
end program helium
