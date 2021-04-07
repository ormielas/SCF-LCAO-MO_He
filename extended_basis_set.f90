program extanded_basis_set
        !#### GENERAL ####
        implicit none
        double precision, dimension(:,:), allocatable :: &
                S, T, Sbar, Sm12, H, F, C, Ds, R, G, Fprime, tmpm1, tmpm2
        double precision, dimension(:), allocatable :: alpha, d, e, LCAO, LCAOnew
        double precision :: pqrs, tmp, e1, eps = 10e-8, IE
        double precision, dimension(2) :: Et
        integer :: M, i, j, Z = 2, p, q, k, it = 1
        logical :: condition = .true., equal
        character :: slater

        !#### OPERATOR'S CHOICE OF THE SIZE LIMIT OF THE BASIS SET ####
        write(*,*) ''
        write(*,*) 'Chose the basis set dimension M :'
        write(*,*) ''
        read(*,*) M
        write(*,*) ''

        !#### ALLOCATING MATRICES ####
        allocate(S(M,M), T(M,M), Sbar(M,M), Sm12(M,M), H(M,M), F(M,M), C(M,M), &
        Ds(M,M), R(M,M), G(M,M), Fprime(M,M), tmpm1(M,M), tmpm2(M,M))
        
        !#### ALLOCATING VECTORS ####
        allocate(alpha(M), LCAO(M), d(M), e(M), LCAOnew(M))

        !#### INITIALIZATION OF THE SLATER COEFFS ####
        do i = 1, M
                if (i .eq. 1) then
                        alpha(i) = 1.d0
                else 
                        alpha(i) = 1.2d0 * alpha(i-1)
                end if
        end do  
                
        !#### INITIALIZATION OF THE LCAO COEFFS ####
        do i = 1, M
                if (i .eq. 1) then
                        LCAO(i) = 1.d0
                else 
                        LCAO(i) = 0.d0
                end if
        end do  

        !#### BUILDING OVERLAP MATRIX S ####
        do i = 1, M
                do j = 1, M
                        if (i .eq. j) then
                                S(i,i) = 1
                        else
                                S(i,j) = (2 * sqrt(alpha(i) * alpha(j))/ (alpha(i) + alpha(j))) ** 3
                        end if
                end do
        end do

        !#### DIAGONALIZATION OF THE OVERLAP MATRIX S ####
        T = S
        d = 0.d0
        e = 0.d0
        call tred2(T,M,M,d,e)
        call tqli(d,e,M,M,T)       

        !#### COMPUTING Sbar**-1/2 AND S**-1/2 (Sm12) ####
        Sbar = 0.d0
        do i = 1, M
                Sbar(i,i) = d(i)
        end do
        do i = 1, M
                Sbar(i,i) = 1 / sqrt(Sbar(i,i))
        end do        
        Sm12 = matmul(matmul(T,Sbar),transpose(T))
        
        !#### BUILDING THE HAMILTONIAN MATRIX H ####
        do i = 1, M
                do j = 1, M
                        H(i,j) = 4.d0 * (alpha(i) * alpha(j) - Z * (alpha(i) + alpha(j))) * &
                        (sqrt(alpha(i) * alpha(j)) / (alpha(i) + alpha(j))) ** 3.d0
                end do
        end do

        !#### FORMATING OUTPUT TABLE ####
        write(*,'(A)') 'Iterative SCF process :'
        write(*,'(A)') '***********************'
        !write(*,*) ''
        write(*,'(9A10)') 'Basis size', 'Iteration', 'c1', 'c2', 'F11', 'F12', 'F22', 'e', 'E'

        !#### ITERATING SCF PROCESS ####
        it = 1
        F = 0.d0
        condition = .true.
        Et = 0.d0
        do while (condition .eqv. .true.)
        
                !#### COMPUTATION OF THE FOCK MATRIX F ####
                G = 0.d0
                do i = 1, M
                        do j = 1, M
                                tmp = 0
                                do p = 1, M
                                        do q = 1, M
                                                if (LCAO(p) * LCAO(q) .ne. 0) then
                                                        call pqrsint(pqrs,i,j,p,q,M,alpha)
                                                        tmp = tmp + LCAO(p) * LCAO(q) * pqrs
                                                end if
                                        end do
                                end do
                                G(i,j) = tmp
                        end do
                end do
                F = H + G
                
                !#### COMPUTATION OF THE TOTAL ENERGY Et ####
                Et(2) = 2.d0 * dot_product(matmul(LCAO,H),LCAO) + dot_product(matmul(LCAO,G),LCAO)
                
                !#### LOWDIN ORTHOGONALIZATION FOR THE FOCK MATRIX F ####
                Fprime = matmul(matmul(Sm12,F),Sm12)
                d = 0.d0
                e = 0.d0
                call tred2(Fprime,M,M,d,e)
                call tqli(d,e,M,M,Fprime)
                C = matmul(Sm12,Fprime)
                Fprime = 0.d0
                do i = 1, M
                        Fprime(i,i) = d(i)
                end do
                
                !#### STORING LOWEST ENERGY EIGENVALUE AND EIGENVECTOR ####
                e1 = minval(d)
                LCAOnew = C(:,minloc(d,DIM=1))
                
                !#### BUILDING THE DENSITY MATRIX R ####
                do i = 1, M
                        do j = 1, M
                                R(i,j) = LCAOnew(i) * LCAOnew(j)
                        end do
                end do
                
                !#### VERIFYING EQ. (8) ####
                tmp = 0.d0
                do i = 1, M
                        do j = 1, M
                                tmp = tmp + 2 * R(i,j)*S(i,j)
                        end do
                end do



                write(*,'(2I11.1,7F10.6)') M, it, LCAO(1), LCAO(2), F(1,1), F(1,2), F(2,2), e1, Et(2)
                condition = (abs(Et(1) - Et(2)) .gt. eps) .and. (nint(tmp) .eq. Z)
                Et(1) = Et(2)
                LCAO = LCAOnew
                it = it + 1
        end do
        write(*,*) ''





        !do i = 1, M
        !write(*,*) (F(i,j), j=1,M)
        !end do
        !write(*,*) ''
        !write(*,*) alpha
        !write(*,*)
        !#### DEALLOCATING TABLES ####
        deallocate(alpha, S, T, Sbar, Sm12, d, e, H, LCAO, F, C, Ds, R, G, Fprime, LCAOnew, tmpm1, tmpm2)
        
        !end do
        
        !write(*,*) e1
        !write(*,*) ''
        !write(*,*) Et(2)
        !write(*,*)
end program
