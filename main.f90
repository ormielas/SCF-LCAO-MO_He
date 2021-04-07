! DEADLINE = 4 weeks from 10/03/21

program helium
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
        write(*,*) 'Note that for M > 2, Slater coefficients are built recursively'
        !following : alpha(n+1) = 1.2 * alpha(n), alpha(1)
        ! = 1.d0'
        write(*,*) ''
        read(*,*) M
        write(*,*) ''

        if (M .eq. 2) then
                !#### OPERATOR'S CHOICE OF SLATER COEFFS ####
                write(*,*) ''
                write(*,*) 'Chose the set of Slater coefficients for the following SCF calulation :'
                write(*,*) '(S)tandard Slater coefficients, (O)ptimized Slater coeffients, (R)ecursive Slater coefficients'
                write(*,*) ''
                read(*,*) slater
                write(*,*) ''
        
                !#### READING INPUT PARAMS ####
                if ((slater .eq. 'S') .or. (slater .eq. 's')) then        
                        open(unit=1, file="input_standard.txt", status="old")
                        read(1,*) M
                        allocate(alpha(M))
                        read(1,*) alpha
                        close(1)
                else if ((slater .eq. 'O') .or. (slater .eq. 'o')) then
                        open(unit=2, file="input_optimized.txt", status="old")
                        read(2,*) M
                        allocate(alpha(M))
                        read(2,*) alpha
                        close(2)
                else if ((slater .eq. 'R') .or. (slater .eq. 'r')) then
                        open(unit=3, file="input_recursive.txt", status="old")
                        read(3,*) M
                        allocate(alpha(M))
                        read(3,*) alpha
                        close(3)
                end if  
        else
                !#### INITIALIZATION OF THE SLATER COEFFS ####
                allocate(alpha(M))
                do i = 1, M
                        if (i .eq. 1) then
                                alpha(i) = 1.d0
                        else 
                                alpha(i) = 1.2d0 * alpha(i-1)
                        end if
                end do  

        end if

        !#### ALLOCATING MATRICES ####
        allocate(S(M,M), T(M,M), Sbar(M,M), Sm12(M,M), H(M,M), F(M,M), C(M,M), &
        Ds(M,M), R(M,M), G(M,M), Fprime(M,M), tmpm1(M,M), tmpm2(M,M))
        
        !#### ALLOCATING VECTORS ####
        allocate(LCAO(M), d(M), e(M), LCAOnew(M))

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
        call tred2(T,M,M,d,e)
        call tqli(d,e,M,M,T)       
        
        !#### COMPUTING Sbar**-1/2 AND S**-1/2 (Sm12) ####
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
                        H(i,j) = 4 * (alpha(i) * alpha(j) - Z * (alpha(i) + alpha(j))) * &
                        (sqrt(alpha(i) * alpha(j)) / (alpha(i) + alpha(j))) ** 3
                end do
        end do
        
        !#### INITIALIZATION OF THE LCAO COEFFS ####
        do i = 1, M
                if (i .eq. 1) then
                        LCAO(i) = 1.d0
                else 
                        LCAO(i) = 0.d0
                end if
        end do  

        !#### FORMATING OUTPUT TABLE ####
        write(*,'(A)') 'Iterative SCF process :'
        write(*,'(A)') '***********************'
        write(*,*) ''
        write(*,'(9A10)') 'Basis size', 'Iteration', 'c1', 'c2', 'F11', 'F12', 'F22', 'e', 'E'

        !#### ITERATING SCF PROCESS ####
        open(unit=4, file="output.txt", status="new")    
        write(4,'(9A10)') 'Basis size', 'Iteration', 'c1', 'c2', 'F11', 'F12', 'F22', 'e', 'E'
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
                write(4,'(2I11.1,7F10.6)') M, it, LCAO(1), LCAO(2), F(1,1), F(1,2), F(2,2), e1, Et(2)
                condition = (abs(Et(1) - Et(2)) .gt. eps) .and. (nint(tmp) .eq. Z)
                Et(1) = Et(2)
                LCAO = LCAOnew
                it = it + 1
        end do
        close(4)
        
        !#### COMPUTING TOTAL ENERGY WITH EQ. (10) ####
        Et(2) = 2 * e1 - dot_product(matmul(LCAO,G),LCAO)
        
        !#### CHECKING IDEMPOTENCY OF DENSITY MATRIX R ####
        tmpm1 = matmul(matmul(R,S),R)
        call frobenius(R,tmpm1,M,0.0001d0,equal)
        write(*,*) ''
        write(*,'(A)') 'After SCF convergence :'
        write(*,'(A)') '***********************'
        write(*,*) ''
        write(*,'(A,L1)') 'The statement "The density matrix R is idempotent" is : ', equal
        write(*,*) ''
        !#### CHECKING THAT F AND R COMMUTE ####
        tmpm1 = matmul(matmul(F,R),S)
        tmpm2 = matmul(matmul(S,R),F)
        call frobenius(tmpm1,tmpm2,M,0.0001d0,equal)
        write(*,'(A,L1)') 'The statement "FRS = SRF" is : ', equal   
        write(*,*) ''

        !#### COMPUTING FIRST IONIZATION ENERGY IE ####
        IE = dot_product(matmul(LCAO,H),LCAO) - Et(1)
        write(*,'(A,1F10.6,A)') 'The first ionization energy IE is IE =', IE, ' au.'
        write(*,*) ''
        
        !#### COMPUTING FIRST IONIZATIN ENERGY IE ACCORDING TO KOOPMAN ####
        write(*,'(A,1F10.6,A)') "According to Koopman's theorem the first ionization energy is IE =", -e1, ' au.'
        write(*,*) ''

        !#### DEALLOCATING TABLES ####
        deallocate(alpha, S, T, Sbar, Sm12, d, e, H, LCAO, F, C, Ds, R, G, Fprime, LCAOnew, tmpm1, tmpm2)

end program helium
