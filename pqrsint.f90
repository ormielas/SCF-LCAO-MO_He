subroutine pqrsint(pqrs,p,q,r,s,M,alpha)
        implicit none
        integer :: p, q, r, s, M
        double precision :: pqrs, a, b, c, d, sumalpha, prodalpha
        double precision, dimension(M) :: alpha

        ! Calculation for the (pq|rs) integral, given p, q, r and s
        
        sumalpha = alpha(p) + alpha(q) + alpha(r) + alpha(s)
        prodalpha = alpha(p) * alpha(q) * alpha(r) * alpha(s)
        a = 32 * (prodalpha ** (3.d0/2.d0)) 
        b = 1 / (((alpha(p) + alpha(q)) ** 3) * ((alpha(r) + alpha(s)) ** 2))
        c = 1 / (((alpha(p) + alpha(q)) ** 3) * (sumalpha ** 2))
        d = 1 / (((alpha(p) + alpha(q)) ** 2) * (sumalpha ** 3))
        pqrs = a * (b - c - d)



end subroutine pqrsint
