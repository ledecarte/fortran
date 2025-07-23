program test
    implicit none
    real(8) :: k , rho, cp , alpha, T0, Ts
    real(8) :: dt, dx
    real(8) :: L
    integer :: n, nt, i, it, time
    
    real(8), dimension(:), allocatable :: T, x, T_old, T_analytical

    k = 116.d0
    rho = 1850.d0
    cp = 719.d0
    L = 1.0d0
    time = 110
    n = 100
    T0 = 300
    Ts = 3500

    dt = 0.001d0
    dx = L / (n-1)
    nt = int(time/dt)
    alpha = k / (rho * cp)
    
    allocate(x(0:n-1), T(0:n-1), T_old(0:n-1), T_analytical(0:n-1))

    ! Intialisation
    T_old(1:n-1) = T0
    T_old(0) = Ts
 
    ! Adding an openacc if parallel runs required 
    do it=1, nt

        do i=1, n-2
            T(i) = T_old(i) + (alpha * (dt/dx**2.d0) *( T_old(i-1) - (2.d0*T_old(i)) + T_old(i+1) ) )
        end do 
        
        T(0) = Ts
        T(n-1) = T(n-2)

        T_old = T
    end do 

    
    open(unit=222, file='anal.dat', status='unknown')
    open(unit=111, file='num.dat', status='unknown')

    do i=0, n-1
        x(i) = i * dx
        T_analytical(i) = Ts + ( (T0 - Ts) * erf( x(i) / (2.d0 * sqrt(alpha * time)) ) )
        write(222,*) x(i) , T_analytical(i)
        write(111,*) x(i) , T(i)
    end do 

end program test
