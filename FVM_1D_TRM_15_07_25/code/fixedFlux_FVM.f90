! Author : Akash Asokan
! Description : 1D FVm Transient Thermal Reponse Model
! Date :  15.07.2025 
! File Name : fvm_v1

program test
    implicit none
    real(8) :: k , rho, cp , alpha, T0, Ts, Fo, hf, T_b
    real(8) :: dt, dx
    real(8) :: L

    real(8) :: c1, c2, c3, h1, h2, h3, pi
    
    integer :: n, nt, i, it, time
    
    real(8), dimension(:), allocatable :: T, x, T_old, T_analytical

    k = 116.d0
    rho = 1850.d0
    cp = 719.d0
    L = 1.0d0
    time = 110
    
    ! Number of cells required
    n = 100

    T0 = 300
    Ts = 3500

    pi = 4 * atan(1.0d0)

    dt = 0.001d0

    ! Cell Size of each Control Volume assuming A ~ 1
    dx = L / (n)
    hf = 2000*10**4.d0

    nt = int(time/dt)
    alpha = k / (rho * cp)
    Fo = dt * alpha / dx**2.d0
    
    allocate(x(0:n-1), T(0:n-1), T_old(0:n-1), T_analytical(0:n-1))

    ! Intialisation
    T_old(0:n-1) = T0

    ! Time Marching
    do it=1, nt

        ! Updating Left Boundary Cell
        T(0) = T_old(0) + (Fo * (T_old(1) - T_old(0) + (hf * dx / k) ) ) 

        ! Internal Cells
        do i=1, n-2
            T(i) = T_old(i) + (Fo *( T_old(i-1) - (2.d0*T_old(i)) + T_old(i+1) ) )
        end do      

        ! Updating Right Boundary Cell
        T(n-1) = T_old(n-1) + (Fo * (T_old(n-2) -T_old(n-1)))

        ! Updating T_old array
        T_old = T
    end do 

    ! Interpolating face centre value from bondary cell
    T_b = T(0) + (hf * dx / 2 /k )

    ! Temperature plotting along the length of the rod
    open(unit=333,file='numerical.dat', status='unknown')
    write(333,*) 0.d0, T_b
    
    do i=0, n-1
        x(i) = ( (2*i) + 1 ) *dx/2
        write(333,*) x(i), T(i)
    end do
    
    write(333,*) L, T(n-1)
    close(333)
    
    ! Analytical Solution
    open(unit=222, file='analytical.dat', status='unknown')
    do i=0, n-1
        x(i) = i * dx
        h1 = 2.d0 * hf * ((alpha *time) / pi)**0.5 / k 
        h2 = exp( -x(i)**2 / (4.d0 * alpha * time) ) 
        h3 = ( hf * x(i) / k) * erfc(x(i) / (2.d0 * sqrt(alpha*time)) )
        
        T_analytical(i) = T0 + ( h1 * h2 ) - h3  
        write(222,*) x(i) , T_analytical(i)
    end do
    close(222)
end program test