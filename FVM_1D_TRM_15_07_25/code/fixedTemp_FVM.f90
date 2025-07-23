! Author : Akash Asokan
! Description : 1D FVm Transient Thermal Reponse Model
! Date :  15.07.2025 
! File Name : fvm_v1

program test
    implicit none
    real(8) :: k , rho, cp , alpha, T0, Ts, Fo
    real(8) :: dt, dx
    real(8) :: L
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

    dt = 0.001d0

    ! Cell Size of each Control Volume assuming A ~ 1
    dx = L / (n)

    nt = int(time/dt)
    alpha = k / (rho * cp)
    Fo = dt * alpha / dx**2.d0
    
    allocate(x(0:n-1), T(0:n-1), T_old(0:n-1), T_analytical(0:n-1))

    ! Intialisation
    T_old(0:n-1) = T0

    ! Time Marching
    do it=1, nt

        ! Updating Left Boundary Cell
        T(0) = T_old(0) + (Fo * (T_old(1) - (3 * T_old(0)) + (2*Ts)) ) 

        ! Internal Cells
        do i=1, n-2
            T(i) = T_old(i) + (Fo *( T_old(i-1) - (2.d0*T_old(i)) + T_old(i+1) ) )
        end do      

        ! Updating Right Boundary Cell
        T(n-1) = T_old(n-1) + (Fo * (T_old(n-2) -T_old(n-1)))

        ! Updating T_old array
        T_old = T
    end do 

    ! Temperature plotting along the length of the rod
    open(unit=333,file='fvm.dat', status='unknown')
    write(333,*) 0.d0, Ts
    
    do i=0, n-1
        x(i) = ( (2*i) + 1 ) *dx/2
        write(333,*) x(i), T(i)
    end do
    
    write(333,*) L, T(n-1)
    close(333)
    
    open(unit=222, file='anal.dat', status='unknown')
    write(111,*) 0.d0, Ts
    do i=0, n-1
        x(i) = i * dx
        T_analytical(i) = Ts + ( (T0 - Ts) * erf( x(i) / (2.d0 * sqrt(alpha * time)) ) )
        write(222,*) x(i) , T_analytical(i)
    end do
    close(111)

end program test