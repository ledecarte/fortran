! Author : Akash Asokan
! Description : 2D FVM Transient Thermal Reponse Model
! Date :  21.07.2025 
! File Name : fvm_v3

program fvm_v3
    implicit none
    
    ! Variable declaration
    real(8) :: k_xx, k_yy, cp, rho, alpha_xx, alpha_yy
    real(8) :: dt, dx_e, dx_w, dy_n, dy_s, dx, dy
    real(8), allocatable :: T(:,:), T_old(:,:), x(:,:)
    real(8) :: T0, Tb, htc, qf
    real(8) :: a_W, a_E, a_N, a_S, a_P0, a_P
    real(8) :: Lx, Ly
    integer :: nx, ny, nt, time, unit
    integer :: i, j, it 
    
    ! Material Properties
    k_xx    = 116.d0    ! Through Plane Conductivity
    k_yy    = 116.d0    ! In-plane Conductivity
    cp      = 719.d0    ! Temperature dependent in later implementation    
    rho     = 1850.d0   ! Density in kg / m3

        ! Secondary Compuation
        alpha_xx   = k_xx / (rho * cp)      ! Computed Thermal Diffusivity Through Plane
        alpha_yy   = k_yy / (rho * cp)      ! Computed Thermal Diffusivity In-Plane

    ! 2D Grid Generation
    Lx      = 1.d0      ! Length of Domain in X Direction
    Ly      = 1.d0      ! Length of Domain in Y Direction
    nx      = 100       ! Number of cells in X Directions
    ny      = 100       ! Number of cells in Y Directions

        ! Spacing Computations
        dx  = Lx / nx   ! Distance between two nodes along X Direction
        dy  = Ly / ny   ! Distance between two nodes along Y Direction

        ! Directional Spacing 
        dx_e    = dx    ! Uniform spacing along x direction
        dx_w    = dx    ! Uniform spacing along x direction
        dy_n    = dy    ! Uniform spacing along y direction
        dy_s    = dy    ! Uniform spacing along y direction
    
    ! Time Step and Number of Timesteps
    time    = 100.0     ! Total time in seconds
    dt      = 0.001d0   ! Time step and this value is directly related to the stabilty 
    nt      = int(time / dt) ! Number of time steps 
    
    ! Standard Boundary Conditions and values 
    T0      = 300.d0        ! K
    Tb      = 3500.d0       ! K
    htc     = 9000.d0       ! W/m2K
    qf      = 2000*10**4    ! W/m2

    ! Array Allcoation and Initialisation
    allocate(T(0:ny,0:nx), T_old(0:ny,0:nx), x(0:ny,0:nx))

    ! Initialisation 
        ! Internal Cells  
        T_old(1:ny-2, 1:nx-2) = T0
        
        ! Edge Cells
        T_old(0, 1:nx-2)    = 300.0      ! North Edge Cell
        T_old(ny-1, 1:nx-2) = 300.0      ! South Edge Cell
        T_old(1:ny-2, 0)    = 300.0      ! West Edge Cell
        T_old(1:ny-2, nx-1) = 300.0      ! East Edge Cell

        ! Corner Cells 
        T_old(0,0)          = 300.0      ! Top West Corner
        T_old(ny-1,0)       = 300.0      ! Bottom West Corner
        T_old(0,nx-1)       = 300.0      ! Top East Corner
        T_old(ny-1,nx-1)    = 300.0      ! Bottom East Corner

    ! Coefficent Calculations for aE, aW, aS, aN and aP
    a_E     = k_xx * dy / dx_e
    a_W     = k_xx * dy / dx_w
    a_N     = k_yy * dx / dy_n
    a_S     = k_yy * dx / dy_s
    a_P0    = rho * cp * dx * dy / dt
    a_P     = a_P0

    ! Time Marching
    do it =1, nt 
        ! Looping inside internal cells
        do i = 1, ny-2
            do j =1, nx-2
                ! Internal Cells
            T(i,j) = (  (a_E * T_old(i, j+1)) + &
                        (a_W * T_old(i, j-1)) + &
                        (a_N * T_old(i-1, j)) + &
                        (a_S * T_old(i+1, j)) + &
                        ( (a_P0 - a_E - a_W - a_N - a_S) * T_old(i,j) ) &
                    ) / a_P
            end do   
        end do

        ! Boundary Calculations
        ! Corner Boundary Cells
            ! Top West Corner
            T(0,0)      = ((a_E * T_old(0,1)) +   & 
                            (a_S * T_old(1,0)) +    &
                            ( (a_P0 - a_E - (2*a_W) - a_S) * T_old(0,0) ) + &
                            ( 2 * a_W * Tb)) / a_P
            ! Bottom West Corner
            T(ny-1,0)   = ((a_E * T_old(ny-1,1)) +   & 
                            (a_N * T_old(ny-2,0)) +    &
                            ( (a_P0 - a_E - (2*a_W) - a_N) * T_old(ny-1,0) ) + &
                            ( 2 * a_W * Tb)) / a_P

            T(0,nx-1)   = ((a_W * T_old(0,nx-2)) +   & 
                            (a_S * T_old(1,nx-1)) +    &
                            ( (a_P0 - a_W - a_S) * T_old(0,nx-1) ) ) / a_P
            ! Bottom Right Corner
            T(ny-1,nx-1)= ((a_W * T_old(ny-1,nx-2)) +   & 
                            (a_N * T_old(ny-2,nx-1)) +    &
                            ( (a_P0 - a_W - a_N) * T_old(ny-1,nx-1) ) ) / a_P
                            
        ! Edge Boundary Cells   
        ! West Side - Constant Temperature
        do i = 1, ny-2
            T(i, 0) = (( a_E * T_old(i, 1) ) + &
                        ( a_N * T_old(i-1, 0) ) + &
                        (a_S * T_old(i+1, 0)) + &
                        ((a_P0 - a_E - (2*a_W) - a_N - a_S ) * T_old(i, 0)) + &
                        (2 * a_W * Tb)) / a_P 
        end do 
        ! East Side - Insulated
        do i = 1, ny-2
            T(i, nx-1) = (( a_W * T_old(i,nx-2) ) + &
                            ( a_N * T_old(i-1, nx-1) ) + &
                            ( a_S * T_old(i+1, nx-1) ) + &
                            ((a_P0 - a_N - a_S - a_W) * T_old(i,nx-1)) ) / a_P                 
        end do
        ! North Side - Insulated
        do j = 1, nx-2
            T(0, j)    = (( a_E * T_old(0, j+1)) + & 
                            ( a_W * T_old(0, j-1)) + &
                            (a_S * T_old(1, j)) + & 
                            ((a_P0 - a_E -a_W -a_S) * T_old(0, j )) ) / a_P
        end do
        ! South Side - Insulated
        do j = 1, nx-2
            T(ny-1, j)  = ((a_E * T_old(ny-1, j+1) ) + &
                            (a_W * T_old(ny-1, j-1) ) + &
                            (a_N * T_old(ny -2, j)) + &
                            ((a_P0 - a_E -a_W -a_N) * T_old(ny-1, j)) ) / a_P
        end do
        ! Update T_old
        T_old = T 
    end do 
    
    ! Extracting data 
    unit = 111
    open(unit=unit, file='fvm.dat', status='unknown')

    ! Face Centre Value Interpolation
    
    ! Post Processing
    ! write(unit,*) 0, ((2*i)+1) * dy/2, 3500.0
    
    do i = 0, ny-1
        write(unit,*) 0, ((2*i)+1) * dy/2, 3500.0 
        do j = 0, nx-1
            write(unit, *) ( ((2*j)+1) * dx / 2), ( ((2*i)+1) * dy / 2), T(i,j)            
        end do
        write(unit,*) Lx, ((2*i)+1) * dy/2, 300.0   ! should be changed later 
        write(unit,*)
    end do

    close(unit)

end program fvm_v3