! Initialize particle positions r(t=0) and velocities v(t=0)
! Input: n number of elements

MODULE init
  IMPLICIT NONE
  REAL ::  L, rho, NUC                                          ! l is the length of cube, rho is the 
  REAL, DIMENSION :: r(N,3), v(N,3)
  INTEGER :: N, a
                                                                ! Calculate cell dimension of cubical side l
  rho = 1                                                       ! rho is the system density and depends on simulation being run
  L = (N/rho)**(1/3)                                            ! Determining length of the cubic cell
  a = 1                                                         ! Particles per unit cell (simple cubic)
  NUC = CEIL((N/A)**(1/3))                                      ! Number of unit cells
  w = L/NUC                                                     ! Unit cell width
  
  CONTAINS
  SUBROUTINE rvinit(N,rho,T,sigma,eps,l)
    IMPLICIT NONE
    INTEGER :: i,j,k                                            ! Length of cell (L), x, y, z coordinates of
                                                                ! cell (i,j,k)
    REAL :: T, sigma, eps                                       ! Total number of particles (N), number of
                                                                ! desnity (rho), temperature (T), Lennard-Jones constants (sigma,
                                                                ! eps)
    REAL :: kB = 1                                              ! Boltzmann constant in atomic units
    REAL :: m = 1                                               ! Mass of particle in atomic units

  ! Assign positions then velocities for particles 1 through N
    DO i=0,L-1                                              ! x coordinate
     DO j=0,L-1                                             ! y coordinate
       DO k=0,L-1                                           ! z coordinate
         r(k+j*l+i*L**2+1,1) = a*k                          ! particle 1 located at 0,0,0
         r(k+j*l+i*L**2+1,2) = a*j                          ! particle 2 located at 1,0,0
         r(k+j*l+i*L**2+1,3) = a*i                          ! particle L+1 located ate 0,1,0 and so on.

         
       END DO
     END DO
    END DO
    
    ! Initialize velocity.
    CALL random_number(v)                              ! Assign random number from uniform
                                                       ! distribution
    v = 2*v-1                                          ! Normalize range from -1 to 1
    v = v * SQRT(3*N*kB*T/m*SUM(SUM(v*v,1))            ! Equipartition theorem
    v = v - SUM(v,1)/DOUBLE(N)                         ! subtract off COM velocity

  END SUBROUTINE rvinit
END MODULE init(N)

  SUBROUTINE rdf(r,N)                                           !The subroutine for creating histogram such as radial distribution function        
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NBINS = 100                             !We initialize the number of bins.
  INTEGER :: i,j,N                                              !We initialize two counter variables i,j and the number of particles, N.
  REAL, DIMENSION :: r(1:N,1:3), g(1:NBINS)                     !We initialize the arrays for distance frequency.
  REAL, DIMENSION :: ri(1:3), dr(1:3),w                         !We initialize the position and distance
  INTEGER, PARAMETER :: NBINS = 100                             ! Number of bins

  SAVE g,krdf                                                   !SAVE the g and krdf  
  IF (switch == 1) THEN                                         !if it is after first  iteration / time step
  krdf = krdf + 1
  ELSE IF( switch == 0) THEN                                    ! First time step switch = 0 
  DO i=1, N-1
          ri = r(i,1:DIMS)                                      !position of particle i 
          DO j=i+1,N                        
                  dr = ri - r(j,1:DIMS)                         ! distance  between particle i&j
                  dr = dr - L*ANINT(dr/L)                       ! Minimum image condition
                  dist2=SUM(dr*dr)                              ! distance2 between i and j    
                  IF( dist2 < (L/2.)**2) THEN                   ! keep the nearest distance if the image particle is actually closer                   
                    ig = NINT( SQRT(dist2)/w)                   !  the bin (distance range) value
                    g(ig) = g(ig) + 1                           ! frequency of the bins (distance range)
                  ENDIF
          END DO

  ELSE IF( switch == 0) THEN
    g = 0;
    w = (L/2.d0)/DOUBLE(NBINS)
  ELSE IF( switch == 2) THEN
    ! Write array g to file for subsequent visualization
  END IF
     
  
  
  
  
