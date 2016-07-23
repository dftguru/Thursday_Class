! Initialize particle positions r(t=0) and velocities v(t=0)
! Input: n number of elements

MODULE init(N)
  IMPLICIT NONE
  REAL ::  l, rho                                               ! l is the length of cube, rho is the 
  REAL, DIMENSION :: r, v
                                                                ! Calculate cell dimension of cubical side l
  rho = 1                                                       ! rho = 1 for a simple cubic lattice and 4 for
                                                                ! a face-centered cubic lattice
  l = (N/rho)**(1/3)                                            ! Determining length of cell
  l = CEIL(l)                                                   ! Rounding length of cell to integer

  SUBROUTINE rvinit(N,rho,T,sigma,eps,l)
    IMPLICIT NONE
    INTEGER :: l,i,j,k,counter                                  ! Length of cell (l), x, y, z coordinates of
                                                                ! cell (i,j,k), counter to keep track of number
                                                                ! of particles (counter) 
    REAL :: N, rho, T, sigma, eps                               ! Total number of particles (N), number of
                                                                ! particles per lattice (rho), temperature
                                                                ! (T), Maxwell-Boltzmann constants (sigma,
                                                                ! eps)
    REAL :: kB = 1                                              ! Boltzmann constant in atomic units
    REAL :: m = 1                                               ! Mass of particle in atomic units
    REAL, DIMENSION(N,3) :: r, v

                                                                ! Assign positions then velocities for particles 1 through N
    counter = 1                                                 ! Start counting number of particles
    DO WHILE (counter <= N)
      DO i=0,l                                                  ! x coordinate
         DO j=0,l                                               ! y coordinate
           DO k=0,l                                             ! z coordinate
             r(i,j,k) = counter                                 ! Assign number of particle to position r(i,j,k)
             call random_number(v)                              ! Assign random number from uniform
                                                                ! distribution
             v = 2*v-1                                          ! Normalize range from -1 to 1
             v = v * SQRT(3*N*kB*T/m*SUM(SUM(v*v))              ! Equipartition theorem
             
             counter = counter + 1                              ! Go to next particle
           END DO
         END DO
      END DO
    END DO
  END SUBROUTINE rvinit
END MODULE init(N)

  SUBROUTINE rdf(r,N)                                           !The subroutine for creating histogram such as radial distribution function        
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NBINS = 100                             !We initialize the number of bins.
  INTEGER :: i,j,N                                              !We initialize two counter variables i,j and the number of particles, N.
  REAL, DIMENSION :: r(1:N,1:3), g(1:N)                         !We initialize the arrays for distance frequency.
  REAL, DIMENSION :: ri(1:3), dr(1:3)                           !We initialize the position and distance

  SAVE g,krdf                                                   !SAVE the g and krdf  
  IF (switch == 1) THEN                                         !if it is after first  iteration / time step
  krdf = krdf + 1
  ELSE IF( switch == 0) THEN                                    ! First time step switch = 0 
  DO i=1, N-1
          ri = r(i,1:DIMS)                                      !position of particle i 
          DO j=i+1,N                        
                  dr = ri - r(j,1:DIMS)                         ! distance  between particle i&j
                  dist2=SUM(dr*dr)                              ! distance2 between i and j    
                  IF( dist2 < (L/2.)**2) THEN                   ! keep the nearest distance if the image particle is actually closer                   
                    ig = NINT( SQRT(dist2)/w)                   !  the bin (distance range) value
                    g(ig) = g(ig) + 1                           ! frequency of the bins (distance range)
                  ENDIF
                  IF (dist2 > (L/2.)**2) THEN
                    dr=L/2-dr                                   ! distance between i & j 
                    dist2=SUM(dr*dr)
                    ig = NINT(SQRT(dist2)/w)
                    g(ig)=g(ig)+1 
          END DO

  ELSE IF( switch == 2) THEN
  
  END IF
     
  
  
  
  
