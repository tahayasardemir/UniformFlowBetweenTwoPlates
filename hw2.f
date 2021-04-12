c..   Taha Yaşar Demir / 1881978
c..   Ce-580 / Homework #2
	  program ParallelPlates
	  parameter(mx=5001)
	  double precision u(mx), y(mx), A(mx), B(mx), C(mx), D(mx)
	  double precision E(mx), F(mx)
	  double precision H, Cp, nu, dy, Umax
	  double precision BC1, BC2
	  double precision ua(mx)
	  double precision ee
c	  real u(mx), y(mx), A(mx), B(mx), C(mx), D(mx), E(mx), F(mx)
c	  real H, Cp, nu, dy, Umax
c	  real BC1, BC2
c	  real ua(mx)
c	  real ee
	  integer N
	  character boundary
c-----------------------------------------------------------------------
c..Variable names:
c     u(x):velocity vevtor ,  y(x):distance from wall
c     A,B,C,D: coefficient vectors
c     E,F: coefficient of thomas algorithm
c     H: Half of the distance between plates, Umax: Analutical Max velocity
c     Cp: Pres. Coef , nu: viscousity ,dy:delta y
c     ua: Analytical velocity distribution vector
c     N: Toral grid points
c     ee: normalized error
c-----------------------------------------------------------------------
c..Create Output Files
	  open(1,file='Error.dat')
	  open(2,file='center_vel.dat')
	  open(3,file='analytic.dat')
c..Get the centerline boundary information
	  print*, "Specify centerline B.C"
	  print*, "d => Drichlet, n => Neumann"
	  read*, boundary

	  H    = 0.1  !m
	  Cp   = -1.  !Pa/m
	  nu   = 1E-3 !N.s/m^2
	  Umax = 5.   !m/s
c..Specify boundary conditions
	  BC1  = 0.0     ! wall B.C => 0 means no slip cond.
	  if (boundary.eq."d") then	! initiate boundary conditions
	  BC2  = Umax	! centerline B.C
	  elseif(boundary.eq."n") then
	  BC2  = 0.0     ! centerline B.C => zero means 0 gradient=symmetry
	  endif
c..Set first thomas coefficients
	  E(1) = 0.0
	  F(1) = BC1

c..Main program loop
	  do N=10,5000,10
	  	dy = H/(N-1)

	  	do i=2,N-1 ! Coefficient vectors
	  		A(i)   = 1.0
	  		B(i)   = 2.0
	  		C(i)   = 1.0
	  		D(i)   =-(Cp*dy**2)/nu
	  		y(i+1) = y(i)+dy
	  	enddo

	  	call thomas(N,A,B,C,D,E,F) ! Returns Thomas Coefficients

	  	if (boundary.eq."d") then	! set Umax based on B.C
	    u(N) = BC2
	  	elseif(boundary.eq."n") then
	  	u(N) = (F(N-1)+dy*BC2)/(1-E(N-1))	
	  	endif 

	  	do k=N-1,1,-1 ! Get the velocity distribution recursively
			u(k)   = E(k)*u(k+1) + F(k) ! Numerical solution
	  	enddo

	  	call Analytic(N,dy,y,ua) ! Get analytic solution
	  	call Error(N,Umax,u,ua,ee) ! Calculate normalized error
	  	call Output(N,ee,u(N),ua,y) ! Output the data

	  enddo

	  close(1)
	  close(2)
	  close(3)

	  stop
	  end

c..Thomas Algorithm
	  subroutine thomas(nn,ca,cb,cc,cd,ce,cf)
	  parameter(mx=1001)
c	  real ca(mx),cb(mx),cc(mx),cd(mx),ce(mx),cf(mx)
	  double precision ca(mx),cb(mx),cc(mx),cd(mx),ce(mx),cf(mx)	
	  integer nn

	  do i=2,nn-1
	  	ce(i) = ca(i)/(cb(i)-cc(i)*ce(i-1))
	  	cf(i) = (cd(i)+cc(i)*cf(i-1))/(cb(i)-cc(i)*ce(i-1))
	  enddo

	  return
	  end
c..Analytic Solution
	  subroutine Analytic(p,deltaY,ydist,uu)
	  parameter(mx=5001)
c	  real uu(mx), ydist(mx), deltaY
      double precision uu(mx), ydist(mx), deltaY
	  integer p

	  ydist(1) = 0.
	  do i=1,p
	  	uu(i) = -500*ydist(i)**2 + 100*ydist(i)
	  	ydist(i+1) = ydist(i) + deltaY
	  enddo

	  return
	  end
c..Error Calculation
	  subroutine Error(gnum,max,nu,ru,err)
	  parameter(mx=5001)
c	  real max,nu(mx),ru(mx),err
	  double precision max,nu(mx),ru(mx),err
	  integer gnum

	  err = 0.
	  do k=2,gnum
	  	err = err + abs(nu(k)-ru(k))
	  enddo
	  err = (err)/((gnum-1)*max)

	  return
	  end	
c..Oupput Data
	  subroutine Output(grid,e_r,center_vel,analytic_u,yy)
	  parameter(mx=5001)
c	  real e_r,center_vel,analytic_u(mx),yy(mx)
	  double precision e_r,center_vel,analytic_u(mx),yy(mx)
	  integer grid

	  write(1,*) grid, e_r
	  write(2,*) grid, center_vel
	  if(grid.eq.100) then
	  	do k=1,grid
	  	write(3,*) yy(k), analytic_u(k)
	   enddo
	  endif

	  return
	  end

