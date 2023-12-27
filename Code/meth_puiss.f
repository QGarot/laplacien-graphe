      subroutine meth_puiss(A,n,niter,v1,vpmax)
      implicit none
      
      double precision, dimension(n,n) :: A
      integer :: n, niter, i
      double precision, dimension(n,1) :: v1
      double precision, intent(out) :: vpmax
      
      double precision, dimension(:,:), allocatable :: x, xbis
      double precision :: lambda, rand_val
      
      allocate(x(n,1))
      allocate(xbis(n,1))
      
      call random_seed()
      
      ! vecteur initial aleatoire
      do i=1,n
      	call random_number(rand_val)
      	x(i, 1) = rand_val
      enddo
      
      ! methode de puissance
      do i=1,niter
      	xbis = matmul(A,x)
      	lambda=norm2(xbis)
      	
      	! eviter le vecteur propre NaN en prenant le dernier plus proche donc x_{i-1}
      	if (lambda == 0.d0) then
      	  exit
      	endif
      	
      	x=xbis/lambda
      enddo
      
      vpmax=lambda
      v1=x
      
      return
      end
