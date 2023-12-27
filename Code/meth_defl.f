      program meth_defl

      implicit none
      
      integer, parameter:: niter=250
      integer :: n, i
      double precision :: vpmax
      double precision, dimension(:,:), allocatable :: A, v1
      
      read(*,*) n
      
      allocate(A(n,n))
      allocate(v1(n,1))
      
      ! lecture de la matrice
      do i=1,n
      	read(*,*) A(i,:)
      enddo
      
      ! methode de deflation
      do i=1,n
      	call meth_puiss(A,n,niter,v1,vpmax)
      	
      	write(*,*) "val propre:", vpmax, "; vec propre:", v1
      	
      	A = A - (vpmax/(norm2(v1)**2))*matmul(v1,transpose(v1))
      enddo
      
      end
