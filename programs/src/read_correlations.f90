!
!  Reads a (undirected) network and puts it into vector format
!
!  The iput file "edge_list.dat" should already be ordered,
!  and whithout repetitions, using:
!  sort -n file.dat | uniq -u > edge_list.dat
!
!
program main
implicit none
integer :: N, E  ! number of nodes, edges
integer, allocatable :: A(:), Pi(:), Pf(:), k(:)  ! Adjacency vector, initial and final pointers, degree of the nodes
integer, allocatable :: P(:), Pc(:) ! Degree distributions
integer :: x,y  ! Indices for reading
integer :: i,j  ! Indices for loops
! ------------------------------------------------------------------

! First read to count the number of edges
E = 0
N = 0
open(10,file='edge_list.dat') ! input file
open(20,file='array.dat')  ! output file

do while(1.eq.1)
   read(10,*, END=200) x,y
   N = max(N,y)
   E = E + 1
enddo
200 close(10)
N = N + 1   ! Counting the zero node

!write(20,*) "N"
write(20,*) N
!write(20,*) "E"
write(20,*) E
allocate (Pi(N),Pf(N),k(N),A(2*E))
A = 0
Pi = 0
Pf = 0

! We still don't know the number of connections for each node
k = 0
! Second read
open(10,file='edge_list.dat')
do i=1,E
   read(10,*) x,y
   k(x+1) = k(x+1) + 1
   k(y+1) = k(y+1) + 1
enddo
close(10)


! Create pointers
Pi(1) = 1
Pf(1) = 0
do i=2,N
   Pi(i) = Pi(i-1) + k(i-1)
   Pf(i) = Pi(i) - 1
enddo

! Third read to generate the vectors
open(10,file='edge_list.dat')
do i=1,E
   read(10,*) x,y
   x = x + 1
   y = y + 1

! Skip "autointeractions", e.i., nodes connected to themselves
   if (x.eq.y) cycle
! Check when reading (2,1) if (1,2) already exists
   if (x.gt.y) then
      do j=Pi(y),Pf(y)
        if (A(j).eq.x) goto 300
      enddo
   endif

   Pf(x) =  Pf(x) + 1
   A(Pf(x)) = y
   Pf(y) =  Pf(y) + 1
   A(Pf(y)) = x
   300 continue
enddo
close(10)

!write(20,*) "Pi"
do i=1,N
   write(20,*) Pi(i)
enddo

!write(20,*) "Pf"
do i=1,N
   write(20,*) Pf(i)
enddo

!write(20,*) "A"
do i=1,2*E
   write(20,*) A(i)
enddo

close(20)


! correlations
allocate (P(0:maxval(k)),Pc(0:maxval(k)))
P = 0
Pc = 0
do i=1,N
    P(k(i)) = P(k(i)) + 1
enddo

Pc(0) = P(0)
do i=1,maxval(k)
    !Pc(i) = sum(P(0:i))
    Pc(i) = Pc(i-1) + P(i)
    print*, i,dble(P(i))/dble(N),dble(Pc(i))/dble(N)
enddo

k_av = dlbe(sum(k))/dble(N)
theta = 0d0
do i=2,hmaxval(k)
    theta = theta + dlbe(i-1)*P(i)*rho
enddo
theta = theta/k_av

phi = 0d0
do i=2,hmaxval(k)
    phi = phi + dlbe(i-1)*P(i)*recovered
enddo
phi =phi/k_av
end program main
