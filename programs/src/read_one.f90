!
!  Reads a (undirected) network and puts it into vector format
!
!  The iput file "edge_list.dat" should already be ordered,
!  and whithout repetitions, using:
!  sort -n file.dat | uniq -u > edge_list.dat
!
program main
implicit none
integer :: N=0, E=0  ! number of nodes, edges
integer, allocatable :: A(:), Pi(:), Pf(:), k(:)  ! Adjacency vector, initial and final pointers, degree of the nodes
integer, allocatable :: P(:), Pc(:) ! Degree distributions
integer, allocatable :: skip(:,:)
integer :: x,y !x_max=0,y_max=0  ! Indices for reading
integer :: gap_n
integer :: i,j  ! Indices for loops
! ------------------------------------------------------------------

! First read to count the number of edges
open(10,file='edge_list.dat') ! input file
open(20,file='array.dat')  ! output file

do while(1.eq.1)
   read(10,*, END=200) x,y
   N = max(N,x,y)
   !x_max = max(x_max,x)
   !y_max = max(y_max,y)
   !if (x.gt.N) N = N + 1
   E = E + 1
enddo
200 close(10)
!if (y_max.gt.x_max) print*, 'bigger node on right column'


! We still don't know the number of connections for each node
allocate (k(N))
k = 0
! Second read
open(10,file='edge_list.dat')
do i=1,E
   read(10,*) x,y
   k(x) = k(x) + 1
   k(y) = k(y) + 1
enddo
close(10)


!  __ _  __ _ _ __  ___
! / _` |/ _` | '_ \/ __|
!| (_| | (_| | |_) \__ \
! \__, |\__,_| .__/|___/
! |___/      |_|
! ------------------------------------------------------------------

! if necessary re-writes the edge list without number gaps
if (count(k.eq.0).eq.0) goto 201
!print*, N,count(k.eq.0)
allocate (skip(N,2))
i = 0
j = 0 ! counts how big the gap between numbers is
gap_n = 0 ! counts number of gaps
x = N - count(k.eq.0)
do while (N.gt.x)
i = i + 1
    if (j.eq.0) then
        if ( k(i).eq.0 ) then
            j=1
            gap_n = gap_n + 1
            skip(gap_n,1) = i
        endif
    else
        if ( k(i).eq.0 ) then
            j = j + 1
        else
            N = N - j
            skip(gap_n,2) = j
            j = 0
        endif
    endif
enddo
!print*, 'skip matrix, of dim',gap_n
!do i=1,gap_n
    !print*, skip(i,1),skip(i,2)
!enddo

open(10,file='edge_list.dat')
open(11,file='edge_list_2.dat')
do i=1,E
    read(10,*) x,y
    do j=gap_n,1,-1
        if(x.ge.skip(j,1)) x = x - skip(j,2)
        if(y.ge.skip(j,1)) y = y - skip(j,2)
    enddo
    write(11,*) x,y
enddo
close(10)
close(11)
call execute_command_line('mv edge_list_2.dat edge_list.dat')

! Re-do second read
deallocate(k,skip)
allocate (k(N))
k = 0
open(10,file='edge_list.dat')
do i=1,E
   read(10,*) x,y
   k(x) = k(x) + 1
   k(y) = k(y) + 1
enddo
close(10)
! ------------------------------------------------------------------
201 continue

! Create pointers
allocate (Pi(N),Pf(N),A(2*E))
Pi(1) = 1
Pf(1) = 0
do i=2,N
    Pi(i) = Pi(i-1) + k(i-1)
    Pf(i) = Pi(i) - 1
enddo

! Third read to generate the vectors
open(10,file='edge_list.dat')
A = 0
do i=1,E
   read(10,*) x,y

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

write(20,*) N
write(20,*) E

do i=1,N
   write(20,*) Pi(i)
enddo

do i=1,N
   write(20,*) Pf(i)
enddo

do i=1,2*E
   write(20,*) A(i)
enddo

close(20)

! degree distribution
open(30,file='degree.dat') ! input file
allocate (P(0:maxval(k)),Pc(0:maxval(k)))
P = 0
Pc = 0
do i=1,N
    P(k(i)) = P(k(i)) + 1
enddo

Pc(0) = P(0)
do i=1,maxval(k)
    !Pc(i) = sum(P(0:i))
    !k_av = i*P(i)
    Pc(i) = Pc(i-1) + P(i)
    write(30,*) i,dble(P(i))/dble(N),dble(Pc(i))/dble(N)
enddo
close(30)

end program main
