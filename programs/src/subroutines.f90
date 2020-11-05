MODULE SUBROUTINES
CONTAINS

! *******************************************************************
! IMMUNITY
! *******************************************************************
subroutine enlarge_immune(immune,N,immune_length)
! increases the size of immune by 5
implicit none
integer, intent(in) :: N
integer, intent(inout) :: immune_length
integer, dimension(:,:), allocatable,intent(inout) :: immune
integer :: immune_length_old
integer :: immune2(N,immune_length)

immune2 = immune
deallocate(immune)
immune_length_old = immune_length
immune_length = immune_length + 5
allocate(immune(N,immune_length))
immune(:,:immune_length_old) = immune2
immune(:,immune_length_old+1:) = 0
end subroutine enlarge_immune

! *******************************************************************

subroutine enlarge_infected_i(n_infected_i,n_recovered_i,mutation_time_prevalence,infected_i_length)
! increases the size of infected_i_length by 10
implicit none
integer, intent(inout) :: infected_i_length
integer, allocatable, intent(inout) :: n_infected_i(:),n_recovered_i(:)
real*8, allocatable, intent(inout) :: mutation_time_prevalence(:,:)
integer :: n_infected_i_2(infected_i_length),n_recovered_i_2(infected_i_length)
real*8 :: mutation_time_prevalence_2(infected_i_length,2)
integer :: infected_i_length_old

n_infected_i_2 = n_infected_i
n_recovered_i_2 = n_recovered_i
mutation_time_prevalence_2 = mutation_time_prevalence
deallocate(n_infected_i,n_recovered_i,mutation_time_prevalence)
infected_i_length_old = infected_i_length
infected_i_length = infected_i_length + 10
allocate(n_infected_i(infected_i_length),n_recovered_i(infected_i_length))
allocate(mutation_time_prevalence(infected_i_length,2))
n_infected_i(:infected_i_length_old) = n_infected_i_2
n_infected_i(infected_i_length_old+1:) = 0
n_recovered_i(:infected_i_length_old) = n_recovered_i_2
n_recovered_i(infected_i_length_old+1:) = 0
mutation_time_prevalence(:infected_i_length_old,:) = mutation_time_prevalence_2
mutation_time_prevalence(infected_i_length_old+1:,:) = 0
end subroutine enlarge_infected_i

! *******************************************************************

subroutine immunize(immune,state,node_n,N,immune_length)
! immunize of current mutation
! if needed calls "enlarge_immune" subroutine to make immune matrix bigger
implicit none
integer, intent(in) :: N
integer, intent(inout) :: immune_length
integer, allocatable,intent(inout) :: immune(:,:)
integer, intent(in) :: node_n
integer, intent(in) :: state(N)
integer :: i
100 continue
do i=1,immune_length
    if(immune(node_n,i).eq.state(node_n)) then
        print*, 'error immunity'
        print*, immune(node_n,:)
        print*, state(node_n)
        stop
    endif
    if(immune(node_n,i).ne.0) cycle
    immune(node_n,i) = state(node_n)
    exit
enddo
if(i.eq.immune_length+1) then
    call enlarge_immune(immune,N,immune_length)
goto 100
endif
end subroutine immunize

! *******************************************************************
! Pointer to active_edges (array B)
! *******************************************************************

subroutine change_B_1(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E)
implicit none
integer, intent(in) :: j,N,E
integer, intent(inout) :: n_active_edges,active_edges(E,2)
integer, intent(in) :: Pi(N),Pf(N),A(2*E)
integer, intent(inout) :: B(2*E)
integer :: k,index_n_B,i,iii
!print*, "change 1"
index_n_B = B(j)
! move last active edge to the position that is deleted
active_edges(index_n_B,:) = active_edges(n_active_edges,:)
! change values in B
i = active_edges(n_active_edges,1)
iii = active_edges(n_active_edges,2)
do k=Pi(i),Pf(i)
    if(A(k).eq.iii) then
        B(k) = index_n_B
        exit
    endif
enddo
do k=Pi(iii),Pf(iii)
    if(A(k).eq.i) then
        B(k) = index_n_B
        exit
    endif
enddo
n_active_edges = n_active_edges - 1
end subroutine change_B_1

! *******************************************************************

subroutine change_B_2(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E,node_n,index_n)
implicit none
integer, intent(in) :: j,N,E,index_n,node_n
integer, intent(inout) :: n_active_edges,active_edges(E,2)
integer, intent(in) :: Pi(N),Pf(N),A(2*E)
integer, intent(inout) :: B(2*E)
integer :: k
!print*, "change 2"
n_active_edges = n_active_edges + 1
active_edges(n_active_edges,:) = (/node_n, index_n /)
! Write number of active edge in B for both nodes
B(j) = n_active_edges
do k=Pi(index_n),Pf(index_n)
    if(A(k).eq.node_n) then
        B(k) = n_active_edges
        exit
    endif
enddo
end subroutine change_B_2

! *******************************************************************

subroutine change_B_3(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E,node_n,index_n)
implicit none
integer, intent(in) :: j,N,E,index_n,node_n
integer, intent(inout) :: n_active_edges,active_edges(E,2)
integer, intent(in) :: Pi(N),Pf(N),A(2*E)
integer, intent(inout) :: B(2*E)
integer :: k
!print*, "change 3"
n_active_edges = n_active_edges + 1
active_edges(n_active_edges,:) = (/index_n, node_n/)
! Write number of active edge in B for both nodes
B(j) = n_active_edges
do k=Pi(index_n),Pf(index_n)
    if(A(k).eq.node_n) then
        B(k) = n_active_edges
        exit
    endif
enddo
end subroutine change_B_3

! *******************************************************************

subroutine change_B_122(Pi,Pf,A,B22,active_edges22,n_active_edges22,j,N,E)
implicit none
integer, intent(in) :: j,N,E
integer, intent(inout) ::n_active_edges22,active_edges22(E,2)
integer, intent(in) :: Pi(N),Pf(N),A(2*E)
integer, intent(inout) ::B22(2*E)
integer :: k,index_n_B,i,iii
!print*, "change 1"
index_n_B =B22(j)
! move last active edge to the position that is deleted
active_edges22(index_n_B,:) =active_edges22(n_active_edges22,:)
! change values in B
i =active_edges22(n_active_edges22,1)
iii =active_edges22(n_active_edges22,2)
do k=Pi(i),Pf(i)
    if(A(k).eq.iii) then
        B22(k) = index_n_B
        exit
    endif
enddo
do k=Pi(iii),Pf(iii)
    if(A(k).eq.i) then
        B22(k) = index_n_B
        exit
    endif
enddo
n_active_edges22 = n_active_edges22 - 1
end subroutine change_B_122


END MODULE SUBROUTINES

! *******************************************************************
! Prevalence minmax
! *******************************************************************

subroutine prevalence(n_infected,N,time,prev,prev_der,modifier,minimum,maximum)
! call every ~ 100 steps
implicit none
real*8,intent(in) :: time
integer,intent(in) :: N,n_infected
integer,intent(inout) :: prev(2),prev_der(10)
integer,intent(inout) :: modifier ! = 0 init, 1 up, 2 down
integer,intent(inout) :: maximum,minimum
integer :: prev_der_sum,i

prev(2) = prev(1)
prev(1) = n_infected

if(modifier.eq.0) then
    i = 0
    do while(i.le.9)
        i = i + 1
        if (prev_der(10-i).lt.(2*N)) cycle
        prev_der(10-i) = prev(1)-prev(2)
        if(i.eq.9) modifier=1
        return
    enddo
endif

! move all by one place
do i=10,2,-1
    prev_der(i)=prev_der(i-1)
enddo
prev_der(1)=prev(1)-prev(2)
prev_der_sum=sum(prev_der)

if(modifier.eq.1) then
    if(prev_der_sum.gt.0) then
        maximum=max(maximum,n_infected)
    else
        write(35,*) maximum
        modifier=2
        minimum=maximum
        return
    endif
else
    if(prev_der_sum.gt.0) then
        write(36,*) minimum
        modifier=1
        maximum=minimum
        return
    else
        minimum=min(minimum,n_infected)
    endif
endif
end subroutine prevalence

! *******************************************************************
! Running variance
! *******************************************************************

subroutine running_variance(filename,mean,sigma)
implicit none
real*8 :: M,M1,S=0d0
real*8, intent(out) :: mean,sigma
character*32, intent(in) :: filename
real*8 :: x
integer :: i
open(10,file=filename) ! input file

i = 0
do while(1.eq.1)
    i = i + 1
    read(10,*, END=200) x
    mean = mean + x
    if(i.eq.1) then
        M=x
        cycle
    endif
    M1 = M
    M = M + (x-M)/dble(i)
    S = S + (x-M1)*(x-M)
enddo
200 close(10)
mean = mean/dble(i)
sigma = sqrt(S/dble(i-1))
end
