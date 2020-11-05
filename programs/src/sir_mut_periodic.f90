! Susceptible (0), Infected(1) model with mutations
! Network in edge list format should first be converted to arrays using
! either read_zero.f90 or read_one.f90
program main
! random numbers from ran2
USE SUBROUTINES
implicit none
integer :: SEED, NSEED, SEED0, idum
integer :: N, E,N2  ! number of nodes, edges
! Adjacency vector, initial and final pointers, state
integer, allocatable :: A(:), Pi(:), Pf(:), state(:)
! pointer to active_edges
integer, allocatable :: B(:)
integer, allocatable :: B11(:), B12(:), B22(:)
integer, allocatable :: infected(:), active_edges(:,:)
integer, allocatable :: infected1(:), infected2(:)
integer, allocatable :: active_edges11(:,:),active_edges12(:,:),active_edges22(:,:)
integer :: i,iii,j,jjj ! indices for loops
integer :: n_infected, n_active_edges
integer :: n_infected1, n_infected2
integer :: n_active_edges11, n_active_edges12, n_active_edges22
integer :: node_n, index_n, index_inf, index_act
integer :: index_inf1, index_inf2
! initial (microscopic) number of nodes to infect
integer :: x, state_old
integer, allocatable :: immune(:,:)
integer :: immune_length=5, n_mutation !,mutability
integer, allocatable :: endemic(:)
! probability of infection
real*8 :: lambda, lambda_sum, prob_heal1, prob_heal2
real*8 :: prob_infect11,prob_infect22,prob_infect12
real*8 :: lambda0,lambda1,lambda_inv
real*8 :: mutation_rate !, prob_mutate
real*8 :: poisson_time, time, time_max, time_write, time_inc
real*8 :: ran2, random, time_sum, time_sum_t
real*8 :: delta,tau,delta_max,delta_value,delta_inv
real*8, parameter :: pi_greco=4d0*atan(1d0)
logical :: exist
!real :: t_start, t_actual
NAMELIST/DADES/x,lambda0,lambda1,mutation_rate,SEED0,NSEED,time_max,time_inc,tau,delta_max
! ---------------------------------------------------------------
! Read input file
open(unit=11,file='../inp/sir_mut_periodic.inp',status='old') ! parallel
read(11,DADES)
close(11)
lambda0 = lambda0/100d0
lambda1 = lambda1/100d0
mutation_rate = mutation_rate/1d5
! ---------------------------------------------------------------
! Read array containing the network
open(20,file='array.dat')
read(20,*) N
N2 = int(dble(N)/2d0)
read(20,*) E

allocate (Pi(N),Pf(N),A(2*E),B(2*E),state(N))
allocate (infected(N),active_edges(E,2))
allocate(immune(N,immune_length))
allocate(endemic(int(time_max/time_inc)))
allocate(infected1(N2),infected2(N2))
allocate(active_edges11(E,2),active_edges12(E,2),active_edges22(E,2))
allocate(B11(2*E),B12(2*E),B22(2*E))

do i=1,N
   read(20,*) Pi(i)
enddo

do i=1,N
   read(20,*) Pf(i)
enddo

do i=1,2*E
   read(20,*) A(i)
enddo
close(20)


endemic = 0
time_sum = 0d0
time_sum_t = 0d0
!mutability = 0
!call cpu_time(t_start)
! ------------------------------------------------------------
! seed loop
loop_seed: do SEED = SEED0, SEED0+NSEED-1
    idum = -SEED
!-------
! Initialize population of first network
state=0
n_infected = 0
n_infected1 = 0
n_infected2 = 0
n_active_edges = 0
n_active_edges11 = 0
n_active_edges12 = 0
n_active_edges22 = 0
n_mutation = 1
immune = 0
! ===================
jjj = 1
loop0: do while (jjj.lt.x)
! ===================
node_n = int(ran2(idum)*N2) + 1
if(state(node_n).eq.1) cycle loop0
state(node_n) = 1
n_infected = n_infected + 1
n_infected1 = n_infected1 + 1
infected(n_infected) = node_n
infected1(n_infected1) = node_n

! Check all neighbors through A
do j=Pi(node_n),Pf(node_n)
    ! neighbor number is the value of A(j)
    index_n = A(j)
    ! if neighbor is active delete active edge
    if (state(index_n).eq.1) then
        call change_B_1(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E)
        if(index_n.le.N2) then
            call change_B_1(Pi,Pf,A,B11,active_edges11,n_active_edges11,j,N,E)
        else
            call change_B_1(Pi,Pf,A,B12,active_edges12,n_active_edges12,j,N,E)
        endif
    ! if neighbor is susceptible create new active edge
    elseif (state(index_n).eq.0) then
        call change_B_2(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E,node_n,index_n)
        if(index_n.le.N2) then
            call change_B_2(Pi,Pf,A,B11,active_edges11,n_active_edges11,j,N,E,node_n,index_n)
        else
            call change_B_2(Pi,Pf,A,B12,active_edges12,n_active_edges12,j,N,E,node_n,index_n)
        endif
    endif
enddo
jjj = jjj + 1
! ===================
enddo loop0
! ===================
!do i=1,n_active_edges
    !print*, active_edges(i,:)
!enddo
!do i=1,n_active_edges11
    !write(1,*) active_edges11(i,:)
!enddo
!do i=1,n_active_edges12
    !write(2,*) active_edges12(i,:)
!enddo
!stop

! *******************************************************************
!  _   _                  _
! | |_(_)_ __ ___   ___  | | ___   ___  _ __
! | __| | '_ ` _ \ / _ \ | |/ _ \ / _ \| '_ \
! | |_| | | | | | |  __/ | | (_) | (_) | |_) |
!  \__|_|_| |_| |_|\___| |_|\___/ \___/| .__/
!                                      |_|
time = 0d0
delta_value=delta(time,tau,delta_max)
delta_inv=1d0+delta_max-delta_value
lambda = lambda0*(1d0+lambda1*cos(2d0*pi_greco*time/tau))
lambda_inv = lambda0*(1d0-lambda1*cos(2d0*pi_greco*time/tau))
time_write = time_inc
jjj = 1
print*, '#  time, prevalence, number of mutations, immune_length'
print*, time, dble(n_infected1)/dble(N2),dble(n_infected2)/dble(N2),delta_value
! ===================
loop_time: do while ( (n_infected.gt.0) .and. (time.lt.time_max) )
! ===================
lambda_sum = dble(n_infected2)*delta_value+dble(n_infected1)*delta_inv+&
    & dble(n_active_edges11)*lambda+dble(n_active_edges22)*lambda_inv+&
    & dble(n_active_edges12)*sqrt(lambda*lambda_inv)+dble(n_infected)*mutation_rate
prob_heal1 = dble(n_infected1)/lambda_sum * delta_inv
prob_heal2 = dble(n_infected2)/lambda_sum * delta_value
prob_infect11 = dble(n_active_edges11)*lambda/lambda_sum
prob_infect12 = dble(n_active_edges12)*sqrt(lambda*lambda_inv)/lambda_sum
prob_infect22 = dble(n_active_edges22)*lambda_inv/lambda_sum
!print*, prob_heal1, prob_heal2, prob_infect, 1d0-prob_heal1-prob_heal2-prob_infect,lambda_sum
random = ran2(idum)
if (random.lt.prob_heal1) then
! heal1
    index_inf1 = int(ran2(idum)*n_infected1) + 1
    node_n = infected1(index_inf1)
    ! find node from infected1 in global list
    do j=1,n_infected
        if(infected(j).eq.node_n) then
            index_inf=j
            exit
        endif
    enddo
    call immunize(immune,state,node_n,N,immune_length)
    state_old = state(node_n)
    state(node_n) = 0
    infected(index_inf) = infected(n_infected)
    infected1(index_inf1) = infected1(n_infected1)
    n_infected = n_infected - 1
    n_infected1 = n_infected1 - 1
    neighbors_h1: do j=Pi(node_n),Pf(node_n)
        index_n = A(j)
        if (state(index_n).eq.0) then
            ! check if S neighbor was immune to the mutation of the
            ! now healed node
            do i=1,immune_length
                if(immune(index_n,i).eq.0) exit
                if(immune(index_n,i).ne.state_old) cycle
                cycle neighbors_h1
            enddo
        ! if neighbor was S remove active edge
            call change_B_1(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E)
            if(index_n.le.N2) then
                call change_B_1(Pi,Pf,A,B11,active_edges11,n_active_edges11,j,N,E)
            else
                call change_B_1(Pi,Pf,A,B12,active_edges12,n_active_edges12,j,N,E)
            endif
        else
        ! if neighbor is I with a mutation that the node is not
        ! immune to create active edge
            do i=1,immune_length
                if(immune(node_n,i).eq.0) exit
                if(immune(node_n,i).ne.state(index_n)) cycle
                cycle neighbors_h1
            enddo
            call change_B_3(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E,node_n,index_n)
            if(index_n.le.N2) then
                call change_B_3(Pi,Pf,A,B11,active_edges11,n_active_edges11,j,N,E,node_n,index_n)
            else
                call change_B_3(Pi,Pf,A,B12,active_edges12,n_active_edges12,j,N,E,node_n,index_n)
            endif
        endif
    enddo neighbors_h1

elseif ( random .lt. (prob_heal1+prob_heal2) ) then
! heal2
    index_inf2 = int(ran2(idum)*n_infected2) + 1
    node_n = infected2(index_inf2)
    ! find node from infected2 in global list
    do j=1,n_infected
        if(infected(j).eq.node_n) then
            index_inf=j
            exit
        endif
    enddo
    call immunize(immune,state,node_n,N,immune_length)
    state_old = state(node_n)
    state(node_n) = 0
    infected(index_inf) = infected(n_infected)
    n_infected = n_infected - 1
    infected2(index_inf2) = infected2(n_infected2)
    n_infected2 = n_infected2 - 1
    neighbors_h2: do j=Pi(node_n),Pf(node_n)
        index_n = A(j)
        if (state(index_n).eq.0) then
            ! check if S neighbor was immune to the mutation of the
            ! now healed node
            do i=1,immune_length
                if(immune(index_n,i).eq.0) exit
                if(immune(index_n,i).ne.state_old) cycle
                cycle neighbors_h2
            enddo
        ! if neighbor was S remove active edge
            call change_B_1(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E)
            if(index_n.gt.N2) then
                call change_B_1(Pi,Pf,A,B22,active_edges22,n_active_edges22,j,N,E)
                    print*, n_active_edges22
            else
                call change_B_1(Pi,Pf,A,B12,active_edges12,n_active_edges12,j,N,E)
            endif
        else
        ! if neighbor is I with a mutation that the node is not
        ! immune to create active edge
            do i=1,immune_length
                if(immune(node_n,i).eq.0) exit
                if(immune(node_n,i).ne.state(index_n)) cycle
                cycle neighbors_h2
            enddo
            call change_B_3(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E,node_n,index_n)
            if(index_n.gt.N2) then
                call change_B_3(Pi,Pf,A,B22,active_edges22,n_active_edges22,j,N,E,node_n,index_n)
            else
                call change_B_3(Pi,Pf,A,B12,active_edges12,n_active_edges12,j,N,E,node_n,index_n)
            endif
        endif
    enddo neighbors_h2

elseif ( random .lt. (prob_heal1+prob_heal2+prob_infect11) ) then
! infect 11
print*, 11
    index_act = int(ran2(idum)*n_active_edges11) + 1
    ! node to be infected
    node_n = active_edges11(index_act,2)
    state(node_n) = state(active_edges11(index_act,1))
    n_infected = n_infected + 1
    infected(n_infected) = node_n
    if(node_n.le.N2) then
        n_infected1 = n_infected1 + 1
        infected1(n_infected1) = node_n
    else
        n_infected2 = n_infected2 + 1
        infected2(n_infected2) = node_n
    endif

    neighbors_i11: do j=Pi(node_n),Pf(node_n)
        index_n = A(j)
        ! if neighbor is active delete active edge
        if (state(index_n).ne.0) then
            ! check if was immune to the mutation of the infected
            ! neighbor node
            do i=1,immune_length
                if(immune(node_n,i).eq.0) exit
                if(immune(node_n,i).ne.state(index_n)) cycle
                cycle neighbors_i11
            enddo
            ! move last active edge to the position that is deleted
            call change_B_1(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E)
            if(index_n.le.N2) then
                call change_B_1(Pi,Pf,A,B11,active_edges11,n_active_edges11,j,N,E)
            else
                call change_B_1(Pi,Pf,A,B12,active_edges12,n_active_edges12,j,N,E)
            endif
        else
            do i=1,immune_length
                if(immune(index_n,i).eq.0) exit
                if(immune(index_n,i).ne.state(node_n)) cycle
                cycle neighbors_i11
            enddo
            call change_B_2(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E,node_n,index_n)
            if(index_n.le.N2) then
                call change_B_2(Pi,Pf,A,B11,active_edges11,n_active_edges11,j,N,E,node_n,index_n)
            else
                call change_B_2(Pi,Pf,A,B12,active_edges12,n_active_edges12,j,N,E,node_n,index_n)
            endif
        endif
    enddo neighbors_i11
elseif ( random .lt. (prob_heal1+prob_heal2+prob_infect11+prob_infect22) ) then
! infect 22
print*, 22
    index_act = int(ran2(idum)*n_active_edges22) + 1
    ! node to be infected
    node_n = active_edges22(index_act,2)
    state(node_n) = state(active_edges22(index_act,1))
    n_infected = n_infected + 1
    infected(n_infected) = node_n
    if(node_n.le.N2) then
        n_infected1 = n_infected1 + 1
        infected1(n_infected1) = node_n
    else
        n_infected2 = n_infected2 + 1
        infected2(n_infected2) = node_n
    endif

    neighbors_i22: do j=Pi(node_n),Pf(node_n)
        index_n = A(j)
        ! if neighbor is active delete active edge
        if (state(index_n).ne.0) then
            ! check if was immune to the mutation of the infected
            ! neighbor node
            do i=1,immune_length
                if(immune(node_n,i).eq.0) exit
                if(immune(node_n,i).ne.state(index_n)) cycle
                cycle neighbors_i22
            enddo
            ! move last active edge to the position that is deleted
            call change_B_1(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E)
            if(index_n.gt.N2) then
                call change_B_1(Pi,Pf,A,B22,active_edges22,n_active_edges22,j,N,E)
                    print*, n_active_edges22
            else
                call change_B_1(Pi,Pf,A,B12,active_edges12,n_active_edges12,j,N,E)
            endif
        else
            do i=1,immune_length
                if(immune(index_n,i).eq.0) exit
                if(immune(index_n,i).ne.state(node_n)) cycle
                cycle neighbors_i22
            enddo
            call change_B_2(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E,node_n,index_n)
            if(index_n.le.N2) then
                call change_B_2(Pi,Pf,A,B22,active_edges22,n_active_edges22,j,N,E,node_n,index_n)
            else
                call change_B_2(Pi,Pf,A,B12,active_edges12,n_active_edges12,j,N,E,node_n,index_n)
            endif
        endif
    enddo neighbors_i22
elseif ( random .lt. (prob_heal1+prob_heal2+prob_infect11+prob_infect22+prob_infect12) ) then
! infect 12
print*, 12
    index_act = int(ran2(idum)*n_active_edges12) + 1
    ! node to be infected
    node_n = active_edges12(index_act,2)
    state(node_n) = state(active_edges12(index_act,1))
    n_infected = n_infected + 1
    infected(n_infected) = node_n
    if(node_n.le.N2) then
        n_infected1 = n_infected1 + 1
        infected1(n_infected1) = node_n
    else
        n_infected2 = n_infected2 + 1
        infected2(n_infected2) = node_n
    endif

    neighbors_i12: do j=Pi(node_n),Pf(node_n)
        index_n = A(j)
        ! if neighbor is active delete active edge
        if (state(index_n).ne.0) then
            ! check if was immune to the mutation of the infected
            ! neighbor node
            do i=1,immune_length
                if(immune(node_n,i).eq.0) exit
                if(immune(node_n,i).ne.state(index_n)) cycle
                cycle neighbors_i12
            enddo
            ! move last active edge to the position that is deleted
            call change_B_1(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E)
            if(node_n.gt.N2) then
                if(index_n.gt.N2) then
                    call change_B_1(Pi,Pf,A,B22,active_edges22,n_active_edges22,j,N,E)
                    print*, n_active_edges22
                else
                    call change_B_1(Pi,Pf,A,B12,active_edges12,n_active_edges12,j,N,E)
                endif
            else
                call change_B_1(Pi,Pf,A,B11,active_edges11,n_active_edges11,j,N,E)
            endif
        else
            do i=1,immune_length
                if(immune(index_n,i).eq.0) exit
                if(immune(index_n,i).ne.state(node_n)) cycle
                cycle neighbors_i12
            enddo
            call change_B_2(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E,node_n,index_n)
            if(node_n.gt.N2) then
                if(index_n.le.N2) then
                    call change_B_2(Pi,Pf,A,B22,active_edges22,n_active_edges22,j,N,E,node_n,index_n)
                else
                    call change_B_2(Pi,Pf,A,B12,active_edges12,n_active_edges12,j,N,E,node_n,index_n)
                endif
            else
                call change_B_2(Pi,Pf,A,B11,active_edges11,n_active_edges11,j,N,E,node_n,index_n)
            endif
        endif
    enddo neighbors_i12
else
!mutate
    !mutability = mutability + 1
    !exit loop_time
    index_inf = int(ran2(idum)*n_infected) + 1
    node_n = infected(index_inf)
    ! immunize to current mutation before mutating
    call immunize(immune,state,node_n,N,immune_length)
    n_mutation = n_mutation + 1
    state_old = state(node_n)
    state(node_n) = n_mutation
    ! neighbors that were immune are not anymore
    ! so active edges must be created
    neighbors_m: do j=Pi(node_n),Pf(node_n)
        index_n = A(j)
        if (state(index_n).eq.0) then
            iii=0
            do i=1,immune_length
                if(immune(index_n,i).eq.0) exit
                if(immune(index_n,i).eq.state_old) iii=1
            enddo
            if(iii.eq.1) then
                call change_B_2(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E,node_n,index_n)
                if(node_n.gt.N2) then
                    if(index_n.le.N2) then
                        call change_B_2(Pi,Pf,A,B22,active_edges22,n_active_edges22,j,N,E,node_n,index_n)
                    else
                        call change_B_2(Pi,Pf,A,B12,active_edges12,n_active_edges12,j,N,E,node_n,index_n)
                    endif
                else
                    call change_B_2(Pi,Pf,A,B11,active_edges11,n_active_edges11,j,N,E,node_n,index_n)
                endif
            endif
        endif
    enddo neighbors_m
endif
print*, time, dble(n_infected1)/dble(N2),dble(n_infected2)/dble(N2),delta_value
print*, n_active_edges11,n_active_edges12,n_active_edges22,n_active_edges,&
    & n_active_edges11+n_active_edges12+n_active_edges22

! endemic array instead of number to do many times with just one run
if (time.gt.time_write) then
    if (jjj.eq.(int(time_max/time_inc))) exit loop_time
    endemic(jjj) = endemic(jjj) + 1
    time_write = time_write + time_inc
    jjj = jjj + 1
endif
time = time + poisson_time(lambda_sum,idum)
delta_value=delta(time,tau,delta_max)
delta_inv=1d0+delta_max-delta_value
lambda = lambda0*(1d0+lambda1*cos(2d0*pi_greco*time/tau))
lambda_inv = lambda0*(1d0-lambda1*cos(2d0*pi_greco*time/tau))
! ===================
enddo loop_time
! ===================
if(time.gt.time_max) then
    endemic(jjj) = endemic(jjj) + 1
    !time_sum_t = time_sum_t + time
!else
    !time_sum = time_sum + time
endif
!call cpu_time(t_actual)
!print*, mutation_rate, int(dble(seed-seed0+1)/dble(nseed)*1d2),'%',t_actual-t_start,&
    !& (t_actual-t_start)/(1d0-dble(seed-seed0+1)/dble(nseed))
! ===================
enddo loop_seed
! ===================
inquire(file='endemic.dat', exist=exist)
if(.not.exist) call execute_command_line('touch endemic.dat')
open(25,file='endemic.dat',access='append',status='old')
write(25,*) lambda, mutation_rate,time_max,time_inc,NSEED, &
    & dble(endemic)/dble(NSEED) !,time_sum/dble(NSEED),time_sum_t/dble(NSEED)
close(25)
!open(26,file='mutability.dat',access='append',status='old')
!write(26,*) lambda, mutation_rate, dble(mutability)/dble(NSEED),time
!close(26)
! *******************************************************************
end program main

real*8 function poisson_time(lambda,idum)
implicit none
real*8 :: lambda, ran2
integer :: idum
poisson_time = -log(1d0-ran2(idum))/lambda
return
end function poisson_time

real*8 function delta(time,tau,delta_max)
implicit none
real*8, intent(in) :: time
real*8 :: time_int
real*8, intent(in)  :: tau,delta_max
real*8 :: pi=4d0*atan(1d0)
time_int=time+tau/4d0
delta=sign(1d0,sin(2d0*pi*time_int/tau))*tanh( (mod(time_int,tau/2d0)-pi/2d0)*pi )
!delta= 10d0*(delta+11d0/9d0)/(20d0/9d0)
!delta=(delta_max-1d0)/2d0 * cos(2d0*pi*time_int/tau) + (delta_max+1d0)/2d0
!delta=sign(1d0,sin(2d0*pi*time_int/tau))
delta=(delta_max-1d0)/2d0 * delta + (delta_max+1d0)/2d0
return
end function delta


FUNCTION ran2(idum)
INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
DOUBLE PRECISION ran2,AM,EPS,RNMX
PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
       &   IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211, &
       &   IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
!Long period (> 2 x 10 18 ) random number generator of L'Ecuyer with Bays-Durham shuffle
!and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
!of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
!alter idum between successive deviates in a sequence. RNMX should approximate the largest
!floating value that is less than 1.
INTEGER idum2,j,k,iv(NTAB),iy
SAVE iv,iy,idum2
DATA idum2/123456789/, iv/NTAB*0/, iy/0/

if (idum.le.0) then               !Initialize.
    idum=max(-idum,1)             !Be sure to prevent idum = 0.
    idum2=idum
    do j=NTAB+8,1,-1           !Load the shuffle table (after 8 warm-ups).
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        if (j.le.NTAB) iv(j)=idum
    enddo
    iy=iv(1)
endif

k=idum/IQ1                        !Start here when not initializing.
idum=IA1*(idum-k*IQ1)-k*IR1       !Compute idum=mod(IA1*idum,IM1) without over-
if (idum.lt.0) idum=idum+IM1      !flows by Schrage's method.
k=idum2/IQ2
idum2=IA2*(idum2-k*IQ2)-k*IR2     !Compute idum2=mod(IA2*idum2,IM2) likewise.
if (idum2.lt.0) idum2=idum2+IM2
j=1+iy/NDIV                       !Will be in the range 1:NTAB.
iy=iv(j)-idum2                    !Here idum is shuffled, idum and idum2 are com-
iv(j)=idum                        !bined to generate output.
if(iy.lt.1)iy=iy+IMM1
ran2=dmin1(AM*dble(iy),RNMX)      !Because users don't expect endpoint values.
return
END
