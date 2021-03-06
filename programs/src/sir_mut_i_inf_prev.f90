! Susceptible (0), Infected(1) model with mutations
! Network in edge list format should first be converted to arrays using
! either read_zero.f90 or read_one.f90
! saves the prevalence at the time of infection
program main
! random numbers from ran2
USE SUBROUTINES
implicit none
integer :: SEED, NSEED, SEED0, idum
integer :: N, E  ! number of nodes, edges
! Adjacency vector, initial and final pointers, state
integer, allocatable :: A(:), Pi(:), Pf(:), state(:)
! pointer to active_edges
integer, allocatable :: B(:)
integer, allocatable :: infected(:), active_edges(:,:)
integer :: i,iii,j,jjj ! indices for loops
integer :: n_infected, n_active_edges
integer, allocatable :: n_infected_i(:), n_recovered_i(:)
integer :: node_n, index_n, index_inf, index_act
! initial (microscopic) number of nodes to infect
integer :: x, state_old
integer, allocatable :: immune(:,:)
integer :: immune_length=5, infected_i_length=26, n_mutation
! probability of infection
real*8 :: lambda, lambda_sum, prob_heal, prob_infect
real*8 :: mutation_rate !, prob_mutate
real*8 :: poisson_time, time, time_max
! time of mutation for each new mutation and prevalence at the time
real*8, allocatable :: mutation_time_prevalence(:,:)
real*8, allocatable :: infection_prevalence(:)
real*8 :: ran2, random
character*100 :: name
NAMELIST/DADES/x,lambda,mutation_rate,SEED0,NSEED,time_max
! ---------------------------------------------------------------
! Read input file
open(unit=11,file='../inp/sir_mut_i.inp',status='old')
read(11,DADES)
close(11)
lambda = lambda/100d0
mutation_rate = 0d0
! ---------------------------------------------------------------
! Read array containing the network
open(20,file='array.dat')
read(20,*) N
read(20,*) E

allocate (Pi(N),Pf(N),A(2*E),B(2*E),state(N))
allocate (infected(N),active_edges(E,2))
! allocate(recovered(N))
allocate(immune(N,immune_length))
allocate(n_infected_i(infected_i_length),n_recovered_i(infected_i_length))
allocate(mutation_time_prevalence(infected_i_length,2))
allocate(infection_prevalence(N))

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

! ------------------------------------------------------------
! seed loop
open(31,file='infection_prevalence.dat')
loop_seed: do SEED = SEED0, SEED0+NSEED-1
    idum = -SEED
!-------
! Initialize population
state=0
n_infected = 0
n_infected_i = 0
n_recovered_i = 0
!mutation_time_prevalence = 0
infection_prevalence= 2d0 ! for the ones that don't get infected
n_active_edges = 0
n_mutation = 1
immune = 0
! ===================
loop0: do jjj=1,x
! ===================
node_n = int(ran2(idum)*N) + 1
if(state(node_n).eq.1) cycle loop0
state(node_n) = 1
n_infected = n_infected + 1
n_infected_i(n_mutation) = n_infected_i(n_mutation) + 1
infected(n_infected) = node_n

! Check all neighbors through A
do j=Pi(node_n),Pf(node_n)
    ! neighbor number is the value of A(j)
    index_n = A(j)
    ! if neighbor is active delete active edge
    if (state(index_n).eq.1) then
        call change_B_1(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E)
    ! if neighbor is susceptible create new active edge
    elseif (state(index_n).eq.0) then
        call change_B_2(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E,node_n,index_n)
    endif
enddo
! ===================
enddo loop0
! ===================

! *******************************************************************
!  _   _                  _
! | |_(_)_ __ ___   ___  | | ___   ___  _ __
! | __| | '_ ` _ \ / _ \ | |/ _ \ / _ \| '_ \
! | |_| | | | | | |  __/ | | (_) | (_) | |_) |
!  \__|_|_| |_| |_|\___| |_|\___/ \___/| .__/
!                                      |_|
time = 0d0
!print*, '#  time, prevalence, number of mutations, immune_length'
!print*, time, dble(n_infected)/dble(N), n_mutation, immune_length, dble(n_infected_i)/dble(N)
! ===================
loop_time: do while ( (n_infected.gt.0) .and. (time.lt.time_max) )
! ===================
lambda_sum = dble(n_infected)+dble(n_active_edges)*lambda+dble(n_infected)*mutation_rate
time = time + poisson_time(lambda_sum,idum)
prob_heal = dble(n_infected)/lambda_sum
!prob_mutate = dble(n_infected)*mutation_rate/lambda_sum
prob_infect = dble(n_active_edges)*lambda/lambda_sum
random = ran2(idum)
if (random.lt. prob_heal) then
! heal
    index_inf = int(ran2(idum)*n_infected) + 1
    node_n = infected(index_inf)
    call immunize(immune,state,node_n,N,immune_length)
    state_old = state(node_n)
    state(node_n) = 0
    infected(index_inf) = infected(n_infected)
    n_infected = n_infected - 1
    n_infected_i(state_old) = n_infected_i(state_old) - 1
    n_recovered_i(state_old) = n_recovered_i(state_old) + 1
    neighbors_h: do j=Pi(node_n),Pf(node_n)
        index_n = A(j)
        if (state(index_n).eq.0) then
            ! check if S neighbor was immune to the mutation of the
            ! now healed node
            do i=1,immune_length
                if(immune(index_n,i).eq.0) exit
                if(immune(index_n,i).ne.state_old) cycle
                cycle neighbors_h
            enddo
        ! if neighbor was S remove active edge
            call change_B_1(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E)
        else
        ! if neighbor is I with a mutation that the node is not
        ! immune to create active edge
            do i=1,immune_length
                if(immune(node_n,i).eq.0) exit
                if(immune(node_n,i).ne.state(index_n)) cycle
                cycle neighbors_h
            enddo
            call change_B_3(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E,node_n,index_n)
        endif
    enddo neighbors_h

elseif ( random .lt. (prob_heal+prob_infect) ) then
! infect
    index_act = int(ran2(idum)*n_active_edges) + 1
    ! node to be infected
    node_n = active_edges(index_act,2)
    state(node_n) = state(active_edges(index_act,1))
    infection_prevalence(node_n) = dble(n_infected)/dble(N)
    n_infected = n_infected + 1
    n_infected_i(state(node_n)) = n_infected_i(state(node_n)) + 1
    infected(n_infected) = node_n

    neighbors_i: do j=Pi(node_n),Pf(node_n)
        index_n = A(j)
        ! if neighbor is active delete active edge
        if (state(index_n).ne.0) then
            ! check if was immune to the mutation of the infected
            ! neighbor node
            do i=1,immune_length
                if(immune(node_n,i).eq.0) exit
                if(immune(node_n,i).ne.state(index_n)) cycle
                cycle neighbors_i
            enddo
            ! move last active edge to the position that is deleted
            call change_B_1(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E)
        else
            do i=1,immune_length
                if(immune(index_n,i).eq.0) exit
                if(immune(index_n,i).ne.state(node_n)) cycle
                cycle neighbors_i
            enddo
            call change_B_2(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E,node_n,index_n)
        endif
    enddo neighbors_i
else
!mutate
    index_inf = int(ran2(idum)*n_infected) + 1
    node_n = infected(index_inf)
    ! immunize to current mutation before mutating
    call immunize(immune,state,node_n,N,immune_length)
    n_mutation = n_mutation + 1
    ! time and prevalence at time of mutation
    if (n_mutation.gt.infected_i_length) then
        !call enlarge_infected_i(n_infected_i,n_recovered_i,mutation_time_prevalence,infected_i_length)
    endif
    state_old = state(node_n)
    state(node_n) = n_mutation
    n_infected_i(state_old) = n_infected_i(state_old) - 1
    n_recovered_i(state_old) = n_recovered_i(state_old) + 1
    n_infected_i(n_mutation) = n_infected_i(n_mutation) + 1
    !mutation_time_prevalence(n_mutation,:) = (/ time, dble(n_infected)/dble(N) /)
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
            endif
        endif
    enddo neighbors_m
endif
!print*, time, dble(n_infected)/dble(N), n_mutation, immune_length, dble(n_infected_i)/dble(N)
! ===================
enddo loop_time
do i=1,N
    write(31,*) i,infection_prevalence(i),seed
enddo
enddo loop_seed
! ===================
close(31)
write(name,'("./histo_inf_prev > histo_inf_prev_",I3.3,".dat")') int(lambda*100d0)
call execute_command_line(name)
!write(name,'("mv infection_prevalence.dat infection_prevalence_",I3.3,".dat")') int(lambda*100d0)
!call execute_command_line(name)
call execute_command_line("rm infection_prevalence.dat")
! *******************************************************************
end program main

real*8 function poisson_time(lambda,idum)
implicit none
real*8 :: lambda, ran2
integer :: idum
poisson_time = -log(1d0-ran2(idum))/lambda
return
end function poisson_time

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
