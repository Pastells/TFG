! Susceptible (0), Infected(n), Removed (-1)  model
! Network in edge list format should first be converted to arrays using read_zero.f90
! Two diferent "modes":
! 1) If x > 1 use different seeds to check prevalence(lambda)
! 2) If x = 1 use different seeds to check probability that virus survives
program main
! random numbers from ran2
USE SUBROUTINES
implicit none
integer :: SEED, NSEED=1, SEED0=1, idum
integer :: N, E  ! number of nodes, edges
! Adjacency vector, initial and final pointers, state
integer, allocatable :: A(:), Pi(:), Pf(:), state(:)
! pointer to active_edges
integer, allocatable :: B(:)
integer, allocatable :: infected(:),active_edges(:,:),removed(:),susceptible(:)
integer :: i,iii,j,jjj,k ! indices for loops
integer :: n_infected,n_active_edges,n_removed,n_sus
integer :: node_n, index_n, index_n_B
! initial (microscopic) number of nodes to infect
integer :: x, state_old
integer, allocatable :: immune(:,:)
integer :: immune_length=5, n_mutation !,mutability
integer, allocatable :: endemic(:)
! number of viruses (seeds) that get past time=1
integer :: survival
! rate of infection
real*8 :: lambda,lambda_sum
! rate of births/deaths
real*8 :: mu
! % of removed population
real*8 :: immune_0
real*8 :: prob_heal,prob_infect,prob_mutate,prob_birth,prob_death_I
real*8 :: poisson_time, time, time_max, time_write, time_inc
real*8 :: M_var, M_var_1, S_var, mean=0d0
real*8 :: M_varp, M_varp_1, S_varp, meanp=0d0
real*8 :: ran2, random, time_sum, time_sum_t
NAMELIST/DADES/x,lambda,mu,immune_0,mutation_rate,SEED0,NSEED,time_max,time_inc
! ------------------------------------------------------------------
! Read input file
open(unit=11,file='../inp/sir_dynamic.inp',status='old')
read(11,DADES)
close(11)
lambda = lambda/100d0
mu = mu/100d0
immune_0 = immune_0/100d0
mutation_rate = mutation_rate/1d5
! ------------------------------------------------------------------
! Read array containing the network
open(20,file='array.dat')
read(20,*) N
read(20,*) E

allocate (Pi(N),Pf(N),A(2*E),B(2*E),state(N))
allocate (infected(N),active_edges(E,2))
!allocate(susceptible(N),removed(int(0.1d0*dble(N))))
allocate(susceptible(N),removed(N))

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

survival = 0
! ------------------------------------------------------------
! seed loop
loop_seed: do SEED = SEED0, SEED0+NSEED-1
    !print*, seed
    idum = -seed
!-------
! Initialize population
state=0
n_infected = 0
n_active_edges = 0
n_mutation = 1
immune = 0
n_removed = 0
!B = 0
!active_edges = 0
! susceptible list
n_sus = N
do jjj=1,N
    susceptible(jjj) = jjj
enddo

! remove 10% of the population
!print*, int(0.1*N), int(0.1d0*N), int(0.1d0*dble(N)),N,n_removed
do jjj=1,int(immune_0*N)
    index_n = int(ran2(idum)*n_sus) + 1
    node_n = susceptible(index_n)
    state(node_n) = -1
    susceptible(index_n) = susceptible(n_sus)
    n_sus = n_sus - 1
    n_removed = n_removed + 1
    removed(n_removed) = node_n
enddo
! infect microscopic number of the population "x"
! ===================
jjj = 1
loop0: do while (jjj.lt.x)
! ===================
node_n = int(ran2(idum)*N) + 1
if(state(node_n).eq.1) cycle loop0
if(state(node_n).eq.-1) cycle loop0
state(node_n) = 1
n_infected = n_infected + 1
infected(n_infected) = node_n

do j=1,n_sus
    if(susceptible(j).eq.node_n) then
        susceptible(j)=susceptible(n_sus)
        n_sus = n_sus - 1
        exit
    endif
enddo

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
jjj = jjj + 1
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
time_write = time_inc
jjj = 1
!print*, time, dble(n_infected)/dble(N*(1d0-immune_0)),n_infected/(1d0-immune_0),&
   !&   n_recovered/(1d0-immune_0), n_sus/(1d0-immune_0), n_removed/(1d0-immune_0)
! ===================
loop_time: do while ( (n_infected.gt.0) .and. (time.lt.time_max) )
! ===================
lambda_sum = dble(n_infected)*(1d0+mutation_rate)+dble(n_active_edges)*lambda+ &
    & mu*( (1d0-immune_0)*dble(N) + dble(n_infected) + dble(n_sus) )
!             births     +      deaths I    +  deaths S
time = time + poisson_time(lambda_sum,idum)
prob_heal = dble(n_infected)/lambda_sum
prob_infect = dble(n_active_edges)*lambda/lambda_sum
prob_mutate = dble(n_infected)*mutation_rate/lambda_sum
prob_birth = mu*(1d0-immune_0)*dble(N)/lambda_sum
prob_death_I = mu*dble(n_infected)/lambda_sum
random = ran2(idum)
if (random.lt.prob_heal) then
! heal
    index_inf = int(ran2(idum)*n_infected) + 1
    node_n = infected(index_inf)
    call immunize(immune,state,node_n,N,immune_length)
    state_old = state(node_n)
    state(node_n) = 0
    n_sus = n_sus + 1
    susceptible(n_sus) = node_n
    infected(index_inf) = infected(n_infected)
    n_infected = n_infected - 1
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
        elseif(state(index_n).ne.-1) then
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
    n_infected = n_infected + 1
    infected(n_infected) = node_n

    do j=1,n_sus
        if(susceptible(j).eq.node_n) then
            susceptible(j)=susceptible(n_sus)
            n_sus = n_sus - 1
            exit
        endif
    enddo

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
        elseif(state(index_n).ne.-1) then
            do i=1,immune_length
                if(immune(index_n,i).eq.0) exit
                if(immune(index_n,i).ne.state(node_n)) cycle
                cycle neighbors_i
            enddo
            call change_B_2(Pi,Pf,A,B,active_edges,n_active_edges,j,N,E,node_n,index_n)
        endif
    enddo neighbors_i

elseif ( random .lt. (prob_heal+prob_infect+prob_mutate) ) then
!mutate
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
            endif
        endif
    enddo neighbors_m
elseif ( random .lt. (prob_heal+prob_infect+prob_mutate+prob_birth) ) then
! birth
    if (n_removed.eq.0) then
        print*, 'REMOVED ERROR'
        STOP
    endif
    index_n = int(ran2(idum)*n_removed) + 1
    node_n = removed(index_n)
    state(node_n) = 0
    removed(index_n) = removed(n_removed)
    n_removed = n_removed - 1
    n_sus = n_sus + 1
    susceptible(n_sus) = node_n
    ! active edges
    do j=Pi(node_n),Pf(node_n)
        index_n = A(j)
        if (state(index_n).eq.1) then
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
        endif
    enddo
elseif ( random .lt. (prob_heal+prob_infect+prob_mutate+prob_birth+prob_death_I) ) then
! death I
    index_n = int(ran2(idum)*n_infected) + 1
    node_n = infected(index_n)
    state(node_n) = -1
    n_removed = n_removed + 1
    removed(n_removed) = node_n
    infected(index_n) = infected(n_infected)
    n_infected = n_infected - 1
    ! active edges
    do j=Pi(node_n),Pf(node_n)
        index_n = A(j)
        if (state(index_n).eq.0) then
            index_n_B = B(j)
            if (index_n_B.gt.n_active_edges) then
                print*, 'ERROR 3'
                print*, '#',index_n_B, n_active_edges,seed
                !stop
            endif
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
        endif
    enddo
else
! death S
    index_n = int(ran2(idum)*n_sus) + 1
    node_n = susceptible(index_n)
    state(node_n) = -1
    n_removed = n_removed + 1
    removed(n_removed) = node_n
    susceptible(index_n) = susceptible(n_sus)
    n_sus = n_sus - 1

    do j=Pi(node_n),Pf(node_n)
        index_n = A(j)
        ! if neighbor is active delete active edge
        if (state(index_n).eq.1) then
            ! move last active edge to the position that is deleted
            index_n_B = B(j)
            if (index_n_B.gt.n_active_edges) then
                print*, 'ERROR 4'
                print*, '#',index_n_B, n_active_edges,seed
                stop
            endif
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
        endif
    enddo
endif
print*, time, dble(n_infected)/dble(N*(1d0-immune_0)), n_infected, n_recovered, n_sus, n_removed
!print*, time, dble(n_infected)/dble(N*(1d0-immune_0)),n_infected/(1d0-immune_0),&
   !&   n_recovered/(1d0-immune_0), n_sus/(1d0-immune_0), n_removed/(1d0-immune_0)
! ===================
enddo loop_time
! ===================
! ===================
enddo loop_seed
! ===================
! end time loop and seed loop
if (x.gt.1) then
    open(25,file='recovered.dat',access='append',status='old')
    write(25,*) lambda,immune_0,mu,mean/dble(NSEED), sqrt( S_var/dble(NSEED-1) ), &
        & meanp/dble(NSEED), sqrt( S_varp/dble(NSEED-1) )
    close(25)
elseif (x.eq.1) then
    open(25,file='survival.dat',access='append',status='old')
    write(25,*) lambda,immune_0,mu, dble(survival)/dble(NSEED)
    close(25)
endif
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
