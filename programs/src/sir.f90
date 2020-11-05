! Susceptible (0), Infected(1), Recovered (2) model
! Network in edge list format should first be converted to arrays using read_zero.f90
! Two diferent "modes":
! 1) If x > 1 use different seeds to check prevalence(lambda)
! 2) If x = 1 use different seeds to check probability that virus survives
program main
! random numbers from ran2
implicit none
integer :: SEED, NSEED=1, SEED0=1, idum
integer :: N, E  ! number of nodes, edges
! Adjacency vector, initial and final pointers, state
integer, allocatable :: A(:), Pi(:), Pf(:), state(:)
! pointer to active_edges
integer, allocatable :: B(:)
integer, allocatable :: infected(:), recovered(:), active_edges(:,:)
integer :: i,iii,j,jjj,k ! indices for loops
integer :: n_infected, n_recovered, n_active_edges
integer :: node_n, index_n, index_n_B
! initial (microscopic) number of nodes to infect
integer :: x
! number of viruses (seeds) that get past time=1
integer :: survival
! probability of infection
real*8 :: lambda
real*8 :: lambda_sum, prob_heal
real*8 :: poisson_time, time
real*8 :: M_var, M_var_1, S_var, mean=0d0
real*8 :: ran2, random
NAMELIST/DADES/x,lambda,SEED0,NSEED
! ------------------------------------------------------------------
! Read input file
open(unit=11,file='../inp/sir.inp',status='old')
read(11,DADES)
close(11)
lambda = lambda/100d0
! ------------------------------------------------------------------
! Read array containing the network
open(20,file='array.dat')
read(20,*) N
read(20,*) E

allocate (Pi(N),Pf(N),A(2*E),B(2*E),state(N))
allocate (infected(N),recovered(N),active_edges(E,2))

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
n_recovered = 0
n_active_edges = 0
!B = 0
!active_edges = 0
! ===================
jjj = 0
loop0: do while (jjj.lt.x)
! ===================
jjj = jjj + 1
node_n = int(ran2(idum)*N) + 1
if(state(node_n).eq.1) cycle loop0
state(node_n) = 1
n_infected = n_infected + 1
infected(n_infected) = node_n

! Check all neighbors through A
do j=Pi(node_n),Pf(node_n)
    ! neighbor number is the value of A(j)
    index_n = A(j)
    ! if neighbor is active delete active edge
    if (state(index_n).eq.1) then
        ! move last active edge to the position that is deleted
        index_n_B = B(j)
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
    ! if neighbor is susceptible create new active edge
    elseif (state(index_n).eq.0) then
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
!print*, time, dble(n_infected)/dble(N), n_infected, n_recovered
! ===================
!loop_time: do while (n_active_edges.gt.0)
loop_time: do while (n_infected.gt.0)
! ===================
lambda_sum = dble(n_infected)+dble(n_active_edges)*lambda
time = time + poisson_time(lambda_sum,idum)
prob_heal = dble(n_infected)/lambda_sum
if (ran2(idum).lt. prob_heal) then
! heal
    index_n = int(ran2(idum)*n_infected) + 1
    node_n = infected(index_n)
    state(node_n) = 2
    n_recovered = n_recovered + 1
    recovered(n_recovered) = node_n
    infected(index_n) = infected(n_infected)
    n_infected = n_infected - 1
    ! active edges
    do j=Pi(node_n),Pf(node_n)
        index_n = A(j)
        if (state(index_n).eq.0) then
            index_n_B = B(j)
            if (index_n_B.gt.n_active_edges) then
                print*, 'ERROR'
                print*, index_n_B, n_active_edges, SEED
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
! infect
    index_n = int(ran2(idum)*n_active_edges) + 1
    node_n = active_edges(index_n,2)
    state(node_n) = 1
    n_infected = n_infected + 1
    infected(n_infected) = node_n

    do j=Pi(node_n),Pf(node_n)
        index_n = A(j)
        ! if neighbor is active delete active edge
        if (state(index_n).eq.1) then
            ! move last active edge to the position that is deleted
            index_n_B = B(j)
            if (index_n_B.gt.n_active_edges) then
                print*, 'ERROR'
                print*, index_n_B, n_active_edges
                stop
            endif
            active_edges(index_n_B,:) = active_edges(n_active_edges,:)
            ! change values in B
            i = active_edges(n_active_edges,1)
            iii = active_edges(n_active_edges,2)
            do k=Pi(i),Pf(i)
                if(A(k).eq.iii) B(k) = index_n_B
            enddo
            do k=Pi(iii),Pf(iii)
                if(A(k).eq.i) B(k) = index_n_B
            enddo
            n_active_edges = n_active_edges - 1
        elseif (state(index_n).eq.0) then
            n_active_edges = n_active_edges + 1
            active_edges(n_active_edges,:) = (/node_n, index_n /)
            ! Write number of active edge in B for both nodes
            B(j) = n_active_edges
            do k=Pi(index_n),Pf(index_n)
                if(A(k).eq.node_n) B(k) = n_active_edges
            enddo
        endif
    enddo
endif
!print*, time, dble(n_infected)/dble(N), n_infected, n_recovered
! ===================
enddo loop_time
! ===================
if (x.gt.1) then
    jjj = 0
    n_recovered = n_recovered + jjj
    mean = mean + dble(n_recovered)/dble(N)
    if (SEED.EQ.SEED0) then
        M_var = dble(n_recovered)/dble(N)
        S_var = 0d0
        cycle loop_seed
    endif
    M_var_1 = M_var
    M_var = M_var_1 + ( dble(n_recovered)/dble(N)-M_var_1 ) / dble(SEED-SEED0)
    S_var = S_var + (dble(n_recovered)/dble(N)-M_var_1)*(dble(n_recovered)/dble(N)-M_var)
elseif (x.eq.1) then
    !if (time.gt.10d0)  survival = survival + 1
    if (n_recovered.gt.100)  survival = survival + 1
endif
! ===================
enddo loop_seed
! ===================
! end time loop and seed loop
if (x.gt.1) then
    open(25,file='recovered.dat',access='append',status='old')
    write(25,*) lambda, mean/dble(NSEED), sqrt( S_var/dble(NSEED-1) )
    close(25)
elseif (x.eq.1) then
    open(25,file='survival.dat',access='append',status='old')
    write(25,*) lambda, dble(survival)/dble(NSEED)
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
