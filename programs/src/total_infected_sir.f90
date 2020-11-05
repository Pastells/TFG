! Find the total number of infected individuals in the sir
! model for both the homogenous and complex network case
! needs input file degree.dat with the degree distribution
real*8 :: lambda, k_av, k2_av
real*8 :: r_homo, r_inf_homo, r_com, r_inf_com, lambda_c
real*8, allocatable :: P(:)
integer :: i,k,n=0
integer :: lambda_0,lambda_f,lambda_step
NAMELIST/DADES/lambda_0,lambda_f,lambda_step
! ------------------------------------------------------------------
! Read input file
open(unit=11,file='../inp/total_infected_sir.inp',status='old')
read(11,DADES)
close(11)
! ------------------------------------------------------------------


open(10,file='degree.dat')
do while(1.eq.1)
   read(10,*, END=200) r,r
   n = n + 1
enddo
200 close(10)

allocate(P(n))
lambda_c = 0d0
k_av = 0d0
k2_av = 0d0
open(10,file='degree.dat')
do i=1,n
    read(10,*) k,P(k)
    k_av = k_av + dble(k)*P(k)
    k2_av = k2_av + dble(k**2)*P(k)
    !lambda_c = lambda_c + dble(k*(k-1))*P(k)
enddo
close(10)
!print*, k_av, k2_av
!print*, lambda_c, k_av/(k2_av-k_av)
!lambda_c = k_av/lambda_c
lambda_c = k_av/(k2_av-k_av)

open(5,file='total_inf.dat')
write(5,*) '#lambda,r_com,eps_c_com,r_homo,eps_homo'
do i=lambda_0,lambda_f,lambda_step
    lambda=dble(i)*0.01
    r_homo = r_inf_homo(lambda,k_av)
    r_com = r_inf_com(lambda,k_av,P,n,lambda_c)
    eps_c = lambda/(lambda-lambda_c)
    write(5,*) lambda,r_com,eps_c/r_com,r_homo,eps_c/r_homo
enddo
end

! #######################################################
! #######################################################

real*8 function r_inf_homo(lambda,k_av)
implicit none
real*8 :: lambda, k_av
real*8 :: diff, r_inf_homo_0

!k_av = k_av - 1d0
if(lambda.lt.1d0/k_av) then
    r_inf_homo = 0
    return
endif
r_inf_homo_0 = 1d0
diff = 1d0
do while(diff.gt.1e-8)
    r_inf_homo = 1d0-exp(-lambda*k_av*r_inf_homo_0)
    !print*, r_inf_homo
    diff=abs(r_inf_homo-r_inf_homo_0)
    r_inf_homo_0 = r_inf_homo
enddo
return
end function r_inf_homo


real*8 function r_inf_com(lambda,k_av,P,n,lambda_c)
implicit none
real*8 lambda, lambda_c, k_av
integer :: k,n
real*8 :: P(n)
real*8 :: phi,phi_inf

if(lambda.lt.lambda_c) then
    r_inf_com = 0
    return
endif
phi = phi_inf(lambda,k_av,P,n)
do k=1,n
    r_inf_com = r_inf_com + P(k)*(1d0-exp(-lambda*k*phi))
enddo
return
end function r_inf_com

real*8 function phi_inf(lambda,k_av,P,n)
implicit none
real*8 :: lambda, k_av
integer :: k,n
real*8 :: P(n)
real*8 :: diff, phi_inf_0,sum

phi_inf_0 = 1d0
diff = 1d0
do while(diff.gt.1e-8)
    sum = 0d0
    do k=1,n
        sum = sum + dble(k-1)*P(k)*exp(-lambda*k*phi_inf_0)
    enddo
    phi_inf = 1d0-1d0/k_av*(1d0+sum)
    !print*, phi_inf
    diff=abs(phi_inf-phi_inf_0)
    phi_inf_0 = phi_inf
enddo
return
end function phi_inf
