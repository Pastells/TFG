integer :: histo(20)=0 ! bin every 0.05
real*8 :: a,b,n_out,N
integer :: i
open(unit=15,file='infection_prevalence.dat')
do while(1.eq.1)
    read(15,*, END=200) a,b
    i = int(b/0.05d0)+1
    if(i.gt.20) then
        n_out = n_out + 1
        cycle
    endif
    histo(i) = histo(i) + 1
enddo
200 close(15)
N=n_out+sum(histo)
do i=1,20
    print*, dble(i-1)*0.05,histo(i)/N
enddo
print*, '#, N,non_inf, R_inf'
print*, '#',N,n_out/N,sum(histo)/N
end
