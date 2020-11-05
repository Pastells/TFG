implicit none
real*8 :: a,b,c,d,d1
real*8 :: integral(16), integral1(16)
integer :: i,j
open(20,file='survival.dat')
do i=1,4
    integral(i) = 0d0
    integral1(i) = 0d0
    d = 0d0
    do j=1,15
        d1 = d
        read(20,*) a,b,c,d
        integral(i) = integral(i) + c*d
        integral1(i) = integral1(i) + c*d1
    enddo
    print*, a,integral(i),integral1(i)
enddo
end
