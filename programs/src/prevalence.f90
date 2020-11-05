real*8 :: histo(16,2) ! bin every 0.05
real*8 :: a,c,d
integer :: i
open(unit=15,file='mutation_time_prevalence.dat')
do while(1.eq.1)
   read(15,*, END=200) a,a,c,d
   if (d.gt.100) then
       histo(int(c/0.025)+1,1) = histo(int(c/0.025)+1,1) + 1
    else
       histo(int(c/0.025)+1,2) = histo(int(c/0.025)+1,2) + 1
    endif
enddo
200 close(15)
!open(unit=16,file='prevalence.dat')
do i=1,16
    !write(16,*) dble(i)*0.025,histo(i,:)
    print*, dble(i)*0.025,histo(i,:)
enddo
end
