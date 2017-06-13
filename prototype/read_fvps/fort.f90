program main

real(kind=4)  :: rpar(1:10)
real(kind=4)  :: vpos(0:7)
integer :: ipar(1:10)
integer :: ip



open(10,file='out00.bin',form='unformatted')
read(10)ipar
read(10)rpar
read(10)vpos
write(*,*)vpos





close(10)


end program main