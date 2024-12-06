
!            ********************************************************************************
!            * "open file .F90" is a little program for teaching the way of open a file     * 
!            *  and save the results of the program in it.                                  * 
!            ********************************************************************************


program openfile

implicit none

!**********************************************************************************************************************
!                                       Determine variables
!**********************************************************************************************************************

integer             i        ,n

real*8              a        ,b       ,x         ,y              ,z          ,p         ,q         ,pi         ,deltax                                    

character(len=30)   aa       ,bb      ,nam       ,filename1      ,filename2  ,filename3
   

!**********************************************************************************************************************
!                                         Zero to variables
!**********************************************************************************************************************
 
                    i=0      ;n=0

                    a=0.     ;b=0.    ;x=0.      ;y=0.           ;z=0.       ;p=0.      ;q=0.     ,deltax=0.
    
!********************************************************************************************************************
!                                             Inputs		  
!********************************************************************************************************************

write(*,*)
write(*,*)

write(*,*)
write(*,'(a,\)') '   Enter  a = '
 read(*,*)a
write(*,'(a,\)') '   Again    = '
 read(*,*)aa

write(*,*)
write(*,'(a,\)') '   Enter  b = '
 read(*,*)b
write(*,'(a,\)') '   Again    = '
 read(*,*)bb

!**********************************************************************************************************************
!                                 Determine  Filenames & Open files
!**********************************************************************************************************************

write(*,*)

nam=' a='//trim(aa)//' b='//trim(bb)

filename1=Trim(name)//' .txt'
filename2=Trim(name)//' .plt'
filename3=Trim(name)//' out.txt'

open (1,file=filename1)
open (2,file=filename2)
open (3,file=filename3)

write(*,*)'  Note : Results,will be saved in these files : '
write(*,*)'                                                ',filename1
write(*,*)'                                                ',filename2
write(*,*)'                                                ',filename3

!**********************************************************************************************************************
!                                        Determine  Constants
!**********************************************************************************************************************

    pi = 4*atan(1.) ! pi = 3.141592654  

     n = 100
deltax = (b-a)/n

!**********************************************************************************************************************
!                                      Main  &  printe results     
!**********************************************************************************************************************

write(1,*)
write(1,'(a,/)')' This File is the output of "open file .F90" program.'
write(1,'(a,/)')' Results Will be Ploted with "Tecplot" software.'
write(1,'(a,/)')'   x           y=sinx       x=Arcsiny '

write(2,*)
write(2,'(a,/)')' This File is the output of "open file .F90" program.'
write(2,'(a,/)')' Results Will be Ploted with "Tecplot" software.'
write(2,'(a,/)')'   x           y=sinx       x=Arcsiny '

do i=1,n
   x=a+i*deltax
   y=sin(x)
   
   write(1,'(f8.3,5x,f8.3)')x,y
   write(2,'(f8.3,5x,f8.3)')x,y

end do !x    
             
!--------------------------- read from file 1  then write  in to file 3
write(3,*)
write(3,'(a,/)')' This File is the output of "open file .F90" program.'
write(3,'(a,/)')' Results Will be Ploted with "Tecplot" software.'
write(3,'(a,/)')'   x           y=sinx       x=Arcsiny '

rewind(1)
read(1,'(6/)')

do i=1,n

   read(1,*)p,q
   z=asin(q)
   
!   backspace(1)
   write(3,'(f8.3,5x,f8.3,5x,f8.3)')p,q,z
!   read(1,*)

end do !x

write(1,*)
write(1,'(a,/)')' Note : We khow that "Arcsin[sin(x)]=sin[arcsin(x)]=x",but for x=2 and y=sin(2)=0.909 ,   '
write(1,'(a,/)')'        We have x=arcsin(0.909)=1.141 while we expected the resut would be 2. Why?'
write(2,*)
write(2,'(a,/)')' Note : We khow that "Arcsin[sin(x)]=sin[arcsin(x)]=x",but for x=2 and y=sin(2)=0.909 ,   '
write(2,'(a,/)')'        We have x=arcsin(0.909)=1.141 while we expected the resut would be 2. Why?'
write(3,*)
write(3,'(a,/)')' Note : We khow that "Arcsin[sin(x)]=sin[arcsin(x)]=x",but for x=2 and y=sin(2)=0.909 ,   '
write(3,'(a,/)')'        We have x=arcsin(0.909)=1.141 while we expected the resut would be 2. Why?'


!**********************************************************************************************************************
!                                      Close files & end program 
!**********************************************************************************************************************
close(1)
close(2)
close(3) 
write(*,'(a,/)')'_____________________________________________________________________'
end program openfile                  
!======================================================================================================================
         
