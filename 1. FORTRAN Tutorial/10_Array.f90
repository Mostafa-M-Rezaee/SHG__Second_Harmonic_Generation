
!            **********************************************************************
!            * "Array .F90" is an example of using arrays in a Fortran program.   * 
!            **********************************************************************

program Array

implicit none

!**********************************************************************************************************************
!                                       Variables Definition
!**********************************************************************************************************************

integer       i          ,j          ,k                                                                                 

real*8        a[allocatable](:)      ,b[allocatable](:,:)   ,c[allocatable](:,:,:)

!**********************************************************************************************************************
!                                        Arrays Allocattion 
!**********************************************************************************************************************

allocate( a(0:100) )     

allocate( b(1:10,0:20) )     

allocate( c(-10:+10,0:10,10:20) )     

!**********************************************************************************************************************
!                                     Giving Zero to Arrays
!********************************************************************************************************************** 

forall (i=0:100)
                                  a(i) = 0.        
end forall 

!------------------
forall (i=1:10,j=0:20)
                                  b(i,j) = 0.        
end forall 

!------------------
forall (i=-10:+10,j=0:10,k=10:20)
                                  c(i,j,k) = 0.        
end forall 

!**********************************************************************************************************************
!                                   Main Block of the Program     
!**********************************************************************************************************************

!------------------
do i=0,+100 

   a(i) = i

end do !i

!------------------
do i=1,+10 

    do j=0,20
           
	   b(i,j) = i+j

    end do !j
end do !i

!------------------
do i=-10,+10 

	do j=0,10
           
	    do k=10,20
      
	       c(i,j,k) = i+j+k
		    
	    end do !k
    end do !j
end do !i

!**********************************************************************************************************************
!                                        Printing Results     
!**********************************************************************************************************************

do i=1,100
   
   write(*,'(2x,i5,5x,f5.1)') i,a(i)

end do !i      						   

!**********************************************************************************************************************
!                                      Closing Files and Ending the Program 
!**********************************************************************************************************************

end program Array                     

!======================================================================================================================
         
 