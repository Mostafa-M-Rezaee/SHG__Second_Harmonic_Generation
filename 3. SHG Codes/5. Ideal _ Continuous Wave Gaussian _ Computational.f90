
				  
!            ********************************************************************************
     
!            * Dr.Mohammad Sabaeian   , Departmant of Physics, Shahid Chamran University    * 

!            * Fateme Sedaghat        , Departmant of Physics, Shahid Chamran University    *

!            * Mostafa Mohammad Rezaee, Departmant of Physics, Shahid Chamran University    *

!            * f.sedaghat2010@yahoo.com                                                     *

!            * mostafa_mohammadrezaee@yahoo.com                                             *

!            * "couplied elec-temp-phase.F90" is aprogram to                                * 

!            *                                                                              * 

!            *             originally Written : 21/10/2011  by FS  ,MM                      *

!            *                   last revised : ---------                                   *

!            ********************************************************************************

program coupliedelectempphase

implicit none

!**********************************************************************************************************************
!                                       Determine variables
!**********************************************************************************************************************

integer        i                   ,j                     ,k                      ,f                      ,nz                 
                                                                                           
real*8         l                   ,p                     ,z                       ,c                      ,d               &                                                      		
              ,pi                  ,n1                    ,n2                      ,n3                     ,r1f             &                    
			  ,r2f                 ,r3f                                                                                     &                    
			  ,gama1               ,gama2                 ,gama3                   ,omegaf                                  &
			  ,deltaz              ,deltak                ,omega1                  ,epsilon                                 &
			  ,landa1              ,lenght                                                                                  &                                                                                     
                                                                                                                        

              ,elecp1[allocatable](:)                             ,elecm1[allocatable](:)                                   &    
			  ,elecp2[allocatable](:)                             ,elecm2[allocatable](:)                                   &    
              ,elecp3[allocatable](:)                             ,elecm3[allocatable](:)                                                                      
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         

    																											    
complex*16    Ii                                                                                                            &                                                                                            
              ,psip1[allocatable](:)                              ,psim1[allocatable](:)                                    &                                                                                                                                                     
              ,psip2[allocatable](:)                              ,psim2[allocatable](:)                                    &                                                                                      
              ,psip3[allocatable](:)                              ,psim3[allocatable](:)                                                                                                                             
                 

character*30       filenamez

!**********************************************************************************************************************
!                                         Zero to variables
!**********************************************************************************************************************
                        
                                                                             
              ;l = 0                   ;p = 0                     ;z = 0                      ;c = 0                    ;d = 0                                                                                                              		
              ;pi = 0                  ;n1 = 0                    ;n2 = 0                     ;n3 = 0                                   
			  ;gama1 = 0               ;gama2 = 0                 ;gama3 = 0                  ;epsilon = 0                             
			  ;deltaz = 0              ;deltak = 0                ;omega1 = 0                 ;omegaf = 0                              
			  ;landa1 = 0 
                
                   

!**********************************************************************************************************************
!                                             Inputs		  
!**********************************************************************************************************************

!write(*,'(/,2x,a,\)') 'Enter the fundemental power of one wave   : '
 !read(*,*)p

!write(*,'(/,2x,a,\)') 'Enter omegaf  : '                
 !read(*,*)omegaf
 
write(*,'(/,2x,a,\)') 'Enter the deltak value  : '
 read(*,*)deltak


!write(*,'(/,2x,a,\)') 'Enter gama1  : '
 !read(*,*)gama1

 
!write(*,'(/,2x,a,\)') 'Enter gama2  : '
 !read(*,*)gama2

 
!write(*,'(/,2x,a,\)') 'Enter gama3  : '
 !read(*,*)gama3
 

!**********************************************************************************************************************
!                                 Determination of Filenames & Open files
!**********************************************************************************************************************

filenamez='elec1p-z.txt'
open(1,file=filenamez)

filenamez='elec2p-z.txt'
open(2,file=filenamez)

filenamez='elec3p-z.txt'
open(3,file=filenamez)

filenamez='elec1,2,3.txt'  
open(4,file=filenamez) 

filenamez='elec1m-z.txt'
open(5,file=filenamez)

filenamez='elec2m-z.txt'
open(6,file=filenamez)

filenamez='elec3m-z.txt'
open(7,file=filenamez)


!**********************************************************************************************************************
!                                        Determination of  Constants
!**********************************************************************************************************************
        p = 40
        c = 3.e8   
       pi = 4*atan(1.)
       n1 = 1.8296
       n2 = 1.7466 
       n3 = 1.7881  
       nz = 1000
       Ii = (0,1) 
       gama1 = 0.5
       gama2 = 0.5
       gama3 = 4.
       landa1 = 1.064     !micro meter 
       lenght = 0.02
	   omegaf = 0.0001
       epsilon = 8.85e-12                        
 	         d = 7.3e-12*epsilon
       omega1 = 2*pi*c*1000000/landa1            ;l = (epsilon**3*c**3*pi*omegaf**2/((4*p)*omega1**2*d**2))**.5     
       deltaz = lenght/nz  
    
      	 
!**********************************************************************************************************************
!                                        printe Constants     
!**********************************************************************************************************************



!**********************************************************************************************************************
!                                         Allocate Arrys
!**********************************************************************************************************************


allocate(elecp1(0:nz))                  ;   allocate(elecm1(0:nz))                  
allocate(elecp2(0:nz))                  ;   allocate(elecm2(0:nz))       
allocate(elecp3(0:nz))                  ;   allocate(elecm3(0:nz))    
                    
allocate(psip1(0:nz))                   ;   allocate(psim1(0:nz))                           
allocate(psip2(0:nz))                   ;   allocate(psim2(0:nz))    
allocate(psip3(0:nz))                   ;   allocate(psim3(0:nz))      

          
   
!**********************************************************************************************************************
!                                         Zero to Arrays
!**********************************************************************************************************************

forall (k=0:nz)
       elecp1(k) = 0                    ;   elecm1(k) = 0  
	   elecp2(k) = 0                    ;   elecm2(k) = 0                
	   elecp3(k) = 0                    ;   elecm3(k) = 0
	   psip1(k)  = 0                    ;   psim1(k)  = 0
	   psip2(k)  = 0                    ;   psim2(k)  = 0
	   psip3(k)  = 0                    ;   psim3(k)  = 0
     
end forall !k     
                                                                                

!**********************************************************************************************************************
!                                              Main     
!**********************************************************************************************************************

!psip1(0)=(p/(2*pi*n1*c*epsilon*(omegaf**2)))**0.5
!psip2(0)=(p/(2*pi*n2*c*epsilon*(omegaf**2)))**0.5
!psip3(0)=0.
psip1(0)=1
psip2(0)=1
psip3(0)=0	 		       		 
!------------------------------------------------------------------------ 
!eq. for  center of crystal  (in centences exp(-2r^2/wf^2) r=0)    
do k=0,nz-1 
   z=deltaz*k
   psip1(k+1) = psip1(k)*(1-gama1*deltaz/2) + deltaz/(l*(n1*n2*n3)**.5)*Ii*psip3(k)*conjg(psip2(k))*exp(-Ii*deltak*z)
    
   psip2(k+1) = psip2(k)*(1-gama2*deltaz/2) + deltaz/(l*(n1*n2*n3)**.5)*Ii*psip3(k)*conjg(psip1(k))*exp(-Ii*deltak*z)
      
   psip3(k+1) = psip3(k)*(1-gama3*deltaz/2) + deltaz/(l*(n1*n2*n3)**.5)*Ii*psip1(k)*psip2(k)*exp(Ii*deltak*z)
   
end do   	
!----------------------
r1f = 0.99
r2f = 0.99
r3f = 0.99

psim1(nz) = r1f*psip1(nz)
psim2(nz) = r2f*psip2(nz)
psim3(nz) = r3f*psip3(nz)

!-----------------------
do k=nz,1,-1
   z=deltaz*k
   psim1(k-1) = psim1(k)*(1-gama1*deltaz/2) + deltaz/(l*(n1*n2*n3)**.5)*Ii*psim3(k)*conjg(psim2(k))*exp(Ii*deltak*z)
   psim2(k-1) = psim2(k)*(1-gama2*deltaz/2) + deltaz/(l*(n1*n2*n3)**.5)*Ii*psim3(k)*conjg(psim1(k))*exp(Ii*deltak*z)
   psim3(k-1) = psim3(k)*(1-gama3*deltaz/2) + deltaz/(l*(n1*n2*n3)**.5)*Ii*psim1(k)*psim2(k)*exp(-Ii*deltak*z)   	
end do	  	  
!**************************************************************************************************************************
!                                          printe results     
!**************************************************************************************************************************
do k=0,nz
   elecp1(k)=psip1(k)*conjg(psip1(k))                     !*exp(0)=1
   elecp2(k)=psip2(k)*conjg(psip2(k))
   elecp3(k)=psip3(k)*conjg(psip3(k)) 
end do
do k=nz,0,-1
   elecm1(k)=psim1(k)*conjg(psim1(k))                     !*exp(0)=1
   elecm2(k)=psim2(k)*conjg(psim2(k))
   elecm3(k)=psim3(k)*conjg(psim3(k)) 
end do!k
!-------------------------------
do k=0,nz
   z=k*deltaz

   write(1,'(f13.8 ,5X, F20.12)' ) z, elecp1(k)*100
   write(2,'(f13.8 ,5x, f20.12)' ) z , elecp2(k)*100 
   write(3,'(f13.8 ,5x, f20.12)' ) z, elecp3(k)*100
   
   write(5,'(F13.8 ,5X ,f20.12)' ) z ,  elecm1(k)*100 
   write(6,'(f13.8 ,5x,f20.12)' )  z ,  elecm2(k)*100 
   write(7,'(f13.8 ,5x,f20.12)' )  z ,  elecm3(k)*100

 !  write(1,'(F20.12)' )  elecp1(k)*100
  ! write(2,'(f20.12)' )  elecp2(k)*100 
 !  write(3,'( f20.12)' )  elecp3(k)*100
   
  ! write(5,'(f20.12)' )  elecm1(k)*100 
  ! write(6,'(f20.12)' )   elecm2(k)*100 
  ! write(7,'(f20.12)' )   elecm3(k)*100


write(4,'(f13.8 ,5x, f20.12 ,5x, f20.12 ,5x, f20.12 ,5x, f20.12 ,5x, f20.12 ,5x, f20.12)')                             &
 
                 z ,  elecp1(k)*100 , elecp2(k)*100 , elecp3(k)*100 , elecm1(k)*100 , elecm2(k)*100 , elecm3(k)*100
 


end do !k


   						   

!**********************************************************************************************************************
!                                      Close files & end program 
!**********************************************************************************************************************

close(1)
end program coupliedelectempphase 
                 
!======================================================================================================================
                                                                                                                                            




