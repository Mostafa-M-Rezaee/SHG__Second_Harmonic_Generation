
!            ********************************************************************************
!            * "Readable Code Structure .F90" is an example of a readable Fortran program.  * 
!            ********************************************************************************

program Temp_G_PW

implicit none

!**********************************************************************************************************************
!                                       Variables Definition
!**********************************************************************************************************************

integer       i          ,j          ,k          ,l                                                       &
             ,nt         ,nr         ,nz         ,Np                                                                       

real*8        t          ,z          ,E          ,h          ,r                                           &                                                                                    
             ,T0         ,ka         ,pi         ,Cp         ,tp                                          &
			 ,roh        ,aa1        ,aa2        ,aa3        ,aa4                                         &
             ,freq       ,gama                                                                            & 
			 ,timet      ,sigma                                                                           &
			 ,omegaf     ,length     ,deltar     ,deltaz     ,deltat                                      &
			 ,radius                                                                                      &
			 ,epsilong   ,tbetween                                                                        &
			 ,stability                                                                                   &
			  
			 ,temperature[allocatable](:,:,:)  

character*30  filenameTt   ,filenameTr   ,filenameTz     ,Npf       ,freqf     ,tpf       ,EE             &
             ,plot_extention

                                                     

!**********************************************************************************************************************
!                                    Giving Zero to variables
!**********************************************************************************************************************

                    i = 0           ;j = 0           ;k = 0          ;l = 0            
                   nt = 0          ;nr = 0          ;nz = 0         ;Np = 0

                    t = 0.          ;z = 0.          ;E = 0.         ;h = 0.                  ;r = 0.                                                                                                          
                   T0 = 0.         ;ka = 0.         ;pi = 0.        ;Cp = 0.                 ;tp = 0.                        
                  roh = 0.        ;aa1 = 0.        ;aa2 = 0.       ;aa3 = 0.                ;aa4 = 0.                                                                        
                 freq = 0.       ;gama = 0.
			    timet = 0.      ;sigma = 0.                                                                       
			   omegaf = 0.     ;length = 0.     ;deltar = 0.    ;deltaz = 0.             ;deltat = 0.                                           
			   radius = 0.
			 epsilong = 0.   ;tbetween = 0.                                   
	        stability = 0.   

!**********************************************************************************************************************
!                                             Inputs		  
!**********************************************************************************************************************

! Note: 
!     This code lets the user enter values twice: once numerically (for calculations) 
!     and once as a string (for filenames or labels).  
!     For example, `E` is number,while `EE` store the same values as strings.  
!     This dual input ensures accurate calculations and meaningful file naming.

write(*,'(/,2x,a,\)') '            Enter the Energy value  : '
read(*,*) E

write(*,'(/,2x,a,\)') '                             Again  : '
read(*,*) EE

write(*,'(/,2x,a,\)') '         Enter the frequency value  : '
read(*,*) freq

write(*,'(/,2x,a,\)') '                             Again  : '
read(*,*) freqf

write(*,'(/,2x,a,\)') '        Enter the Number of Pulses  : '
read(*,*) Np

write(*,'(/,2x,a,\)') '                             Again  : '
read(*,*) Npf

write(*,'(/,2x,a,\)') '                      Enter the tp  : '
read(*,*) tp

write(*,'(/,2x,a,\)') '                             Again  : '
read(*,*) tpf

!**********************************************************************************************************************
!                          Determination of Filenames and Opening files
!**********************************************************************************************************************

! Note:
!      To achieve both efficiency and clarity in managing output data,
!      below, we generate filenames based on input information.

plot_extention = '.plt'

filenameTt = 'E_'//trim(EE)//'_f_'//trim(freqf)//'_Np_'//trim(Npf)//'_tp_'//trim(tpf)//'T_t'//plot_extention 
open(1,file=filenameTt)
!write(1,'(/,a,/)')    '                     t                               temperature'

filenameTr = 'E_'//trim(EE)//'_f_'//trim(freqf)//'_Np_'//trim(Npf)//'_tp_'//trim(tpf)//'_T_r'//plot_extention 
open(2,file=filenameTr)
!write(2,'(/,a,/)')    '                     r                               temperature'

filenameTz = 'E_'//trim(EE)//'_f_'//trim(freqf)//'_Np_'//trim(Npf)//'_tp_'//trim(tpf)//'_T_z'//plot_extention 
open(3,file=filenameTz)
!write(3,'(/,a,/)')    '                     z                               temperature'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenameTt, filenameTr, filenameTz
! read(*,*)


!**********************************************************************************************************************
!                                           Constants
!**********************************************************************************************************************

        h = 10.                 !heat transfer coefficient (convection - cylinder)        W/(m^2.K)
       T0 = 300.                !initial temperature                                      K
       pi = 4*atan(1.)                                                                   !dimensionless
       ka = 2.75                !the average or thermal conductivity of KTP crystal       W/(m.K)
                                !k1=2 k2=3 k3=3.3  
	   Cp = 728.016             !heat capacity at constant pressure                       J/(kg.K)

      roh = 2945.               !mass density                                             kg/m^3
     gama = 60.                 !absorption coefficient                                   1/m
    sigma = 5.669e-8            !Stephan-Bultzman constant                                W/(m^2.K^4) 

 tbetween = 1/freq                                                                       !s
    timet = Np*tbetween                                                                  !s
   deltat = tp / 10                                                                      !s     
	   nt = int(tbetween/deltat)                                                         !dimensionless

       nz = 150                                                                          !dimensionless
   length = 0.02                !length of crystal                                        m 
   deltaz = length/nz                                                                    !m

   radius = 0.002               !radius of crystal                                       !m
   omegaf = 0.001               !spot size                                               !m
   deltar = omegaf/10                                                                    !m
       nr = int(radius/deltar)                                                           !dimensionless 

 epsilong = 0.9                 !surface emissivity                                      !dimensionless

stability = ( (2*ka*deltat)/(roh*Cp) ) * ( (deltar**2+deltaz**2)/(deltar**2*deltaz**2) ) !stability coefficient  
                                                                                         !dimensionless

!**********************************************************************************************************************
!                                        Arrays Allocattion 
!**********************************************************************************************************************

allocate(temperature(1:2,0:nr,0:nz))     

!**********************************************************************************************************************
!                                     Giving Zero to Arrays
!********************************************************************************************************************** 

forall (i=1:2,j=0:nr,k=0:nz)
                                 temperature(i,j,k)=0.        
end forall !i

!**********************************************************************************************************************
!                                       Printing Constants     
!**********************************************************************************************************************

write(*,*)
write(*,*)'------- Heat Equation Constants --------------------------------------------'
write(*,*)
write(*,'(A13,I9       )') '        tp = ',tp
write(*,'(A13,I9       )') '        Nt = ',Nt
write(*,'(A13,I9       )') '        Nr = ',Nr
write(*,'(A13,I9       )') '        Np = ',Np
write(*,'(A13,I9     ,/)') '        Nz = ',Nz

		  
write(*,'(A13,F15.10   )') '         h = ',h
write(*,'(A13,F15.10   )') '         E = ',E

write(*,'(A13,F15.10   )') '        T0 = ',T0
write(*,'(A13,F15.10   )') '        pi = ',pi
write(*,'(A13,F15.10 ,/)') '        Cp = ',Cp

write(*,'(A13,F15.10   )') '       kT0 = ',kT0
write(*,'(A13,F15.10 ,/)') '       roh = ',roh

write(*,'(A13,F15.10   )') '      freq = ',freq
write(*,'(A13,F15.10   )') '      gama = ',gama
write(*,'(A13,F15.10   )') '      Tinf = ',Tinf 
write(*,'(A13,F15.10 ,/)') '      Tamb = ',Tamb

write(*,'(A13,F15.10   )') '     timet = ',timet
write(*,'(A13,F15.10 ,/)') '     sigma = ',sigma

write(*,'(A13,F15.10   )') '    omegaf = ',omegaf
write(*,'(A13,F15.10 ,/)') '    length = ',length

write(*,'(A13,F15.10   )') '    deltat = ',deltat
write(*,'(A13,F15.10   )') '    deltar = ',deltar
write(*,'(A13,F15.10 ,/)') '    deltaz = ',deltaz

write(*,'(A13,F15.10 ,/)') '    radius = ',radius

write(*,'(A13,F15.10   )') '  epsilong = ',epsilong
write(*,'(A13,F15.10 ,/)') '  tbetween = ',tbetween

write(*,'(A13,F15.10 ,/)') ' stability = ',stability   
                                                                 
write(*,*)'----------------------------------------------------------------------------'
write(*,'(A,\)')' Please press any key to continue '
read(*,*)


!**********************************************************************************************************************
!                                   Main Block of the Program     
!**********************************************************************************************************************

do j=0,nr
    do k=0,nz

       temperature(i,j,k) = T0
	   
    end do !k
end do !j	     

!------------------ for save time 
aa1 = 1 - (h*deltaz)/(ka)
aa2 = (epsilong*sigma*deltaz)/(ka)
aa3 = (deltat*ka)/(roh*Cp)
aa4 = ( deltat/(roh*Cp) )  * ( (2*gama*E)/(tp*pi*omegaf**2) )
!------------------

do l=1,Np !Run program for Np pulses 
   
   !--------------------------------------------- Run program for one pulse 
   do i=0,nt
      t=i*deltat
      
	  do j=1,nr-1
         r=j*deltar  
      
	     do k=1,nz-1
	        z=k*deltaz 

            !------------------------------------ Boundary conditions
		    temperature(1,0 ,k)  = temperature(1,1,k)                                            !Thermal insulation condition for crystal axis
            temperature(1,nr,k)  = T0                                                            !Temperature-fixed condition for lateral surface
            temperature(1,j ,0)  = aa1*temperature(1,j,1)    - aa2*( temperature(1,j,1)    )**4  !Convection & Radiation condition for input  surface
            temperature(1,j,nz)  = aa1*temperature(1,j,nz-1) - aa2*( temperature(1,j,nz-1) )**4  !Convection & Radiation condition for output surface
            !---------------------
	        temperature(1,0 ,0 ) = aa1*temperature(1,0,1)    - aa2*( temperature(1,0,1)    )**4  !Convection & Radiation condition for ( 0,0 )
		    temperature(1,0 ,nz) = aa1*temperature(1,0,nz-1) - aa2*( temperature(1,0,nz-1) )**4  !Convection & Radiation condition for ( 0,nz)
            temperature(1,nr,0 ) = T0                                                            !Temperature-fixed condition for (nr,0 )
		    temperature(1,nr,nz) = T0                                                            !Temperature-fixed condition for (nr,nz)
			!------------------------------------ End of Boundary conditions

            !------------------------------------ Heat Equation
			temperature(2,j,k) = temperature(1,j,k)                                                               &
		                      
        			 		     + aa3 * (   temperature(1,j+1,k) * (2*r+deltar)/(r * 2 * deltar**2)              &
								 
								 - temperature(1,j  ,k) * (    4*r   )/(r * 2 * deltar**2)                        & 
 
 							     + temperature(1,j-1,k) * (2*r-deltar)/(r * 2 * deltar**2)                        & 

							     +(temperature(1,j,k-1) - 2*temperature(1,j,k) + temperature(1,j,k+1))/deltaz**2) &                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

	        				     + aa4 * exp( (-2*r**2)/(omegaf**2) ) * exp(-gama*z) * exp(-( (t-2*tp)/tp )**2 ) 
              
          end do !k
      end do !j
   !--------------------------------------------- End of run for each deltat 

   !============================================= Print Results for each deltat
   t=(l-1)*nt*deltat + i*deltat 
   write(1,'(2x,f25.10,5x,f25.10)') t , temperature(1,0,0)
   !=============================================

   !--------------------------------------------- End-temprature of each deltat  ==> Initial temperature for next deltat
   do j=1,nr-1
      do k=1,nz-1
 	  
         temperature(1,j,k) = temperature(2,j,k)

      end do !k
   end do !j	     
   !---------------------------------------------

   end do !i
end do !l


!**********************************************************************************************************************
!                                        Printing Results     
!**********************************************************************************************************************

!------------------------------------------------ 
do j=0,nr
   r=j*deltar 
   write(2,'(2x,f25.10,5x,f25.10)') r , temperature(1,j,0)
end do !j      						   

!------------------------------------------------
do k=0,nz
   z=k*deltaz 
   write(3,'(2x,f25.10,5x,f25.10)') z , temperature(1,0,k)
end do !k      						   


!**********************************************************************************************************************
!                                      Closing Files and Ending the Program 
!**********************************************************************************************************************

close(1)
close(2)
close(3)

end program Temp_G_PW                     

!======================================================================================================================
         
 
