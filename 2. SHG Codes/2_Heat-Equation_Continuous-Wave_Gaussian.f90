! Note: 
!     We have aligned the code with VS Code's formatting standards.
!     But it may appears cluttered on GitHub due to formatting differences.  



!            **************************************************************************************
!            *                                                                                    *
!            * File name:                                                                         *
!            *     2. Heat Equation _ Continuous Wave Gaussian _ Computational.F90                *
!            *                                                                                    *
!            * This Fortran code is developed specifically for the article titled:                *
!            *     Temperature Distribution in a Gaussian End-Pumped Nonlinear KTP Crystal:       *
!            *     the Temperature Dependence of Thermal Conductivity and Radiation Boundary      *
!            *     Condition                                                                      *
!            *                                                                                    *
!            * Cite Us:                                                                           *
!            *     Sabaeian, M., Jalil-Abadi, F.S., Rezaee, M.M., Motazedian, A. and Shahzadeh,   *
!            *     M., 2015. Temperature distribution in a Gaussian end-pumped nonlinear KTP      *
!            *     crystal: the temperature dependence of thermal conductivity and radiation      *
!            *     boundary condition. Brazilian Journal of Physics, 45, pp.1-9.                  *
!            *                                                                                    *
!            **************************************************************************************

program Temp_G_CW

implicit none

!**********************************************************************************************************************
!                                       Determine variables
!**********************************************************************************************************************

integer       i          ,j          ,k          ,f                                 &
             ,nt         ,nr         ,nz                                                                               

real*8        t          ,z          ,p          ,h           ,r           ,G       &
             ,T0         ,pi 	   	 ,Cp          					                     &
			    ,roh        ,aa1        ,aa2        ,aa3        ,aa4         ,aa5      &
			    ,KT0        ,Tinf       ,Tamb       ,gama                              &
			    ,timet      ,sigma                                                     &
			    ,omegaf     ,length     ,deltar     ,deltaz      ,deltat               &
			    ,radius     ,epsilong   ,Fidegree   ,Firadian                          &
			    ,stability                                                             &
			    ,r_integral ,z_integral                                                &
    		    ,temperature[allocatable](:,:,:)    , KT[allocatable] (:,:)
			 
complex*16    Ii  

character*30  pp                                                                    &
             ,timetf     ,omegafch                                                  &
		   	 ,filenamet  ,filenamer  ,filenamez                                     &
             ,stabilityf                                                            &
			    ,plot_extension                                                          
			 
!**********************************************************************************************************************
!                                         Zero to variables
!**********************************************************************************************************************
 
                 i = 0.          ;j = 0.            ;k = 0.             ;f = 0.
                nt = 0.         ;nr = 0.           ;nz = 0.

                 t = 0.          ;z = 0.            ;p = 0.             ;h = 0.           ;r = 0.           ;G = 0.  
                T0 = 0.         ;pi = 0.           ;Cp = 0.         
			      roh = 0.        ;aa1 = 0.          ;aa2 = 0.           ;aa3 = 0.         ;aa4 = 0.         ;aa5 = 0.          
				   KT0 = 0.
              Tinf = 0.       ;Tamb = 0.         ;gama = 0.
			    timet = 0.      ;sigma = 0.
			   omegaf = 0.     ;length = 0.       ;deltar = 0.        ;deltaz = 0.      ;deltat = 0.
			   radius = 0.   ;epsilong = 0.     ;Fidegree = 0.      ;Firadian = 0.
	      stability = 0.  
        r_integral = 0. ;z_integral = 0.    

!**********************************************************************************************************************
!                                             Inputs		  
!**********************************************************************************************************************

! Note: 
!     This code lets the user enter values twice: once numerically (for calculations) 
!     and once as a string (for filenames or labels).  
!     For example, `stability` is number,while `stabilityf` store the same values as strings.  
!     This dual input ensures accurate calculations and meaningful file naming.

!write(*,'(/,2x,a,\)') '                     Enter the stability value  : '
 !read(*,*)stability
!write(*,'(/,2x,a,\)') 'Enter the stability value without decimal point : '
 !read(*,*)stabilityf

!write(*,'(/,2x,a,\)') '                      Enter the total time value : '
 !read(*,*)timet
!write(*,'(/,2x,a,\)') 'Enter the total time value without decimal point : '
 !read(*,*)timetf 

! For Calculation        
timet = 1
stability = 0.85

! For Generating Filenames based on the values above
timetf = '1'
stabilityf = '085'

!**********************************************************************************************************************
!                                 Determine  Filenames & Open files
!**********************************************************************************************************************

! Note:
!      To achieve both efficiency and clarity in managing output data,
!      below, we generate filenames based on input information.

plot_extension = '.plt'

filenamet = 'ST_'//trim(stabilityf)//'_time_'//trim(timetf)//'_T_t'//plot_extension
open(1,file=filenamet)

filenamer = 'ST_'//trim(stabilityf)//'_time_'//trim(timetf)//'_T_r'//plot_extension
open(2,file=filenamer)

filenamez = 'ST_'//trim(stabilityf)//'_time_'//trim(timetf)//'_T_z'//plot_extension
open(3,file=filenamez)

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :', filenamet, filenamer, filenamez

!------------------------------------------------ contour Temp Step 1
!open(11,file= 'contour_temp_2D.plt')
!write(11,* ) 'TITLE = "EXAMPLE : MULTI-ZONE 2D PLOT"'  
!write(11,* ) 'variables = "r" , "z" , "temp"'

!open(21,file= 'contour_temp_3D.plt')
!write(21,* ) 'TITLE = "EXAMPLE : MULTI-ZONE 3D PLOT"'  
!write(21,* ) 'variables = "r" , "z" , "fi" ,  "temp"'

!**********************************************************************************************************************
!                                        Determine  Constants
!*********************************************************************************************************************
      
	     p = 80.                 !power of laser                                          W
	     h = 10.                 !heat transfer coefficient (convection - cylinder)       W/(m^2.K) 
       pi = 4*atan(1.)                                                                  !dimensionless
      KT0 = 13.                 !thermal conductivity of KTP crystal                     W/(m.K)
	  	                            ! k1=2 , k2=3 , k3=3.3 
       Cp = 728.016             !heat capacity at constant pressure                      J/(kg.K)
       T0 = 300.                !initial temperature                                     K
      roh = 2945.               !mass density                                            kg/m^3
     Tamb = 300.                !K
	  Tinf = 300.                !K              
     gama = 1.                  !absorption coefficient                                  1/m
    sigma = 5.669e-8            !Stephan-Bultzman constant                               W/(m^2.K^4) 

   radius = 0.002               !radius of crystal                                      !m
   omegaf = 0.00005             !spot size                                              !m 
   deltar = omegaf/10                                                                   !m
       nr = int(radius/deltar)                                                          !dimensionless
  
       nz = 150.                                                                        !dimensionless
   length = 0.02                !length of crystal                                       m 
   deltaz = length/nz                                                                   !m 

   deltat = (stability*roh*Cp*0.5/KT0)*(deltar**2*deltaz**2/(deltar**2+deltaz**2))      !s     
   	 nt = int(timet/deltat)                                                           !dimensionless 

 epsilong = 0.9                 !surface emissivity                                     !dimensionless

       Ii = (0,1)

! To Calculate normalization constant (G) - Formula (3) in the Article
do k = 0, nz
   z = k * deltaz 
 
   r_integral = 0.
   do j = 0, nr
      r = j * deltar
      r_integral = r_integral + exp(-2 * r**2 / omegaf**2) * r * deltar
   end do

   z_integral= z_integral+ exp(-gama * z) * r_integral * deltaz
end do
G = z_integral

      
!**********************************************************************************************************************                                                
!                                         Allocate Arrys
!**********************************************************************************************************************

allocate(temperature(1:2,0:nr,0:nz))
allocate(KT(0:nr,0:nz)) 

!**********************************************************************************************************************
!                                         Zero to Arrays
!**********************************************************************************************************************
 forall (i=1:2,j=0:nr,k=0:nz)
                                             
                                temperature(i,j,k) = 0.                                     
 end forall 

!--------------------------

forall (j=0:nr,k=0:nz )
                               KT(j,k) = 0. 
 end forall

!**********************************************************************************************************************
!                                        print Constants     
!**********************************************************************************************************************

write(*,*)'------- Constants ----------------------------------------------------------'
write(*,*)
write(*,'(A13,I9       )') '        Nt = ',Nt
write(*,'(A13,I9       )') '        Nr = ',Nr
write(*,'(A13,I9     ,/)') '        Nz = ',Nz

write(*,'(A13,F15.10   )') '         h = ',h
write(*,'(A13,F15.10   )') '         p = ',p
write(*,'(A13,F15.10 ,/)') '         G = ',G

write(*,'(A13,F15.10   )') '        T0 = ',T0
write(*,'(A13,F15.10   )') '       KT0 = ',KT0
write(*,'(A13,F15.10   )') '        pi = ',pi
write(*,'(A13,F15.10   )') '        Cp = ',Cp
write(*,'(A13,2F15.10,/)') '        Ii = ',Ii   

write(*,'(A13,F15.10 ,/)') '       roh = ',roh

write(*,'(A13,F15.10   )') '     timet = ',timet
write(*,'(A13,F15.10 ,/)') '      gama = ',gama     

write(*,'(A13,F15.10   )') '    omegaf = ',omegaf
write(*,'(A13,F15.10   )') '    length = ',length
write(*,'(A13,F15.10   )') '    deltat = ',deltat
write(*,'(A13,F15.10   )') '    deltar = ',deltar
write(*,'(A13,F15.10 ,/)') '    deltaz = ',deltaz

write(*,'(A13,F15.10 ,/)') '    radius = ',radius   

write(*,'(A13,F15.10 ,/)') ' stability = ',stability   

                                                                        
write(*,*)'----------------------------------------------------------------------------'
   write(*,'(A,\)')' Please press any key to continue '
   read(*,*)

!**********************************************************************************************************************
!                                   The Main Block of the Program     
!**********************************************************************************************************************

! Display estimated execution time information
write(*,*)
write(*,*) '--- This code takes approximately 20 minutes to execute on &
	        a medium-performance    laptop. Execution time may vary depending on &
			the system''s CPU, RAM, and        background tasks. ---!'	

write(*,*) 	


do j=0,nr
  do k=0,nz
  
  temperature(1,j,k) = T0
			 KT(j,k) = KT0
   
  end do !k
end do !j	     
!-------------------------- 

do i = 0,nt
   t = i*deltat
   
   do j = 1,nr-1   
      r = j * deltar
	   
	   do k = 1,nz-1
	      z = k * deltaz 
         aa1 = (h*deltaz)/(kT(j,k))
         
         aa2 = (epsilong*sigma*deltaz)/(kT(j,k))
         
         aa3 = ( deltat/(roh*Cp) ) * kT(j,k)           

         aa4 = ( (deltat)/(roh * Cp) ) * ( ( p)/(2 * pi * G) )

         aa5 = deltat /(4 * roh * Cp  )

         !------------------------------------ Boundary conditions		 

         temperature(1,0 ,k)  = temperature(1,1,k)                !Thermal insulation condition for crystal axis
         temperature(1,nr,k)  = T0                                !Temperature-fixed condition for lateral surface
         temperature(1,j ,0)  = temperature(1,j,1) - aa1*( temperature(1,j,1) - Tinf )             &
                                                   - aa2*( temperature(1,j,1)**4 - Tamb**4 )
                                                                  !Convection & Radiation condition for input  surface
         temperature(1,j,nz)  = temperature(1,j,nz-1) - aa1*( temperature(1,j,nz-1) - Tinf )       &
                                                      - aa2*( temperature(1,j,nz-1)**4 - Tamb**4 )
                                                                  !Convection & Radiation condition for output surface
         !---------------------
         temperature(1,0 ,0 ) = temperature(1,0,1) - aa1*( temperature(1,0,1) - Tinf )             &
                                                   - aa2*( temperature(1,0,1)**4 - Tamb**4 ) 
                                                                  !Convection & Radiation condition for ( 0,0 )
         temperature(1,0 ,nz) = temperature(1,0,nz-1) - aa1*( temperature(1,0,nz-1) - Tinf )       &
                                                      - aa2*( temperature(1,0,nz-1)**4 - Tamb**4 ) 
                                                                  !Convection & Radiation condition for ( 0,nz)
         temperature(1,nr,0 ) = T0                                !Temperature-fixed condition for (nr,0 )
         temperature(1,nr,nz) = T0                                !Temperature-fixed condition for (nr,nz)
			!------------------------------------ Ending of Boundary conditions
  

         !------------------------------------ Heat Equation 
		    
							  
      temperature(2,j,k) = temperature(1,j,k)                                                                                   &
                           
                           + aa3 * (   (temperature(1,j-1,k) - 2 * temperature(1,j,k) + temperature(1,j+1,k))/(deltar ** 2)     &
                        
                           + (temperature(1,j+1,k) - temperature(1,j-1,k)) /(r * 2 * deltar)                                    &

                           + (temperature(1,j,k-1) - 2 * temperature(1,j,k) + temperature(1,j,k+1))/deltaz**2   )               &

                           + aa4 * exp( (-2 * r ** 2)/(omegaf ** 2) ) * exp(- gama * z)                                         &

                           + aa5 * ( (temperature(1,j+1,k) - temperature(1,j-1,k)) * ( KT(j+1,k) - KT(j-1,k)) / deltar**2  )    &
                              
                           + aa5 * ( (temperature(1,j,k+1) - temperature(1,j,k-1)) * ( KT(j,k+1) - KT(j,k-1)) / deltaz**2  )  
                                          
	   end do !k
 
   end do !j

!--------------------------------------------- End of run for each deltat 

!============================================= Print Results for each deltat
 t = deltat * i
  write(1,'(2x,f25.10,5x,f27.12)') t , temperature(1,0,0)
!=================================================
  
!--------------------------------------------- End-temprature of each deltat  ==> Initial temperature for next deltat

   do j=1,nr-1
      do k=1,nz-1
 	  	
         temperature(1,j,k) = temperature(2,j,k)
		    
         ! kT(j,k) = kT0 * T0 / temperature(2,j,k) 
      end do !k
   end do !j
	 

   !----------------------------------------------
   do j=0,nr
      do k=0,nz
 	  	
		 kT(j,k) = kT0   * T0 / temperature(1,j,k)
         
      end do !k
   end do !j

    
   !----------------------------------------------
end do !i 
	
!**********************************************************************************************************************
!                                          printe results     
!**********************************************************************************************************************      						   

do j=0,nr
   r=j*deltar 
   write(2,'(2x,f25.10,5x,f27.12)')   r , temperature(1,j,0)
end do !j     						   
!-------------------

do k=0,nz
   z=k*deltaz 
   write(3,'(2x,f25.10,5x,f27.12)')  z , temperature(1,0,k)
end do !k     
 						   
!**********************************************************************************************************************
!                                      Close files & end program 
!**********************************************************************************************************************
close(1)
close(2)
close(3)
close(11)
close(21)

write(*,*) 
write(*,*) '---- The results are stored in `.plt` format.                                  &
	         If a different format is required, users can set the desried extension in      &
			   "Determine Filenames & Open files" section of the code or rename the file      & 
			   manually and open it with their preferred software. ----!'	

write(*,*) 	
write(*,*) '---- Program Completed ----!'

end program Temp_G_CW
                    
				
!======================================================================================================================