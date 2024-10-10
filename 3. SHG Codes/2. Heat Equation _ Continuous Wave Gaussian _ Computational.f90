

!            ********************************************************************************
!            * File: "2. Heat Equation _ Continuous Wave Gaussian _ Computational.F90"      *
!            *                                                                              *
!            * Note: This Fortran code is developed specifically for the article titled:    *
!            * "Temperature Distribution in a Gaussian End-Pumped Nonlinear KTP Crystal:    *
!            * the Temperature Dependence of Thermal Conductivity and Radiation Boundary    *
!            * Condition"                                                                   *
!            *                                                                              *
!            * Authors: Sabaeian, M., Jalil-Abadi, F.S., Rezaee, M.M., Motazedian, A.       *
!            * and Shahzadeh, M.                                                            *
!            *                                                                              *
!            * Harvard style:                                                               *
!            * Sabaeian, M., Jalil-Abadi, F.S., Rezaee, M.M., Motazedian, A. and Shahzadeh, *
!            * M., 2015. Temperature distribution in a Gaussian end-pumped nonlinear KTP    *
!            * crystal: the temperature dependence of thermal conductivity and radiation    *
!            * boundary condition. Brazilian Journal of Physics, 45, pp.1-9.                *
!            *                                                                              *
!            ********************************************************************************

program Temp_G_CW

implicit none

!**********************************************************************************************************************
!                                       Determine variables
!**********************************************************************************************************************

integer       i          ,j          ,k          ,f                                                                    &
             ,nt         ,nr         ,nz                                                                               

real*8        t          ,z          ,p          ,h           ,r                                                       &                                                                                    
             ,G          ,T0         ,pi         ,Cp          ,roh                                                     &
             ,aa1        ,aa2        ,aa3        ,aa4         ,aa5                                                     &
			 ,KT0        ,Tinf       ,Tamb                                                                             & 
			 ,gama       ,timet      ,sigma                                                                            &
			 ,omegaf     ,length     ,deltar     ,deltaz      ,deltat                                                  &
			 ,radius                                                                                                   &
			 ,epsilong   ,stability  ,Fidegree   ,Firadian                                                             &
			  
    		 ,temperature[allocatable](:,:,:)    , KT[allocatable] (:,:)
			 
complex*16    Ii  

character*30  filenamet  ,filenamer  ,filenamez  ,stabilityf ,timetf ,pp   ,omegafch

!**********************************************************************************************************************
!                                         Zero to variables
!**********************************************************************************************************************
 
                    i = 0.          ;j = 0.            ;k = 0.             ;f = 0.
                   nt = 0.         ;nr = 0.           ;nz = 0.

                    t = 0.          ;z = 0.            ;p = 0.             ;h = 0.           ;r = 0.                                                                                                  
                    G = 0.         ;T0 = 0.           ;pi = 0.            ;Cp = 0.         ;roh = 0.                       
			      aa1 = 0.        ;aa2 = 0.          ;aa3 = 0.           ;aa4 = 0.         ;aa5 = 0. 
				  KT0 = 0.       ;Tinf = 0.         ;Tamb = 0.                                                     
			     gama = 0.      ;timet = 0.        ;sigma = 0.                                                                           
			   omegaf = 0.     ;length = 0.       ;deltat = 0.        ;deltar = 0.      ;deltaz = 0.                               
			   radius = 0.   ;epsilong = 0.     ;Fidegree = 0.      ;Firadian = 0.                                                                                    
	        stability = 0.    

!**********************************************************************************************************************
!                                             Inputs		  
!**********************************************************************************************************************

write(*,'(/,2x,a,\)') 'Enter the stability value  : '
 !read(*,*)stability
write(*,'(/,2x,a,\)') '                    Again  : '
 !read(*,*)stabilityf

write(*,'(/,2x,a,\)') 'Enter the total time value : '
 !read(*,*)timet
write(*,'(/,2x,a,\)') '                     Again : '
 !read(*,*)timetf 

stability = 0.85
timet = 1

stabilityf = '85'
timetf = '1'

!**********************************************************************************************************************
!                                 Determine  Filenames & Open files
!**********************************************************************************************************************

filenamet = 'ST '//trim(stabilityf)//' time '//trim(timetf)//' T t .plt'
open(1,file=filenamet)

filenamer = 'ST '//trim(stabilityf)//' Tim '//trim(timetf)//' T r .plt'
open(2,file=filenamer)

filenamez = 'ST '//trim(stabilityf)//' Tim '//trim(timetf)//' T z .plt'
open(3,file=filenamez)

write(*,*) filenamet  ,filenamer  ,filenamez

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
      
	    p = 80.                 ! power of laser                                         W
	    h = 10.                 !heat transfer coefficient (convection - cylinder)       W/(m^2.K) 
  	!	G = 2.910714020e-9      ! normalization constant(from maple program)
        G = 4.805e-11
       pi = 4*atan(1.)                                                                  !dimensionless
      KT0 = 13.                 !thermal conductivity of KTP crystal                     W/(m.K)
	  	                            ! k1=2 , k2=3 , k3=3.3 
       Cp = 728.016             !heat capacity at constant pressure                      J/(kg.K)
       T0 = 300.                !initial temperature                                     K
      roh = 2945.               !mass density                                            kg/m^3
     Tamb = 300.                !K
	 Tinf = 300.                !K              
     gama = 4.                  !absorption coefficient                                  1/m
    sigma = 5.669e-8            !Stephan-Bultzman constant                               W/(m^2.K^4) 
   radius = 0.002               !radius of crystal                                      !m
   
   omegaf = 0.0001              !spot size                                              !m 
   deltar = omegaf/10                                                                   !m
       nr = int(radius/deltar)                                                          !dimensionless
  
       nz = 150.                                                                        !dimensionless
   length = 0.02                !length of crystal                                       m 
   deltaz = length/nz                                                                   !m 
   deltat = (stability*roh*Cp*0.5/KT0)*(deltar**2*deltaz**2/(deltar**2+deltaz**2))      !s     
   ! deltat = .00001 
	   nt = int(timet/deltat)                                                           !dimensionless 
 epsilong = 0.9                 !surface emissivity                                     !dimensionless

     Ii = (0,1)
      
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
!                                        printe Constants     
!**********************************************************************************************************************

write(*,*)
write(*,*)'------- Constants ----------------------------------------------------------'
write(*,*)
write(*,'(A13,I5    ,/,      &
          A13,I5    ,/,      &
		  A13,I5    ,//,     &
		  
		  A13,F15.10,/,      &
		  A13,F15.10,//,     &

		  A13,F15.10,/,      &
		  A13,F15.10,/,      &
		  A13,F15.10,/,      &
		  A13,F15.10,/,      &
		  A13,2F15.10,//,    &

		  A13,F15.10,//,     &

		  A13,F15.10,/,      &
		  A13,F15.10,//,     &

		  A13,F15.10,/,      &
		  A13,F15.10,/,      &
		  A13,F15.10,/,      &
		  A13,F15.10,/,      &
		  A13,F15.10,//,     &

		  A13,F15.10,//)')   &
          
'        Nt = ',Nt     ,     &
'        Nr = ',Nr     ,     &
'        Nz = ',Nz     ,     &

'         h = ',h      ,     &
'         p = ',p      ,     & 

'        T0 = ',T0     ,     &
'        KT0 = ',KT0   ,     &
'        pi = ',pi     ,     &
'        Cp = ',Cp     ,     &
'        Ii = ',Ii     ,     &       

'       roh = ',roh    ,     &

'     timet = ',timet  ,     &
'      gama = ',gama  ,      &

'    omegaf = ',omegaf ,     &
'    length = ',length ,     &
'    deltat = ',deltat ,     &
'    deltar = ',deltar ,     &
'    deltaz = ',deltaz ,     &       

'   radius = ',radius                                                                     

write(*,*)'----------------------------------------------------------------------------'
   write(*,'(A,\)')' Press Enter to continue '
   read(*,*)

!------------------------------------------------ Control Condition

if ((omegaf .NE. 0.001) .OR. (radius .NE. 0.002) .OR. (length .NE. 0.02)) Then
    write(*,*) 'Error-------please enter new value in maple file'
end if
  
!**********************************************************************************************************************
!                                              Main     
!**********************************************************************************************************************

   do j=0,nr
      do k=0,nz
	  
      temperature(1,j,k) = T0
	             KT(j,k) = KT0
	   
      end do !k
  end do !j	     
!-------------------------- 


do i = 0,nt
   t = i*deltat
   
 !  if (mod(i,100)==0) then
    
 !     f=f+1
	  !------------------------------------ contour Temp step 2   
!      write(11,'(a,i5,a)')  'ZONE T= "',f,'"'   
      !------------------------------------
!      write(21,'(a,i5,a)')  'ZONE T= "',f,'"' 
      !------------------------------------ contor Temp step 3
!	  do j = 0,nr
!         r = j * deltar
	   
!         do k= 0,nz
!	        z = k * deltaz
	     
!		    write(11,'(2x,f25.4,2x,f25.4,2x,f25.4)')  z,r,temperature(1,j,k) 

!			do Fidegree = 6,360,6
!			   Firadian = Fidegree * pi/180

!			   write(21,'(2x,f25.4,2x,f25.4,2x,f25.4)') z,r,Firadian,temperature(1,j,k)
         
!		   end do !Fi 

!         end do !k
!      end do !j
   
!   end if	    	  
   !------------------------------------
 
 !------------------
     do j = 1,nr-1   
        r = j * deltar
	   
	    do k = 1,nz-1
	       z = k * deltaz 
            aa1 = (h*deltaz)/(kT(j,k))
            
			aa2 = (epsilong*sigma*deltaz)/(kT(j,k))
            
			aa3 = ( deltat/(roh*Cp) ) * kT(j,k)           !@@@@@@@@@@ in khat va sharayete marzi avaz shodana  2_06_91
 
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
		    
							  
			temperature(2,j,k) = temperature(1,j,k)                                                                                                   &
		                       
        			 		     + aa3 * (   (temperature(1,j-1,k) - 2 * temperature(1,j,k) + temperature(1,j+1,k))/(deltar ** 2)                     &
								 
								           + (temperature(1,j+1,k) - temperature(1,j-1,k)) /(r * 2 * deltar)                                          & 					            

							               + (temperature(1,j,k-1) - 2 * temperature(1,j,k) + temperature(1,j,k+1))/deltaz**2   )                     &                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

	        				               + aa4 * exp( (-2 * r ** 2)/(omegaf ** 2) ) * exp(- gama * z)                                               &

                                           + aa5 * ( (temperature(1,j+1,k) - temperature(1,j-1,k)) * ( KT(j+1,k) - KT(j-1,k)) / deltar**2  )          &
										   
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
end program Temp_G_CW
                    
!======================================================================================================================
         
