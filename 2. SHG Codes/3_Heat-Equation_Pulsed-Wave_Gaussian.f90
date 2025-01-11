! Note: 
!     We have aligned the code with VS Code's formatting standards.
!     But it may appears cluttered on GitHub due to formatting differences.


!            **********************************************************************************
!            *                                                                                *
!            * File name:                                                                     *
!            *     3. Heat Equation _ Pulsed Wave Gaussian _ Computational.F90                *
!            *                                                                                *
!            * This Fortran code is developed specifically for the article titled:            *
!            *     Complete Anisotropic Time-Dependent Heat Equation in KTP Crystal under     *
!            *     Repetitively Pulsed Gaussian Beams: A Numerical Approach                   *
!            *                                                                                *
!            * Cite Us:                                                                       *
!            *     Rezaee, M.M., Sabaeian, M., Motazedian, A., Jalil-Abadi, F.S. and Khaldi-  *
!            *     Nasab, A., 2015. Complete anisotropic time-dependent heat equation in KTP  *
!            *     crystal under repetitively pulsed Gaussian beams: a numerical approach.    *
!            *     Applied Optics, 54(6), pp.1241-1249.                                       *
!            *                                                                                *
!            **********************************************************************************

program Temp_G_PW

implicit none

!**********************************************************************************************************************
!                                       Variables Definition
!**********************************************************************************************************************

integer       i          ,j          ,k          ,l                                                  &
             ,nt         ,nr         ,nz         ,Np                                                                       

real          t          ,z          ,E          ,h          ,r          ,G           ,P             &                                                                                    
             ,T0         ,pi         ,Cp         ,tp         ,Q0                                     &
	          ,roh        ,kT0        ,aa1        ,aa2        ,aa3        ,aa4         ,aa5           &
             ,freq       ,gama       ,Tinf       ,Tamb                                               & 
	          ,timet      ,sigma                                                                      &
	          ,omegaf     ,length     ,deltar     ,deltaz     ,deltat     ,radius                     &
	          ,epsilong   ,tbetween                                                                   &
	          ,stability                                                                              &                                                                             
			    ,r_integral ,z_integral                                                                 &
             ,temperature[allocatable](:,:,:)    ,kT[allocatable](:,:)

character*35  EE                                                                                     &
             ,Npf            ,tpf                                                                    &
             ,freqf                                                                                  &           
             ,filenameTt     ,filenameTr   ,filenameTz   ,filenameKt   ,filenameKr   ,filenameKz     &
			    ,plot_extension
                                                  

!**********************************************************************************************************************
!                                    Giving Zero to variables
!**********************************************************************************************************************

                    i = 0           ;j = 0           ;k = 0          ;l = 0            
                   nt = 0          ;nr = 0          ;nz = 0         ;Np = 0

                    t = 0.          ;z = 0.          ;E = 0.         ;h = 0.       ;r = 0.       ;G = 0.    ;p = 0.                                                                                                             
                   T0 = 0.         ;pi = 0.         ;Cp = 0.        ;tp = 0.      ;Q0 = 0.                              
		            roh = 0.        ;kT0 = 0.        ;aa1 = 0.       ;aa2 = 0.     ;aa3 = 0.     ;aa4 = 0.  ;aa5 = 0.                                                                    
                 freq = 0.       ;gama = 0.       ;Tinf = 0.      ;Tamb = 0.
			   	 timet = 0.      ;sigma = 0.                                                                       
			      omegaf = 0.     ;length = 0.     ;deltar = 0.    ;deltaz = 0.  ;deltat = 0.  ;radius = 0.                                        			      
			    epsilong = 0.   ;tbetween = 0.                                   
			   stability = 0.
           r_integral = 0. ;z_integral = 0.   

!**********************************************************************************************************************
!                                             Inputs		  
!**********************************************************************************************************************

! Note: 
!     This code lets the user enter values twice: once numerically (for calculations) 
!     and once as a string (for filenames or labels).  
!     For example, `E` is number,while `EE` store the same values as strings.  
!     This dual input ensures accurate calculations and meaningful file naming.

!write(*,'(/,2x,a,\)') '                      Enter the Energy value : '
!read(*,*) E          
!write(*,'(/,2x,a,\)') 'Enter the Energy value without decimal point : '
!read(*,*) EE           

!write(*,'(/,2x,a,\)') '                      Enter the frequency value : '
!read(*,*) freq
!write(*,'(/,2x,a,\)') 'Enter the frequency value without decimal point : '
!read(*,*) freqf

!write(*,'(/,2x,a,\)') '                   Enter the Number of Pulses : '
!read(*,*) Np
!write(*,'(/,2x,a,\)') 'Enter the Pulses' value without decimal point : '
!read(*,*) Npf

!write(*,'(/,2x,a,\)') '                            Enter the tp : '
!read(*,*) tp
!write(*,'(/,2x,a,\)') 'Enter the tp value without decimal point : '
!read(*,*) tpf

! For Calculation
    E = 0.09
   Np = 10
   tp = 50e-6
 freq = 500

! For Generating Filenames based on the values above
   EE = '009' 
  Npf = '10'
  tpf = '50'
freqf = '500'


!**********************************************************************************************************************
!                          Determination of Filenames and Opening files
!**********************************************************************************************************************

! Note:
!      To achieve both efficiency and clarity in managing output data,
!      below, we generate filenames based on input information.

plot_extension = '.plt'

filenameTt = 'E_'//trim(EE)//'_f_'//trim(freqf)//'_Np_'//trim(Npf)//'_tp_'//trim(tpf)//'_T_t'//plot_extension 
open(1,file=filenameTt)
!write(1,'(/,a,/)')    '                     t                               temperature'

filenameTr = 'E_'//trim(EE)//'_f_'//trim(freqf)//'_Np_'//trim(Npf)//'_tp_'//trim(tpf)//'_T_r'//plot_extension 
open(2,file=filenameTr)
!write(2,'(/,a,/)')    '                     r                               temperature'

filenameTz = 'E_'//trim(EE)//'_f_'//trim(freqf)//'_Np_'//trim(Npf)//'_tp_'//trim(tpf)//'_T_z'//plot_extension 
open(3,file=filenameTz)
!write(3,'(/,a,/)')    '                     z                               temperature'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :', filenameTt, filenameTr, filenameTz
! read(*,*)

!-----------------------
filenameKt = 'E_'//trim(EE)//'_f_'//trim(freqf)//'_Np_'//trim(Npf)//'_tp_'//trim(tpf)//'_K_t'//plot_extension 
open(4,file=filenameKt)
!write(1,'(/,a,/)')    '                     t                               K(T)'

filenameKr = 'E_'//trim(EE)//'_f_'//trim(freqf)//'_Np_'//trim(Npf)//'_tp_'//trim(tpf)//'_K_r'//plot_extension 
open(5,file=filenameKr)
!write(2,'(/,a,/)')    '                     r                               K(T)'

filenameKz = 'E_'//trim(EE)//'_f_'//trim(freqf)//'_Np_'//trim(Npf)//'_tp_'//trim(tpf)//'_K_z'//plot_extension 
open(6,file=filenameKz)
!write(3,'(/,a,/)')    '                     z                               K(T)'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :', filenameKt, filenameKr, filenameKz
! read(*,*)
!**********************************************************************************************************************
!                                           Constants
!**********************************************************************************************************************

        h = 10.                 !heat transfer coefficient (convection - cylinder)        W/(m^2.K)
       T0 = 300.                !initial temperature                                      K
       pi = 4*atan(1.)                                                                   !dimensionless
       Cp = 728.016             !heat capacity at constant pressure                       J/(kg.K)

      kT0 = 13.                 !the average or thermal conductivity of KTP crystal       W/(m.K)
                                !k1=2 k2=3 k3=3.3  
      roh = 2945.               !mass density                                             kg/m^3
     gama = 0.1 ! ??4.                  !absorption coefficient                                   1/m
     Tinf = 300.                                                                         !K
     Tamb = 300.                                                                         !K   
    sigma = 5.669e-8            !Stephan-Bultzman constant                                W/(m^2.K^4) 

 tbetween = 1/freq                                                                       !s
    timet = Np*tbetween                                                                  !s
   deltat = tp / 10                                                                      !s     
       nt = int(tbetween/deltat)                                                         !dimensionless

       nz = 58                                                                           !dimensionless
   length = 0.02                !length of crystal                                        m 
   deltaz = length/nz                                                                    !m

   radius = 5e-3                !radius of crystal                                       !m
   omegaf = 100.e-6             !spot size                                               !m
   deltar = omegaf/10                                                                    !m
       nr = int(radius/deltar)                                                           !dimensionless 

 epsilong = 0.9                 !surface emissivity                                      !dimensionless

stability = ( (2*kT0*deltat)/(roh*Cp) ) * ( (deltar**2+deltaz**2)/(deltar**2*deltaz**2) )!stability coefficient  
                                                                                         !dimensionless

!**********************************************************************************************************************
!                                        Arrays Allocattion 
!**********************************************************************************************************************
allocate(temperature(1:2,0:nr,0:nz))     

allocate( kT(0:nr,0:nz) )     

!**********************************************************************************************************************
!                                     Giving Zero to Arrays
!********************************************************************************************************************** 

forall (i=1:2,j=0:nr,k=0:nz)
                      temperature(i,j,k)=0.        
end forall !i

!------------------
forall (j=0:nr,k=0:nz)
                                 kT(j,k)=0.       
end forall 

!**********************************************************************************************************************
!                                       Print Constants     
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
!                                   The Main Block of the Program     
!**********************************************************************************************************************

! Display estimated execution time information
write(*,*)
write(*,*) '--- This code takes approximately 1 minute to execute on &
	        a medium-performance      laptop. Execution time may vary depending on &
			the system''s CPU, RAM, and        background tasks. ---!'	

write(*,*) 

!-------- 
! To Calculate normalization constant (G) - Formula (3) in the Article
do k = 0, nz
   z = k * deltaz 
 
   r_integral = 0.
   do j = 0, nr
      r = j * deltar
      r_integral = r_integral + exp(-2 * r**2 / omegaf**2) * r * deltar
   end do !j

   z_integral= z_integral + exp(-gama * z) * r_integral * deltaz
end do !k

Q0 = 1 / (sqrt(pi) * z_integral)    !Normalization	m^-3	-	Formula (14) in the Article

 p = E / (sqrt(pi) * tp)            !total power of the pulse

!--------------------------------------------------------
do j=0,nr
   do k=0,nz
	  
      temperature(1,j,k) = T0
	  
	   kT(j,k) = kT0 
!if (j==1 .and. k==1) write(*,*)kT(j,k)  
   end do !k
end do !j	     


do l=1,Np !Running the program for Np pulses 
   
   !--------------------------------------------- Running the program for one pulse 
   do i=0,nt
      t=i*deltat
      
		do j=1,nr-1
		   r=j*deltar  
      
			do k=1,nz-1
			   z=k*deltaz 

            !------------------  
            aa1 = (h*deltaz)/(kT(j,k))

            aa2 = (epsilong*sigma*deltaz)/(kT(j,k))

            aa3 = ( deltat/(roh*Cp) ) * kT(j,k)          

            aa4 = ( deltat/(roh*Cp) ) * p * Q0

            aa5 = ( deltat/(roh*Cp) ) * (1/4)

            !------------------
            !------------------------------------ Boundary conditions
            temperature(1,0 ,k)  = temperature(1,1,k)                    !Thermal insulation condition for crystal axis

            temperature(1,nr,k)  = T0                                    !Temperature-fixed condition for lateral surface

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

            temperature(1,nr,0 ) = T0                       !Temperature-fixed condition for (nr,0 )

            temperature(1,nr,nz) = T0                       !Temperature-fixed condition for (nr,nz)

            !------------------------------------ Ending of Boundary conditions

            !------------------------------------ Heat Equation
            temperature(2,j,k) = temperature(1,j,k)                                                                             &
                        
                               + aa3 * ( (temperature(1,j+1,k) -  temperature(1,j-1,k))/(2*r*deltar)                            &
                           
                               +(temperature(1,j+1,k) -2*temperature(1,j,k) + temperature(1,j-1,k))/(deltar**2) )               & 

                               + aa3 * ( (temperature(1,j,k-1) -2*temperature(1,j,k) + temperature(1,j,k+1))/(deltaz**2) )      &                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

                               + aa4 * exp( (-2*r**2)/(omegaf**2) ) * exp(-gama*z) * exp(-( (t-2*tp)/tp )**2 )                  &

                               + aa5 * ( ( kT(j+1,k)-kT(j-1,k) ) * ( temperature(1,j+1,k)-temperature(1,j-1,k) ) / (deltar**2)  &

                               +( kT(j,k+1)-kT(j,k-1) ) * ( temperature(1,j,k+1)-temperature(1,j,k-1) ) / (deltaz**2) )
                  
			end do !k
		end do !j
      !--------------------------------------------- End of running for each deltat 
      !============================================= Printing Results for each deltat
      t=(l-1)*nt*deltat + i*deltat 
      write(1,'(2x,f25.10,5x,f25.10)') t , temperature(1,0,0)
      
      write(4,'(2x,f25.10,5x,f25.10)') t , KT(0,0)
      !=============================================

      !--------------------------------------------- End-temprature of each deltat  ==> Initial temperature for next deltat
      do j=1,nr-1
         do k=1,nz-1
      
            temperature(1,j,k) = temperature(2,j,k)

         end do !k
      end do !j	    
      
      !---------------
      do j=0,nr
         do k=0,nz
      
            kT(j,k) = kT0 * T0 / temperature(1,j,k)

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

!=============================================== 
do j=0,nr
   r=j*deltar 
   write(5,'(2x,f25.10,5x,f25.10)') r , KT(j,0)
end do !j      						   

!------------------------------------------------
do k=0,nz
   z=k*deltaz 
   write(6,'(2x,f25.10,5x,f25.10)') z , KT(0,k)
end do !k      						   

!**********************************************************************************************************************
!                                      Closing Files and Ending the Program 
!**********************************************************************************************************************

close(1)
close(2)
close(3)

close(4)
close(5)
close(6)

write(*,*) 
write(*,*) '---- The results are stored in `.plt` format.                                  &
	         If a different format is required, users can set the desried extension in      &
			   "Determine Filenames & Open files" section of the code or rename the file      & 
			   manually and open it with their preferred software. ----!'	
			
write(*,*) 	
write(*,*) '---- Program Completed ----!'

end program Temp_G_PW                     

!======================================================================================================================


