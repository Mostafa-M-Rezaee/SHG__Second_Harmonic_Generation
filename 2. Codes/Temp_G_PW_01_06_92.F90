

!            ********************************************************************************
     
!            * Dr.Mohammad Sabaeian    , Department of Physics, Shahid Chamran University   * 

!            * Mostafa Mohammad-Rezaee , Department of Physics, Shahid Chamran University   *

!            * Alireza Motazedian      , Department of Physics, Shahid Chamran University   *

!            * Fatemeh Sedaghat        , Department of Physics, Shahid Chamran University   *

!            *                                                                              * 

!            * m_sabaeian@yahoo.com                                                         *

!            * mostafa_mohammadrezaee@yahoo.com                                             *

!            * alireza.motazedian@yahoo.com                                                 *

!            * f.sedaghat2010@yahoo.com                                                     *

!            *                                                                              * 

!            * The "Temp_G_PW.F90" is a program solving time-dependent Heat-equation with a *

!            *                     repatatively pulsed-source by  FDM.                      * 

!            *                                                                              * 

!            *             originally Written : 30-Des-2011   by MM  &  AM  &  FS           *

!            *                   last revised : 23-Aug-2013   by MM  &  AM  &  FS           *

!            ********************************************************************************

program Temp_G_PW

implicit none

!**********************************************************************************************************************
!                                       Variables Definition
!**********************************************************************************************************************

integer       i          ,j          ,k          ,l                                                                   &
             ,nt         ,nr         ,nz         ,Np                                                                       

real*8        t          ,z          ,E          ,h          ,r          ,G           ,P                              &                                                                                    
             ,T0         ,pi         ,Cp         ,tp                     ,Q0                                          &
	     ,roh        ,kT0        ,aa1        ,aa2        ,aa3        ,aa4         ,aa5                            &
             ,freq       ,gama       ,Tinf       ,Tamb                                                                & 
	     ,timet      ,sigma                                                                                       &
	     ,omegaf     ,length     ,deltar     ,deltaz     ,deltat                                                  &
	     ,radius                                                                                                  &
	     ,epsilong   ,tbetween                                                                                    &
	     ,stability                                                                                               &
			  
	     ,temperature[allocatable](:,:,:)    ,kT[allocatable](:,:)

character*30  filenameTt   ,filenameTr   ,filenameTz     ,Npf       ,freqf     ,tpf       ,EE                         &
             ,filenameKt   ,filenameKr   ,filenameKz 
                                                     

!**********************************************************************************************************************
!                                    Giving Zero to variables
!**********************************************************************************************************************

                    i = 0           ;j = 0           ;k = 0          ;l = 0            
                   nt = 0          ;nr = 0          ;nz = 0         ;Np = 0

                    t = 0.          ;z = 0.          ;E = 0.         ;h = 0.       ;r = 0.    ;G = 0.    ;p = 0.                                                                                                             
                   T0 = 0.         ;pi = 0.         ;Cp = 0.        ;tp = 0.                 ;Q0 = 0.                              
		  roh = 0.        ;kT0 = 0.        ;aa1 = 0.       ;aa2 = 0.     ;aa3 = 0.  ;aa4 = 0.  ;aa5 = 0.                                                                    
                 freq = 0.       ;gama = 0.       ;Tinf = 0.      ;Tamb = 0.
		timet = 0.      ;sigma = 0.                                                                       
	       omegaf = 0.     ;length = 0.     ;deltar = 0.    ;deltaz = 0.  ;deltat = 0.                                           
	       radius = 0.
	     epsilong = 0.   ;tbetween = 0.                                   
	    stability = 0.   

!**********************************************************************************************************************
!                                             Inputs		  
!**********************************************************************************************************************

write(*,'(/,2x,a,\)') '            Enter the Energy value  : '
!read(*,*) E
E = 0.09           
write(*,'(/,2x,a,\)') '                             Again  : '
!read(*,*) EE
EE = '009'            

write(*,'(/,2x,a,\)') '         Enter the frequency value  : '
!read(*,*) freq
freq = 500
write(*,'(/,2x,a,\)') '                             Again  : '
!read(*,*) freqf
freqf = '500'

write(*,'(/,2x,a,\)') '        Enter the Number of Pulses  : '
!read(*,*) Np
Np=1
write(*,'(/,2x,a,\)') '                             Again  : '
!read(*,*) Npf
Npf='1'

write(*,'(/,2x,a,\)') '                      Enter the tp  : '
!read(*,*) tp
tp = 50e-6
write(*,'(/,2x,a,\)') '                             Again  : '
!read(*,*) tpf
tpf = '50'

!**********************************************************************************************************************
!                          Determination of Filenames and Opening files
!**********************************************************************************************************************

filenameTt = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Tt.plt'
open(1,file=filenameTt)
!write(1,'(/,a,/)')    '                     t                               temperature'

filenameTr = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Tr.plt'
open(2,file=filenameTr)
!write(2,'(/,a,/)')    '                     r                               temperature'

filenameTz = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Tz.plt'
open(3,file=filenameTz)
!write(3,'(/,a,/)')    '                     z                               temperature'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenameTt,filenameTr,filenameTz
! read(*,*)

!-----------------------
filenameKt = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Kt.plt'
open(4,file=filenameKt)
!write(1,'(/,a,/)')    '                     t                               K(T)'

filenameKr = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Kr.plt'
open(5,file=filenameKr)
!write(2,'(/,a,/)')    '                     r                               K(T)'

filenameKz = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Kz.plt'
open(6,file=filenameKz)
!write(3,'(/,a,/)')    '                     z                               K(T)'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenameKt,filenameKr,filenameKz
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
     gama = 4.                  !absorption coefficient                                   1/m
     Tinf = 300.                                                                         !K
     Tamb = 300.                                                                         !K   
    sigma = 5.669e-8            !Stephan-Bultzman constant                                W/(m^2.K^4) 

 tbetween = 1/freq                                                                       !s
    timet = Np*tbetween                                                                  !s
   deltat = tp / 10                                                                      !s     
       nt = int(tbetween/deltat)                                                         !dimensionless

       nz = 150                                                                          !dimensionless
   length = 0.02                !length of crystal                                        m 
   deltaz = length/nz                                                                    !m

   radius = 0.005               !radius of crystal                                       !m
   omegaf = 100.e-6              !spot size                                               !m
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
!                                       Printing Constants     
!**********************************************************************************************************************

write(*,*)
write(*,*)'------- Heat Equation Constants --------------------------------------------'
write(*,*)
write(*,'(A13,I5    ,/,      &
          A13,I5    ,/,      &
		  A13,I5    ,//,     &
		  
		  A13,F15.10,/,      &
		  A13,F15.10,/,      &

		  A13,F15.10,/,      &
		  A13,F15.10,/,      &
		  A13,F15.10,//,     &

		  A13,F15.10,/,      &
		  A13,F15.10,//,     &

		  A13,F15.10,/,      &
		  A13,F15.10,/,      &
		  A13,F15.10,//,     &

		  A13,F15.10,/,      &
		  A13,F15.10,/,      &

		  A13,F15.10,/,      &
		  A13,F15.10,/,      &
		  A13,F15.10,/,      &
		  A13,F15.10,/,      &
		  A13,F15.10,//,     &

		  A13,F15.10,//,     &

		  A13,F15.10,/,      &
		  A13,F15.10,//,     &

		  A13,F15.10,//)')   &
          
'        Nt = ',Nt         ,     &
'        Nr = ',Nr         ,     &
'        Nz = ',Nz         ,     &

'         h = ',h          ,     &
'         E = ',E          ,     & 

'        T0 = ',T0         ,     &
'        pi = ',pi         ,     &
'        Cp = ',Cp         ,     &

'       kT0 = ',kT0        ,     &
'       roh = ',roh        ,     &

'      gama = ',gama       ,     &
'      Tinf = ',Tinf       ,     &
'      Tamb = ',Tamb       ,     &

'     timet = ',timet      ,     &
'     sigma = ',sigma      ,     &
  
'    omegaf = ',omegaf     ,     &
'    length = ',length     ,     &
'    deltat = ',deltat     ,     &
'    deltar = ',deltar     ,     &
'    deltaz = ',deltaz     ,     &       

'    radius = ',radius     ,     &

'  epsilong = ',epsilong   ,     &
'  tbetween = ',tbetween   ,     &

' stability = ',stability                                                                

write(*,*)'----------------------------------------------------------------------------'
write(*,'(A,\)')' Press Enter to continue '
!read(*,*)


!**********************************************************************************************************************
!                                   The Main Block of the Program     
!**********************************************************************************************************************

!-------------------------------------- control condition
 if ( (omegaf .NE. 100.e-6) .OR. (length .NE. 2e-2) .OR. (radius .NE. 5e-3) .OR. (gama .NE. 4.) ) then
      
	write(*,'(a/)')'Error --- please Enter the new variables in maple program and then obtain normalization(G)'
  
  else
 

!-------- 
 G = 1.703409099e-10 * pi * tp                         

Q0 = sqrt(pi) * tp / G              !Normalization      m^-3

 p = E / (tp * sqrt(pi))            !total power of the pulse
!--------------------------------------------------------

do j=0,nr
   do k=0,nz
	  
      temperature(1,j,k) = T0
	  
	         kT(j,k) = kT0 
   
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
	    temperature(2,j,k) = temperature(1,j,k)                                                                                  &
		                      
        		         + aa3 * ( (temperature(1,j+1,k) -  temperature(1,j-1,k))/(2*r*deltar)                               &
								 
				          +(temperature(1,j+1,k) -2*temperature(1,j,k) + temperature(1,j-1,k))/(deltar**2) )         & 
 
                                 + aa3 * ( (temperature(1,j,k-1) -2*temperature(1,j,k) + temperature(1,j,k+1))/(deltaz**2) )         &                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

	        		 + aa4 * exp( (-2*r**2)/(omegaf**2) ) * exp(-gama*z) * exp(-( (t-2*tp)/tp )**2 )                     &

				 + aa5 * ( ( kT(j+1,k)-kT(j-1,k) ) * ( temperature(1,j+1,k)-temperature(1,j-1,k) ) / (deltar**2)     &

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
end if

close(1)
close(2)
close(3)

close(4)
close(5)
close(6)

end program Temp_G_PW                     

!======================================================================================================================


