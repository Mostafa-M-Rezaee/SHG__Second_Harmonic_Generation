! Note: 
!     We have aligned the code with VS Code's formatting standards.
!     But it may appears cluttered on GitHub due to formatting differences.


!            **********************************************************************************
!            * File name:                                                                     *
!            *     7. Coupled _ Continuous Wave Gaussian _ Computational.F90                  *
!            *                                                                                *
!            * This Fortran code is developed specifically for the article titled:            *
!            *     Heat coupled Gaussian continuous-wave double-pass type-II second harmonic  *
!            *     generation: inclusion of thermally induced phase mismatching and thermal   *
!            *     lensing                                                                    *
!            *                                                                                *
!            * Cite Us:                                                                       *
!            *     Sabaeian, M., Jalil-Abadi, F.S., Rezaee, M.M. and Motazedian, A., 2014.    *
!            *     Heat coupled Gaussian continuous-wave double-pass type-II second harmonic  *
!            *     generation: inclusion of thermally induced phase mismatching and thermal   *
!            *     lensing. Optics Express, 22(21), pp.25615-25628.                           *
!            *                                                                                *
!            **********************************************************************************

program Elec_temp_phase_CW
implicit none

!*********************************************************************************************************************
!                                       Variables Definition
!*********************************************************************************************************************

!--------------------------------- common variables

integer       i            ,j            ,k            ,f            ,m                                              &
             ,nr           ,nt           ,nz                                                                         &
       	    ,fmax

real*8        t            ,r            ,z            ,p            ,A            ,B            ,S                  &                                              
		       ,pi                                                                                                     &
             ,gama1        ,gama2        ,gama3        ,timet                                                        &
             ,omegaf       ,length       ,deltar       ,deltaz       ,deltat       ,radius                           &		                                                                
		       ,lambda1      ,lambda2 
			 		    
!--------------------------------- fields variables
  
real*8        c            ,L                                                                                        &
             ,r1f          ,r2f         ,r3f           ,cc3          ,dd3          ,ee3                              &	 	                                           
             ,deff                                                                                                   & 
	 	       ,omega                                                                                                  & 
	 	       ,psi12pr      ,psi22pr      ,psi32pr      ,psi12pz      ,psi22pz      ,psi32pz                          &	 	                                    
		       ,psi12mr      ,psi22mr      ,psi32mr      ,psi12mz      ,psi22mz      ,psi32mz                          &		                                         
		       ,epsilon0                                                  
		


complex*16    Ii                                                                                                     &

             ,cc2          ,dd2          ,ee2          ,cc4          ,dd4          ,ee4                              &	                                          
		       ,cc5          ,dd5          ,ee5                                                                        &  			
		       ,psi1p[allocatable](:,:)    ,psi1m[allocatable](:,:)                                                    &
		       ,psi2p[allocatable](:,:)    ,psi2m[allocatable](:,:)                                                    &
		       ,psi3p[allocatable](:,:)    ,psi3m[allocatable](:,:) 			                                                                       
			                                                                                                    
character*30  filenamepsi12pr            ,filenamepsi12pz                                                            & 
             ,filenamepsi22pr            ,filenamepsi22pz                                                            &
			    ,filenamepsi32pr            ,filenamepsi32pz                                                            &
	          ,filenamepsi12mr            ,filenamepsi12mz                                                            & 
             ,filenamepsi22mr            ,filenamepsi22mz                                                            &
			    ,filenamepsi32mr            ,filenamepsi32mz  		    
			
!--------------------------------- temperature variables

real*8        h                                                                                                      &
			    ,T0           ,Cp                                                                                       &
	          ,aa1          ,aa2          ,aa3          ,aa4          ,aa5          ,roh                              &			                                  
			    ,KT0                                                                                                    &      
			    ,Tinf         ,Tamb                                                                                     &        
		       ,sigma                                                                                                  &
			    ,epsilong                                                                                               &
			    ,fidegree     ,firadian                                                                                 & 
		       ,stability                                                                                              &
			
	    	    ,Temperature[allocatable](:,:,:)          ,KT[allocatable](:,:)                         
			
character*30  timetch      ,filenamet    ,filenamer    ,filenamez                                                    &
             ,stabilitych  
			   
!--------------------------------- phase variables

real*8        phi                                                                                                    &
             ,B1T0         ,B2T0         ,C1T0         ,C2T0         ,B1rT         ,B2rT         ,C1rT               &
	          ,C2rT                                                                                                   &		 

			    ,aa1T0        ,bb1T0        ,cc1T0        ,nx1T0        ,ny1T0        ,nz1T0        ,no1T0              &
			    ,aa2T0        ,bb2T0        ,cc2T0        ,nx2T0        ,ny2T0        ,nz2T0        ,ne1T0              &
			    ,aa1rT        ,bb1rT        ,cc1rT        ,nx1rT        ,ny1rT        ,nz1rT        ,ne2T0              &
			    ,aa2rT        ,bb2rT        ,cc2rT        ,nx2rT        ,ny2rT        ,nz2rT        ,no1rT              &
			    ,B1r0T        ,B2r0T        ,C1r0T        ,C2r0T        ,theta        ,ne2rT        ,ne1rT              &                                
             ,Term1        ,Term2        ,Term3                                                                      &            

			    ,aa1r0T       ,bb1r0T       ,cc1r0T       ,nx1r0T       ,ny1r0T       ,nz1r0T       ,dnx1dT             &
			    ,aa2r0T       ,bb2r0T       ,cc2r0T       ,nx2r0T       ,ny2r0T       ,nz2r0T       ,dnx2dT             &
             ,dny1dT       ,dny2dT       ,dnz1dT       ,dnz2dT       ,no1r0T       ,ne1r0T       ,ne2r0T             &   
			    ,deltano1rT   ,deltane1rT   ,deltane2rT                                                                 &                                                           
			    ,deltano1r0T  ,deltane1r0T  ,deltane2r0T                                                     

 complex*8    deltaphase[allocatable](:,:)

character*30  filenamePt   ,filenamePr  ,filenamePz                                                                  &
             ,plot_extension                                                                    
     			              			                                      

!**********************************************************************************************************************
!                                    Giving Zero to variables
!**********************************************************************************************************************

!--------------------------------- Giving Zero to common variables

              i = 0             ;j = 0             ;k = 0         ;f = 0         ;m = 0              
             nr = 0            ;nt = 0            ;nz = 0                                       
           fmax = 0

              p = 0.            ;A = 0.            ;B = 0.        ;S = 0.        ;t = 0.        ;r = 0.        ;z = 0.                                                                               
		       pi = 0.                                                                
          gama1 = 0.        ;gama2 = 0.        ;gama3 = 0.    ;timet = 0.                                   
         omegaf = 0.       ;length = 0.       ;deltar = 0.   ;deltaz = 0.   ;deltat = 0.   ;radius = 0.                                                                	   
	     lambda1 = 0.      ;lambda2 = 0.

!--------------------------------- Giving Zero to field equation variable
		    
              c = 0.            ;L = 0.                                                       
            r1f = 0.          ;r2f = 0.          ;r3f = 0.      ;cc3 = 0.      ;dd3 = 0.      ;ee3 = 0.                                              
           deff = 0.                                                             
	   	 omega = 0.                                                            
        psi12pr = 0.      ;psi22pr = 0       ;psi32pr = 0.  ;psi12pz = 0.  ;psi22pz = 0.  ;psi32pz = 0.                                                                          
        psi12mr = 0.      ;psi22mr = 0.      ;psi32mr = 0.  ;psi12mz = 0.  ;psi22mz = 0.  ;psi32mz = 0.                                                                        
       epsilon0 = 0.  	
	   		                                                
		       Ii = (0.,0.)                                                               
            cc2 = (0.,0.)     ;dd2 = (0.,0.)     ;ee2 = (0.,0.) ;cc4 = (0.,0.) ;dd4 = (0.,0.) ;ee4 = (0.,0.)                                                                                
	                                                  
		                                                                                                    
!-------------------------------- Giving Zero to temperature variable

              h = 0.                                                         
	          T0 = 0.           ;Cp = 0.                                            
	         aa1 = 0.          ;aa2 = 0.          ;aa3 = 0.      ;aa4 = 0.      ;aa5 = 0.      ;roh = 0.                                                                   
            KT0 = 0.
           Tinf = 0.         ;Tamb = 0.                                                 
	       sigma = 0.                                                     
       epsilong = 0.     ;fidegree = 0.     ;firadian = 0.                                                                                    
      stability = 0.                           

                                                  
!------------------------------- Giving Zero to phase variable
            phi = 0.                                                                                                    
           B1T0 = 0.         ;B2T0 = 0.         ;C1T0 = 0.     ;C2T0 = 0.     ;B1rT = 0.     ;B2rT = 0.             
		     C1rT = 0.         ;C2rT = 0.    
			 
		    aa1T0 = 0.        ;bb1T0 = 0.        ;cc1T0 = 0.    ;nx1T0 = 0.    ;ny1T0 = 0.    ;nz1T0 = 0.    ;no1T0 = 0.         
		    aa2T0 = 0.        ;bb2T0 = 0.        ;cc2T0 = 0.    ;nx2T0 = 0.    ;ny2T0 = 0.    ;nz2T0 = 0.    ;ne2T0 = 0.         
		    aa1rT = 0.        ;bb1rT = 0.        ;cc1rT = 0.    ;nx1rT = 0.    ;ny1rT = 0.    ;nz1rT = 0.    ;ne1T0 = 0.               
		    aa2rT = 0.        ;bb2rT = 0.        ;cc2rT = 0.    ;nx2rT = 0.    ;ny2rT = 0.    ;nz2rT = 0.    ;no1rT = 0.        
	       B1r0T = 0.        ;B2r0T = 0.        ;C1r0T = 0.    ;C2r0T = 0.    ;theta = 0.    ;ne2rT = 0.    ;ne1rT = 0.                                                         
          Term1 = 0.        ;Term2 = 0.        ;Term3 = 0.
            
	      aa1r0T = 0.       ;bb1r0T = 0.       ;cc1r0T = 0.   ;nx1r0T = 0.   ;ny1r0T = 0.   ;nz1r0T = 0.   ;dnx1dT = 0.       
		   aa2r0T = 0.       ;bb2r0T = 0.       ;cc2r0T = 0.   ;nx2r0T = 0.   ;ny2r0T = 0.   ;nz2r0T = 0.   ;dnx2dT = 0.     
         dny1dT = 0.       ;dny2dT = 0.       ;dnz1dT = 0.   ;dnz2dT = 0.   ;no1r0T = 0.   ;ne1r0T = 0.   ;ne2r0T = 0.                                             
		     
	  deltano1rT = 0.   ;deltane1rT = 0.   ;deltane2rT = 0.                                                                                                                           
    deltano1r0T = 0.  ;deltane1r0T = 0.  ;deltane2r0T = 0.                                                    

!**********************************************************************************************************************
!                                             Inputs		  
!**********************************************************************************************************************

! Note: 
!     This code lets the user enter values twice: once numerically (for calculations) 
!     and once as a string (for filenames or labels).  
!     For example, `stability` is number,while `stabilitych` store the same values as strings.  
!     This dual input ensures accurate calculations and meaningful file naming.

!write(*,'(/,2x,a,\)') '                      Enter the stability value : '
!read(*,*)stability
!write(*,'(/,2x,a,\)') 'Enter the stability value without decimal point : '
!read(*,*)stabilitych

!write(*,'(/,2x,a,\)') '                      Enter the total time value : '
!read(*,*)timet
!write(*,'(/,2x,a,\)') 'Enter the total time value without decimal point : '
!read(*,*)timetch

! For Calculation
stability = 0.85
timet = 0.0001

! For Generating Filenames based on the values above
stabilitych = '85'
timetch = '01'

!**********************************************************************************************************************
!                          Determination of Filenames and Opening files
!**********************************************************************************************************************

! Note:
!      To achieve both efficiency and clarity in managing output data,
!      below, we generate filenames based on input information.

plot_extension = '.plt'

!------------------------------------  temperature 

filenamet = 'ST_'//trim(stabilitych)//'_time_'//trim(timetch)//'_T_t'//plot_extension
open(1,file=filenamet)

filenamer = 'ST_'//trim(stabilitych)//'_time_'//trim(timetch)//'_T_r'//plot_extension
open(2,file=filenamer)

filenamez = 'ST_'//trim(stabilitych)//'_time_'//trim(timetch)//'_T_z'//plot_extension
open(3,file=filenamez)

!-----------------------------------  phase equation

filenamet = 'ST_'//trim(stabilitych)//'_time_'//trim(timetch)//'_p_t'//plot_extension
open(4,file=filenamet)

filenamer = 'ST_'//trim(stabilitych)//'_time_'//trim(timetch)//'_p_r'//plot_extension
open(5,file=filenamer)

filenamez = 'ST_'//trim(stabilitych)//'_time_'//trim(timetch)//'_p_z'//plot_extension
open(6,file=filenamez)

write(*,*) filenamet  ,filenamer  ,filenamez

!-----------------------------------  field
 
filenamePsi12pr = 'Psi_12_p_r'//plot_extension
open(7,file=filenamePsi12pr)
!write(7,'(/,a,/)')    ' variables =       "r"                               "Psi1p ** 2"'

filenamePsi12pz = 'Psi_12_p_z'//plot_extension
open(8,file=filenamePsi12pz)
!write(8,'(/,a,/)')    ' variables =       "z"                               "Psi1p ** 2"'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenamePsi12pr ,filenamePsi12pz
 write(*,'(A,\)')' Please press any key to continue '
 read(*,*)

!------------------

filenamePsi22pr = ' Psi_22_p_r'//plot_extension
open(9,file=filenamePsi22pr)
!write(9,'(/,a,/)')    ' variables =        "r"                               "Psi2p ** 2"'

filenamePsi22pz =' Psi_22_p_z'//plot_extension
open(10,file=filenamePsi22pz)
!write(10,'(/,a,/)')    ' variables =       "z"                               "Psi2p ** 2"'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenamePsi22pr ,filenamePsi22pz
 write(*,'(A,\)')' Please press any key to continue '
 read(*,*)

!------------------

filenamePsi32pr = ' Psi_32_p_r'//plot_extension
open(11,file=filenamePsi32pr)
!write(11,'(/,a,/)')    ' variables =       "r"                               "Psi3p ** 2"'

filenamePsi32pz = ' Psi_32_p_z'//plot_extension
open(12,file=filenamePsi32pz)
!write(12,'(/,a,/)')    ' variables =       "z"                               "Psi3p** 2"'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenamePsi32pr ,filenamePsi32pz
 write(*,'(A,\)')' Please press any key to continue '
 read(*,*)

!--------------------------------------------------------------------------------------------

filenamePsi12mr = 'Psi_12_m_r'//plot_extension
open(13,file=filenamePsi12mr)
!write(13,'(/,a,/)')    ' variables =       "r"                               "Psi1m ** 2"'

filenamePsi12mz = 'Psi_12_m_z'//plot_extension
open(14,file=filenamePsi12mz)
!write(14,'(/,a,/)')    ' variables =       "z"                               "Psi1m ** 2"'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenamePsi12mr ,filenamePsi12mz
 write(*,'(A,\)')' Please press any key to continue '
 read(*,*)

!------------------

filenamePsi22mr = ' Psi_22_m_r'//plot_extension
open(15,file=filenamePsi22mr)
!write(15,'(/,a,/)')    ' variables =       "r"                               "Psi2m ** 2"'

filenamePsi22mz =' Psi_22_m_z'//plot_extension
open(16,file=filenamePsi22mz)
!write(16,'(/,a,/)')    ' variables =       "z"                               "Psi2m ** 2"'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenamePsi22mr ,filenamePsi22mz
 write(*,'(A,\)')' Please press any key to continue '
 read(*,*)

!------------------

filenamePsi32mr = ' Psi_32_m_r'//plot_extension
open(17,file=filenamePsi32mr)
!write(17,'(/,a,/)')    ' variables =       "r"                               "Psi3m ** 2"'

filenamePsi32mz = ' Psi_32_m_z'//plot_extension
open(18,file=filenamePsi32mz)
!write(18,'(/,a,/)')    ' variables =       "z"                               "Psi3m** 2"'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenamePsi32mr ,filenamePsi32mz
 write(*,'(A,\)')' Please press any key to continue '
 read(*,*)
!--------------------------------------------------------------------------


!**********************************************************************************************************************
!                                           Constants
!**********************************************************************************************************************

!----------------- The constants of Common

        p = 40.                 !power of laser                                          W
       nz = 200                                                                         !dimensionless
       pi = 4.*atan(1.)                                                                 !dimensionless	   	   	   	   	
    gama1 = 0.5                 !the absorption coefficient of fundomental wave          1/m
    gama2 = 0.5                 !the absorption coefficient of fundomental wave          1/m
    gama3 = 4.                  !the absorption coefficient of SHW                       1/m 
   radius = 0.002               !radius of crystal                                      !m
   omegaf = 200e-6              !spot size                                              !m
   length = 0.01                !length of crystal                                       m       
   deltaz = length/nz                                                                   !m
   deltar = omegaf/10.                                                                  !m 
       nr = int(radius/deltar)                                                          !dimensionless

!----------------- The constants of Temperature 

        h = 10.                 !heat transfer coefficient (convection - cylinder)       W/(m^2.K)
       T0 = 300.                !initial temperature                                     K
       Cp = 728.016             !heat capacity at constant pressure                      J/(kg.K)                                                                    !dimensionless 
      KT0 = 13.                 !thermal conductivity of KTP crystal                     W/(m.K)                 
      roh = 2945.               !mass density                                            kg/m^3 	     
     Tamb = 300.
     Tinf = 300.	  
    sigma = 5.669e-8            !Stephan-Bultzman constant                               W/(m^2.K^4) 
 epsilong = 0.9                 !surface emissivity                                     !dimensionless
       Ii = (0.,1.)  
   deltat = (stability*roh*Cp*0.5/KT0)*(deltar**2.*deltaz**2./(deltar**2.+deltaz**2.))  !s       
       nt = int(timet/deltat)                                                           !dimensionless
 	   	   
!------------------ The constants of phase equation
 
      phi = 23.5*pi/180.
    theta = 90.*pi/180.
  lambda1 = 1064e-9  
  lambda2 =  532e-9

   dnx1dT = (0.1323*(lambda1*1e6)**(-3.) - 0.4385*(lambda1*1e6)**(-2.) + 1.2307*(lambda1*1e6)**(-1.) + 0.7709)*1e-5
   dny1dT = (0.5014*(lambda1*1e6)**(-3.) - 2.0030*(lambda1*1e6)**(-2.) + 3.3016*(lambda1*1e6)**(-1.) + 0.7498)*1e-5
   dnz1dT = (0.3896*(lambda1*1e6)**(-3.) - 1.3332*(lambda1*1e6)**(-2.) + 2.2762*(lambda1*1e6)**(-1.) + 2.1151)*1e-5

   dnx2dT = (0.1323*(lambda2*1e6)**(-3.) - 0.4385*(lambda2*1e6)**(-2.) + 1.2307*(lambda2*1e6)**(-1.) + 0.7709)*1e-5
   dny2dT = (0.5014*(lambda2*1e6)**(-3.) - 2.0030*(lambda2*1e6)**(-2.) + 3.3016*(lambda2*1e6)**(-1.) + 0.7498)*1e-5
   dnz2dT = (0.3896*(lambda2*1e6)**(-3.) - 1.3332*(lambda2*1e6)**(-2.) + 2.2762*(lambda2*1e6)**(-1.) + 2.1151)*1e-5

    nx1T0 = sqrt(3.0065+0.03901/((lambda1*1e6)**2.-0.04251)-0.01327*(lambda1*1e6)**2.) 
    ny1T0 = sqrt(3.0333+0.04154/((lambda1*1e6)**2.-0.04547)-0.01408*(lambda1*1e6)**2.) 
    nz1T0 = sqrt(3.3134+0.05694/((lambda1*1e6)**2.-0.05658)-0.01682*(lambda1*1e6)**2.) 

    nx2T0 = sqrt(3.0065+0.03901/((lambda2*1e6)**2.-0.04251)-0.01327*(lambda2*1e6)**2.) 
    ny2T0 = sqrt(3.0333+0.04154/((lambda2*1e6)**2.-0.04547)-0.01408*(lambda2*1e6)**2.) 
    nz2T0 = sqrt(3.3134+0.05694/((lambda2*1e6)**2.-0.05658)-0.01682*(lambda2*1e6)**2.) 

    aa1T0 = 1. / nx1T0 ** 2.
    bb1T0 = 1. / ny1T0 ** 2.
    cc1T0 = 1. / nz1T0 ** 2.

    aa2T0 = 1. / nx2T0 ** 2.
    bb2T0 = 1. / ny2T0 ** 2. 
    cc2T0 = 1. / nz2T0 ** 2.

    Term1 = sin(theta)**2. * cos(phi)**2.
    Term2 = sin(theta)**2. * sin(phi)**2.
    Term3 = cos(theta)**2.

     B1T0 = -Term1 * ( bb1T0 + cc1T0 )              &
		      -Term2 * ( aa1T0 + cc1T0 )              &
		      -Term3 * ( aa1T0 + bb1T0 ) 
     		 
     C1T0 =  Term1 * bb1T0 * cc1T0                  &
	 	      +Term2 * aa1T0 * cc1T0                  &
		      +Term3 * aa1T0 * bb1T0 

     B2T0 = -Term1 * ( bb2T0 + cc2T0 )              &
		      -Term2 * ( aa2T0 + cc2T0 )              &
		      -Term3 * ( aa2T0 + bb2T0 )
             
     C2T0 =  Term1 * bb2T0 * cc2T0                  &
		      +Term2 * aa2T0 * cc2T0                  &
		      +Term3 * aa2T0 * bb2T0 

    no1T0 = (2**0.5) / sqrt( -B1T0 - sqrt( B1T0 ** 2. - 4. * C1T0 ) )  
    ne1T0 = (2**0.5) / sqrt( -B1T0 + sqrt( B1T0 ** 2. - 4. * C1T0 ) ) 
    ne2T0 = (2**0.5) / sqrt( -B2T0 + sqrt( B2T0 ** 2. - 4. * C2T0 ) ) 

!----------------- The field equations constants

        c = 3.e8                                                                          !m/s
      r1f = 0.99                !radiation constant (back mirror)
      r2f = 0.99                !radiation constant (back mirror)  
      r3f = 0.99                !radiation constant (back mirror)
    omega = 2.*pi*c/lambda1                                                               !rad/s   
 epsilon0 = 8.85e-12                                                                      !C**2/N*m**2
     deff = 7.3e-12                                                                       !C**2/N*m**2 
        L = sqrt( ( epsilon0 * c**3. * pi * omegaf**2. ) / ( 4. * omega**2. * deff**2. * p ) ) 

!**********************************************************************************************************************
!                                        Arrays Allocattion 
!**********************************************************************************************************************
 
!------------------- Allocate Array temperature

allocate(temperature(1:2,0:nr,0:nz))
allocate(KT(0:nr,0:nz))

!------------------- Allocate Array phase

allocate(deltaphase(0:nr,0:nz))

!------------------- Allocate Array field equations

allocate(psi1p(0:nr,1:2))      ;allocate(psi1m(0:nr,1:2))               
allocate(psi2p(0:nr,1:2))      ;allocate(psi2m(0:nr,1:2))          
allocate(psi3p(0:nr,1:2))      ;allocate(psi3m(0:nr,1:2))        
           
!**********************************************************************************************************************
!                                     Giving Zero to Arrays
!********************************************************************************************************************** 
  
!------------------- Zero to Array temperature

forall (i=1:2,j=0:nr,k=0:nz)
                            temperature(i,j,k) = 0.
end forall

!---------
forall (j=0:nr,k=0:nz)
                            KT(j,k) = 0.
end forall

!-------------------- Zero to Array phase equation

forall (j=0:nr,k=0:nz)
                            deltaphase(j,k)=(0.,0.)        
end forall 

!--------------------- Zero to Array field equations

forall (j=0:nr,k=1:2)
                            psi1p(j,k)=(0.,0.)     ;psi1m(j,k)=(0.,0.)      
				                psi2p(j,k)=(0.,0.)     ;psi2m(j,k)=(0.,0.)     
                            psi3p(j,k)=(0.,0.)     ;psi3m(j,k)=(0.,0.)        			                               						                                                      
end forall !i

!**********************************************************************************************************************
!                                       Printing Constants     
!**********************************************************************************************************************

!------------------------------------------------ Common 
write(*,*)
write(*,*)'------- Common Constants ---------------------------------------------------'
write(*,*)
write(*,'(A13,F15.10 ,/)')  '      Power = ',P                

write(*,'(A13,2F15.10  )')  '        Ii = ',Ii                  
write(*,'(A13,I4       )')  '        Nt = ',Nt              
write(*,'(A13,I6       )')  '        Nr = ',Nr                
write(*,'(A13,I6       )')  '        Nz = ',Nz
write(*,'(A13,F15.10 ,/)')  '        pi = ',pi

write(*,'(A13,F15.10 ,/)')  '     timet = ',timet  

write(*,'(A13,f15.10   )')  '     gama1 = ',gama1               
write(*,'(A13,f15.10   )')  '     gama2 = ',gama2               
write(*,'(A13,f15.10 ,/)')  '     gama3 = ',gama3                      

write(*,'(A13,F15.10   )')  '    omegaf = ',omegaf          
write(*,'(A13,F15.10   )')  '    length = ',length          
write(*,'(A13,F15.10   )')  '    deltat = ',deltat          
write(*,'(A13,F15.10   )')  '    deltar = ',deltar          
write(*,'(A13,F15.10 ,/)')  '    deltaz = ',deltaz          

write(*,'(A13,F15.10 ,/)')  '    radius = ',radius           

write(*,'(A13,F15.10   )')  '   lambda1 = ',lambda1         
write(*,'(A13,f15.10 ,/)')  '   lambda2 = ',lambda2        
                                   

!------------------------------------------------ For Heat Equation 
write(*,*)
write(*,*)'------- Heat Equation Constants --------------------------------------------'
write(*,*)
write(*,'(A13,F15.10   )')  '         h = ',h               

write(*,'(A13,F15.10   )')  '        T0 = ',T0              
write(*,'(A13,F15.10 ,/)')  '        Cp = ',Cp              

write(*,'(A13,F15.10   )')  '       KT0 = ',KT0             
write(*,'(A13,F15.10 ,/)')  '       roh = ',roh              

write(*,'(A13,F15.10 ,/)')  '     sigma = ',sigma           

write(*,'(A13,F15.10 ,/)')  '  ePsilong = ',ePsilong        

write(*,'(A13,F15.10 ,/)')  ' stability = ',stability                                                                

write(*,*)'----------------------------------------------------------------------------'
write(*,'(A,\)')' Please press any key to continue '
read(*,*)

!------------------------------------------------ For Phase Equation 
write(*,*)
write(*,*)'------- Phase Equation Constants -------------------------------------------'
write(*,*)
write(*,'(A13,F15.10 ,/)')  '       phi = ',phi             

write(*,'(A13,f15.10   )')  '      B1T0 = ',B1T0                  
write(*,'(A13,F15.10   )')  '      B2T0 = ',B2T0            
write(*,'(A13,F15.10   )')  '      C1T0 = ',C1T0            
write(*,'(A13,F15.10 ,/)')  '      C2T0 = ',C2T0            

write(*,'(A13,F15.10 ,/)')  '     theta = ',theta           

write(*,'(A13,F15.10   )')  '     aa1T0 = ',aa1T0           
write(*,'(A13,F15.10   )')  '     bb1T0 = ',bb1T0           
write(*,'(A13,F15.10   )')  '     cc1T0 = ',cc1T0           
write(*,'(A13,F15.10   )')  '     aa2T0 = ',aa2T0           
write(*,'(A13,F15.10   )')  '     bb2T0 = ',bb2T0           
write(*,'(A13,F15.10   )')  '     cc2T0 = ',cc2T0           
write(*,'(A13,F15.10   )')  '     nx1T0 = ',nx1T0           
write(*,'(A13,F15.10   )')  '     ny1T0 = ',ny1T0           
write(*,'(A13,F15.10   )')  '     nz1T0 = ',nz1T0           
write(*,'(A13,F15.10   )')  '     nx2T0 = ',nx2T0           
write(*,'(A13,F15.10   )')  '     ny2T0 = ',ny2T0           
write(*,'(A13,F15.10 ,/)')  '     nz2T0 = ',nz2T0           

write(*,'(A13,F15.10   )')  '     no1T0 = ',no1T0           
write(*,'(A13,F15.10   )')  '     ne1T0 = ',ne1T0           
write(*,'(A13,F15.10 ,/)')  '     ne2T0 = ',ne2T0           

write(*,'(A13,F15.10   )')  '     Term1 = ',Term1           
write(*,'(A13,F15.10   )')  '     Term2 = ',Term2           
write(*,'(A13,F15.10 ,/)')  '     Term3 = ',Term3           

write(*,'(A13,F15.10   )')  '    dnx1dT = ',dnx1dT          
write(*,'(A13,F15.10   )')  '    dny1dT = ',dny1dT          
write(*,'(A13,F15.10   )')  '    dnz1dT = ',dnz1dT          
write(*,'(A13,F15.10   )')  '    dnx2dT = ',dnx2dT          
write(*,'(A13,F15.10   )')  '    dny2dT = ',dny2dT          
write(*,'(A13,F15.10 ,/)')  '    dnz2dT = ',dnz2dT          

write(*,*)'----------------------------------------------------------------------------'
write(*,'(A,\)')' Please press any key to continue '
read(*,*)

!------------------------------------------------ For fields Equation 
write(*,*)
write(*,*)'------- Field Equations Constants ------------------------------------------'
write(*,*)
write(*,'(A13,f15.3    )')  '         c = ',c 
write(*,'(A13,f35.30 ,/)')  '         L = ',L                 

write(*,'(A13,f35.30   )')  '      deff = ',deff                
write(*,'(A13,f25.5  ,/)')  '     omega = ',omega                            

write(*,'(A13,f25.20 ,/)')  '  ePsilon0 = ',ePsilon0           

  
!**********************************************************************************************************************
!                                   Main Block of the Program     
!**********************************************************************************************************************

! Display estimated execution time information
write(*,*)
write(*,*) '--- This code takes approximately 1 minute to execute on &
	        a medium-performance      laptop. Execution time may vary depending on &
			the system''s CPU, RAM, and        background tasks. ---!'	

write(*,*) 

!------------------------------------------------ Field equations
    no1rT = no1T0 
    ne1rT = ne1T0 
    ne2rT = ne2T0 
!--------------------------------------------- Run program for double passed 

  fmax = 1
 ! S = 100000000.
 ! A = 0.
 !do while (S>0.0001)

do f=1,fmax

write(*,*) 'step =',f                     ! Number of iteration

!---------------------------------------------elecp

write(*,*)   'start of field equation' ,f
   
   do k=0,nz
      z=k*deltaz 

      do j=1,nr-1
         r=j*deltar  
	     
!------------------------------------ Bounday conditions
!------------------------------------ For Field Equations

!----------------------------------------- psi1p

             if(k==0 )    psi1p(j ,1 ) = exp(-r**2. / omegaf**2.)     !for input  surface

	     if(k==0 )    psi1p(nr,1 ) = (0.,0.)                           !for (nr,0)

	     if(k==0 )    psi1p(0 ,1 ) = (1.,0.)                           !for (0,0)			 

             if(k==nz)    psi1p(nr,1 ) = (0.,0.)                      !for (nr,nz)

	     if(k==nz)    psi1p(0 ,1 ) = psi1p(1,1 )                       !for (0,nz)
		
	                  psi1p(0 ,1 ) = psi1p(1 ,1 )                      !for crystal axis
					      
	                  psi1p(nr,1 ) = (0.,0.)                           !for lateral surface 
 
!---------------------------------------- psi2p

              if(k==0 )    psi2p(j ,1 ) = exp(-r**2. / omegaf**2.)    !for input  surface

	      if(k==0 )    psi2p(nr,1 ) = (0.,0.)                          !for (nr,0)

	      if(k==0 )    psi2p(0 ,1 ) = (1.,0.)                          !for (0,0)			 

              if(k==nz)    psi2p(nr,1 ) = (0.,0.)                     !for (nr,nz)

	      if(k==nz)    psi2p(0 ,1 ) = psi2p(1,1 )                      !for (0,nz)
		
	                  psi2p(0 ,1 ) = psi2p(1 ,1 )                      !for crystal axis
					      
	                  psi2p(nr,1 ) = (0.,0.)                           !for lateral surface                       

!----------------------------------------- psi3p
		    
             if(k==0 )    psi3p(j ,1 ) = (0.,0.)                      !for input  surface
           
	     if(k==0 )    psi3p(0 ,1 ) = (0.,0.)                           !for (0 ,0 )
           
	     if(k==0 )    psi3p(nr,1 ) = (0.,0.)                           !for (nr,0 ) 

	     if(k==nz)    psi3p(nr,1 ) = (0.,0.)                           !for (nr,nr) 

	     if(k==nz)    psi3p(0 ,1 ) = psi3p(1,1)                        !for (0 ,nz) 

	                  psi3p(0 ,1 ) = psi3p(1 ,1 )                      !for crystal axis
					      
	                  psi3p(nr ,1) = (0.,0.)                           !for lateral surface 

!------------------------------------ End of Bounday conditions
         	
!------------------------------------ Field Equations
   
           !      cc2 = (Ii*c*deltaz)/(2.*no1rT*omega)  
			
		   	     cc3 = (deltaz*gama1)/2. 
  
                 cc4 = (Ii/L) * deltaz * sqrt( (1.)/(no1rT*ne1rT*ne2rT) )  

!-------------  
   
          !       dd2 = (Ii*c*deltaz)/(2.*ne1rT*omega) 

			     dd3 = (deltaz*gama2)/2.
   
                 dd4 = cc4  

!--------------
   
          !      ee2 = (Ii*c*deltaz)/(4.*ne2rT*omega) 
			
			    ee3 = (deltaz*gama3)/2. 
   
                ee4 = cc4     
!--------------			 
	     psi1p(j,2) =  psi1p(j  ,1)                                                                   &
					      
				        + cc2 * ( (psi1p(j+1,1) - 2 * psi1p(j,1) + psi1p(j-1,1))/(deltar**2) )           &

				        + cc2 * ( (psi1p(j+1,1) - psi1p(j-1,1))/(2*r*deltar) )                           &

			           - cc3 * psi1p(j,1)                                                               &

				        + cc4 * conjg(psi2p(j,1)) * psi3p(j,1) * exp(-Ii*deltaphase(j,k) )                                                                                              
           
		   !---------------------------------
	     psi2p(j,2) =  psi2p(j ,1)                                                                    &

				        + dd2 * ( (psi2p(j+1,1) - 2 * psi2p(j,1) + psi2p(j-1,1))/(deltar**2) )           &

				        + dd2 * ( (psi2p(j+1,1) - psi2p(j-1,1))/(2*r*deltar) )                           &

				        - dd3 * psi2p(j,1)                                                               &

				        + dd4 * conjg(psi1p(j,1)) * psi3p(j,1) * exp(-Ii*deltaphase(j,k) )                                                                                               
		   
		   !---------------------------------			
	     psi3p(j,2 ) =  psi3p(j  ,1)                                                                  &
					      
			            + ee2 * ( (psi3p(j+1,1) - 2 * psi3p(j,1) + psi3p(j-1,1))/(deltar**2) )          &

			            + ee2 * ( (psi3p(j+1,1) - psi3p(j-1,1))/(2*r*deltar) )                          &

			 	         - ee3 * psi3p(j,1)                                                              &

			            + ee4 * psi1p(j,1) * psi2p(j,1) * exp(Ii*deltaphase(j,k) )                                                                                               

			 
 

   
      end do !j

           !----------------------------------- print results
        psi12pz = psi1p(0,1) * conjg(psi1p(0,1))
        write(8,'(2x,f25.10,5x,f25.10)') z    , psi12pz * 100

	    psi22pz = psi2p(0,1) * conjg(psi2p(0,1))
        write(10,'(2x,f25.10,5x,f25.10)') z   , psi22pz * 100

	    psi32pz = psi3p(0,1) * conjg(psi3p(0,1))
        write(12,'(2x,f25.10,5x,f25.10)') z   , psi32pz * 100
!-----------------------------------------------------------

     do j=0,nr
	        psi1p(j,1)=psi1p(j,2) 
		psi2p(j,1)=psi2p(j,2)
		psi3p(j,1)=psi3p(j,2)
	end do !j
   end do !k
   !---------------------------------------------------------------------------------------------------------------- 
   !-------------------------------------------elecm (backward  )
   do k=nz,0,-1
      z=k*deltaz 
   
	  do j=1,nr-1
         r=j*deltar 
	     
!------------------------------------ Bounday conditions 
 
!------------------------------------psi1m
         if(k==0)     psi1m( 0 ,2) = psi1m(1,2)
	
         if(k==0)     psi1m(nr ,2) = (0.,0.)

         if(k==nz)    psi1m(j ,2 ) = r1f * psi1p(j,1)

         if(k==nz)    psi1m(0 ,2 ) = psi1m(1 ,2)
	
         if(k==nz)    psi1m(nr,2 ) = (0.,0.)

	              psi1m(0,2  ) = psi1m(1 ,2 )

	              psi1m(nr,2 ) = (0.,0.)


!------------------------------------psi2m
         if(k==0)     psi2m(0 ,2 ) = psi2m(1,2)
	
         if(k==0)     psi2m(nr,2)  = (0.,0.)

         if(k==nz)    psi2m(j ,2 ) = r2f * psi2p(j,1)

         if(k==nz)     psi2m(0 ,2) = psi2m(1 ,2)
	
         if(k==nz)     psi2m(nr,2) = (0.,0.)

				      psi2m(0 ,2 ) = psi2m(1 ,2 )

	                  psi2m(nr,2 ) = (0.,0.)


!-------------------------------------psi3m
         if(k==0)     psi3m(0 ,2 ) = psi3m(1,2)

         if(k==0)     psi3m(nr ,2) = (0.,0.) 

         if(k==nz)    psi3m(j,2  ) = r3f * psi3p(j,1)

         if(k==nz)    psi3m(nr,2 ) = (0.,0.)

         if(k==nz)    psi3m(0,2  ) =  psi3m(1,2 )

                      psi3m(0 ,2 ) = psi3m(1,2)

	              psi3m(nr,2 ) = (0.,0.)

!------------------------------------ End of Bounday conditions   

         psi1m(j,1 ) =  psi1m(j  ,2)                                                                   &
					       
			             + cc2 * ( (psi1m(j+1,2) - 2 * psi1m(j,2) + psi1m(j-1,2))/(deltar**2) )           &

			             + cc2 * ( (psi1m(j+1,2) - psi1m(j-1,2))/(2*r*deltar) )                           &

			             - cc3 * psi1m(j,2)                                                               & 

		                + cc4 * conjg(psi2m(j,2)) * psi3m(j,2) * exp(Ii*deltaphase(j,k))                                                                                              
    
!---------------------------------
         psi2m(j,1 )=  psi2m(j  ,2)                                                                    &

			            + dd2 * ( (psi2m(j+1,2) - 2 * psi2m(j,2) + psi2m(j-1,2))/(deltar**2) )            &

			            + dd2 * ( (psi2m(j+1,2) - psi2m(j-1,2))/(2*r*deltar) )                            &

		    	         - dd3 * psi2m(j,2)                                                                &

			            + dd4 * conjg(psi1m(j,2)) * psi3m(j,2) * exp(Ii*deltaphase(j,k) )                                                                                               
		   
!---------------------------------			
         psi3m(j,1 ) =  psi3m(j  ,2)                                                                   &
					      
			             + ee2 * ( (psi3m(j+1,2) - 2 * psi3m(j,2) + psi3m(j-1,2))/(deltar**2) )           &

			             + ee2 * ( (psi3m(j+1,2) - psi3m(j-1,2))/(2*r*deltar) )                           &

		  	             - ee3 * psi3m(j,2)                                                               &

	                   + ee4 * psi1m(j,2) * psi2m(j,2) * exp(-Ii*deltaphase(j,k) )          
				                                                  

			 
      end do !j
!----------------------------------print results

	  psi12mz = psi1m(0,2) * conjg(psi1m(0,2))
	  write(14,'(2x,f25.10,5x,f25.10)')  z     ,psi12mz*100
	  
	  psi22mz = psi2m(0,2) * conjg(psi2m(0,2))
	  write(16,'(2x,f25.10,5x,f25.10)')  z     ,psi22mz*100
	  
	  psi32mz = psi3m(0,2) * conjg(psi3m(0,2))
	  write(18,'(2x,f25.10,5x,f25.10)')  z     ,psi32mz*100

!------------------------------------
      do j=0,nr

	     psi1m(j,2)=psi1m(j,1)
	     psi2m(j,2)=psi2m(j,1)
	     psi3m(j,2)=psi3m(j,1)

	  end do !j
	   
   end do !k
!-------------------------------
         
   write(*,*) 'end of field equation'  ,f       						   
!============================== heat equation 
    do j=0,nr
       do k=0,nz
	  
          temperature(1,j,k) = T0
	      KT(j,k) = KT0
	   
      end do !k
   end do !j
	     
!----------------
write(*,*) 'start of temperature ' , f
!-----------------
!-----------------------------
   do i=0,nt
      t = deltat*i  
    	  
!--------------------------------------
 
         do j=1,nr-1
            r = j*deltar
			   
	        do k=1,nz-1
	           z=k*deltaz 
!---------------------------------------for save time

               aa1 = (h*deltaz)/(KT(j,k))
            
	           aa2 = (epsilong*sigma*deltaz)/(KT(j,k))
            
	           aa3 = ( deltat/(roh*Cp) ) * KT(j,k)           

               aa4 = ( (deltat)/(roh * Cp) ) * ( ( p)/( pi * omegaf**2.) )

	           aa5 = deltat /(4. * roh * Cp  )

!----------------------------------------------------------boundry condition

!-------------------------------- boundry condition for head equation
		 
	           temperature(1,0 ,k)  = temperature(1,1,k)                !Thermal insulation condition for crystal axis

               temperature(1,nr,k)  = T0                               !Temperature-fixed condition for lateral surface

               temperature(1,j ,0)  = temperature(1,j,1) - aa1*( temperature(1,j,1) - Tinf )                         &
			                                                - aa2*( temperature(1,j,1)**4. - Tamb**4. )
			                                                              !Convection & Radiation condition for input  surface

               temperature(1,j,nz)  = temperature(1,j,nz-1) - aa1*( temperature(1,j,nz-1) - Tinf )                   &
			                                                   - aa2*( temperature(1,j,nz-1)**4. - Tamb**4. )
			                                                              !Convection & Radiation condition for output surface
            
               temperature(1,0 ,0 ) = temperature(1,0,1) - aa1*( temperature(1,0,1) - Tinf )                         &
			                                                - aa2*( temperature(1,0,1)**4. - Tamb**4. ) 
			                                                              !Convection & Radiation condition for ( 0,0 )

               temperature(1,0 ,nz) = temperature(1,0,nz-1) - aa1*( temperature(1,0,nz-1) - Tinf )                   &
			                                                   - aa2*( temperature(1,0,nz-1)**4. - Tamb**4. ) 
			                                                              !Convection & Radiation condition for ( 0,nz)

               temperature(1,nr,0 ) = T0                               !Temperature-fixed condition for (nr,0 )

               temperature(1,nr,nz) = T0                               !Temperature-fixed condition for (nr,nz)

!------------------------------------ Ending of Boundary conditions


!------------------------------------------------Heat equation
									 
									
		       temperature(2,j,k) = temperature(1,j,k)                                                                                                                                                                                                                                                                                                                                                                                                                                                 &
		                      
        			 	          + aa3 * ( (temperature(1,j+1,k) -  temperature(1,j-1,k))/(2.*r*deltar)                            &
								  
					             + (temperature(1,j+1,k) -2.*temperature(1,j,k) + temperature(1,j-1,k))/(deltar**2.) )             & 
 
                            + aa3 * ( (temperature(1,j,k-1) -2.*temperature(1,j,k) + temperature(1,j,k+1))/(deltaz**2.) )     &                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

					  	          + aa4 * gama1 * psi1p(j,1)*conjg(psi1p(j,1))                                                      &                     
								 
						          + aa4 * gama1 * psi1m(j,2)*conjg(psi1m(j,2))                                                      &                                                
							                                     
					             + aa4 * gama2 * psi2p(j,1)*conjg(psi2p(j,1))                                                      &
							  										  
						          + aa4 * gama2 * psi2m(j,2)*conjg(psi2m(j,2))                                                      &            
																	       
						          + aa4 * 2. * gama3 * psi3p(j,1)*conjg(psi3p(j,1))                                                 &
								 
                            + aa4 * 2. * gama3 * psi3m(j,2)*conjg(psi3m(j,2))                                                 &							   
								 
						          + aa5 * ( (temperature(1,j+1,k) - temperature(1,j-1,k)) * ( KT(j+1,k) - KT(j-1,k)) / deltar**2.)  &
										   
	    				          + aa5 * ( (temperature(1,j,k+1) - temperature(1,j,k-1)) * ( KT(j,k+1) - KT(j,k-1)) / deltaz**2.) 
                 		
		                                                                                                                                                                                                                                                                                                                                                 						 
	         !----------------------------------- Phase Equation constants

   		      nx1r0T = nx1T0 + dnx1dT * ( temperature(1,0,k) - T0 )
	            ny1r0T = ny1T0 + dny1dT * ( temperature(1,0,k) - T0 )
	            nz1r0T = nz1T0 + dnz1dT * ( temperature(1,0,k) - T0 )				
    	 
	            nx2r0T = nx2T0 + dnx2dT * ( temperature(1,0,k) - T0 )
	            ny2r0T = ny2T0 + dny2dT * ( temperature(1,0,k) - T0 )
	  	         nz2r0T = nz2T0 + dnz2dT * ( temperature(1,0,k) - T0 )
        		   
	   	      aa1r0T  = 1. / ( nx1r0T )**2. 
               bb1r0T  = 1. / ( ny1r0T )**2. 
               cc1r0T  = 1. / ( nz1r0T )**2. 

               aa2r0T  = 1. / ( nx2r0T )**2. 
               bb2r0T  = 1. / ( ny2r0T )**2. 
               cc2r0T  = 1. / ( nz2r0T )**2.       
  
    	         B1r0T  =  -Term1 * ( bb1r0T + cc1r0T )                &
			                -Term2 * ( aa1r0T + cc1r0T )                &
			                -Term3 * ( aa1r0T + bb1r0T ) 
     		 
	   	      C1r0T  =  Term1 * bb1r0T * cc1r0T                     &
			               +Term2 * aa1r0T * cc1r0T                     &
			               +Term3 * aa1r0T * bb1r0T 

               B2r0T  = -Term1 * ( bb2r0T + cc2r0T )                 &
			               -Term2 * ( aa2r0T + cc2r0T )                 &
			               -Term3 * ( aa2r0T + bb2r0T )
             
	 	         C2r0T  =  Term1 * bb2r0T * cc2r0T                     &
			               +Term2 * aa2r0T * cc2r0T                     &
			               +Term3 * aa2r0T * bb2r0T 

             
 	   	      no1r0T = (2.**0.5) / sqrt( -B1r0T  - sqrt( B1r0T**2. - 4.*C1r0T ) )  
               ne1r0T = (2.**0.5) / sqrt( -B1r0T  + sqrt( B1r0T**2. - 4.*C1r0T ) )
               ne2r0T = (2.**0.5) / sqrt( -B2r0T  + sqrt( B2r0T**2. - 4.*C2r0T ) ) 
			 
	            deltano1r0T = no1r0T - no1T0
               deltane1r0T = ne1r0T - ne1T0
               deltane2r0T = ne2r0T - ne2T0
               
!----------------------------------- Phase Equation constants

  		       nx1rT = nx1T0 + dnx1dT * ( temperature(1,j,k) - T0 )
	 	       ny1rT = ny1T0 + dny1dT * ( temperature(1,j,k) - T0 )
	 	       nz1rT = nz1T0 + dnz1dT * ( temperature(1,j,k) - T0 )
			 
		       nx2rT = nx2T0 + dnx2dT * ( temperature(1,j,k) - T0 )
		       ny2rT = ny2T0 + dny2dT * ( temperature(1,j,k) - T0 )
	 	       nz2rT = nz2T0 + dnz2dT * ( temperature(1,j,k) - T0 )
          		   
		       aa1rT  = 1. / ( nx1rT )**2. 
             bb1rT  = 1. / ( ny1rT )**2. 
             cc1rT  = 1. / ( nz1rT )**2. 

             aa2rT  = 1. / ( nx2rT )**2. 
             bb2rT  = 1. / ( ny2rT )**2. 
             cc2rT  = 1. / ( nz2rT )**2.       

	          B1rT  =  -Term1 * ( bb1rT + cc1rT )                 &
			             -Term2 * ( aa1rT + cc1rT )                 &
			             -Term3 * ( aa1rT + bb1rT ) 
     		 
		       C1rT  =  Term1 * bb1rT * cc1rT                      &
			            +Term2 * aa1rT * cc1rT                      &
			            +Term3 * aa1rT * bb1rT 

             B2rT  = -Term1 * ( bb2rT + cc2rT )                  &
			              -Term2 * ( aa2rT + cc2rT )                &
			              -Term3 * ( aa2rT + bb2rT )
             
		       C2rT  =   Term1 * bb2rT * cc2rT                     &
			             +Term2 * aa2rT * cc2rT                     &
			             +Term3 * aa2rT * bb2rT 

             
 	           no1rT = (2**0.5) / sqrt( -B1rT  - sqrt( B1rT**2. - 4.*C1rT ) )  
              ne1rT = (2**0.5) / sqrt( -B1rT  + sqrt( B1rT**2. - 4.*C1rT ) )
              ne2rT = (2**0.5) / sqrt( -B2rT  + sqrt( B2rT**2. - 4.*C2rT ) ) 
			 
              deltano1rT = no1rT - no1T0
              deltane1rT = ne1rT - ne1T0
              deltane2rT = ne2rT - ne2T0


       	!------------------------------------ For Phase Equation
               deltaphase(j ,0  ) = (0.,0.)                               !for input surface
            
               deltaphase(nr,k  ) = (0.,0.)                               !for lateral surface

       	      deltaphase(j,nz  ) = deltaphase(j,nz-1)                                             &          			                        
				                      + ( 2.*pi*deltaz / lambda1 )                                     &
					                   * ( deltano1rT + deltane1rT - 2.*deltane2rT )    
		                                                                  !for output surface
	

              deltaphase(nr,0  ) = (0.,0.)                               !for (nr,0 )

		        deltaphase(nr,nz ) = (0.,0.)                               !for (nr,nz) 

	           deltaphase(0 ,0  ) = (0.,0.)                               !for ( 0,0 )

			   deltaphase(0 ,k  ) = deltaphase(0,k-1)                                                 &           
				                   + ( 2.*pi*deltaz / lambda1 )                                        &
							          * ( deltano1r0T  + deltane1r0T  - 2.*deltane2r0T  )    
								                                          !for crystal axis    
								   
			   deltaphase(0 ,nz  ) = deltaphase(0,nz-1)                                               &           
				                   + ( 2.*pi*deltaz / lambda1 )                                        &
						    	   * ( deltano1r0T  + deltane1r0T  - 2.*deltane2r0T  ) 	
								                                          !for crystal axis				                            
 !----------------------------------- Phase Equation
			  
		       deltaphase(j,k-1+1) = deltaphase(j,k-1)                                               &
				                     + ( 2.*pi*deltaz / lambda1 )                                      & 
							            * ( deltano1rT + deltane1rT - 2.*deltane2rT )    
	        end do !j						 
	 					                              
         end do !k	
	
   !============================================= Print Results for each deltat
   m= 0. 
   !if (S<0.0001 ) then   
   !--------------------------------------------- For Heat Equation
      t= i*deltat 
      write(1,'(2x,f25.10,5x,f25.10)')  t , temperature(1,0,0)

   !--------------------------------------------- For Phase Equation
      t= i*deltat 
      write(4,'(2x,f25.10,5x,2f25.10)') t , deltaphase(0,0) 
      !m= m+1
   
   !=============================================
   !end if 
   !--------------------------------------------- End-temprature of each deltat  ==> Initial temperature for next deltat
   do j=1,nr-1
      do k=1,nz-1
 	  
         temperature(1,j,k) = temperature(2,j,k)

      end do !k
   end do !j	 
   
   !---------------------------------------------
   
   do j=0,nr
      do k=0,nz
	  
	  KT(j,k) = KT0 * T0 / Temperature(1,j,k) 
	  
	  end do ! k
  end do !j 	      
   !--------------------------------------------- End-phase-changes of each pulse ==> Initial-phase-changes for next pulse
   do j=0,nr
      do k=0,nz
 	  
         deltaphase(j,k) = deltaphase(j,k)
	   
      end do !k
   end do !j	     
   !---------------------------------------------

   end do !i
   write(*,*) 'end of temperature ' , f
 !  f=f+1
   write(*,*) 'm=',m
   !--------------------------------------------- End of run 

end do !f, do while


!**********************************************************************************************************************
!                                        Printing Results     
!**********************************************************************************************************************
write(*,*)  'print result'
!------------------------------------------------For heat equation
do j=1,nr
   r=j*deltar 
   write(2,'(2x,f25.10,5x,f25.10)')  r , temperature(1,j,1)
 end do !j      						   

!------------------------------------------------
 do k=0,nz
    z=k*deltaz 
    write(3,'(2x,f25.10,5x,f25.10)')  z , temperature(1,1,k)
 end do !k      						   


!------------------------------------------------ For Phase Equation
 do j=1,nr
    r=j*deltar 
    write(5,'(2x,f25.10,5x,2f25.10)') r , deltaphase(j,1)
 end do !j      						   

!------------------------------------------------
 do k=0,nz
    z=k*deltaz 
    write(6,'(2x,f25.10,5x,2f25.10)') z , deltaphase(1,k)
 end do !k 

!------------------------------------------------ For feild equations    	
do j=0,nr
   r=j*deltar 
   psi12pr = psi1p(j,1)* conjg(psi1p(j,1))
   write(7,'(2x,f25.10,5x,f25.10)') r    , psi12pr*100
end do !j      						   

!------------------

do j=1,nr
   r=j*deltar 
   psi22pr = psi2p(j,1) * conjg(psi2p(j,1))
   write(9,'(2x,f25.10,5x,f25.10)') r   , psi22pr*100
end do !j      						   

!------------------

do j=1,nr
   r=j*deltar 
   psi32pr = psi3p(j,1) * conjg(psi3p(j,1))
   write(11,'(2x,f25.10,5x,f25.10)') r   , psi32pr*100
end do !j      						   

!------------------------------------------------------------------

do j=1,nr
   r=j*deltar 
   psi12mr = psi1m(j,2) * conjg(psi1m(j,2))
   write(13,'(2x,f25.10,5x,f25.10)') r , psi12mr*100
end do !j      						   

!------------------

do j=1,nr
   r=j*deltar 
   psi22mr = psi2m(j,2) * conjg(psi2m(j,2))
   write(15,'(2x,f25.10,5x,f25.10)') r , psi22mr*100
end do !j      						   

!------------------

do j=1,nr
   r=j*deltar 
   psi32mr = psi3m(j,2) * conjg(psi3m(j,2))
   write(17,'(2x,f25.10,5x,f25.10)') r  , psi32mr*100
end do !j      						  


!**********************************************************************************************************************
!                              Closing Files and Ending the Program 
!**********************************************************************************************************************

!----------------- for temperature 

close(1)
close(2)
close(3)

!----------------- for phase equation

close(4)
close(5)
close(6)

!----------------- for field equations

close(7)
close(8)

close(9) 
close(10)

close(11)
close(12)

close(13)
close(14)

close(15)
close(16)

close(17)
close(18)

!----------------------

write(*,*) 
write(*,*) '---- The results are stored in `.plt` format.                                  &
	         If a different format is required, users can set the desried extension in      &
			   "Determine Filenames & Open files" section of the code or rename the file      & 
			   manually and open it with their preferred software. ----!'	

			
write(*,*) 	
write(*,*) '---- Program Completed ----!'

end program Elec_temp_phase_CW                     

