! Note: 
!     We have aligned the code with VS Code's formatting standards.
!     But it may appears cluttered on GitHub due to formatting differences.


!            ************************************************************************************
!            *                                                                                  *
!            * File name:                                                                       *
!            *     4. Phase Mismatch _ Pulsed Wave Gaussian _ Computational.F90                 *
!            *                                                                                  *
!            * This Fortran code is developed specifically for the article titled:              *
!            *     Thermally Induced Phase Mismatching in a Repetitively Gaussian Pulsed        *
!            *     Pumping KTP Crystal: A Spatio-Temporal Treatment                             *
!            *                                                                                  *
!            * Cite Us:                                                                         *
!            *     Rezaee, M.M., Sabaeian, M., Motazedian, A., Jalil-Abadi, F.S., Askari, H.    *
!            *     and Khazrk, I., 2015. Thermally induced phase mismatching in a repetitively  *
!            *     Gaussian pulsed pumping KTP crystal: a spatiotemporal treatment. Applied     *
!            *     Optics, 54(15), pp.4781-4788.                                                *
!            *                                                                                  *
!            ************************************************************************************

program Temp_Phase_PW

implicit none

!**********************************************************************************************************************
!                                       Variables Definition
!**********************************************************************************************************************

!------------------------------------------------ Thermal Variables
integer       i          ,j          ,k          ,l                                                                  &
             ,nt         ,nr         ,nz         ,Np                                                                       


real*8        t          ,z          ,E          ,h          ,r          ,G          ,P                              &
             ,T0         ,pi         ,Cp         ,tp         ,Q0                                                     &
             ,roh        ,aa1        ,aa2        ,aa3        ,aa4        ,aa5                                        &
             ,KT0        ,freq       ,gama       ,Tinf       ,Tamb                                                   & 
             ,timet      ,sigma                                                                                      &
			 ,omegaf     ,length     ,deltar     ,deltaz     ,deltat     ,radius                                     &
			 ,epsilong   ,tbetween                                                                                   &
			 ,stability                                                                                              &
			 ,r_integral ,z_integral                                                                                 &
             ,temperature[allocatable](:,:,:)            ,KT[allocatable](:,:)  


character*30  EE                                                                                                     &
             ,Npf          ,tpf                                                                                      &
			 ,freqf                                                                                                  &
             ,filenameTt   ,filenameTr   ,filenameTz     

!-------------------------------------- Phase Variables
real*8        phi                                                                                                    &
             ,B1T0         ,B2T0         ,C1T0          ,C2T0         ,B1rT         ,B2rT         ,C1rT              &
			 ,C2rT                                                                                                   &		 
             ,Phase                                                                                                  &
			 ,aa1T0        ,bb1T0        ,cc1T0         ,nx1T0        ,ny1T0        ,nz1T0        ,Term1             &
			 ,aa2T0        ,bb2T0        ,cc2T0         ,nx2T0        ,ny2T0        ,nz2T0        ,Term2             &
			 ,aa1rT        ,bb1rT        ,cc1rT         ,nx1rT        ,ny1rT        ,nz1rT        ,Term3             &
			 ,aa2rT        ,bb2rT        ,cc2rT         ,nx2rT        ,ny2rT        ,nz2rT        ,theta             &
             ,B1r0T        ,B2r0T        ,C1r0T         ,C2r0T        ,no1T0        ,ne1T0        ,ne2T0             &
             ,no1rT        ,ne1rT        ,ne2rT                                                                      &
                                                          
			 ,aa1r0T       ,bb1r0T       ,cc1r0T        ,nx1r0T       ,ny1r0T       ,nz1r0T       ,dnx1dT            &
			 ,aa2r0T       ,bb2r0T       ,cc2r0T        ,nx2r0T       ,ny2r0T       ,nz2r0T       ,dnx2dT            &
             ,dny1dT       ,dny2dT       ,dnz1dT        ,dnz2dT       ,no1r0T       ,ne1r0T       ,ne2r0T            &   

             ,lambda1      ,lambda2                                                                                  &
             ,Phasemin                                                                                               &

			 ,deltano1rT   ,deltane1rT   ,deltane2rT                                                                 &   
			 ,deltano1r0T  ,deltane1r0T  ,deltane2r0T                                                                 	                                                                                      		 

complex*8     deltaphase[allocatable](:,:)                                                                           

character*35  filenamePt   ,filenamePr   ,filenamePz                                                                 &
			 ,plot_extension                                                              

!**********************************************************************************************************************
!                                    Giving Zero to variables
!**********************************************************************************************************************

!------------------------------------------------ Giving Zero to Thermal Variables
            i = 0           ;j = 0          ;k = 0          ;l = 0  
           nt = 0          ;nr = 0         ;nz = 0         ;Np = 0                                   

            t = 0.          ;z = 0.         ;E = 0.         ;h = 0.         ;r = 0.        ;G = 0.       ;P = 0.
           T0 = 0.         ;pi = 0.        ;cp = 0.        ;tp = 0.        ;Q0 = 0.
          roh = 0.        ;KT0 = 0.       ;aa1 = 0.       ;aa2 = 0.       ;aa3 = 0.      ;aa4 = 0.     ;aa5 = 0.
         freq = 0.       ;gama = 0.      ;Tinf = 0.      ;Tamb = 0.       
        timet = 0.      ;sigma = 0.
       omegaf = 0.     ;length = 0.    ;deltar = 0.    ;deltaz = 0.    ;deltat = 0.   ;radius = 0.
     epsilong = 0.   ;tbetween = 0.
    stability = 0.                                                             
   r_integral = 0. ;z_integral = 0.   
!------------------------------------------------ Giving Zero to Phase Variables
          phi = 0.                                                                                         

         B1T0 = 0.        ;B2T0 = 0.        ;C1T0 = 0.    ;C2T0 = 0.     ;B1rT = 0.    ;B2rT = 0.    ;C1rT = 0.
		 C2rT = 0.                                                                                        

        Phase = 0.                                                                                       
		aa1T0 = 0.       ;bb1T0 = 0.       ;cc1T0 = 0.   ;nx1T0 = 0.    ;ny1T0 = 0.   ;nz1T0 = 0.   ;Term1 = 0.
		aa2T0 = 0.       ;bb2T0 = 0.       ;cc2T0 = 0.   ;nx2T0 = 0.    ;ny2T0 = 0.   ;nz2T0 = 0.   ;Term2 = 0.
        aa1rT = 0.       ;bb1rT = 0.       ;cc1rT = 0.   ;nx1rT = 0.    ;ny1rT = 0.   ;nz1rT = 0.   ;Term3 = 0.
		aa2rT = 0.       ;bb2rT = 0.       ;cc2rT = 0.   ;nx2rT = 0.    ;ny2rT = 0.   ;nz2rT = 0.   ;theta = 0.   
		B1r0T = 0.       ;B2r0T = 0.       ;C1r0T = 0.   ;C2r0T = 0.    ;no1T0 = 0    ;ne1T0 = 0    ;ne2T0 = 0.
        no1rT = 0.       ;ne1rT = 0.       ;ne2rT = 0.                                                              

       aa1r0T = 0.      ;bb1r0T = 0.      ;cc1r0T = 0.   ;nx1r0T = 0.  ;ny1r0T = 0.  ;nz1r0T = 0.  ;dnx1dT = 0.  
       aa2r0T = 0.      ;bb2r0T = 0.      ;cc2r0T = 0.   ;nx2r0T = 0.  ;ny2r0T = 0.  ;nz2r0T = 0.  ;dnx2dT = 0.  
       dny1dT = 0.      ;dny2dT = 0.      ;dnz1dT = 0.   ;dnz2dT = 0.  ;no1r0T = 0.  ;ne1r0T = 0.  ;ne2r0T = 0.   

      lambda1 = 0.     ;lambda2 = 0.             
     Phasemin = 0.                                                                                     

   deltano1rT = 0.  ;deltane1rT = 0.  ;deltane2rT = 0.                                                 
  deltano1r0T = 0. ;deltane1r0T = 0. ;deltane2r0T = 0.                                               

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

!write(*,'(/,2x,a,\)') '                   Enter the Number of Pulses : '
 !read(*,*) Np
!write(*,'(/,2x,a,\)') 'Enter the Pulses' value without decimal point : '
 !read(*,*) Npf

!write(*,'(/,2x,a,\)') '                            Enter the tp : '
 !read(*,*) tp
!write(*,'(/,2x,a,\)') 'Enter the tp value without decimal point : '
 !read(*,*) tpf

!write(*,'(/,2x,a,\)') '                      Enter the frequency value : '
 !read(*,*) freq
!write(*,'(/,2x,a,\)') 'Enter the frequency value without decimal point : '
 !read(*,*) freqf
 
! For Calculation
   E = 0.09  
  Np = 1
  tp = 50e-6
freq = 500

! For Generating Filenames based on the values above
   EE = '009'   
  Npf = '1'
  tpf = '50'
freqf = '500'


!**********************************************************************************************************************
!                          Determination of Filenames and Opening files
!**********************************************************************************************************************

! Note:
!      To achieve both efficiency and clarity in managing output data,
!      below, we generate filenames based on input information.

plot_extension = '.plt'

!------------------------------------------------ Heat Equation Files
filenameTt = 'E_'//trim(EE)//'_f_'//trim(freqf)//'_Np_'//trim(Npf)//'_tp_'//trim(tpf)//'_T_t'//plot_extension
open(1,file=filenameTt)
write(1,'(/,a,/)')   ! ' variables=         "t"                             "temperature"'

filenameTr = 'E_'//trim(EE)//'_f_'//trim(freqf)//'_Np_'//trim(Npf)//'_tp_'//trim(tpf)//'_T_r'//plot_extension
open(2,file=filenameTr)
write(2,'(/,a,/)')    !' variables=         "r"                             "temperature"'

filenameTz = 'E_'//trim(EE)//'_f_'//trim(freqf)//'_Np_'//trim(Npf)//'_tp_'//trim(tpf)//'_T_z'//plot_extension
open(3,file=filenameTz)
write(3,'(/,a,/)')    !' variables=         "z"                             "temperature"' 

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :', filenameTt, filenameTr, filenameTz
 write(*,'(A,\)')' Please press any key to continue '
 read(*,*)

!------------------------------------------------ Phase Equation Files
filenamePt = 'E_'//trim(EE)//'_f_'//trim(freqf)//'_Np_'//trim(Npf)//'_tp_'//trim(tpf)//'_P_t'//plot_extension
open(4,file=filenamePt)
write(4,'(/,a,/)')    !' variables=         "t"          "deltaphase_real"      "deltaphase_imaginary"'

filenamePr = 'E_'//trim(EE)//'_f_'//trim(freqf)//'_Np_'//trim(Npf)//'_tp_'//trim(tpf)//'_P_r'//plot_extension
open(5,file=filenamePr)
write(5,'(/,a,/)')    !' variables=         "r"          "deltaphase_real"      "deltaphase_imaginary"'

filenamePz = 'E_'//trim(EE)//'_f_'//trim(freqf)//'_Np_'//trim(Npf)//'_tp_'//trim(tpf)//'_P_z'//plot_extension
open(6,file=filenamePz)
write(6,'(/,a,/)')    !' variables=         "z"           "deltaphase_real"      "deltaphase_imaginary"'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :', filenamePt, filenamePr, filenamePz
 write(*,'(A,\)')' Please press any key to continue '
 read(*,*)

!**********************************************************************************************************************
!                                           Constants
!**********************************************************************************************************************

!------------------------------------------------ Thermal properties
        h = 10                  !heat transfer coefficient (convection - cylinder)          W/(m^2.K)
       T0 = 300.                !initial temperature                                        K
       pi = 4*atan(1.)                                                                     !dimensionless
       Cp = 728.016             !heat capacity at constant pressure                         J/(kg.K)

      KT0 = 13.                 !thermal conductivity of KTP crystal                        W/(m.K)
      roh = 2945.               !mass density                                               kg/m^3

     gama = 4.                  !absorption coefficient                                     1/m
     Tinf = 300.
     Tamb = 300.
    sigma = 5.669e-8            !Stephan-Bultzman constant                                  W/(m^2.K^4) 

 tbetween = 1/freq !1/(5*freq)                                                             !s
    timet = Np*tbetween                                                                    !s
   deltat = tp / 10                                                                        !s     
       nt = int(tbetween/deltat)                                                           !dimensionless

       nz = 58 !150 !200                                                                   !dimensionless
   length = 0.02                !length of crystal                                          m 
   deltaz = length/nz                                                                      !m

   radius = 0.005               !radius of crystal                                         !m
   omegaf = 100.e-6             !spot size                                                 !m
   deltar = omegaf/10                                                                      !m
       nr = int(radius/deltar)                                                             !dimensionless 

 epsilong = 0.9                 !surface emissivity                                        !dimensionless

stability = ( (2*KT0*deltat)/(roh*Cp) ) * ( (deltar**2+deltaz**2)/(deltar**2*deltaz**2) )  !stability coefficient  

!------------------------------------------------ Phase properties
       phi = 24.77*pi/180
     theta =    90*pi/180

   lambda1 = 1064e-9  
   lambda2 =  532e-9
    dnx1dT = (0.1323*(lambda1*1e6)**(-3) - 0.4385*(lambda1*1e6)**(-2) + 1.2307*(lambda1*1e6)**(-1) + 0.7709)*1e-5
    dny1dT = (0.5014*(lambda1*1e6)**(-3) - 2.0030*(lambda1*1e6)**(-2) + 3.3016*(lambda1*1e6)**(-1) + 0.7498)*1e-5
    dnz1dT = (0.3896*(lambda1*1e6)**(-3) - 1.3332*(lambda1*1e6)**(-2) + 2.2762*(lambda1*1e6)**(-1) + 2.1151)*1e-5

    dnx2dT = (0.1323*(lambda2*1e6)**(-3) - 0.4385*(lambda2*1e6)**(-2) + 1.2307*(lambda2*1e6)**(-1) + 0.7709)*1e-5
    dny2dT = (0.5014*(lambda2*1e6)**(-3) - 2.0030*(lambda2*1e6)**(-2) + 3.3016*(lambda2*1e6)**(-1) + 0.7498)*1e-5
    dnz2dT = (0.3896*(lambda2*1e6)**(-3) - 1.3332*(lambda2*1e6)**(-2) + 2.2762*(lambda2*1e6)**(-1) + 2.1151)*1e-5

     nx1T0 = sqrt(3.0065+0.03901/((lambda1*1e6)**2-0.04251)-0.01327*(lambda1*1e6)**2) 
     ny1T0 = sqrt(3.0333+0.04154/((lambda1*1e6)**2-0.04547)-0.01408*(lambda1*1e6)**2) 
     nz1T0 = sqrt(3.3134+0.05694/((lambda1*1e6)**2-0.05658)-0.01682*(lambda1*1e6)**2) 

     nx2T0 = sqrt(3.0065+0.03901/((lambda2*1e6)**2-0.04251)-0.01327*(lambda2*1e6)**2) 
     ny2T0 = sqrt(3.0333+0.04154/((lambda2*1e6)**2-0.04547)-0.01408*(lambda2*1e6)**2) 
     nz2T0 = sqrt(3.3134+0.05694/((lambda2*1e6)**2-0.05658)-0.01682*(lambda2*1e6)**2) 

     aa1T0 = 1 / nx1T0 ** 2
     bb1T0 = 1 / ny1T0 ** 2
     cc1T0 = 1 / nz1T0 ** 2

     aa2T0 = 1 / nx2T0 ** 2
     bb2T0 = 1 / ny2T0 ** 2 
     cc2T0 = 1 / nz2T0 ** 2

     Term1 = sin(theta)**2 * cos(phi)**2
     Term2 = sin(theta)**2 * sin(phi)**2
     Term3 = cos(theta)**2

      B1T0 = -Term1 * ( bb1T0 + cc1T0 )            &
	         -Term2 * ( aa1T0 + cc1T0 )            &
	         -Term3 * ( aa1T0 + bb1T0 ) 
     		 
      C1T0 =  Term1 * bb1T0 * cc1T0                &
	         +Term2 * aa1T0 * cc1T0                &
	         +Term3 * aa1T0 * bb1T0 


      B2T0 = -Term1 * ( bb2T0 + cc2T0 )            &
	     	 -Term2 * ( aa2T0 + cc2T0 )            &
	     	 -Term3 * ( aa2T0 + bb2T0 )
             
      C2T0 =  Term1 * bb2T0 * cc2T0                &
	     	 +Term2 * aa2T0 * cc2T0                &
	     	 +Term3 * aa2T0 * bb2T0 


     no1T0 = (2**0.5) / sqrt( -B1T0 - sqrt( B1T0 ** 2 - 4 * C1T0 ) )  
     ne1T0 = (2**0.5) / sqrt( -B1T0 + sqrt( B1T0 ** 2 - 4 * C1T0 ) ) 
     ne2T0 = (2**0.5) / sqrt( -B2T0 + sqrt( B2T0 ** 2 - 4 * C2T0 ) ) 
   
	          
!**********************************************************************************************************************
!                                        Arrays Allocattion 
!**********************************************************************************************************************
!----------------------------------- Allocate Arrays Thermal
allocate(temperature(1:2,0:nr,0:nz))     
allocate(KT(0:nr,0:nz))

!----------------------------------- Allocate Arrays phase
allocate(deltaphase(0:nr,0:nz))

!**********************************************************************************************************************
!                                     Giving Zero to Arrays
!********************************************************************************************************************** 

!----------------------------------- Giving Zero to Arrys Thermal
forall (i=1:2,j=0:nr,k=0:nz)
                     temperature(i,j,k)=0.        
end forall !i

!----------------------------------

forall (j=0:nr,k=0:nz)
                                KT(j,k)=0.
end forall 								 

!----------------------------------- Giving Zero to Arrys phase
forall (j=0:nr,k=0:nz)
                        deltaphase(j,k)=(0.,0.)        
end forall !i

!**********************************************************************************************************************
!                                       Printing Constants     
!**********************************************************************************************************************

!------------------------------------------------ For Heat Equation 
write(*,*)
write(*,*)'------- Heat Equation Constants --------------------------------------------'
write(*,*)
write(*,'(A13,I9       )') '        Nt = ',Nt          
write(*,'(A13,I9       )') '        Nr = ',Nr          
write(*,'(A13,I9       )') '        Np = ',Np
write(*,'(A13,I9     ,/)') '        Nz = ',Nz          
		  
write(*,'(A13,F15.10   )') '         h = ',h           
write(*,'(A13,F15.10   )') '         E = ',E           

write(*,'(A13,F15.10   )') '        tp = ',tp
write(*,'(A13,F15.10   )') '        T0 = ',T0          
write(*,'(A13,F15.10   )') '       KT0 = ',KT0         
write(*,'(A13,F15.10   )') '        pi = ',pi          
write(*,'(A13,F15.10 ,/)') '        Cp = ',Cp          

write(*,'(A13,F15.10 ,/)') '       roh = ',roh         

write(*,'(A13,F15.10   )') '      freq = ',freq
write(*,'(A13,F15.10 ,/)') '      gama = ',gama        

write(*,'(A13,F15.10   )') '     timet = ',timet       
write(*,'(A13,F15.10 ,/)') '     sigma = ',sigma       

write(*,'(A13,F15.10   )') '    omegaf = ',omegaf      
write(*,'(A13,F15.10   )') '    length = ',length      
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

!------------------------------------------------ For Phase Equation 
write(*,*)
write(*,*)'------- Phase Equation Constants -------------------------------------------'
write(*,*)
write(*,'(A13,f15.10 ,/)') '       phi = ',phi       

write(*,'(A13,f15.10   )') '      B1T0 = ',B1T0      
write(*,'(A13,f15.10   )') '      B2T0 = ',B2T0      
write(*,'(A13,f15.10   )') '      C1T0 = ',C1T0      
write(*,'(A13,f15.10 ,/)') '      C2T0 = ',C2T0      

write(*,'(A13,f15.10   )') '     aa1T0 = ',aa1T0     
write(*,'(A13,f15.10   )') '     bb1T0 = ',bb1T0     
write(*,'(A13,f15.10   )') '     cc1T0 = ',cc1T0     
write(*,'(A13,f15.10   )') '     aa2T0 = ',aa2T0     
write(*,'(A13,f15.10   )') '     bb2T0 = ',bb2T0     
write(*,'(A13,f15.10 ,/)') '     cc2T0 = ',cc2T0     

write(*,'(A13,f15.10   )') '     nx1T0 = ',nx1T0     
write(*,'(A13,f15.10   )') '     ny1T0 = ',ny1T0     
write(*,'(A13,f15.10   )') '     nz1T0 = ',nz1T0     
write(*,'(A13,f15.10   )') '     nx2T0 = ',nx2T0     
write(*,'(A13,f15.10   )') '     ny2T0 = ',ny2T0     
write(*,'(A13,f15.10 ,/)') '     nz2T0 = ',nz2T0     

write(*,'(A13,f15.10   )') '     no1T0 = ',no1T0     
write(*,'(A13,f15.10   )') '     ne1T0 = ',ne1T0     
write(*,'(A13,f15.10 ,/)') '     ne2T0 = ',ne2T0     

write(*,'(A13,f15.10   )') '    dnx1dT = ',dnx1dT    
write(*,'(A13,f15.10   )') '    dny1dT = ',dny1dT    
write(*,'(A13,f15.10   )') '    dnz1dT = ',dnz1dT    
write(*,'(A13,f15.10   )') '    dnx2dT = ',dnx2dT    
write(*,'(A13,f15.10   )') '    dny2dT = ',dny2dT    
write(*,'(A13,f15.10 ,/)') '    dnz2dT = ',dnz2dT    

write(*,'(A13,f15.10   )') '     Term1 = ',Term1     
write(*,'(A13,f15.10   )') '     Term2 = ',Term2     
write(*,'(A13,f15.10 ,/)') '     Term3 = ',Term3     

write(*,'(A13,f15.10 ,/)') '     theta = ',theta     

write(*,'(A13,f15.10   )') '   lambda1 = ',lambda1   
write(*,'(A13,f15.10 ,/)') '   lambda2 = ',lambda2   
       

write(*,*)'----------------------------------------------------------------------------'
write(*,'(A,\)')' Please press any key to continue '
read(*,*)

!**********************************************************************************************************************
!                                   Main Block of the Program     
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
 
 p = E / (tp * sqrt(pi))            !total power of the pulse
!--------------------------------------------------------
 
do  j=0,nr
    do k=0,nz
	  
       temperature(1,j,k) = T0
	   KT(j,k) = KT0
	   
    end do !k
end do !j	     


do  l=1,Np     !Runing program for Np pulses 
         
!--------------------------------------------- Run program for one pulse 
    do  i=0,nt
	    t=deltat*i
	  
	    do  j=1,nr-1
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
			   temperature(1,0 ,k)  = temperature(1,1,k)            !Thermal insulation condition for crystal axis

			   temperature(1,nr,k)  = T0                            !Temperature-fixed condition for lateral surface

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
			   temperature(1,nr,0 ) = T0                            !Temperature-fixed condition for (nr,0 )

			   temperature(1,nr,nz) = T0                            !Temperature-fixed condition for (nr,nz)
			   !------------------------------------ Ending of Boundary conditions

			   !------------------------------------ Heat Equation
			   temperature(2,j,k) = temperature(1,j,k)                                                                               &
							
						          + aa3 * ( (temperature(1,j+1,k) -  temperature(1,j-1,k))/(2*r*deltar)                              &
								
						          + (temperature(1,j+1,k) - 2*temperature(1,j,k) + temperature(1,j-1,k))/(deltar**2) )               & 

						          + aa3 * ( (temperature(1,j,k-1) - 2*temperature(1,j,k) + temperature(1,j,k+1))/(deltaz**2) )       &                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

						          + aa4 * exp( (-2*r**2)/(omegaf**2) ) * exp(-gama*z) * exp(-( (t-2*tp)/tp )**2 )                    &

						          + aa5 * ( ( kT(j+1,k)-kT(j-1,k) ) * ( temperature(1,j+1,k)-temperature(1,j-1,k) ) / (deltar**2)    &       

						          + ( kT(j,k+1)-kT(j,k-1) ) * ( temperature(1,j,k+1)-temperature(1,j,k-1) ) / (deltaz**2) )
										

			   !----------------------------------- Phase Equation constants
			   nx1r0T = nx1T0 + dnx1dT * ( temperature(1,0,k) - t0 )
			   ny1r0T = ny1T0 + dny1dT * ( temperature(1,0,k) - t0 )
			   nz1r0T = nz1T0 + dnz1dT * ( temperature(1,0,k) - t0 )
			
			   nx2r0T = nx2T0 + dnx2dT * ( temperature(1,0,k) - t0 )
			   ny2r0T = ny2T0 + dny2dT * ( temperature(1,0,k) - t0 )
			   nz2r0T = nz2T0 + dnz2dT * ( temperature(1,0,k) - t0 )
				
			   aa1r0T = 1 / ( nx1r0T )**2 
			   bb1r0T = 1 / ( ny1r0T )**2 
			   cc1r0T = 1 / ( nz1r0T )**2 

			   aa2r0T = 1 / ( nx2r0T )**2 
			   bb2r0T = 1 / ( ny2r0T )**2 
			   cc2r0T = 1 / ( nz2r0T )**2       

			   B1r0T = -Term1 * ( bb1r0T + cc1r0T )            &
			           -Term2 * ( aa1r0T + cc1r0T )            &
			           -Term3 * ( aa1r0T + bb1r0T ) 
			
			   C1r0T =  Term1 * bb1r0T * cc1r0T                &
			           +Term2 * aa1r0T * cc1r0T                &
			           +Term3 * aa1r0T * bb1r0T 

			   B2r0T = -Term1 * ( bb2r0T + cc2r0T )            &
			           -Term2 * ( aa2r0T + cc2r0T )            &
			           -Term3 * ( aa2r0T + bb2r0T )
			
			   C2r0T = Term1 * bb2r0T * cc2r0T                 &
			          +Term2 * aa2r0T * cc2r0T                 &
			          +Term3 * aa2r0T * bb2r0T 

			
			   no1r0T = (2**0.5) / sqrt( -B1r0T  - sqrt( B1r0T**2 - 4*C1r0T ) )  
			   ne1r0T = (2**0.5) / sqrt( -B1r0T  + sqrt( B1r0T**2 - 4*C1r0T ) )
			   ne2r0T = (2**0.5) / sqrt( -B2r0T  + sqrt( B2r0T**2 - 4*C2r0T ) ) 
			
			   deltano1r0T = no1r0T - no1T0
			   deltane1r0T = ne1r0T - ne1T0
			   deltane2r0T = ne2r0T - ne2T0

			   !----------------------------------- Phase Equation constants
			   nx1rT = nx1T0 + dnx1dT * ( temperature(1,j,k) - t0 )
			   ny1rT = ny1T0 + dny1dT * ( temperature(1,j,k) - t0 )
			   nz1rT = nz1T0 + dnz1dT * ( temperature(1,j,k) - t0 )
			
			   nx2rT = nx2T0 + dnx2dT * ( temperature(1,j,k) - t0 )
			   ny2rT = ny2T0 + dny2dT * ( temperature(1,j,k) - t0 )
			   nz2rT = nz2T0 + dnz2dT * ( temperature(1,j,k) - t0 )
				
			   aa1rT = 1 / ( nx1rT )**2 
			   bb1rT = 1 / ( ny1rT )**2 
			   cc1rT = 1 / ( nz1rT )**2 

			   aa2rT = 1 / ( nx2rT )**2 
			   bb2rT = 1 / ( ny2rT )**2 
			   cc2rT = 1 / ( nz2rT )**2       

			   B1rT = -Term1 * ( bb1rT + cc1rT )                &
			          -Term2 * ( aa1rT + cc1rT )                &
			          -Term3 * ( aa1rT + bb1rT ) 
			
			   C1rT =  Term1 * bb1rT * cc1rT                    &
			          +Term2 * aa1rT * cc1rT                    &
			          +Term3 * aa1rT * bb1rT 

			   B2rT = -Term1 * ( bb2rT + cc2rT )                &
			          -Term2 * ( aa2rT + cc2rT )                &
			          -Term3 * ( aa2rT + bb2rT )
			
			   C2rT =  Term1 * bb2rT * cc2rT                    &
			          +Term2 * aa2rT * cc2rT                    &
			          +Term3 * aa2rT * bb2rT 

			
			   no1rT = (2**0.5) / sqrt( -B1rT  - sqrt( B1rT**2 - 4*C1rT ) )  
			   ne1rT = (2**0.5) / sqrt( -B1rT  + sqrt( B1rT**2 - 4*C1rT ) )
			   ne2rT = (2**0.5) / sqrt( -B2rT  + sqrt( B2rT**2 - 4*C2rT ) ) 
			
			   deltano1rT = no1rT - no1T0
			   deltane1rT = ne1rT - ne1T0
			   deltane2rT = ne2rT - ne2T0

			   !------------------------------------ For Phase Equation
			   deltaphase(j ,0  ) = (0.,0.)                                      	          !for input surface

			   deltaphase(nr,k  ) = (0.,0.)                                       	          !for lateral surface

			   deltaphase(j,nz  ) = deltaphase(j,nz-1)                           	   &      !for output surface
				                  + ( 2*pi*deltaz / lambda1 )                    	   &
				                  * ( deltano1rT + deltane1rT - 2*deltane2rT )                     

			   deltaphase(nr,0  ) = (0.,0.)                                        	          !for (nr,0 )

			   deltaphase(nr,nz ) = (0.,0.)                                     	          !for (nr,nz)

			   deltaphase(0 ,0  ) = (0.,0.)                                       	          !for ( 0,0 )

			   !------
			   deltaphase(0 ,k  ) = deltaphase(0,k-1)                                  &      !for crystal axis
				                  + ( 2*pi*deltaz / lambda1 )                          &
				   	              * ( deltano1r0T  + deltane1r0T  - 2*deltane2r0T  )                                            


			   deltaphase(0 ,nz ) = deltaphase(0,nz-1)                                 &      !for ( 0,nz)
							      + ( 2*pi*deltaz / lambda1 )                          &
					              * ( deltano1rT + deltane1rT - 2*deltane2rT )                     


			   !----------------------------------- Phase Equation
																				
			   deltaphase(j,k) = deltaphase(j,k-1)                                     &
				               + ( 2*pi*deltaz / lambda1 )                             &  
				               * ( deltano1rT + deltane1rT  - 2*deltane2rT  )                                 

			   !-----------------------------------
			end do !k
		end do !j


		!--------------------------------------------- End of run for each deltat 

		!============================================= Print Results for each deltat
		t=(l-1)*nt*deltat + i*deltat  
		!--------------------------------------------- For Heat Equation
		write(1,'(2x,f25.10,5x,f25.10)')  t , temperature(1,0,0)

		!--------------------------------------------- For Phase Equation
		write(4,'(2x,f25.10,5x,2f25.10)') t , deltaphase(1,1) 

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

			KT(j,k) = KT0  * T0 / temperature(1,j,k)

			end do !k
		end do !j	     

		!---------------------------------------------

	end do !i

end do !l

!**********************************************************************************************************************
!                                        Printing Results     
!**********************************************************************************************************************

!------------------------------------------------ For Heat Equation
do j=0,nr
   r=j*deltar 
   write(2,'(2x,f25.10,5x,f25.10)')  r , temperature(1,j,1)
end do !j      						   

!------------------------------------------------
do k=0,nz
   z=k*deltaz 
   write(3,'(2x,f25.10,5x,f25.10)')  z , temperature(1,0,k)
end do !k      						   

!------------------------------------------------ For Phase Equation
do j=0,nr
   r=j*deltar 
   write(5,'(2x,f25.10,5x,2f25.10)') r , deltaphase(j,1)
end do !j      						   

!------------------------------------------------
do k=0,nz
   z=k*deltaz 
   write(6,'(2x,f25.10,5x,2f25.10)') z , deltaphase(0,k)
end do !k      						   

!**********************************************************************************************************************
!                                      Closing Files and Ending the Program 
!**********************************************************************************************************************

!------------------------------------------------ For Heat Equation
close(1)
close(2)
close(3)
!------------------------------------------------ For Phase Equation
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

end program Temp_Phase_pw                     

!======================================================================================================================
         
 
