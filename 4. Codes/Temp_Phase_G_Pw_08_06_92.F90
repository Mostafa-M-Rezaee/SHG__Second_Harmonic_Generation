
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

!            * "Temp_Phase_PW.F90" is a program to solve Heat-equation & Phase by FDM.      * 

!            *                                                                              * 

!            *             originally Written : 30-DES-2011   by MM  &  AM  &  FS           *

!            *                   last revised : 31-Aug-2012   by MM  &  AM  &  FS           *

!            ********************************************************************************

program Temp_Phase_PW

implicit none

!**********************************************************************************************************************
!                                       Variables Definition
!**********************************************************************************************************************

!------------------------------------------------ Thermal Variables
integer       i          ,j          ,k          ,l                                                                   &
             ,nt         ,nr         ,nz         ,Np                                                                       


real*8        t          ,z          ,E          ,h        ,r            ,G          ,P                               &
             ,T0         ,pi         ,Cp         ,tp         ,Q0                                                      &
             ,roh        ,aa1        ,aa2        ,aa3        ,aa4        ,aa5                                         &
             ,KT0        ,freq       ,gama       ,Tinf       ,Tamb                                                    & 
             ,timet      ,sigma                                                                                       &
	     ,omegaf     ,length     ,deltar     ,deltaz     ,deltat                                                  &
	     ,radius                                                                                                  &
	     ,epsilong   ,tbetween                                                                                    &
	     ,stability                                                                                               &
			  
	    ,temperature[allocatable](:,:,:)          ,KT[allocatable](:,:)  


character*30  filenameTt   ,filenameTr   ,filenameTz     ,Npf       ,freqf     ,tpf       ,EE

!------------------------------------------------ Phase Variables
real*8        phi                    ,aa1T0       ,aa2T0       ,aa1r0T       ,aa2r0T                        &
             ,B1T0       ,B2T0       ,bb1T0       ,bb2T0       ,bb1r0T       ,bb2r0T                        &
	     ,C1T0       ,C2T0       ,cc1T0       ,cc2T0       ,cc1r0T       ,cc2r0T                        &                                

	     ,B1rT       ,B2rT       ,B1r0T       ,B2r0T       ,nx1r0T       ,nx2r0T                        &
             ,C1rT       ,C2rT       ,C1r0T       ,C2r0T       ,ny1r0T       ,ny2r0T                        &
	     ,theta      ,nz1r0T     ,nz2r0T                                                                &

             ,nx1T0      ,nx2T0      ,no1T0       ,nx1rT       ,dnx1dT       ,dnx2dT                        &
	     ,ny1T0      ,ny2T0      ,ne1T0       ,ny1rT       ,dny1dT       ,dny2dT                        &
	     ,nz1T0      ,nz2T0      ,ne2T0       ,nz1rT       ,dnz1dT       ,dnz2dT                        &

	     ,Term1      ,aa1rT      ,aa2rT       ,no1rT       ,no1r0T                                      &
	     ,Term2      ,bb1rT      ,bb2rT       ,ne1rT       ,ne1r0T       ,lambda1                       &
             ,Term3      ,cc1rT      ,cc2rT       ,ne2rT       ,ne2r0T       ,lambda2                       &                           
             
	     ,nx2rT      ,deltano1r0T             ,deltano1rT                                               &
             ,ny2rT      ,deltane1r0T             ,deltane1rT                                               &
	     ,nz2rT      ,deltane2r0T             ,deltane2rT                                 
 
 complex*8    deltaphase[allocatable](:,:)

character*30  filenamePt   ,filenamePr  ,filenamePz                                                                   

!**********************************************************************************************************************
!                                    Giving Zero to variables
!**********************************************************************************************************************

!------------------------------------------------ Giving Zero to Thermal Variables
 
                    i = 0           ;j = 0           ;k = 0          ;l = 0            
                   nt = 0          ;nr = 0          ;nz = 0         ;Np = 0

                    t = 0.          ;z = 0.          ;E = 0.         ;h = 0.         ;r = 0.         ;G = 0.        ;P = 0.                                                                                                               
                   T0 = 0.         ;pi = 0.         ;Cp = 0.        ;tp = 0.        ;Q0 = 0.                                     
	          roh = 0.        ;aa1 = 0.        ;aa2 = 0.       ;aa3 = 0.       ;aa4 = 0.       ;aa5 = 0.                                                                        
                 Tinf = 0.       ;Tamb = 0.       ;freq = 0.      ;gama = 0.
	          KT0 = 0.      ;timet = 0.      ;sigma = 0.                                                                       
	       omegaf = 0.     ;length = 0.     ;deltar = 0.    ;deltaz = 0.    ;deltat = 0.                                           
	       radius = 0.
 	     epsilong = 0.   ;tbetween = 0.                                   
	    stability = 0.   

!------------------------------------------------ Giving Zero to Phase Variables
             phi = 0.                        ;aa1T0 = 0.       ;aa2T0 = 0.       ;aa1r0T = 0.       ;aa2r0T = 0.                        
            B1T0 = 0.       ;B2T0 = 0.       ;bb1T0 = 0.       ;bb2T0 = 0.       ;bb1r0T = 0.       ;bb2r0T = 0.                        
	    C1T0 = 0.       ;C2T0 = 0.       ;cc1T0 = 0.       ;cc2T0 = 0.       ;cc1r0T = 0.       ;cc2r0T = 0.                                                        

	    B1rT = 0.       ;B2rT = 0.       ;B1r0T = 0.       ;B2r0T = 0.       ;nx1r0T = 0.       ;nx2r0T = 0.                        
            C1rT = 0.       ;C2rT = 0.       ;C1r0T = 0.       ;C2r0T = 0.       ;ny1r0T = 0.       ;ny2r0T = 0.                        
	  ;theta = 0.     ;nz1r0T = 0.      ;nz2r0T = 0.                        

           nx1T0 = 0.      ;nx2T0 = 0.       ;no1T0 = 0.       ;nx1rT = 0.       ;dnx1dT = 0.       ;dnx2dT = 0.                        
	   ny1T0 = 0.      ;ny2T0 = 0.       ;ne1T0 = 0.       ;ny1rT = 0.       ;dny1dT = 0.       ;dny2dT = 0.                        
	   nz1T0 = 0.      ;nz2T0 = 0.       ;ne2T0 = 0.       ;nz1rT = 0.       ;dnz1dT = 0.       ;dnz2dT = 0.                        

	   Term1 = 0.      ;aa1rT = 0.       ;aa2rT = 0.       ;no1rT = 0.       ;no1r0T = 0.                                      
	   Term2 = 0.      ;bb1rT = 0.       ;bb2rT = 0.       ;ne1rT = 0.       ;ne1r0T = 0.       ;lambda1 = 0.                       
           Term3 = 0.      ;cc1rT = 0.       ;cc2rT = 0.       ;ne2rT = 0.       ;ne2r0T = 0.       ;lambda2 = 0.                                                 
             
           nx2rT = 0.      ;deltano1r0T = 0.              ;deltano1rT = 0.                                               
           ny2rT = 0.      ;deltane1r0T = 0.              ;deltane1rT = 0.                                               
	   nz2rT = 0.      ;deltane2r0T = 0.              ;deltane2rT = 0.                                 

!**********************************************************************************************************************
!                                             Inputs		  
!**********************************************************************************************************************

write(*,'(/,2x,a,\)') '            Enter the Energy value  : '
!read(*,*) E
E = 0.09          !????????????????????????????????????????????????????  
write(*,'(/,2x,a,\)') '                             Again  : '
!read(*,*) EE
EE = '0.09'            

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

!------------------------------------------------ Heat Equation Files
filenameTt = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Tt.plt'
open(1,file=filenameTt)
write(1,'(/,a,/)')   ! ' variables=         "t"                             "temperature"'

filenameTr = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Tr.plt'
open(2,file=filenameTr)
write(2,'(/,a,/)')    !' variables=         "r"                             "temperature"'

filenameTz = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Tz.plt'
open(3,file=filenameTz)
write(3,'(/,a,/)')    !' variables=         "z"                             "temperature"' 

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenameTt,filenameTr,filenameTz
 read(*,*)

!------------------------------------------------ Phase Equation Files
filenamePt = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Pt.plt'
open(4,file=filenamePt)
write(4,'(/,a,/)')    !' variables=         "t"          "deltaphase_real"      "deltaphase_imaginary"'

filenamePr = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Pr.plt'
open(5,file=filenamePr)
write(5,'(/,a,/)')    !' variables=         "r"          "deltaphase_real"      "deltaphase_imaginary"'

filenamePz = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Pz.plt'
open(6,file=filenamePz)
write(6,'(/,a,/)')    !' variables=         "z"           "deltaphase_real"      "deltaphase_imaginary"'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenamePt,filenamePr,filenamePz
 read(*,*)

!**********************************************************************************************************************
!                                           Constants
!**********************************************************************************************************************

!------------------------------------------------ Thermal properties
        h = 10                  !heat transfer coefficient (convection - cylinder)        W/(m^2.K)
       T0 = 300.                !initial temperature                                      K
       pi = 4*atan(1.)                                                                   !dimensionless
       Cp = 728.016             !heat capacity at constant pressure                       J/(kg.K)
    !  Cp = 590.             !heat capacity at constant pressure                       J/(kg.K)

!      KT0 = 2.75                !thermal conductivity of KTP crystal                      W/(m.K)

      KT0 = 13.                !thermal conductivity of KTP crystal                      W/(m.K)
      roh = 2945.               !mass density                                             kg/m^3
!     roh = 560.               !mass density                                             kg/m^3
     gama = 4.                 !absorption coefficient                                   1/m
     Tinf = 300.
     Tamb = 300.
    sigma = 5.669e-8            !Stephan-Bultzman constant                                W/(m^2.K^4) 

 tbetween = 1/freq !1/(5*freq)                                                                       !s
    timet = Np*tbetween                                                                  !s
   deltat = tp / 10                                                                      !s     
       nt = int(tbetween/deltat)                                                         !dimensionless

       nz = 150 !200                                                                          !dimensionless
   length = 0.02                !length of crystal                                        m 
   deltaz = length/nz                                                                    !m

   radius = 0.005               !radius of crystal                                       !m
   omegaf = 100.e-6               !spot size                                               !m
   deltar = omegaf/10                                                                    !m
       nr = int(radius/deltar)                                                           !dimensionless 

 epsilong = 0.9                 !surface emissivity                                      !dimensionless

stability = ( (2*KT0*deltat)/(roh*Cp) ) * ( (deltar**2+deltaz**2)/(deltar**2*deltaz**2) ) !stability coefficient  

!------------------------------------------------ Phase properties
    !  phi = 0.4324341714      !????????????????????????????????????????????????????????????????????????????????????????
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

      B1T0 = -Term1 * ( bb1T0 + cc1T0 )                                   &
	     -Term2 * ( aa1T0 + cc1T0 )                                   &
	     -Term3 * ( aa1T0 + bb1T0 ) 
     		 
      C1T0 =  Term1 * bb1T0 * cc1T0                                       &
	     +Term2 * aa1T0 * cc1T0                                       &
	     +Term3 * aa1T0 * bb1T0 


      B2T0 = -Term1 * ( bb2T0 + cc2T0 )                                   &
	     -Term2 * ( aa2T0 + cc2T0 )                                   &
	     -Term3 * ( aa2T0 + bb2T0 )
             
      C2T0 =  Term1 * bb2T0 * cc2T0                                       &
	     +Term2 * aa2T0 * cc2T0                                       &
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
write(*,'(A13,I5    ,/,      &
          A13,I5    ,/,      &
		  A13,I5    ,//,     &
		  
		  A13,F15.10,/,      &
		  A13,F15.10,/,      &

		  A13,F15.10,/,      &
		  A13,F15.10,/,      &
		  A13,F15.10,/,      &
		  A13,F15.10,//,     &

		  A13,F15.10,//,     &

		  A13,F15.10,//,     &

		  A13,F15.10,/,      &
		  A13,F15.10,//,     &

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
'       KT0 = ',KT0        ,     &
'        pi = ',pi         ,     &
'        Cp = ',Cp         ,     &

'       roh = ',roh        ,     &

'      gama = ',gama       ,     &

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
read(*,*)

!------------------------------------------------ For Phase Equation 
write(*,*)
write(*,*)'------- Phase Equation Constants -------------------------------------------'
write(*,*)
write(*,'(A13,f15.10,//,     &

		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,//,     &

		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,//,     &

		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,//,     &

		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,//,     &

		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,//,     &

		  A13,f15.10,/,      &
		  A13,f15.10,/,      &
		  A13,f15.10,//,     &

		  A13,f15.10,//,     &

		  A13,f15.10,/,      &
		  A13,f15.10,// )')  &

'       phi = ',phi        ,     &

'      B1T0 = ',B1T0       ,     &
'      B2T0 = ',B2T0       ,     &
'      C1T0 = ',C1T0       ,     &
'      C2T0 = ',C2T0       ,     &

'     aa1T0 = ',aa1T0      ,     &
'     bb1T0 = ',bb1T0      ,     &
'     cc1T0 = ',cc1T0      ,     &
'     aa2T0 = ',aa2T0      ,     &
'     bb2T0 = ',bb2T0      ,     &
'     cc2T0 = ',cc2T0      ,     &

'     nx1T0 = ',nx1T0      ,     &
'     ny1T0 = ',ny1T0      ,     &
'     nz1T0 = ',nz1T0      ,     &
'     nx2T0 = ',nx2T0      ,     &
'     ny2T0 = ',ny2T0      ,     &
'     nz2T0 = ',nz2T0      ,     &

'     no1T0 = ',no1T0      ,     &
'     ne1T0 = ',ne1T0      ,     &
'     ne2T0 = ',ne2T0      ,     &

'    dnx1dT = ',dnx1dT     ,     &
'    dny1dT = ',dny1dT     ,     &
'    dnz1dT = ',dnz1dT     ,     &
'    dnx2dT = ',dnx2dT     ,     &
'    dny2dT = ',dny2dT     ,     &
'    dnz2dT = ',dnz2dT     ,     &

'     Term1 = ',Term1      ,     &
'     Term2 = ',Term2      ,     &
'     Term3 = ',Term3      ,     &

'     theta = ',theta      ,     &

'   lambda1 = ',lambda1    ,     &
'   lambda2 = ',lambda2        

write(*,*)'----------------------------------------------------------------------------'
write(*,'(A,\)')' Press Enter to continue '
read(*,*)

!**********************************************************************************************************************
!                                   Main Block of the Program     
!**********************************************************************************************************************
!-------------------------------------- control condition
 if ( (omegaf .NE. 100.e-6) .OR. (length .NE. 2e-2) .OR. (radius .NE. 5e-3) .OR. (gama .NE. 4.) ) then
      
	write(*,'(a/)')'Error --- please Enter the new variables in maple program and then obtain normalization(G)'
  
end if

!-------- 
 G = 1.703409099e-10 * pi * tp                         

Q0 = sqrt(pi) * tp / G              !Normalization      m^-3

 p = E / (tp * sqrt(pi))            !total power of the pulse
!--------------------------------------------------------
 
do j=0,nr
   do k=0,nz
	  
      temperature(1,j,k) = T0
	         KT(j,k) = KT0
	   
   end do !k
end do !j	     


do l=1,Np !Runing program for Np pulses 
         
!--------------------------------------------- Run program for one pulse 
   do i=0,nt
      t=deltat*i
      
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

	     B1r0T = -Term1 * ( bb1r0T + cc1r0T )                                   &
		     -Term2 * ( aa1r0T + cc1r0T )                                   &
		     -Term3 * ( aa1r0T + bb1r0T ) 
     		 
	     C1r0T =  Term1 * bb1r0T * cc1r0T                                       &
		     +Term2 * aa1r0T * cc1r0T                                       &
		     +Term3 * aa1r0T * bb1r0T 

             B2r0T = -Term1 * ( bb2r0T + cc2r0T )                                   &
		     -Term2 * ( aa2r0T + cc2r0T )                                   &
		     -Term3 * ( aa2r0T + bb2r0T )
             
	     C2r0T = Term1 * bb2r0T * cc2r0T                                       &
		    +Term2 * aa2r0T * cc2r0T                                       &
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

	     B1rT = -Term1 * ( bb1rT + cc1rT )                                   &
		    -Term2 * ( aa1rT + cc1rT )                                   &
		    -Term3 * ( aa1rT + bb1rT ) 
     		 
	     C1rT =  Term1 * bb1rT * cc1rT                                       &
		    +Term2 * aa1rT * cc1rT                                       &
		    +Term3 * aa1rT * bb1rT 

             B2rT = -Term1 * ( bb2rT + cc2rT )                                   &
		    -Term2 * ( aa2rT + cc2rT )                                   &
		    -Term3 * ( aa2rT + bb2rT )
             
	     C2rT =  Term1 * bb2rT * cc2rT                                       &
		    +Term2 * aa2rT * cc2rT                                       &
		    +Term3 * aa2rT * bb2rT 

             
 	    no1rT = (2**0.5) / sqrt( -B1rT  - sqrt( B1rT**2 - 4*C1rT ) )  
            ne1rT = (2**0.5) / sqrt( -B1rT  + sqrt( B1rT**2 - 4*C1rT ) )
            ne2rT = (2**0.5) / sqrt( -B2rT  + sqrt( B2rT**2 - 4*C2rT ) ) 
			 
       deltano1rT = no1rT - no1T0
       deltane1rT = ne1rT - ne1T0
      deltane2rT = ne2rT - ne2T0
			
            !------------------------------------ For Phase Equation
            deltaphase(j ,0  ) = (0.,0.)                                                                            !for input surface
            
            deltaphase(nr,k  ) = (0.,0.)                                                                            !for lateral surface

	    deltaphase(j,nz  ) = deltaphase(j,nz-1)                                                   &           !for output surface
				+ ( 2*pi*deltaz / lambda1 )                                           &
				* ( deltano1rT + deltane1rT - 2*deltane2rT )                     

            deltaphase(nr,0  ) = (0.,0.)                                                                            !for (nr,0 )

	    deltaphase(nr,nz ) = (0.,0.)                                                                            !for (nr,nz)

	    deltaphase(0 ,0  ) = (0.,0.)                                                                            !for ( 0,0 )

	    !------
	    deltaphase(0 ,k  ) = deltaphase(0,k-1)                                                    &           !for crystal axis
				+ ( 2*pi*deltaz / lambda1 )                                           &
			        * ( deltano1r0T  + deltane1r0T  - 2*deltane2r0T  )                                            


            deltaphase(0 ,nz ) = deltaphase(0,nz-1)                                                   &           !for ( 0,nz)
	                        + ( 2*pi*deltaz / lambda1 )                                           &
			        * ( deltano1rT + deltane1rT - 2*deltane2rT )                     

         
             !----------------------------------- Phase Equation
			                                                                       
	     deltaphase(j,k) = deltaphase(j,k-1)                                                       &
			      + ( 2*pi*deltaz / lambda1 )                                              &  
			      * ( deltano1rT + deltane1rT  - 2*deltane2rT  )                                 

	    !-----------------------------------
         end do !k
      end do !j


   !--------------------------------------------- End of run for each deltat 

   !============================================= Print Results for each deltat
   
   !--------------------------------------------- For Heat Equation
   t=(l-1)*nt*deltat + i*deltat 
   write(1,'(2x,f25.10,5x,f25.10)')  t , temperature(1,0,0)

   !--------------------------------------------- For Phase Equation
   t=(l-1)*nt*deltat + i*deltat 
   write(4,'(2x,f25.10,5x,2f25.10)') t , deltaphase(1,1) 
   
   !=============================================
    
	 !write(*,*)'j,deltaphase(1,j,2)',j,deltaphase(1,j,2) 
     !read(*,*)

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

end program Temp_Phase_pw                     

!======================================================================================================================
         
 
