

!            ********************************************************************************
!            * File: 6. Ideal _ Pulsed Wave Bessel Gaussian _ Computational.F90             *
!            *                                                                              *
!            * Note: This Fortran code is developed specifically for the article titled:    *
!            * Pulsed Bessel-Gauss beams: A depleted wave model for type II second harmonic *
!            * generation                                                                   *
!            *                                                                              *
!            * Authors: Sabaeian, M., Motazedian, A., Rezaee, M.M. and Jalil-Abadi, F.S.    *
!            *                                                                              *
!            * Harvard style:                                                               *
!            * Sabaeian, M., Motazedian, A., Rezaee, M.M. and Jalil-Abadi, F.S., 2014.      *
!            * Pulsed Bessel–Gauss beams: a depleted wave model for type II second-harmonic *
!            * generation. Applied Optics, 53(32), pp.7691-7696.                            *
!            *                                                                              *
!            ********************************************************************************

program Elec_BG_PW

implicit none

!**********************************************************************************************************************
!                                       Variables Definition
!**********************************************************************************************************************
!-------------------------------------- Common Variables
integer       i            ,j          ,k          ,l          ,m                                                    &
             ,nt           ,nr         ,nz         ,Np         ,inn         ,kn                                      &
			 ,nt1          ,run                                                                                      &
	         ,nomegaf                                                          

real*8        E            ,t          ,z          ,r          ,x          ,y                                        &                                 
             ,pi           ,tp         ,Cp         ,s1         ,s2         ,y1         ,y2                           &
             ,roh          ,KT0                                                                                      &
             ,freq         ,beta                                                                                     &
			 ,sigma        ,timet      ,alpha      ,gama1      ,gama2      ,gama3      ,power      ,nnrom            &   
			 ,no1T0        ,ne1T0      ,ne2T0      ,no1rT      ,ne2rT      ,ne1rT                                    &
			 ,omegaf       ,length     ,deltar     ,deltaz     ,deltat     ,radius                                   &  
			 ,lambda1      ,lambda2    ,deltar1    ,deltar2                                                          & 
			 ,tbetween     ,Fidegree   ,Firadian   ,J0alphar                                                         &
			 ,stability                                                                                              &
             ,teta_bessel                                                                                            &
	         ,dteta_bessel                                                                                          


complex*16    Ii                                                                                                       

character*35  freqf      ,Npf        ,tpf        ,EE

!-------------------------------------- Fields Variables
integer       f            ,ibest

real*8        c                                                                                                      &
             ,fi                                                                                                     &
			 ,deff                                                                                                   &
             ,omega        ,Psi22        ,Psi32                                                                      &
             ,Lscale       ,Elec12       ,Elec22        ,Elec32                                                      &
			 ,epsilon0     ,Psi2max      ,Psi3max                                                                    &                       
             ,Lscalemax

complex*16    cc1          ,cc2          ,cc3           ,cc4          ,cc5                                           &
             ,dd1          ,dd2          ,dd3           ,dd4          ,dd5                                           &
             ,ee1          ,ee2          ,ee3           ,ee4          ,ee5                                           &

             ,Psi1[allocatable](:,:,:)   ,Elec1[allocatable](:,:,:)                                                  &
			 ,Psi2[allocatable](:,:,:)   ,Elec2[allocatable](:,:,:)                                                  &  
			 ,Psi3[allocatable](:,:,:)   ,Elec3[allocatable](:,:,:) 
                     

character*35  filenameibestl                                                                                         &
             ,filenameElec12t      ,filenameElec12r       ,filenameElec12z                                           &
             ,filenameElec22t      ,filenameElec22r       ,filenameElec22z                                           &         
			 ,filenameElec32t      ,filenameElec32r       ,filenameElec32z                                           &

			 ,filenamePsi3picksl   ,filenamePsi2picksl    ,filenameLscalemaxl

!**********************************************************************************************************************
!                                       Giving Zero to variables
!**********************************************************************************************************************
!-------------------------------------- Giving Zero to Common Variables
              i = 0             ;j = 0             ;k = 0         ;l = 0         ;m = 0          
             nt = 0            ;nr = 0            ;nz = 0        ;Np = 0       ;inn = 0         ;kn = 0.
			nt1 = 0           ;run = 0                                                 
        nomegaf = 0

              E = 0.            ;t = 0.            ;z = 0.        ;r = 0.        ;x = 0.         ;y = 0.         
             pi = 0.           ;tp = 0.           ;Cp = 0.       ;s1 = 0.       ;s2 = 0.        ;y1 = 0.       ;y2 = 0.           
            roh = 0.          ;KT0 = 0. 
           freq = 0.         ;beta = 0.                                                                                     
		  sigma = 0.        ;timet = 0.        ;alpha = 0.    ;gama1 = 0.    ;gama2 = 0.     ;gama3 = 0.    ;power = 0.   
          nnrom = 0.        ;no1T0 = 0.        ;ne1T0 = 0.    ;ne2T0 = 0.    ;no1rT = 0.     ;ne2rT = 0.    ;ne1rT = 0.
		 omegaf = 0.       ;length = 0.       ;deltar = 0.   ;deltaz = 0.   ;deltat = 0.    ;radius = 0.                                               
		lambda1 = 0.      ;lambda2 = 0.      ;deltar1 = 0.  ;deltar2 = 0.
       tbetween = 0.     ;Fidegree = 0.     ;Firadian = 0. ;J0alphar = 0.   
      stability = 0.
   dteta_bessel = 0. 

             Ii = (0.,0.)  

!------------------------------------------------ Giving Zero to Fields Variables
              f = 0.        ;ibest = 0.
              c = 0.                                                    
             fi = 0.         

            cc1 = 0.          ;cc2 = 0.          ;cc3 = 0.      ;cc4 = 0.      ;cc5 = 0.                                           
            dd1 = 0.          ;dd2 = 0.          ;dd3 = 0.      ;dd4 = 0.      ;dd5 = 0.                                          
            ee1 = 0.          ;ee2 = 0.          ;ee3 = 0.      ;ee4 = 0.      ;ee5 = 0.                                          
 
           deff = 0.                                        
		  omega = 0.        ;Psi22 = 0.        ;Psi32 = 0.
		 Lscale = 0.       ;Elec12 = 0.       ;Elec22 = 0.   ;Elec32 = 0.  
       epsilon0 = 0.      ;Psi2max = 0.      ;Psi3max = 0.
      Lscalemax = 0.

!**********************************************************************************************************************
!                                             Inputs		  
!**********************************************************************************************************************
write(*,'(/,2x,a,\)') '            Enter the Energy value  : '
!read(*,*) E
E = 0.8        
write(*,'(/,2x,a,\)') '                             Again  : '
!read(*,*) EE
EE = '0.08'            

write(*,'(/,2x,a,\)') '         Enter the frequency value  : '
!read(*,*) freq
freq = 4000
write(*,'(/,2x,a,\)') '                             Again  : '
!read(*,*) freqf
freqf = '4000'

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

!------------------------------------------------ Field Equations Files
filenameElec12t = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Elec12t.plt'
open(7,file=filenameElec12t)
!write(7,'(/,a,/)')    ' variables =     "t"                              "Elec1 ** 2"'

filenameElec12r = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Elec12r.plt'
open(8,file=filenameElec12r)
!write(8,'(/,a,/)')    ' variables =     "r"                              "Elec1 ** 2"'

filenameElec12z = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Elec12z.plt'
open(9,file=filenameElec12z)
!write(9,'(/,a,/)')    ' variables =     "z"                              "Elec1 ** 2"'


write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenameElec12t  &
                                                                                    ,filenameElec12r  &
                                                                                    ,filenameElec12z
 read(*,*)

!------------------
filenameElec22t = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Elec22t.plt'
open(10,file=filenameElec22t)
!write(10,'(/,a,/)')    ' variables =     "t"                              "Elec2 ** 2"'

filenameElec22r = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Elec22r.plt'
open(11,file=filenameElec22r)
!write(11,'(/,a,/)')    ' variables =     "r"                              "Elec2 ** 2"'

filenameElec22z = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Elec22z.plt'
open(12,file=filenameElec22z)
!write(12,'(/,a,/)')    ' variables =     "z"                              "Elec2 ** 2"'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenameElec22t  &
                                                                                    ,filenameElec22r  &
										    ,filenameElec22z
 read(*,*)

!------------------
filenameElec32t = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Elec32t.plt'
open(13,file=filenameElec32t)
!write(13,'(/,a,/)')    ' variables =     "t"                              "Elec3 ** 2"'

filenameElec32r = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Elec32r.plt'
open(14,file=filenameElec32r)
!write(14,'(/,a,/)')    ' variables =     "r"                              "Elec3 ** 2"'

filenameElec32z = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Elec32z.plt'
open(15,file=filenameElec32z)
!write(15,'(/,a,/)')    ' variables =     "z"                             "Elec3 ** 2"'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenameElec32t  &
                                                                                    ,filenameElec32r  &
										    ,filenameElec32z
 read(*,*)

!------------------
filenamePsi3picksl = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Psi3picks l.plt'
open(16,file=filenamePsi3picksl)
!write(16,'(/,a,/)')    ' variables =     "l"                              "Psi3 picks ** 2"'

!------------------
filenamePsi2picksl = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Psi2picks l.plt'
open(17,file=filenamePsi2picksl)
!write(16,'(/,a,/)')    ' variables =     "l"                              "Psi2 picks ** 2"'

!------------------
filenameLscalemaxl= 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' Lscalemax l.plt'
open(20,file=filenameLscalemaxl)
!write(19,'(/,a,/)')    ' variables =     "l"                              "Lscalemax"'

!------------------
filenameibestl = 'E'//trim(EE)//' f'//trim(freqf)//' Np'//trim(Npf)//' tp'//trim(tpf)//' ibest l.plt'
open(21,file=filenameibestl)
!write(16,'(/,a,/)')    ' variables =     "l"                              "ibest"'

!**********************************************************************************************************************
!                                           Constants
!**********************************************************************************************************************
!-------------------------Optimization
!do E=0.1,0.5,0.1
!do freq=4000.,6000.,1000.
!do length=0.010,0.025,0.001 
!do omegaf=190.e-6,200.e-6,10.e-6

!E=2.
!freq=100.
!length=0.02	
!omegaf=100.e-6

!------------------------------------------------ Common
           pi = 4*atan(1.)                                                              !dimensionless
           Ii = (0.,1.)

     tbetween = 1./freq              !4.*tp                                             !s
        timet = Np*tbetween                                                             !s
          !nt = 1600                 !int(1600./20.) !80   int(tbetween/deltat)         !dimensionless 
      !deltat = tbetween / nt                                                           !s     
          !inn = 20

           nz = 5000 !20000@@@@@@@@@@                                                   !dimensionless
       length = 20.e-3                !length of crystal                                !m 
       deltaz = length/nz                                                               !m
           kn = 100 !10@@@@@@@@@@@@@@@

           nr = 50 !300@@@@@@@@
       radius = 0.002                !radius of crystal                                 !m
       omegaf = 80.e-6              !spot size                                          !m
      nomegaf = 5 
        nnrom = 9./10.
      !deltar = omegaf/10.                                                              !m
      deltar1 = (nomegaf*omegaf)/(int(nnrom*nr))
      deltar2 = (radius-(nomegaf*omegaf))/(int((1-nnrom)*nr))
     !deltar1 = deltar
     !deltar2 = deltar1
          !nr = int(radius/deltar)                                                      !dimensionless 
 
           Cp = 728.016              !heat capacity at constant pressure                !J/(kg.K)
 
          KT0 = 13.                  !thermal conductivity of KTP crystal               !W/(m.K)
          roh = 2945.                !mass density                                      !kg/m^3

    stability = 0.5
       deltat = stability * ( (roh*Cp)/(2.*KT0) ) * ( (deltar1**2.*deltaz**2.)/(deltar1**2.+deltaz**2.) ) !s   
          nt1 = int(tbetween/deltat)                                                    !dimensionless 
          inn = int(nt1/80)
           nt = nt1 - mod(nt1,inn)

        power = E/(sqrt(pi)*tp)                                                         !wat 

        gama1 = 0. !0.5                  !the absorption coefficient of fundomental wave    !1/m
        gama2 = 0. !0.5                  !the absorption coefficient of fundomental wave    !1/m
        gama3 = 0. !4.                   !the absorption coefficient of fundomental SHW     !1/m        

        no1T0 = 1.8296
        ne1T0 = 1.7466   
        ne2T0 = 1.7881   

      lambda1 = 1064.e-9             !wavelength fundamental                            !m
      lambda2 =  532.e-9             !wavelength second harmonic                        !m

         beta = 1.*pi/180.               !Axicon Angle
        alpha = 2.*pi/lambda1 *sin(beta) !wave vector                                   !1/m 

 dteta_bessel = pi/200.

!------------------------------------------------ Fields properties
       c = 3.e8                                                                         !m/s

   omega = 2.*pi*c/lambda1                                                              !rad/s 

epsilon0 = 8.85e-12                                                                     !C**2/N*m**2  or F(farad)/m

    deff = 7.3e-12           !nonlinear  effective coefficient                          !m/v       
	           
!**********************************************************************************************************************
!                                        Arrays Allocattion 
!**********************************************************************************************************************
!----------------------------------- Allocate Arrys Fields
allocate(Psi1(0:nt/inn,0:nr,1:2))                               
allocate(Psi2(0:nt/inn,0:nr,1:2))             
allocate(Psi3(0:nt/inn,0:nr,1:2))             

allocate(Elec1(0:nt/inn,0:nr,0:nz/kn))
allocate(Elec2(0:nt/inn,0:nr,0:nz/kn)) 
allocate(Elec3(0:nt/inn,0:nr,0:nz/kn))  

!**********************************************************************************************************************
!                                     Giving Zero to Arrays
!********************************************************************************************************************** 
!----------------------------------- Giving Zero to Arrys Fields
forall (i=0:nt/inn,j=0:nr,k=1:2)
                           Psi1(i,j,k) = (0.,0.)       
	  	           Psi2(i,j,k) = (0.,0.)       
                           Psi3(i,j,k) = (0.,0.)       
end forall !i

forall (i=0:nt/inn,j=0:nr,k=0:nz/kn)
                          Elec1(i,j,k) = (0.,0.)     
		          Elec2(i,j,k) = (0.,0.)      
                          Elec3(i,j,k) = (0.,0.)      
end forall !i

!**********************************************************************************************************************
!                                       Printing Constants     
!**********************************************************************************************************************
!------------------------------------------------ Common 
write(*,*)
write(*,*)'------- Common Constants ---------------------------------------------------'
write(*,*) 
write(*,'(A15,F15.10,//)')  '            E = ',E                

write(*,'(A15,2F15.10,/)')  '           Ii = ',Ii                  
write(*,'(A15,I5    ,/ )')  '           Nt = ',Nt
write(*,'(A15,I5    ,/ )')  '          inn = ',inn              
write(*,'(A15,I5    ,/ )')  '           Nr = ',Nr                
write(*,'(A15,I5    ,/ )')  '           Nz = ',Nz
write(*,'(A15,I5    ,/ )')  '           kn = ',kn 

write(*,'(A15,F15.10,/ )')  '           pi = ',pi
write(*,'(A15,F15.10,//)')  '           Cp = ',Cp              

write(*,'(A15,F15.10,/ )')  '          KT0 = ',KT0             
write(*,'(A15,F15.10,//)')  '          roh = ',roh              

write(*,'(A15,f15.10,/ )')  '        gama1 = ',gama1               
write(*,'(A15,f15.10,/ )')  '        gama2 = ',gama2               
write(*,'(A15,f15.10,//)')  '        gama3 = ',gama3              

write(*,'(A15,F15.10,/ )')  '        timet = ',timet           
write(*,'(A15,f15.10,/ )')  '        power = ',power              
write(*,'(A15,f15.5 ,//)')  '        alpha = ',alpha              

write(*,'(A15,F15.10,/ )')  '       omegaf = ',omegaf          
write(*,'(A15,F15.10,/ )')  '       length = ',length          
write(*,'(A15,F15.10,/ )')  '       deltat = ',deltat          
write(*,'(A15,F15.10,/ )')  '       deltaz = ',deltaz          
write(*,'(A15,F15.10,/ )')  '      deltar1 = ',deltar1          
write(*,'(A15,F15.10,//)')  '      deltar2 = ',deltar2          

write(*,'(A15,F15.10,//)')  '       radius = ',radius           

write(*,'(A15,F15.10,/ )')  '      lambda1 = ',lambda1         
write(*,'(A15,f15.10,//)')  '      lambda2 = ',lambda2        

write(*,'(A15,F15.10,//)')  '     tbetween = ',tbetween                                        

write(*,'(A15,F15.10,//)')  '    stability = ',stability                                                                

write(*,*)'----------------------------------------------------------------------------'
write(*,'(A,\)')' Press Enter to continue '
!read(*,*)

!------------------------------------------------ For fields Equation 
write(*,*)
write(*,*)'------- Field Equations Constants ------------------------------------------'
write(*,*)
write(*,'(A15,f15.3 ,/ )')  '            c = ',c                  

write(*,'(A15,f35.30,/ )')  '         deff = ',deff                

write(*,'(A15,f25.5 ,/ )')  '        omega = ',omega               

write(*,'(A15,f25.20,//)')  '     epsilon0 = ',epsilon0           

write(*,*)'----------------------------------------------------------------------------'
write(*,'(A,\)')' Press Enter to continue '
!read(*,*)
  
!**********************************************************************************************************************
!                                   Main Block of the Program     
!**********************************************************************************************************************
!----------- Optimization
!write(*,*)

!run=run+1

!write(*,*) 'Run = ',run

no1rT = no1T0 
ne1rT = ne1T0 
ne2rT = ne2T0

!----------------------------------------------------------------- Run program for NP pulse 
do l=1,Np 
         
   !----------------------------------------------- Field
   !--------
   do k=0,nz
      z=k*deltaz 
       
      !--------
      do i=1,nt              
         t=(i-1)*deltat

         if (mod(i,inn)==0) then 

	    !--------	        
   	    do j=1,nr-1
                  
               if (j<=int(nnrom*nr)) then
                  r=j*deltar1
                  deltar=deltar1
                 else
                  r=int(nnrom*nr)*deltar1+(j-int(nnrom*nr))*deltar2
                  deltar=deltar2
               end if  

	       !--------------------------------- Solve of bessel function 
               J0alphar=0.

               s1=0.
               s2=0.

               !------------- odd sentences
               do m=1,199,2

                  teta_bessel = m * dteta_bessel  
      
	          s1=s1+dcos( - (alpha*r)*dsin(teta_bessel) )

               end do !m

               !------------- even sentences
               do m=2,198,2

                  teta_bessel = m * dteta_bessel
   	 
                  s2=s2+dcos( - (alpha*r)*dsin(teta_bessel) )

               end do !m

               J0alphar = (1./pi) * (dteta_bessel/3.) * ( 1. + 4.*S1 + 2.*S2 + 1. ) 
               !--------------------------------- 
            
               !--------------------------------- interaction length
               Lscale =  sqrt(no1rT*ne1rT*ne2rT)                                              &
			   
		       * sqrt( (epsilon0*c**3.*pi*omegaf**2.) / (4.*omega**2.*deff**2.*power) ) 

               !------------- Maximom 
               if (Lscale .GE. Lscalemax) then

                  Lscalemax = Lscale 
                     
               end if

               !--------------------------------- Constants
	       cc1 = deltaz * no1rT / c
               !cc2 = deltaz *  Ii*c / (2.*no1rT*omega)
               cc3 = deltaz * gama1 / 2.
	       cc4 = deltaz *    Ii / Lscale

	       dd1 = deltaz * ne1rT / c
               !dd2 = deltaz *  Ii*c / (2.*ne1rT*omega)
	       dd3 = deltaz * gama2 / 2.
	       dd4 = deltaz *    Ii / Lscale

               ee1 = deltaz * ne2rT / c
               !ee2 = deltaz *  Ii*c / (4.*ne2rT*omega)
	       ee3 = deltaz * gama3 / 2.
	       ee4 = deltaz *    Ii / Lscale
	       !--------------------------------- Bounday conditions For Field Equations			

               !------------- Psi1
               if (k==0) Psi1(i/inn,j,1) = exp( (-(t-2.*tp)**2.)/(tp**2.) ) &

                                         * J0alphar * exp(-r**2./omegaf**2.)              !for input surface

	       Psi1(i/inn,0 ,1 ) = Psi1(i/inn,1,1 )                                           !for crystal axis
           
               Psi1(i/inn,nr,1 ) = (0.,0.)                                                !for lateral surface 
			

	       if (k==0 ) Psi1(i/inn,0 ,1 ) = Psi1(i/inn,1,1 )                                !for (0 ,0 )
           
   	       if (k==0 ) Psi1(i/inn,nr,1 ) = (0.,0.)                                         !for (nr,0 ) 

	       if (k==nz) Psi1(i/inn,nr,1 ) = (0.,0.)                                         !for (nr,nz) 

	       if (k==nz) Psi1(i/inn,0 ,1 ) = Psi1(i/inn,1,1)                                 !for (0 ,nz) 

	       !------------- Psi2
               if (k==0) Psi2(i/inn,j,1) = exp( (-(t-2.*tp)**2.)/(tp**2.) ) &        

                                         * J0alphar * exp(-r**2./omegaf**2.)	          !for input  surface       

               Psi2(i/inn,0 ,1 ) = Psi2(i/inn,1,1 )                                       !for crystal axis
           
               Psi2(i/inn,nr,1 ) = (0.,0.)                                                !for lateral surface 
          
	       if (k==0 ) Psi2(i/inn,0 ,1 ) = Psi2(i/inn,1,1 )                                !for (0 ,0 )
           
	       if (k==0 ) Psi2(i/inn,nr,1 ) = (0.,0.)                                         !for (nr,0 ) 

               if (k==nz) Psi2(i/inn,nr,1) = (0.,0.)                                      !for (nr,nz) 

	       if (k==nz) Psi2(i/inn,0 ,1) = Psi2(i/inn,1 ,1)                                 !for (0 ,nz) 

 	       !------------- Psi3 
               if (k==0)Psi3(i/inn,j,1) = (0.,0.)                                         !for input  surface
            
	       Psi3(i/inn,0  ,1) = Psi3(i/inn,1,1)                                            !for crystal axis
					      
	       Psi3(i/inn,nr ,1) = (0.,0.)                                                    !for lateral surface 

	       if (k==0 ) Psi3(i/inn,0 ,1 ) = (0.,0.)                                         !for (0 ,0 )
           
	       if (k==0 ) Psi3(i/inn,nr,1 ) = (0.,0.)                                         !for (nr,0 ) 

	       if (k==nz) Psi3(i/inn,nr,1 ) = (0.,0.)                                         !for (nr,nz) 

	       if (k==nz) Psi3(i/inn,0 ,1 ) = Psi3(i/inn,1,1)                                 !for (0 ,nz) 
  
	       !--------------------------------- End of Bounday conditions
           
               !--------------------------------- Field Equations		   
	       !-------------
	       Psi1(i/inn,j,2) =  Psi1(i/inn,j,1)                                                                 &
			
		               - cc1  * ( Psi1(i/inn,j,1) - Psi1(i/inn-1,j,1) ) / deltat                              &
			
		               + cc2  * ( Psi1(i/inn,j+1,1) - Psi1(i/inn,j-1,1) ) / (2*r*deltar)                      &
														
		    	       + cc2  * ( Psi1(i/inn,j+1,1) - 2*Psi1(i/inn,j,1) + Psi1(i/inn,j-1,1) ) / deltar**2     &
														
			       - cc3  *   Psi1(i/inn,j,1)                                                                 &

			       + cc4  * conjg(Psi2(i/inn,j,1)) * Psi3(i/inn,j,1) !* exp(-Ii*phasechange(i/inn,j,k) )   

 	       !-------------
	       Psi2(i/inn,j,2) =  Psi2(i/inn,j,1)                                                                 &
		
		               - dd1 * ( Psi2(i/inn,j,1) - Psi2(i/inn-1,j,1) ) / deltat                               &
			
		               + dd2 * ( Psi2(i/inn,j+1,1) - Psi2(i/inn,j-1,1) ) / (2*r*deltar)                       &
													
		               + dd2 * ( Psi2(i/inn,j+1,1) - 2*Psi2(i/inn,j,1) + Psi2(i/inn,j-1,1) ) / deltar**2      &
														
		               - dd3 *   Psi2(i/inn,j,1)                                                              &

		               + dd4 * conjg(Psi1(i/inn,j,1)) * Psi3(i/inn,j,1) !* exp(-Ii*phasechange(i/inn,j,k) )     
 
               !-------------			
	       Psi3(i/inn,j,2) =  Psi3(i/inn,j,1)                                                                 &
		
		               - ee1 * ( Psi3(i/inn,j,1) - Psi3(i/inn-1,j,1) ) / deltat                               &
			
		               + ee2 * ( Psi3(i/inn,j+1,1) -   Psi3(i/inn,j-1,1) ) / (2*r*deltar)                     &
														
		               + ee2 * ( Psi3(i/inn,j+1,1) - 2*Psi3(i/inn,j,1) + Psi3(i/inn,j-1,1) ) / deltar**2      &
														
			       - ee3 *   Psi3(i/inn,j,1)                                                                  &
														
			       + ee4 *   Psi1(i/inn,j,1) * Psi2(i/inn,j,1) !* exp(Ii*phasechange(i/inn,j,k) )          

               !------------- Maximom 
	       if (k==nz) then
			      
	          Psi22 = Psi2(i/inn,j,2) * conjg(Psi2(i/inn,j,2)) * 100
                  Psi32 = Psi3(i/inn,j,2) * conjg(Psi3(i/inn,j,2)) * 100
		  
                  !------
                  if (Psi22 >= Psi2max) then

                     Psi2max = Psi22 
                     ibest = i

                  end if

                  !------
                  if (Psi32 >= Psi3max) then

                     Psi3max = Psi32 
                     
                  end if
				    
               end if

               !------------------------------                 
	       if (mod(k,kn)==0) then
                     
	          !------------- Elec1
                  if (k==0) Elec1(i/inn,j,0) = Psi1(i/inn,j,1)       !for input  surface
 
	          Elec1(i/inn,0 ,k/kn ) = Psi1(i/inn,0,1)                !for crystal axis
           
	          Elec1(i/inn,nr,k/kn ) = (0.,0.)                        !for lateral surface 

                  if (k==0) Elec1(i/inn,0,0) = Psi1(i/inn,0,1)       !for (0 ,0 )
           
                  Elec1(i/inn,nr,0) = (0.,0.)                        !for (nr,0 ) 

                  Elec1(i/inn,nr,nz/kn) = (0.,0.)                    !for (nr,nz) 

                  if (k==nz)Elec1(i/inn,0 ,nz/kn) = Psi1(i/inn,0,1)  !for (0 ,nz) 

	          !------------- Elec2
                  if (k==0) Elec2(i/inn,j,0) = Psi2(i/inn,j,1)       !for input  surface
 
	          Elec2(i/inn,0 ,k/kn ) = Psi2(i/inn,0,1)                !for crystal axis
           
	          Elec2(i/inn,nr,k/kn ) = (0.,0.)                   	 !for lateral surface 

                  if (k==0) Elec2(i/inn,0,0) = Psi2(i/inn,0,1)       !for (0 ,0 )
           
                  Elec2(i/inn,nr,0) = (0.,0.)                        !for (nr,0 ) 

                  Elec2(i/inn,nr,nz/kn) = (0.,0.)                    !for (nr,nz) 

                  if (k==nz) Elec2(i/inn,0 ,nz/kn) = Psi2(i/inn,0,1) !for (0 ,nz) 

	          !------------- Elec3 
                  Elec3(i/inn,j,0) = (0.,0.)                         !for input  surface
 
	          Elec3(i/inn,0 ,k/kn ) = Psi3(i/inn,1,1)           	 !for crystal axis
           
	          Elec3(i/inn,nr,k/kn ) = (0.,0.)                  	     !for lateral surface 

                  Elec3(i/inn,0,0) = (0.,0.)                         !for (0 ,0 )
           
                  Elec3(i/inn,nr,0) = (0.,0.)                        !for (nr,0 ) 

                  Elec3(i/inn,nr,nz/kn) = (0.,0.)                    !for (nr,nz) 

                  if (k==nz) Elec3(i/inn,0,nz/kn) = Psi3(i/inn,1,2)  !for (0 ,nz) 

                  !---------------------
	          Elec1(i/inn,j,k/kn) = Psi1(i/inn,j,2)
		  Elec2(i/inn,j,k/kn) = Psi2(i/inn,j,2)
		  Elec3(i/inn,j,k/kn) = Psi3(i/inn,j,2)
              
	       end if

            end do !j
	 end if
      end do !i

      !-------------- End-Psi of each deltaz  ==> Initial Psi for next deltaz
      do i=0,nt   
         do j=1,nr-1
      	
            Psi1(i/inn,j,1) = Psi1(i/inn,j,2)
            Psi2(i/inn,j,1) = Psi2(i/inn,j,2)
            Psi3(i/inn,j,1) = Psi3(i/inn,j,2)

	     end do !j
       end do !i	     
       !-------------

   end do !k
   !----------------------------------------------- End of Field

   !============================================= Print Results for each deltat
   do i=0,nt
      t=(l-1)*nt*deltat + i*deltat 
    
      if (mod(i,inn)==0) then
	  
          !----- For Field Equations
 	  Elec12 = Elec1(i/inn,0,nz/kn) * conjg(Elec1(i/inn,0,0    )) * 100
	  write( 7,'(2x,f25.10,5x,f25.10)') t , Elec12 

	  Elec22 = Elec2(i/inn,0,nz/kn) * conjg(Elec2(i/inn,0,nz/kn)) * 100      
	  write(10,'(2x,f25.10,5x,f25.10)') t , Elec22 

	  Elec32 = Elec3(i/inn,0,nz/kn) * conjg(Elec3(i/inn,0,nz/kn)) * 100      
	  write(13,'(2x,f25.10,5x,f25.10)') t , Elec32
      
      end if
 
   end do !i 
   !--------------------------------------------- End of run for each deltat     

   !------------- for max & min
   write(16,'(2x,I5,2x,f8.2)') l   ,Psi3max
   psi3max = 0.

   !-----
   write(17,'(2x,I5,2x,f8.2)') l   ,Psi2max  
   psi2max = 0.

   !-----
   write(20,'(2x,I5,2x,f8.4)') l   ,Lscalemax  
   Lscalemax = 0.

   !-----
   write(21,'(2x,I5,2x,I5)') l  ,ibest 

end do !l
!-------------------------------------------------------------------- End of run for each Pulse

!**********************************************************************************************************************
!                                        Arrays Deallocattion 
!**********************************************************************************************************************

!----------------------------------- Deallocate Arrys Fields
!deallocate(Psi1)                               
!deallocate(Psi2)             
!deallocate(Psi3)             

!deallocate(Elec1)
!deallocate(Elec2) 
!deallocate(Elec3)  

!---------------
!end do !omegaf
!end do !Length
!end do !freq
!end do !E

!**********************************************************************************************************************
!                                        Printing Results     
!**********************************************************************************************************************

!------------------ For Field Equations
do i=0,nt 
   
   if (mod(i,4*inn)==0) then
      
      do j=0,nr/3.
   
         if (j<=int(nnrom*nr)) then
            r=j*deltar1
            deltar=deltar1
           else
            r=int(nnrom*nr)*deltar1+(j-int(nnrom*nr))*deltar2
            deltar=deltar2
         end if  

         Elec12 = Elec1(i/inn,j,0) * conjg(Elec1(i/inn,j,0)) * 100      
         write( 8,'(2x,f25.20,5x,f25.20)') r , Elec12
   
      end do !j      						   

      write(8,*) !'ZONE I =',i   
   
   end if !i  

end do !i

!------------------
do i=0,nt

   if (mod(i,4*inn)==0) then     
      
      do k=0,nz
         z=k*deltaz 
         
         if (mod(k,kn)==0) then
         
	    Elec12 = Elec1(i/inn,0,k/kn) * conjg(Elec1(i/inn,0,k/kn)) * 100      
            write( 9,'(2x,f25.20,5x,f25.20)') z , Elec12
   
 	 end if !k  
		    
      end do !k      						   

      write(9,*) !'ZONE I =',i   
   
   end if !i

end do !i

!==================
do i=0,nt 
 
   if (mod(i,4*inn)==0) then

      do j=0,nr/3.

         if (j<=int(nnrom*nr)) then
            r=j*deltar1
            deltar=deltar1
           else
            r=int(nnrom*nr)*deltar1+(j-int(nnrom*nr))*deltar2
            deltar=deltar2
         end if  

         Elec22 = Elec2(i/inn,j,nz/kn) * conjg(Elec2(i/inn,j,nz/kn)) * 100      
         write(11,'(2x,f25.20,5x,f25.20)') r , Elec22
   
      end do !j      						   

      write(11,*) !'ZONE I =',i       
   
   end if !i 

end do !i

!------------------
do i=0,nt 
   
   if (mod(i,4*inn)==0) then
      
      do k=0,nz
         z=k*deltaz 
         
         if (mod(k,kn)==0) then  
      
	    Elec22 = Elec2(i/inn,0,k/kn) * conjg(Elec2(i/inn,0,k/kn)) * 100      
            write(12,'(2x,f25.20,5x,f25.20)') z , Elec22
         
	 end if !k
      end do !k      						   

      write(12,*) !'ZONE I = ',i    
   
   end if !i

end do !i

!==================
do i=0,nt 
   
   if (mod(i,4*inn)==0) then
   
      do j=0,nr/3.

         if (j<=int(nnrom*nr)) then
            r=j*deltar1
            deltar=deltar1
           else
            r=int(nnrom*nr)*deltar1+(j-int(nnrom*nr))*deltar2
            deltar=deltar2
         end if  

         Elec32 = Elec3(i/inn,j,nz/kn) * conjg(Elec3(i/inn,j,nz/kn)) * 100      
         write(14,'(2x,f25.20,5x,f25.20)') r , Elec32
   
      end do !j      						   

      write(14,*) !'ZONE I =',i   
   
   end if !i  

end do !i      						   

!================== Transverse profile
!do j=0,nr

   !if (j<=int(4.*nr/5.)) then
      !r=j*deltar1
      !deltar=deltar1
     !else
      !r=5*omegaf+(j-int(4.*nr/5.))*deltar2
      !deltar=deltar2
   !end if  

   !do f=0,360,6
      !fi=f*pi/180.
 
      !Psi32 = Psi3(nt,j,nz) * conjg(Psi3(nt,j,nz))      
      !write(100,'(2x,f25.20,5x,f25.20,5x,f25.20)') r , fi , Psi32
   
  !end do !f

!end do !j      						   

!================== Transverse profile
  
!do j=0,nr

!   if (j<=int(4.*nr/5.)) then
!      r=j*deltar1
!      deltar=deltar1
!     else
!      r=5*omegaf+(j-int(4.*nr/5.))*deltar2
!      deltar=deltar2
!   end if  

   
!   do f=-200,+200
!      x=f*(r/200.)

!      y=sqrt(r**2-x**2)
!      y1 = y  ; y2 = -y

!      Psi32 = Psi3(nt,j,nz) * conjg(Psi3(nt,j,nz))      
!      write(200,'(2x,f25.20,5x,f25.20,5x,f25.20)') x , y1 , Psi32
!      write(200,'(2x,f25.20,5x,f25.20,5x,f25.20)') x , y2 , Psi32
 
!   end do !i
!end do !j

!------------------
do i=0,nt 
   
   if (mod(i,4*inn)==0) then
      
      do k=0,nz
         z=k*deltaz 
        
	 if (mod(k,kn)==0) then
      
	    Elec32 = Elec3(i/inn,0,k/kn) * conjg(Elec3(i/inn,0,k/kn)) * 100      
            write(15,'(2x,f25.20,5x,f25.20)') z , Elec32
		 
	 end if !k
		    
      end do !k    

      write(15,*) !'ZONE I = ',i    
   
   end if !i

end do !i

!**********************************************************************************************************************
!                                      Closing Files and Ending the Program 
!**********************************************************************************************************************

!------------------ For Field Equations
close( 7)
close( 8)
close( 9)

close(10)
close(11)
close(12)

close(13)
close(14)
close(15)

!----------------- For Max & Min
close(16)
close(17)
close(18)
close(19)
close(20)
close(21)

end program Elec_BG_PW                     

!======================================================================================================================
        
 
