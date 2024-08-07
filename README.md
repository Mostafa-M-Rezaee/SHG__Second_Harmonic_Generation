# Thermal Effect in Second Harmonic Generation


# Table of contents
[1. About this Repository](#1-about-this-repository)       
[2. Second Harmonic Generation ](#2-second-harmonic-generation)       
[3. The Challenge of Second Harmonic Generation](#3-the-challenge-of-second-harmonic-generation)          
[4. Thermal Effects in Second Harmonic Generation](#4-thermal-effects-in-second-harmonic-generation)        
[5. Our Contribution](#5-our-contribution)        
[6. Methodology](#6-methodology)        
&nbsp;&nbsp;&nbsp;&nbsp;[6.1. Computational Approach](##61-computational-approach)        
&nbsp;&nbsp;&nbsp;&nbsp;[6.2. Finite Difference Method](##62-finite-difference-method)        
[7. Research Opportunity](#7-research-opportunity)        
[8. How to Cite Us](#8-how-to-cite-us)        
[9. For Additional Question](#9for-additional-question)        


# 1. About this Repository
This GitHub repository offers comprehensive guidance, from basic to advanced levels, for computationally addressing thermal effects in Second Harmonic Generation (SHG). As an educational resource, this repository starts with covering fundamental aspects of Fortran, including how to install it and master its essential commands. Also, we demonstrate techniques for computationally solving a nonlinear optics phenomenon using the Finite Difference Method (FDM), provide access to the codes utilized in our studies, and explain our research findings clearly. Also, we outline potential research opportunities for future exploration. Our ongoing efforts involve expanding the repository to incorporate further advancements in the field. 


# 2. Second Harmonic Generation 
Second Harmonic Generation (SHG) employs a nonlinear crystal like KTP to convert red light into green light. This conversion is essential because of green light's necessity and difficulty in direct production. During the SHG process, a powerful laser beam interacts with the crystal, causing it to emit light at exactly half the wavelength of the incoming beam, effectively doubling the light's frequency. 

<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image01.png" alt="Image 1">
</p>

<p align="center">Figure 1. Second Harmonic Generation</p>

The crystal absorbs a portion of the laser's energy as heat which potentially damages it and reduces its ability to produce the desired light. To address this issue, the crystal is equipped with a cooling system depicted in the figure 2. A coolant circulates around the crystal, absorbing the heat generated during SHG, thereby maintaining an optimal temperature for a more efficient SHG. The crystal's lateral surface is maintained at a constant temperature through cooling. Typically, a double layer of copper covers these surfaces, with either water or liquid nitrogen passing through it. This ensures a constant temperature condition at the crystal's side surface. Additionally, the input and output surfaces of the crystal are cooled through both radiation and convection. Heat reaches these surfaces through conduction and is then transferred away by convection and radiation processes.

<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image02.png" alt="Image 2" width="65%">
</p>

<p align="center">Figure 2. Schematic of the nonlinear crystal under laser beam</p>


Considering the lateral symmetry allows us to tackle the issue within a half-plane of a cylindrical crystal, as depicted in Figure 4. Visualize rotating this half-plane around the horizontal axis, which effectively encompasses the entire cylindrical crystal. With the axis of the crystal exhibiting the highest temperature and the side surface the lowest, a temperature gradient naturally forms from the axis towards the side surface. On the other hand, the maximum temperature of the crystal axis causes heat to always move from the axis to the surface; like the top of a hill that is higher than all the points around it. Therefore, the axis of the crystal acts like an insulator. In other words, if we want to solve the problem on the half plane as shown in Figure 4, we must apply the isolation boundary condition for the crystal axis.

<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image03.png" alt="Image 3" width="65%">
</p>

<p align="center">Firgure 3. Schematic of the upper half plane of the cross section of the crystal in the longitudinal direction, the vector perpendicular to the entrance and exit surfaces of the crystal, as well as the temperature gradient vector</p>


# 3. The Challenge of Second Harmonic Generation
In the process of Second Harmonic Generation (SHG), a crystal is subjected to laser radiation. A portion of the laser's energy is absorbed by the crystal, converting it into thermal energy. This absorbed thermal energy elevates the crystal's temperature, thereby inducing optical characteristic alterations, significantly reducing its efficiency. To achieve optimal performance, it is crucial to fully understand the underlying physics. The absorbed heat can be characterized through the Heat Equation, which states:

$$
+\rho c \frac{\partial T}{\partial t}-\vec{\nabla} \cdot K(T) \vec{\nabla} T=S
$$

where mass density ($ρ$) in terms of $kg {m}^{-3}$, heat capacity ($c$) in terms of $J k g^{-1} K^{-1}$ and thermal conductivity ($K$) in terms of $W m^{-1} K^{-1}$ at $K(T)=K_0 \times T_0 / T$, which depends on temperature. In this formula, the conductivity coefficient $K_0$ is in temperature $T_0=300 \mathrm{~K}$. Also, ($S$) is the pulsed source that produces heat in terms of $Wm^{-3}$. 


# 4. Thermal Effects in Second Harmonic Generation
The thermal distribution within a crystal subjected to laser radiation is shown in Figure 1. Since the axis of the crystal is under the laser beam and its lateral surface is in heat exchange with the environment, the peak of the temperature gradient is at the center of the crystal, where the laser beam is focused, and gradually decreases towards the surface. Also, the left side of the crystal, where the laser initially contacts, is the hottest, with temperature decreasing towards the right side. This spatial temperature variation is crucial in understanding the thermal behaviour of crystals under laser irradiation.


<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image04.png" alt="Image 4" width="75%">
</p>

<p align="center">Figure 4. Crystal Temprature Gradient</p>


When a crystal is placed under the radiation of a laser, the temperature of different points of the crystal becomes a function of space and time. As the temperature changes, the refractive index of the crystal also changes. In this way, the refractive index of the crystal will also become a function of space and time. Since the speed of light is dependent on the refractive index of the medium, the speed of light passing through different points of the crystal will also be a function of space and time. As a matter of fact, the speed of light changes in the radial directions due to the temperature gradient. Hence,  the different points of the wavefront experience different speeds and the shape of the wavefront is disturbed. This results in a phase mismatch between the primary and second harmonic waves, as shown in Figure 5.


<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image05.png" alt="Image 5" width="75%">
</p>

<p align="center">Figure 5. Schematic of the phase mismatch due to temprature gradient within the crystal.</p>


Phase is a function of temperature and it is clear that due to the presence of the phase difference in the field equations, this quantity is very effective in determining the efficiency of the SHG. In this way, heat and electromagnetic waves are related indirectly. By delineating the correlation between phase ($φ$) and temperature ($T$), as shown below:

$$
\Delta \varphi=\int_0^z \Delta k(T) d z^{\prime}
$$

we can effectively integrate heat considerations into electromagnetic equations, thereby advancing our comprehension of how thermal effects impact the efficiency of nonlinear optical phenomena.


# 5. Our Contribution
Our study introduces a novel computational model to examine the dynamics of heat absorption within a crystal and its consequential impacts during the Second Harmonic Generation. Our model surpasses the accuracy of prior methodologies, offering a more comprehensive understanding of this phenomenon.


# 6. Methodology

## 6.1. Computational Approach
When delving into field, heat and phase equations, we encounter a web of coupled equations, as shown below: 

$$
\begin{aligned}
& \frac{n_1}{c} \frac{d \psi_1}{d t}+\frac{d \psi_1}{d z}-\frac{i c}{2 n_1 \omega} \frac{1}{r} \frac{d \psi_1}{d r}-\frac{i c}{2 n_1 \omega} \frac{d^2 \psi_1}{d^2 r}+\frac{\gamma_1}{2} \psi_1=\frac{i}{L_T} \psi_2^* \psi_3 e^{-i \Delta \phi} \\
& \frac{n_2}{c} \frac{d \psi_2}{d t}+\frac{d \psi_2}{d z}-\frac{i c}{2 n_2 \omega} \frac{1}{r} \frac{d \psi_2}{d r}-\frac{i c}{2 n_2 \omega} \frac{d^2 \psi_2}{d^2 r}+\frac{\gamma_2}{2} \psi_2=\frac{i}{L_T} \psi_1^* \psi_3 e^{-i \Delta \phi} \\
& \frac{n_3}{c} \frac{d \psi_2}{d t}+\frac{d \psi_3}{d z}-\frac{i c}{4 n_3 \omega} \frac{1}{r} \frac{d \psi_3}{d r}-\frac{i c}{4 n_3 \omega} \frac{d^2 \psi_3}{d^2 r}+\frac{\gamma_3}{2} \psi_3=\frac{i}{L_T} \psi_1 \psi_3 e^{i \Delta \phi} \\
& +\rho c \frac{\partial T}{\partial t}-\nabla K_T \cdot \nabla T-K_T \nabla^2 T=\gamma_1 \psi_1+\gamma_2 \psi_2+\gamma_3 \psi_3 \\
& \frac{d \Delta \phi}{d z}=\frac{2 \pi}{\lambda_1}\left[\Delta n_1(T)+\Delta n_2(T)-2 \Delta n_3(T)\right]
\end{aligned}
$$

Analytical solution of these equations requires simplifying assumptions that deviate the model from reality. For example, even the fundamental heat equation which plays a crucial role in this domain, relies on such simplifications. However, through computational approaches, we've pushed the boundaries, avoiding any simplifying assumptions to offer a more precise model. For instance, we no longer assume the thermal conductivity coefficient to be constant; instead, it dynamically varies with temperature throughout time. This shift from traditional analytical models, which rely on simplifying assumptions, enables a more accurate study of nonlinear optics phenomena.

## 6.2. Finite Difference Method
We use the Finite Difference Method (FDM) to model thermal effects in SHG due to its low computational cost and user-friendly nature. FDM offers simplicity in both learning and application. Since heat operates on a macroscopic scale and doesn't vary drastically, FDM provides accurate results without the need for using other complex methods. Its straightforward approach efficiently captures the thermal dynamics involved in SHG without unnecessary complexity. 

The Finite Difference Method (FDM) approximates derivatives in differential equations through discretization of the domain into a grid and replacing derivatives with finite difference expressions. This transforms the equation into a system of algebraic equations solvable with numerical techniques. FDM's accuracy and stability rely on discretization, approximation schemes, and solution methods chosen, offering a versatile and efficient approach for solving complex differential equations when analytical solutions are impractical. By employing FDM, we achieve cost-effective and accurate simulations, making it an ideal choice for modelling thermal effects in SHG processes.


# 7. Research Opportunity


# 8. How to Cite Us
Please refer to the [0. Cite Us](https://github.com/mohammad-ghadri/SHG__Second_Harmonic_Generation/tree/main/0.%20Cite%20Us) folder for accurate citations. It contains essential guidelines for accurate referencing, ensuring accurate acknowledgement of our work.


# 9. For Additional Question





```
Folder PATH listing
. 
|
+---0. Cite Us
|       1. Heat Equation Analytical.pdf
|       2. SHG G_CW Computational.pdf
|       3. Coupled G_CW.pdf
|       4. Heat G_CW.pdf
|       5. Ideal BG_PW.pdf
|       6. Ideal G_CW.pdf
|       7. Heat G_PW.pdf
|       8. Phase Mismatch G_PW.pdf
|       README.md
|       
+---1. ifort Installation Guide
|       README.md
|       
+---2. FORTRAN Tutorial
|   |   1. FORTRAN main commands.md
|   |   2. FORTRAN coding template.md
|   |   3. General Structure .F90
|   |   4. Write & Read & type variables .F90
|   |   5. Legible code.F90
|   |   6. do loop .F90
|   |   7. If then else .F90
|   |   8. open file  .F90
|   |   9. Array .F90
|   |   README.md
|   |   
|   \---Persian materials
|           1. Recommended  Textbooks .ppsx
|           2. Scan - Book - Safari .pdf
|           3. Fortran.pptx
|           Fortran - Framework.pdf
|           
+---3. Literature Review
|       1. SHG _ Gaussian _ Continious wave _ Heat _ Brazilian Journal of Physics.pdf
|       2. SHG _ Gaussian _ Pulsed wave _ Heat _ Applied Optics.pdf
|       21 - BG-PW-ideal fields - Applied Optics.pdf
|       23 - G-CW-ideal fields-Applied Optics.pdf
|       25 - G-PW-Phase-Applied Optics.pdf
|       Applied Optics ___ Phase mismatch __ Source - Gaussian pulse wave.pdf
|       SHG _ Gaussian _ Continious wave - Coupled-Optics Express.pdf
|       
+---4. Codes
|       Temp_G_PW_01_06_92.F90
|       Temp_Phase_G_Pw_08_06_92.F90
```         


