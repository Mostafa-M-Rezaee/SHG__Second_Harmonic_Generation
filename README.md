# Thermal Effect in Second Harmonic Generation (SHG)

# Table of contents
[1. About this Repository](#1-about-this-repository)       
[2. Second Harmonic Generation (SHG)](#2-second-harmonic-generation-(shg))       
[3. The Challenge of SHG](#3-the-challenge-of-shg)          
[4. Thermal Gradient in a Crystal during SHG](#4-thermal-gradient-in-a-crystal-during-shg)  
[5. Phase Mismatch](#5-phase-mismatch)        
[6. Our Contribution](#6-our-contribution)        
[7. Methodology](#7-methodology)        
&nbsp;&nbsp;&nbsp;&nbsp;[7.1. Computational Approach](#71-computational-approach)        
&nbsp;&nbsp;&nbsp;&nbsp;[7.2. Finite Difference Method (FDM)](#72-finite-difference-method-(fdm))        
[8. Research Opportunities](#8-research-opportunities)        
[9. How to Cite Us](#9-how-to-cite-us)        
[10. For Additional Questions](#10-for-additional-questions)        

# 1. About this Repository
This GitHub repository offers comprehensive guidance, from basic to advanced levels, for computationally addressing thermal effects in Second Harmonic Generation (SHG). As an educational resource, this repository starts with covering fundamental aspects of Fortran, including how to install it and master its essential commands. Also, we demonstrate techniques for computationally solving a nonlinear optics phenomenon using the Finite Difference Method (FDM), provide access to the codes utilized in our studies, and explain our research findings clearly. Also, we outline potential research opportunities for future exploration. Our ongoing efforts involve expanding the repository to incorporate further advancements in the field. 

**Who Is This Tutorial For?**   
This tutorial is designed for anyone interested in computational physics, nonlinear optics, or scientific computing, regardless of their prior experience. Whether you're a student, researcher, or professional, this resource will guide you through the process of solving Thermal Effects in SHG using FDM.

**What Will You Learn?**    
By the end of this tutorial, you will:           
- Gain proficiency in Fortran, from installation to mastering essential commands.
- Understand the basic to advanced principles of Thermal Effects in SHG.
- Learn FDM and how to apply it for solving nonlinear optics problems computationally.
- Access and utilize real codes used in cutting-edge research.
- Clearly comprehend the research findings and understand the potential for future exploration in this area.

**Prerequisites**:  
This tutorial assumes no prior knowledge. It starts from the basics and gradually progresses to more advanced topics. A willingness to learn and a basic understanding of physics and mathematics will be helpful, but not strictly required.


```
Folder Path Listing
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



# 2. Second Harmonic Generation (SHG)
SHG employs a nonlinear crystal like Potassium Titanyl Phosphate (KTP) to convert red laser (1064 nm) into green laser (532 nm). This conversion is essential because of green light's necessity and difficulty in direct production. During the SHG process, a powerful laser beam interacts with the crystal, causing it to emit light at exactly half the wavelength of the incoming beam, effectively doubling the light's frequency. 

<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image01.png" alt="Image 1">
</p>

<p align="center">Figure 1. During SHG, a KTP nonlinear crystal converts a 1064 nm red laser (Fundamental Wave) into 532 nm green laser (Second Harmonic Wave), effectively doubling the frequency of the original beam through a nonlinear optical process.</p>
 
A portion of the absorebed energy during SHG is dissipated as heat within the crystal which potentially damages it and reduces its ability to produce the desired laser. To address this issue, the crystal is equipped with a cooling system depicted in the figure 2. A coolant circulates around the crystal, absorbing the heat thereby maintaining an optimal temperature for a more efficient SHG. The crystal's lateral surface is maintained at a constant temperature through cooling. Typically, a double layer of copper covers these surfaces, with either water or liquid nitrogen passing through it. This ensures a constant temperature condition at the crystal's side surface. Additionally, the input and output surfaces of the crystal are cooled through both radiation and convection. Heat reaches these surfaces through conduction and is then transferred away by convection and radiation processes.

<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image02.png" alt="Image 2" width="65%">
</p>

<p align="center">Figure 2. A cooling system uses a double layer of copper and circulating coolant (water or liquid nitrogen) to manage heat dissipation, keeping the crystal at an optimal temperature for efficient laser performance during SHG.</p>

Considering the lateral symmetry allows us to tackle the issue within a half-plane of a cylindrical crystal, as depicted in Figure 3. Visualize rotating this half-plane around the horizontal axis, which effectively encompasses the entire cylindrical crystal. With the axis of the crystal exhibiting the highest temperature and the side surface the lowest, a temperature gradient naturally forms from the axis towards the side surface. On the other hand, the maximum temperature of the crystal axis causes heat to always move from the axis to the surface; like the top of a hill that is higher than all the points around it. Therefore, the axis of the crystal acts like an insulator. In other words, if we want to solve the problem on the half plane as shown in Figure 3, we must apply the isolation boundary condition for the crystal axis.

<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image03.png" alt="Image 3" width="65%">
</p>

<p align="center">Firgure 3. Schematic of the upper half plane of the cross section of the crystal in the longitudinal direction, the vector perpendicular to the entrance and exit surfaces of the crystal, as well as the temperature gradient vector</p>


# 3. Thermal Challenge in SHG
As the SHG process occurs, some of the input energy is not perfectly converted into the desired higher-frequency photons. Instead, a portion of this energy is lost as heat within the nonlinear crystal or medium. The dissipated heat within the nonlinear crystal reduces efficiency by causing thermal dephasing, which disrupts phase matching. This temperature increase can also lead to crystal damage, further lowering the conversion efficiency and output power.


# 4. Thermal Gradient in a Crystal during SHG
The thermal gradient within a crystal subjected to laser radiation is shown in Figure 4. Since the axis of the crystal is under the laser beam and its lateral surface is in heat exchange with the environment, the peak of the temperature gradient is at the center of the crystal, where the laser beam is focused, and gradually decreases towards the surface. Also, the left side of the crystal, where the laser initially contacts, is the hottest, with temperature decreasing towards the right side. This spatial temperature variation is crucial in understanding the thermal behaviour of crystals under laser irradiation.


<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image04.png" alt="Image 4" width="75%">
</p>

<p align="center">Figure 4. Visualization of the thermal gradient in a crystal exposed to laser radiation. The hottest point is at the center where the laser is focused, with temperature decreasing outward toward the edges. The left side, where the laser first hits, is the hottest, with the temperature gradually cooling as it moves to the right.</p>

To achieve optimal performance, it is crucial to fully understand the underlying physics.  The heat conduction equation describes the temporal and spatial evolution of temperature \( T \) within a crystal exposed to a heat source \( S \):

$$
+\rho c \frac{\partial T}{\partial t}-\vec{\nabla} \cdot K(T) \vec{\nabla} T=S
$$

where mass density ($ρ$) in terms of $kg {m}^{-3}$, heat capacity ($c$) in terms of $J k g^{-1} K^{-1}$ and thermal conductivity ($K$) in terms of $W m^{-1} K^{-1}$ at $K(T)=K_0 \times T_0 / T$, which depends on temperature. In this formula, the conductivity coefficient $K_0$ is in temperature $T_0=300 \mathrm{~K}$. Also, ($S$) is the pulsed source that produces heat in terms of $Wm^{-3}$. 


### Variation of the Heat Gradient:

1. **Time Evolution ($\frac{\partial T}{\partial t}$):**
   - The term $\frac{\partial T}{\partial t}$ describes how temperature changes over time at a specific point in the crystal.
   - Initially, the laser heats the crystal rapidly, leading to a steep temperature gradient. Over time, the rate of temperature increase slows as the crystal absorbs more heat.

2. **Spatial Distribution ($\nabla T$):**
   - The gradient $\nabla T$ indicates how temperature varies across the crystal. Near the laser's focal point, the gradient is steep, reflecting a sharp temperature rise. Away from this region, the gradient flattens as heat spreads.

3. **Heat Diffusion:**
   - $\nabla \cdot (K(T) \nabla T)$ describes the diffusion of heat, with $K(T)$ being temperature-dependent thermal conductivity. Heat diffuses more efficiently in regions with higher $K(T)$, leading to a smoother gradient over time.


### Derivative and Gradient Explanation:
- **Derivative**:  
  The derivative is a fundamental concept in calculus that quantifies the rate at which a function changes with respect to its independent variable. For a given function \( f(x) \), the derivative, denoted as \( f'(x) \) or \($\frac{df}{dx}$\), represents the instantaneous rate of change of the function with respect to \( x \). It provides information about the slope of the tangent line to the function at any point \( x \). Mathematically, the derivative is defined as:

  $$
  f'(x) = \lim_{\Delta x \to 0} \frac{f(x + \Delta x) - f(x)}{\Delta x}
  $$

  This limit describes how the function \( f(x) \) changes as the input \( x \) is varied infinitesimally.

- **Gradient**:   
  The gradient is an extension of the derivative concept to functions of multiple variables. For a scalar function \( f(x, y, z) \), which depends on several independent variables, the gradient \($ \nabla f $\) is a vector that points in the direction of the steepest rate of increase of the function. The magnitude of this vector represents the rate of change in that direction. The gradient is composed of the partial derivatives of the function with respect to each independent variable:

  $$
  \nabla f = \left( \frac{\partial f}{\partial x}, \frac{\partial f}{\partial y}, \frac{\partial f}{\partial z} \right)
  $$

  Each component $ \frac{\partial f}{\partial x} $, $ \frac{\partial f}{\partial y} $, and $ \frac{\partial f}{\partial z} $ measures the rate of change of $ f $ with respect to the corresponding variable while holding the other variables constant.


# 5. Phase Mismatch
During SHG, a crystal is subjected to laser radiation, the temperature at various points within the crystal becomes spatially and temporally dependent. This variation in temperature causes corresponding changes in the crystal's refractive index, making the refractive index also a function of position and time. Since the speed of light in a medium is dependent on its refractive index, the speed of light traveling through different regions of the crystal will similarly be a function of position and time. Specifically, the temperature gradient within the crystal causes the speed of light to vary radially. Consequently, different regions of the wavefront experience different speeds, leading to distortions in the wavefront shape. This results in a phase mismatch between the fundamental and second harmonic waves. In different crystals, the wavefronts may be convex or concave. Figure 5 shows a concave wavefront.

<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image05.png" alt="Image 5" width="75%">
</p>

<p align="center">Figure 5. Schematic of the phase mismatch due to temprature gradient within the crystal. In different crystals the wavefronts may be convex or concave. This figure shows a concave wavefront.</p>


Phase is a function of temperature and it is clear that due to the presence of the phase difference in the field equations, this quantity is very effective in determining the efficiency of the SHG. In this way, heat and electromagnetic waves are related indirectly. The below formula is the correlation between phase ($φ$) and temperature ($T$):

$$
\Delta \varphi=\int_0^z \Delta k(T) d z^{\prime}
$$

we can effectively integrate heat considerations into electromagnetic equations, thereby advancing our comprehension of how thermal effects impact the efficiency of nonlinear optical phenomena.


# 6. Our Contribution
Our study introduces a novel computational model utilizing the FDM to examine heat dissipation dynamics within a crystal and their consequential impacts during SHG. This model is more accurate, easier to use, and provides a more comprehensive and detailed understanding of the thermal effects in SHG.

# 7. Methodology

## 7.1. Computational Approach
When delving into field, heat and phase equations, we encounter coupled equations, as shown below: 

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

## 7.2. Finite Difference Method (FDM)
We use the Finite Difference Method (FDM) to model thermal effects in SHG due to its low computational cost and user-friendly nature. FDM offers simplicity in both learning and application. Since heat operates on a macroscopic scale and doesn't vary drastically, FDM provides accurate results without the need for using other complex methods. Its straightforward approach that efficiently captures the thermal dynamics involved in SHG without unnecessary complexity. 

The FDM approximates derivatives in differential equations through discretization of the domain into a grid and replacing derivatives with finite difference expressions. This transforms the equation into a system of algebraic equations solvable with numerical techniques. FDM's accuracy and stability rely on discretization, approximation schemes, and solution methods chosen, offering a versatile and efficient approach for solving complex differential equations when analytical solutions are impractical. By employing FDM, we achieve cost-effective and accurate simulations, making it an ideal choice for modelling thermal effects in SHG processes.


# 8. Research Opportunities


# 9. How to Cite Us
Please refer to the [0. Cite Us](https://github.com/mohammad-ghadri/SHG__Second_Harmonic_Generation/tree/main/0.%20Cite%20Us) folder for accurate citations. It contains essential guidelines for accurate referencing, ensuring accurate acknowledgement of our work.


# 10. For Additional Questions
If you have questions that are not covered in the resources above, the best way to reach [Mostafa M. Rezaee](https://www.linkedin.com/in/mostafa-m-rezaee/).    
- Gmail: mostafa.mohammadrezaee@gmail.com       
- [Linkedin](https://www.linkedin.com/in/mostafa-m-rezaee/)           
- [Personal Website](https://mostafa-mr.com/)       




