# Thermal Effects in Second Harmonic Generation (SHG)

**Contents**  
[1. About this Repository](#1-about-this-repository)    
&nbsp;&nbsp;&nbsp;&nbsp;[1.1. Who Is This Tutorial For?](#11-who-is-this-tutorial-for)           
&nbsp;&nbsp;&nbsp;&nbsp;[1.2. What Will You Learn?](#12-what-will-you-learn)  
&nbsp;&nbsp;&nbsp;&nbsp;[1.3. Prerequisites](#13-prerequisites)   
&nbsp;&nbsp;&nbsp;&nbsp;[1.4. Contents of this Repository](#14-contents-of-this-repository)         
[2. Second Harmonic Generation (SHG)](#2-second-harmonic-generation-shg)       
[3. Thermal Challenge in SHG](#3-thermal-challenge-in-shg)          
[4. Thermal Gradient in a Crystal during SHG](#4-thermal-gradient-in-a-crystal-during-shg)  
[5. Reducing Computational Cost](#5-reducing-computational-cost)  
[6. Boundry Conditions](#6-boundry-conditions)   
&nbsp;&nbsp;&nbsp;&nbsp;[6.1. Our Contribution (1)](#61-our-contribution-1)     
[7. Phase Mismatch](#7-phase-mismatch)      
&nbsp;&nbsp;&nbsp;&nbsp;[7.1. Our Contribution (2)](#71-our-contribution-2)        
[8. Field](#8-field)        
&nbsp;&nbsp;&nbsp;&nbsp;[8.1. Our Contribution (3)](#81-our-contribution-2)        
[9. Coupling Heat, Phase Mismatch, and Field](#9-coupling-heat-phase-mismatch-and-field)        
&nbsp;&nbsp;&nbsp;&nbsp;[9.1. Our Contribution (4)](#91-our-contribution-4)        
[10. Our Contribution](#10-our-contribution)        
&nbsp;&nbsp;&nbsp;&nbsp;[10.1. Methodology](#101-methodology)        
&nbsp;&nbsp;&nbsp;&nbsp;[10.2. Computational Approach using Finite Difference Method (FDM)](#102-computational-approach-using-finite-difference-method-fdm)        
[11. Research Opportunities](#11-research-opportunities)        
[12. How to Cite Us](#12-how-to-cite-us)        
[13. For Additional Questions](#13-for-additional-questions)        


# 1. About this Repository
This GitHub repository offers comprehensive guidance, from basic to advanced levels, for computationally addressing thermal effects in Second Harmonic Generation (SHG). As an educational resource, this repository starts with covering fundamental aspects of Fortran, including how to install it and master its essential commands. Also, we demonstrate techniques for computationally solving a nonlinear optics phenomenon using the Finite Difference Method (FDM), provide access to the codes utilized in our studies, and explain our research findings clearly. Also, we outline potential research opportunities for future exploration. Our ongoing efforts involve expanding the repository to incorporate further advancements in the field. 

## 1.1. Who Is This Tutorial For?   
This tutorial is designed for anyone interested in computational physics, nonlinear optics, or scientific computing, regardless of their prior experience. Whether you're a student, researcher, or professional, this resource will guide you through the process of solving Thermal Effects in SHG using FDM.

## 1.2. What Will You Learn?    
By the end of this tutorial, you will:           
- Gain proficiency in Fortran, from installation to mastering essential commands.
- Understand the basic to advanced principles of Thermal Effects in SHG.
- Learn FDM and how to apply it for solving nonlinear optics problems computationally.
- Access and utilize real codes used in this cutting-edge research.
- Clearly comprehend the research findings and understand the potential for future exploration in this area.

## 1.3. Prerequisites  
This tutorial is designed for three types of researcher:

1. **For those who are familiar with SHG and Fortran**: You can dive straight into the research phase. The codes and topics provided in this repository are meant to deepen your knowledge and assist in further studies.

2. **For those who know Fortran but are new to SHG**: This tutorial will introduce you to the fundamentals of SHG, guiding you step-by-step through the key concepts. By the end, you'll be ready to tackle complex problems in this field.

3. **For beginners with no prior knowledge of SHG or Fortran**: This repository is built with you in mind. We’ll start with the basics, teaching you Fortran from the ground up, followed by an introduction to SHG. Our goal is to help you progress from understanding the basics to solving advanced physics and engineering problems.


## 1.4. Contents of this Repository 

```
Folder PATH listing
. 
|
+---0. Cite Us
|       1. Heat Equation _ Analytical.pdf
|       2. Heat G_CW.pdf
|       3. Heat Equation _ Computational _ Source G_PW.pdf
|       4. Phase Mismatch G_PW.pdf
|       5. Ideal G_CW.pdf
|       6. Ideal BG_PW.pdf
|       7. Coupled G_CW.pdf
|       8. SHG G_CW _ Computational _ Approximate model.pdf
|       README.md
|       
+---1. ifort Installation Guide
|       README.md
|       
+---2. FORTRAN Tutorial
|   |   1. FORTRAN_Main Commands Tutorial.md
|   |   2. FORTRAN_Coding Template Tutorial.md
|   |   3. FORTRAN_Coding Template.F90
|   |   4. Write_Read_Variables Types.F90
|   |   5. Readable Code Structure.F90
|   |   6. do loop.F90
|   |   7. If then else.F90
|   |   8. open file.F90
|   |   9. Array.F90
|   |   README.md
|   |   
|   \---Persian materials
|           1. Recommended  Textbooks .ppsx
|           2. Scan - Book - Safari .pdf
|           3. Fortran.pptx
|           4. Fortran - Framework.pdf
|           
+---3. SHG Codes
|       Temp_G_PW_01_06_92.F90
|       Temp_Phase_G_Pw_08_06_92.F90
```         


# 2. Second Harmonic Generation (SHG)
SHG employs a nonlinear crystal like Potassium Titanyl Phosphate (KTP) to convert red laser (1064 nm) into green laser (532 nm). This conversion is essential because of green light's necessity and difficulty in direct production. During the SHG process, a powerful laser beam interacts with the crystal, causing it to emit light at exactly half the wavelength of the incoming beam, effectively doubling the light's frequency. 

<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image01.png" alt="Image 1">
</p>

<p align="center"> <strong>Figure 1.</strong> During SHG, a nonlinear crystal like KTP converts a 1064 nm red laser (Fundamental Wave) into 532 nm green laser (Second Harmonic Wave), effectively doubling the frequency of the original beam through a nonlinear optical process.</p>

# 3. Thermal Challenge in SHG
As the SHG process occurs, some of the input energy is not perfectly converted into the desired higher-frequency photons. Instead, a portion of this energy is lost as heat within the nonlinear crystal or medium. The dissipated heat within the nonlinear crystal reduces efficiency by causing thermal dephasing, which disrupts phase matching. This temperature increase can also lead to crystal damage, further lowering the conversion efficiency and output power.

To address this issue, the crystal is equipped with a cooling system depicted in the Figure 2. A coolant circulates around the crystal, absorbing the heat thereby maintaining an optimal temperature for a more efficient SHG. The crystal's lateral surface is maintained at a constant temperature through cooling. Typically, a double layer of copper covers these surfaces, with either water or liquid nitrogen passing through it. This ensures a constant temperature condition at the crystal's side surface. Additionally, the input and output surfaces of the crystal are cooled through both radiation and convection. Heat reaches these surfaces through conduction and is then transferred away by convection and radiation processes.

<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image02.png" alt="Image 2" width="65%">
</p>

<p align="center"> <strong>Figure 2.</strong> A cooling system uses a double layer of copper and circulating coolant (water or liquid nitrogen) to manage heat dissipation, keeping the crystal at an optimal temperature for efficient laser performance during SHG.</p>


# 4. Thermal Gradient in a Crystal during SHG
The thermal gradient within a crystal subjected to laser radiation is shown in Figure 3. Since the axis of the crystal is under the laser beam and its lateral surface is in heat exchange with the environment, the peak of the temperature gradient is at the center of the crystal, where the laser beam is focused, and gradually decreases towards the surface. Also, the left side of the crystal, where the laser initially contacts, is the hottest, with temperature decreasing towards the right side. This spatial temperature variation is crucial in understanding the thermal behaviour of crystals under laser irradiation.


<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image03.png" alt="Image 3" width="100%">
</p>

<p align="center"> <strong>Figure 3.</strong> Visualization of the thermal gradient in a crystal exposed to laser radiation. The hottest point is at the center where the laser is focused, with temperature decreasing outward toward the edges. The left side, where the laser first hits, is the hottest, with the temperature gradually cooling as it moves to the right.</p>

To achieve optimal performance, it is crucial to fully understand the underlying physics.  The heat conduction equation describes the temporal and spatial evolution of temperature \( T \) within a crystal exposed to a heat source \( S \):

$$
+\rho c \frac{\partial T}{\partial t}-\vec{\nabla} \cdot (K(T) \vec{\nabla} T) = S
$$

where mass density ($ρ$) in terms of $kg {m}^{-3}$, heat capacity ($c$) in terms of $J k g^{-1} K^{-1}$ and thermal conductivity ($K$) in terms of $W m^{-1} K^{-1}$ at $K(T)=K_0 \times T_0 / T$, which depends on temperature. In this formula, the conductivity coefficient $K_0$ is in temperature $T_0=300 \mathrm{~K}$. Also, ($S$) is the pulsed source that produces heat in terms of $Wm^{-3}$. 


### 4.1. Variation of the Heat Gradient

- **Time Evolution ($\frac{\partial T}{\partial t}$):**
The term $\frac{\partial T}{\partial t}$ describes how temperature changes over time at a specific point in the crystal. Initially, the laser heats the crystal rapidly, leading to a steep temperature gradient. Over time, the rate of temperature increase slows as the crystal absorbs more heat.

- **Spatial Distribution ($\nabla T$):**
The gradient $\nabla T$ indicates how temperature varies across the crystal. Near the laser's focal point, the gradient is steep, reflecting a sharp temperature rise. Away from this region, the gradient flattens as heat spreads.

- **Heat Diffusion:**
$\nabla \cdot (K(T) \nabla T)$ describes the diffusion of heat, with $K(T)$ being temperature-dependent thermal conductivity. Heat diffuses more efficiently in regions with higher $K(T)$, leading to a smoother gradient over time.


### 4.2. Recap
Below is an explanation of the definitions of the Derivative and Gradient:

- **Derivative**:  
  The derivative is a fundamental concept in calculus that quantifies the rate at which a function changes with respect to its independent variable. For a given function \( f(x) \), the derivative, denoted as \( f'(x) \) or \($\frac{df}{dx}$\), represents the instantaneous rate of change of the function with respect to \( x \). It provides information about the slope of the tangent line to the function at any point \( x \). Mathematically, the derivative is defined as:

$$
f^{\prime}(x)=\lim _{\Delta x \rightarrow 0} \frac{f(x+\Delta x)-f(x)}{\Delta x}
$$

  This limit describes how the function \( f(x) \) changes as the input \( x \) is varied infinitesimally.

- **Gradient**:   
  The gradient is an extension of the derivative concept to functions of multiple variables. For a scalar function \( f(x, y, z) \), which depends on several independent variables, the gradient \($\nabla f$\) is a vector that points in the direction of the steepest rate of increase of the function. The magnitude of this vector represents the rate of change in that direction. The gradient is composed of the partial derivatives of the function with respect to each independent variable:

$$
\nabla f=\left(\frac{\partial f}{\partial x} , \frac{\partial f}{\partial y} , \frac{\partial f}{\partial z}\right)
$$

  Each component $\frac{\partial f}{\partial x}$ , $\frac{\partial f}{\partial y}$ , and $\frac{\partial f}{\partial z}$ measures the rate of change of $f$ with respect to the corresponding variable while holding the other variables constant.     


# 5. Reducing Computational Cost  
To efficiently analyze Second Harmonic Generation (SHG) in a KTP crystal modeled as a cylinder, we can reduce the workload by taking advantage of the crystal’s symmetry. Instead of examining the entire cylindrical shape, we focus on a simpler two-dimensional half-plane, a rectangular section that represents one side of the cylinder. This works because of the crystal’s symmetry along its lateral axis; studying this smaller section captures the behavior of the full cylinder. By visualizing this half-plane rotating around the horizontal axis, we effectively account for the whole cylindrical structure. This approach significantly reduces the size of the problem, decreasing the number of calculations needed. It allows for high-resolution analysis with greater efficiency, maintaining the accuracy of the SHG study while saving substantial time and resources.

<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image04.png" alt="Image 4" width="45%">
</p>

<p align="center"> <strong>Figure 4.</strong> Efficient SHG Analysis in KTP Crystals Using Symmetry: The image shows how a two-dimensional half-plane leverages the crystal’s symmetry to represent the entire cylinder, reducing computational complexity while maintaining accurate SHG analysis. 
</p>


# 6. Boundry Conditions  
The boundary conditions for heat transfer within the nonlinear crystal during SHG are critical for accurately modeling thermal behavior. The lateral surfaces of the crystal are maintained at a constant temperature through a cooling system, ensuring effective heat dissipation. Meanwhile, the input and output faces are cooled by radiation and convection, allowing heat to escape efficiently. The crystal axis, which experiences the highest temperatures, is treated as an insulated boundary due to its role as a thermal peak where heat naturally flows outward toward the cooler lateral surfaces. This boundary setup is essential for precisely modeling the heat distribution, which is crucial for optimizing the crystal’s SHG performance and preventing thermal disruptions.

<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image05.png" alt="Image 5" width="65%">
</p>

<p align="center"> <strong>Figure 5.</strong> 
Boundary conditions for heat transfer in SHG crystals: The lateral surfaces are cooled to a constant temperature, facilitating heat dissipation, while input and output faces are cooled by radiation and convection. The crystal axis, experiencing peak temperatures, is treated as an insulated boundary, ensuring accurate heat distribution modeling critical for optimizing SHG performance.
</p>

## 6.1. Our Contribution (1)
Our contributions to understanding these heat transfer dynamics are detailed in the following publications:

- **Heat Equation _ Continuous Wave Gaussian _ Analytical** [(Link)](https://opg.optica.org/ao/abstract.cfm?uri=ao-47-13-2317)   
This work focuses on predicting temperature distributions in laser crystals using a Continuous Wave Gaussian source. The analytical model provided insights into the basic thermal behavior in solid-state lasers, a critical step toward designing more efficient systems by accurately modeling heat within complex crystal structures.

- **Heat Equation _ Continuous Wave Gaussian _ Computational** [(Link)](https://link.springer.com/article/10.1007/s13538-014-0291-x)   
Building upon the analytical work, this computational study incorporated more realistic factors, such as temperature-dependent thermal conductivity and radiation effects. The model demonstrated the significant impact of these often-overlooked factors on heat distribution in KTP crystals, enhancing the thermal modeling of laser systems.

- **Heat Equation _ Pulsed Wave Gaussian _ Computational** [(Link)](https://opg.optica.org/ao/abstract.cfm?uri=ao-54-6-1241)      
This study developed a numerical model for heat distribution under Pulsed Gaussian conditions, highlighting the critical role of variable thermal conductivity, especially when radiation effects are minimal. The findings improved the accuracy of predicting heat behavior in pulsed laser systems, contributing to more effective thermal management strategies.


# 7. Phase Mismatch
During SHG, a crystal is subjected to laser radiation, the temperature at various points within the crystal becomes spatially and temporally dependent. This variation in temperature causes corresponding changes in the crystal's refractive index, making the refractive index also a function of position and time. Since the speed of light in a medium is dependent on its refractive index, the speed of light traveling through different regions of the crystal will similarly be a function of position and time. Specifically, the temperature gradient within the crystal causes the speed of light to vary radially. Consequently, different regions of the wavefront experience different speeds, leading to distortions in the wavefront shape. This results in a phase mismatch between the fundamental and second harmonic waves. In different crystals, the wavefronts may be convex or concave. Figure 6. shows a concave wavefront.

<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image06.png" alt="Image 6" width="75%">
</p>

<p align="center"> <strong>Figure 6.</strong> Schematic of the phase mismatch due to temprature gradient within the crystal. In different crystals the wavefronts may be convex or concave. This figure shows a concave wavefront.</p>


Phase is a function of temperature and it is clear that due to the presence of the phase difference in the field equations, this quantity is very effective in determining the efficiency of the SHG. In this way, heat and electromagnetic waves are related indirectly. The below formula is the correlation between phase ($φ$) and temperature ($T$):

$$
\Delta \varphi=\int_0^z \Delta k(T) d z^{\prime}
$$

we can effectively integrate heat considerations into electromagnetic equations, thereby advancing our comprehension of how thermal effects impact the efficiency of nonlinear optical phenomena.

## 7.1. Our Contribution (2)
Our study on this phenomenon is detailed in the publication:

- **Phase Mismatch _ Pulsed Wave Gaussian _ Computational** [(Link)](https://www.researchgate.net/publication/267926440_Thermally_induced_phase_mismatching_in_a_repetitively_Gaussian_pulsed_pumping_KTP_crystal_A_spatiotemporal_treatment)    
This work addresses the issue of Thermally Induced Phase Mismatching (TIPM) in KTP crystals under Pulsed Wave Gaussian conditions. The study developed a spatiotemporal model to examine how temperature rise influences nonlinear conversion efficiency, highlighting the critical need to manage TIPM to optimize SHG performance in pulsed laser applications. The findings emphasize the importance of precise thermal management strategies to reduce phase mismatches and improve overall system efficiency.

# 8. Field
In ideal conditions where there is no heat dissipation or phase mismatch, all of the fundamental wave is perfectly converted into the Second Harmonic Wave, demonstrating the maximum possible efficiency in nonlinear crystals. 

<p align="center">
  <img src="./Archive/images/3.%20Readme_images/image01.png" alt="Image 1">
</p>

<p align="center"> <strong>Figure 1.</strong> Maximum SHG Efficiency: Achieving full conversion of the Fundamental Wave (FW) to the Second Harmonic Wave (SHW) under ideal thermal and phase conditions.
</p>

## 8.1. Our Contribution (3)
Our research contributions exploring these ideal conditions are detailed in the following publications:

- **Ideal _ Continuous Wave Gaussian _ Computational** [(Link)](https://opg.optica.org/ao/abstract.cfm?uri=ao-54-4-869)    
This study explored SHG efficiency under Continuous Wave Gaussian conditions, highlighting how temperature fluctuations can prevent achieving ideal conversion efficiency. The study found that even minor temperature increases could drastically reduce SHG efficiency due to beam depletion and refractive index changes, highlighting the importance of temperature control in optimizing SHG processes.
 
- **Ideal _ Pulsed Wave Bessel Gaussian _ Computational** [(Link)](https://opg.optica.org/ao/abstract.cfm?uri=ao-53-32-7691)   
This research introduced a model using Pulsed Bessel-Gauss beams, challenging traditional assumptions like the nondepleted wave approximation. This study provided a more accurate framework for SHG by considering wave depletion effects, demonstrating the impact of beam profile on heat and SHG efficiency under pulsed conditions.


# 9. Coupling Heat, Phase Mismatch, and Field
To accurately model SHG, it is essential to recognize that heat, phase mismatch, and the field are interconnected and must be treated as a coupled system. 

The **heat** governs the temperature distribution within the crystal, which directly influences the refractive index through temperature-dependent material properties. This, in turn, affects the phase mismatch between the interacting waves. The **phase mismatch** modifies the interaction conditions for the fundamental and second harmonic waves, directly impacting their efficiency and the field dynamics. The **field** itself defines how the energy is transferred between the waves, how the phase evolves, and how the power distribution affects local heating within the crystal. This heating further alters the temperature profile, creating a feedback loop that perpetuates continuously.

Because of this tightly coupled nature, solving these equations independently would fail to capture the dynamic interactions and feedback mechanisms that occur in real-world conditions. A proper approach is necessary to account for these interdependencies, ensuring a comprehensive and accurate representation of SHG performance under varying thermal and optical conditions.

## 9.1. Our Contribution (4)
Our research contributions for coupling Heat, Phase Mismatch, and Field are detailed in the following publications:

- **Coupled _ Continuous Wave Gaussian _ Computational** [(Link)](https://opg.optica.org/oe/fulltext.cfm?uri=oe-22-21-25615&id=302163)    
This study advanced our understanding by incorporating both Thermally Induced Phase Mismatching (TIPM) and thermal lensing into SHG models using a Continuous Wave Gaussian source. By coupling eight different equations, the model captured the dynamic interactions between heat and SHG efficiency over time. This comprehensive approach provided a realistic simulation that closely matched experimental results, significantly enhancing our understanding of thermal effects in Continuous Wave Gaussian SHG systems.

- **SHG _ Continuous Wave Gaussian _ Computational _ Approximate Model** [(Link)](https://opg.optica.org/oe/fulltext.cfm?uri=oe-18-18-18732&id=205211)    
This article introduced a simplified theoretical model to study TIPM in Continuous Wave Gaussian SHG systems. By coupling heat dissipation and SHG equations, it demonstrated how temperature gradients and thermal dispersion negatively impact conversion efficiency. The study highlighted the importance of managing thermal effects, especially at higher power levels, to maintain optimal SHG performance.


# 10. Our Contribution
Our study introduces a novel computational model utilizing the FDM to examine heat dissipation dynamics within a crystal and their consequential impacts during SHG. This model is more accurate, easier to use, and provides a more comprehensive and detailed understanding of the thermal effects in SHG.

## 10.1. Methodology
Direct measurement of internal temperature in a nonlinear crystal during Second Harmonic Generation (SHG) is not feasible due to the lack of instrumentation capable of probing the localized, transient thermal effects within the crystal. Conventional techniques are insufficient to resolve the rapid, spatially confined temperature changes, making experimental assessment impractical.

Analytical approaches are similarly challenging due to the complexity of the coupled field, heat, and phase equations governing SHG, as shown below: 

$$
\begin{aligned}
& \frac{n_1}{c} \frac{d \psi_1}{d t}+\frac{d \psi_1}{d z}-\frac{i c}{2 n_1 \omega} \frac{1}{r} \frac{d \psi_1}{d r}-\frac{i c}{2 n_1 \omega} \frac{d^2 \psi_1}{d^2 r}+\frac{\gamma_1}{2} \psi_1=\frac{i}{L_T} \psi_2^* \psi_3 e^{-i \Delta \phi} \\
& \frac{n_2}{c} \frac{d \psi_2}{d t}+\frac{d \psi_2}{d z}-\frac{i c}{2 n_2 \omega} \frac{1}{r} \frac{d \psi_2}{d r}-\frac{i c}{2 n_2 \omega} \frac{d^2 \psi_2}{d^2 r}+\frac{\gamma_2}{2} \psi_2=\frac{i}{L_T} \psi_1^* \psi_3 e^{-i \Delta \phi} \\
& \frac{n_3}{c} \frac{d \psi_2}{d t}+\frac{d \psi_3}{d z}-\frac{i c}{4 n_3 \omega} \frac{1}{r} \frac{d \psi_3}{d r}-\frac{i c}{4 n_3 \omega} \frac{d^2 \psi_3}{d^2 r}+\frac{\gamma_3}{2} \psi_3=\frac{i}{L_T} \psi_1 \psi_3 e^{i \Delta \phi} \\
& +\rho c \frac{\partial T}{\partial t}-\vec{\nabla} \cdot (K(T) \vec{\nabla} T)=S \\
& \frac{d \Delta \phi}{d z}=\frac{2 \pi}{\lambda_1}\left[\Delta n_1(T)+\Delta n_2(T)-2 \Delta n_3(T)\right]
\end{aligned}
$$

Analytical solution of these equations requires simplifying assumptions that deviate the model from reality. For example, even the fundamental heat equation which plays a crucial role in this domain, relies on such simplifications. 

However, through computational approaches, we've pushed the boundaries, avoiding any simplifying assumptions to offer a more precise model. For instance, we no longer assume the thermal conductivity coefficient to be constant; instead, it dynamically varies with temperature throughout time. This shift from traditional analytical models, which rely on simplifying assumptions, enables a more accurate study of nonlinear optics phenomena.

## 10.2. Computational Approach using Finite Difference Method (FDM)
We use the FDM as the computational method to model thermal effects in SHG due to its low computational cost and user-friendly nature. FDM offers simplicity in both learning and application. Since heat operates on a macroscopic scale and doesn't vary drastically, FDM provides accurate results without the need for using other complex methods. Its straightforward approach that efficiently captures the thermal dynamics involved in SHG without unnecessary complexity. 

The FDM approximates derivatives in differential equations through discretization of the domain into a grid and replacing derivatives with finite difference expressions. This transforms the equation into a system of algebraic equations solvable with numerical techniques. FDM's accuracy and stability rely on discretization, approximation schemes, and solution methods chosen, offering a versatile and efficient approach for solving complex differential equations when analytical solutions are impractical. By employing FDM, we achieve cost-effective and accurate simulations, making it an ideal choice for modelling thermal effects in SHG processes.


**Finite Difference Approximation**: In FDM, we approximate the derivative by choosing a small but finite value of $( \Delta x $). This turns the derivative into a simple difference equation:

$$
f'(x) \approx \frac{f(x + \Delta x) - f(x)}{\Delta x}
$$

This is called the **Forward Difference** approximation.

**A Simple Example**: Approximating the Derivative of $( f(x) = x^2 $)

Let's use FDM to approximate the derivative of the function $( f(x) = x^2 $) at $( x = 2 $).

1. Function Value at $( x = 2 $):
   
$$
f(2) = 2^2 = 4
$$

2. Function Value at $( x + \Delta x $) with $( \Delta x = 0.1 $):
   
$$
f(2 + 0.1) = (2.1)^2 = 4.41
$$

3. Approximate Derivative:
   
$$
f'(2) \approx \frac{f(2 + 0.1) - f(2)}{0.1} = \frac{4.41 - 4}{0.1} = 4.1
$$

4. Exact Derivative for Comparison:
   
   The exact derivative of $( f(x) = x^2 $) is $( f'(x) = 2x $). At $( x = 2 $), this gives:
   
$$
f'(2) = 2 \times 2 = 4
$$

**Note**: The approximation $( 4.1 $) is close to the exact derivative $( 4 $). As $( \Delta x $) becomes smaller, the approximation improves, demonstrating how FDM essentially breaks down the derivative to its fundamental definition by using a discrete, manageable calculation.

This method helps solve differential equations where exact analytical solutions are not feasible, highlighting how FDM is powerful in practical computational scenarios.

# 11. Research Opportunities
While our research has focused on modeling the thermal effects in KTP crystals under specific wave conditions—namely Gaussian Continuous Wave (GCW), Gaussian Pulsed Wave (GPW), and Bessel-Gauss Pulsed Wave (BGPW)—significant opportunities remain for further exploration in this field. A natural extension of this work involves investigating different types of wave sources and alternative nonlinear crystals. 

The key distinction between different wave types and crystal materials lies in the unique heat transfer behavior each combination exhibits. These variations can fundamentally alter the thermal dynamics within the crystal, leading to different impacts on SHG efficiency and system performance. For instance, the use of other waveforms, such as Hermite-Gaussian or Laguerre-Gaussian beams, or employing different crystal materials like Lithium Niobate or Beta Barium Borate, would result in diverse heat distribution patterns and necessitate distinct thermal management strategies. 

Moreover, while our studies primarily focused on developing detailed models to understand these thermal behaviors, future research can take the next step by conducting comprehensive simulations and experimental validations. These simulations can provide deeper insights into how different thermal properties, such as anisotropic conductivity and varying boundary conditions, influence heat distribution and phase mismatching across various laser and crystal configurations. Such explorations would not only enhance the theoretical understanding but also offer practical guidelines for optimizing SHG systems under varying thermal conditions.

Ultimately, exploring these new avenues will open up entirely new research trajectories, each with its own set of challenges and opportunities. This will enable a more complete understanding of thermal effects in nonlinear optical systems, contributing to the design of more efficient and adaptable laser technologies. To facilitate further research, other researchers can use our GitHub repository as a tutorial, utilizing the provided source code as a foundation for conducting simulations and extending the study in this field. By building on our work, future studies can deepen insights into heat transfer dynamics and refine thermal management strategies across various nonlinear optical configurations.

We are currently working on these topics and are dedicated to pushing this research forward. In the near future, we will share the results of our ongoing studies along with new articles and code updates in this repository. Our goal is to expand thes repository into a more comprehensive resource that will support other researchers in exploring the complexities of thermal effects in nonlinear optics. Stay tuned for these updates, as they will provide even deeper insights and tools to enhance future research in this field.

# 12. How to Cite Us
Please refer to the [0. Cite Us](https://github.com/mohammad-ghadri/SHG__Second_Harmonic_Generation/tree/main/0.%20Cite%20Us) folder for accurate citations. It contains essential guidelines for accurate referencing, ensuring accurate acknowledgement of our work.


# 13. For Additional Questions
If you have questions that are not covered in the resources above, the best way to reach [Mostafa M. Rezaee](https://www.linkedin.com/in/mostafa-m-rezaee/).    
- Gmail: mostafa.mohammadrezaee@gmail.com       
- [Linkedin](https://www.linkedin.com/in/mostafa-m-rezaee/)           
- [Personal Website](https://mostafa-mr.com/)       


# Archive

## 1.3. Prerequisites  
This tutorial assumes no prior knowledge. It starts from the basics and gradually progresses to more advanced topics. A willingness to learn and a basic understanding of physics and mathematics will be helpful.

## Here's a basic example of FDM:   
We will use the following term from the heat equation:

$$
+\rho c \frac{\partial T}{\partial t}-\vec{\nabla} \cdot (K(T) \vec{\nabla} T) = S
$$

where:
- $\rho$ is the density of the material.
- $c$ is the specific heat capacity.
- $T$ is the temperature.
- $t$ is time.
- $K(T)$ is the thermal conductivity, which is temperature-dependent.
- $\vec{\nabla}$ denotes the gradient operator.
- $S$ represents a heat source term.


**Expanding the Heat Equations**:   
The original heat conduction equation given is:

$$
+\rho c \frac{\partial T}{\partial t} - \vec{\nabla} \cdot (K(T) \vec{\nabla} T) = S
$$


First, we to focus on is the divergence of the heat flux:

$$
\vec{\nabla} \cdot (K(T) \vec{\nabla} T)
$$

This term represents the flow of heat within the material.
Using the product rule of the divergence operator, this term can be expanded as:

$$
\vec{\nabla} \cdot (K(T) \vec{\nabla} T) = \nabla K_T \cdot \nabla T + K_T \nabla^2 T
$$

Here's how this expansion works:
1. Divergence of a Product: The divergence of the product of a scalar field $K(T)$ and a vector field $\vec{\nabla} T$ can be expressed as:

$$
\vec{\nabla} \cdot (K(T) \vec{\nabla} T) = (\vec{\nabla} K(T)) \cdot (\vec{\nabla} T) + K(T) \vec{\nabla} \cdot (\vec{\nabla} T)
$$

2. Breaking it Down:
   - $\nabla K_T \cdot \nabla T$: This term represents the dot product of the gradient of the thermal conductivity (which depends on temperature) and the gradient of the temperature.
   - $K_T \nabla^2 T$: This term represents the thermal conductivity multiplied by the Laplacian of the temperature ($\nabla^2 T$), which is essentially the divergence of the temperature gradient.

Substituting this expanded form into the original equation:

$$
+\rho c \frac{\partial T}{\partial t} - (\nabla K_T \cdot \nabla T + K_T \nabla^2 T) = S
$$

Rearrange the terms:

$$
\rho c \frac{\partial T}{\partial t} - \nabla K_T \cdot \nabla T - K_T \nabla^2 T = S
$$

Also, the source term $S$ is then expressed as a sum of separate contributions, each representing different physical sources or effects:

$$
S = \gamma_1 \psi_1 + \gamma_2 \psi_2 + \gamma_3 \psi_3
$$

Here:
- $\gamma_1$, $\gamma_2$, and $\gamma_3$ are coefficients representing the weights or intensities of different source terms.
- $\psi_1$, $\psi_2$, and $\psi_3$ are specific source functions related to different physical processes contributing to heat generation.

Thus, the final form of the equation becomes:

$$
\rho c \frac{\partial T}{\partial t} - \nabla K_T \cdot \nabla T - K_T \nabla^2 T = \gamma_1 \psi_1 + \gamma_2 \psi_2 + \gamma_3 \psi_3
$$

This equation models heat transfer in a medium where the thermal conductivity depends on temperature, and the heat source is composed of multiple contributing factors.



Where:
- $T$ is the temperature.
- $\rho$ is the density.
- $c$ is the specific heat capacity.
- $K_T$ is the thermal conductivity.
- $\gamma_1$, $\gamma_2$, $\gamma_3$ are coefficients related to the fields $\psi_1$, $\psi_2$, $\psi_3$.

### Applying FDM:

To solve this equation using FDM, follow these steps:

#### 1. Discretize Space and Time

- **Time:** Divide time into small intervals with spacing $\Delta t$.
- **Space:** Divide the spatial domain into small intervals (grid points) with spacing $\Delta x$.

#### 2. Approximate Derivatives

- **Time Derivative:**

$$
\frac{\partial T}{\partial t} \approx \frac{T_i^{n+1}-T_i^n}{\Delta t}
$$

  Where $T_i^n$ is the temperature at grid point $x_i$ at time step $t^n$, and $T_i^{n+1}$ is the temperature at the next time step.

- **Spatial Derivatives:**

  - **First Derivative (Gradient):**

$$
\nabla T \approx \frac{T_{i+1}^n-T_{i-1}^n}{2 \Delta x}
$$

-  
  - **Second Derivative (Laplacian):**

$$
\nabla^2 T \approx \frac{T_{i+1}^n-2 T_i^n+T_{i-1}^n}{(\Delta x)^2}
$$

#### 3. Formulate the Finite Difference Equation

Substitute the finite difference approximations into the heat equation:

$$
\rho c \frac{T_i^{n+1} - T_i^n}{\Delta t} - K_T \frac{T_{i+1}^n - 2 T_i^n + T_{i-1}^n}{(\Delta x)^2} = \gamma_1 \psi_1^n + \gamma_2 \psi_2^n + \gamma_3 \psi_3^n
$$

Rearrange to solve for $T_i^{n+1}$:

$$
T_i^{n+1} = T_i^n + \frac{\Delta t}{\rho c} \left( \gamma_1 \psi_1^n + \gamma_2 \psi_2^n + \gamma_3 \psi_3^n + K_T \frac{T_{i+1}^n - 2 T_i^n + T_{i-1}^n}{(\Delta x)^2} \right)
$$

#### 4. Repalcing Numbers

Assume:
- Grid spacing $\Delta x = 0.1$ meter
- Time step $\Delta t = 0.01$ seconds
- Thermal conductivity $K_T = 0.5$
- Density $\rho = 1.0$
- Specific heat capacity $c = 1.0$

For a grid point $x_i$ in the middle of the rod:

$$
T_i^{n+1} = T_i^n + \frac{0.01}{1.0 \times 1.0} \left( \gamma_1 \psi_1^n + \gamma_2 \psi_2^n + \gamma_3 \psi_3^n + 0.5 \frac{T_{i+1}^n - 2 T_i^n + T_{i-1}^n}{(0.1)^2} \right)
$$

This equation allows you to compute the temperature at the next time step, and updating the entire temperature distribution over time.


## 8.3. Our Contribution
The series of our studies presented here aim to address the critical issue of thermal effects in Second Harmonic Generation (SHG) processes, particularly in KTP crystals. Each article contributes a piece to solving this problem by exploring various conditions and wave types, building a comprehensive approach to managing heat in nonlinear optical systems. (Note: To make it easier to understand our achievements, the articles are sorted by thematic order, not by chronological order.)

1. **Heat Equation _ Analytical** [(Link)](https://opg.optica.org/ao/abstract.cfm?uri=ao-47-13-2317)   
We began with "Heat Equation _ Analytical," focusing on the problem of predicting temperature distributions in laser crystals using a Continuous Wave Gaussian source. This analytical model helped us understand the basic thermal behavior in solid-state lasers, crucial for designing better systems by modeling heat within complex crystal structures.

2. **Heat G_CW** [(Link)](https://link.springer.com/article/10.1007/s13538-014-0291-x)   
The next step, "Heat G_CW," expanded on this by incorporating more realistic factors such as temperature-dependent thermal conductivity and radiation effects using a Continuous Wave Gaussian as the source. This computational study showed how these factors, often ignored, significantly impact heat distribution in KTP crystals, enhancing the thermal modeling of lasers.

3. **Heat Equation _ Computational _ Source G_PW** [(Link)](https://opg.optica.org/ao/abstract.cfm?uri=ao-54-6-1241)      
Shifting to Pulsed Wave Gaussian conditions, "Heat Equation _ Computational _ Source G_PW" developed a numerical model for heat distribution in KTP crystals under Pulsed Gaussian beams. This computational study highlighted that while radiation effects are minimal, variable thermal conductivity plays a crucial role, improving accuracy in predicting heat behavior in pulsed laser systems.

4. **Phase Mismatch G_PW** [(Link)](https://www.researchgate.net/publication/267926440_Thermally_induced_phase_mismatching_in_a_repetitively_Gaussian_pulsed_pumping_KTP_crystal_A_spatiotemporal_treatment)    
"Phase Mismatch G_PW" tackled the problem of Thermally Induced Phase Mismatching (TIPM) in KTP under Pulsed Wave Gaussian conditions. It modeled how temperature rise affects nonlinear conversion efficiency, showing the need to manage TIPM to optimize SHG performance in pulsed applications.

5. **Ideal G_CW** [(Link)](https://opg.optica.org/ao/abstract.cfm?uri=ao-54-4-869)    
In "Ideal G_CW," we addressed how temperature affects SHG efficiency in double-pass cavities under Continuous Wave Gaussian conditions. The study found that even minor temperature increases could drastically reduce SHG efficiency due to beam depletion and refractive index changes, highlighting the importance of temperature control in optimizing SHG processes.
 
6. **Ideal BG_PW** [(Link)](https://opg.optica.org/ao/abstract.cfm?uri=ao-53-32-7691)   
"Ideal BG_PW" introduced a model using Pulsed Bessel-Gauss beams, challenging traditional assumptions like the nondepleted wave approximation. This study provided a more accurate framework for SHG by considering wave depletion effects, demonstrating the impact of beam profile on heat and SHG efficiency under pulsed conditions.

7. **Coupled G_CW** [(Link)](https://opg.optica.org/oe/fulltext.cfm?uri=oe-22-21-25615&id=302163)    
"Coupled G_CW" advanced our understanding by integrating TIPM and thermal lensing in SHG with a Continuous Wave Gaussian source. By coupling eight equations, this model captured the dynamic effects of heat on SHG efficiency over time, offering a realistic simulation aligned with experimental data and enhancing our grasp of thermal influences in Continuous Wave Gaussian SHG.

8. **SHG G_CW _ Computational _ Approximate model** [(Link)](https://opg.optica.org/oe/fulltext.cfm?uri=oe-18-18-18732&id=205211)    
Finally, "SHG G_CW _ Computational _ Approximate Model" provided a simpler theoretical approach to understanding TIPM in Continuous Wave Gaussian SHG systems. By coupling heat and SHG equations, it showed how temperature gradients and thermal dispersion reduce conversion efficiency, emphasizing the need for effective thermal management, especially at higher powers.

Each study tackled specific thermal challenges under various conditions—GCW, GPW, and BGPW—showing how heat affects SHG efficiency and offering solutions to reduce these effects. This cohesive approach not only addresses the theoretical aspects of heat distribution and phase mismatching but also provides practical insights for optimizing nonlinear optical systems. Together, these works create a robust framework for understanding and managing thermal impacts in SHG, paving the way for more efficient laser technologies.
