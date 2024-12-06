# FORTRAN Coding Template

## Abstract
Programming efficient code to solve numerical challenges involves various complexities. Therefore, significant coding experience is essential to resolve syntax errors and achieve accurate results. This article introduces a structured format designed to help coders minimize common errors, promoting accurate and straightforward solutions.

## Introduction
Although this article focuses on coding in FORTRAN, the approach presented is applicable to other programming languages as well. A professional coder understands that the complexity of any code lies in its algorithm design, which remains fundamentally language-agnostic. This methodology is organized into eleven steps, each accompanied by explanatory comments marked with an exclamation sign (`!`). The initial comments divide the program into distinct segments. Throughout this article, we briefly outline the purpose of each step. It is recommended to replicate these explanatory comments, similar to guideposts, within the program and then proceed with coding step-by-step, following the provided guidelines. This approach has been tested by students at various stages of their programming projects. Its proven effectiveness inspired us to share it to help researchers improve coding speed and accuracy.

## Program Steps

### 1. Variables Definition:
At first glance, coding in MATLAB, which does not require variable definition, may seem more straightforward. However, the key advantage of low-level languages lies in their ability to define variables explicitly and use the `implicit none` command, effectively overriding FORTRAN’s default settings. This capability enables optimal definition of array dimensions, helping to prevent issues related to insufficient virtual memory, while also enhancing the overall speed of program execution.

```
!******************************************************************************************
!                                       Variables Definition
!******************************************************************************************
```

### 2. Initializing Variables with Zero:
Using uninitialized variables can lead to significant issues in your program. It is essential to initialize all variables before use; otherwise, the compiler may retrieve residual data from a previous program that occupied the same memory segment. This can result in unpredictable behavior and difficult-to-debug errors in your program.

```
!******************************************************************************************
!                                        Giving Zero to variables
!******************************************************************************************
```

### 3. Retrieving Input Parameters from the User:
For variables expected to change values between consecutive runs—especially when the program's goal is to assess variations in these parameters—it is important that they be defined by the user. This approach ensures that the program can be thoroughly tested and evaluated under different scenarios and conditions.

```
!******************************************************************************************
!                                        Inputs
!******************************************************************************************
```

A helpful approach is to input each variable in both numerical and character formats. This allows the numerical format to be used for calculations, while the character format verifies variable names. For example, `E` and `freq` are processed as real numbers, whereas `EE` and `freqf` represent the same values in character format. This distinction is recognized and managed by the compiler.

```
Real E , freq
Character*30 EE, freqf, filename1
write(*,'(/,2x,a,\)') ' Enter the Energy value : '
read(*,*)E
write(*,'(/,2x,a,\)') ' Again : '
read(*,*)EE
write(*,'(/,2x,a,\)') ' Enter the frequency value : '
read(*,*)freq
write(*,'(/,2x,a,\)') ' Again : '
read(*,*)freqf
```

### 4. Assigning File Names Based on Input or Problem Specifications:
This section highlights an important consideration: optimizing the use of output files. Effective utilization can, firstly, improve the processing speed of these files and, secondly, ensure that output information is systematically organized. To achieve this, it is essential to generate file names based on the input information. This straightforward practice can lead to significant improvements in both the efficiency and organization of output data.

```
!******************************************************************************************
!                                        Determination of Filenames and Opening Files
!******************************************************************************************
```

Additionally, when results are represented graphically, it is preferable to save the file with a `.plt` extension, compatible with `Tecplot` software. Tecplot is feature-rich and allows for easy graph generation with a single click on the file icon. 

The `//` operator is used for concatenating two character expressions, while the `trim()` function removes unused characters from a character field. For example, these techniques could be applied to a graph illustrating temperature variations over time at a specified pulse energy and frequency.

```
filenameTt = 'E_'//trim(EE)//'f_'//trim(freqf)//'T_t.plt'
open(1,file=filenameTt)
```

### 5. Entering Program Constants:
It is crucial to centralize all of the program’s values in this section. Additionally, organizing them in a logical and readable sequence is essential for easy management and control.
```
!******************************************************************************************
!                                        Constants
!******************************************************************************************
```

### 6. Optimizing Array Dimensions Based on Inputs and Constants:
Minimizing the number of entries in an array does not necessarily reduce program execution time or decrease memory errors. Therefore, optimizing array allocation is essential. In general, it is preferable to avoid using arrays for problem-solving unless absolutely necessary.

```
!******************************************************************************************
!                                        Arrays Allocation
!******************************************************************************************
```

### 7. Initializing All Array Elements to Zero:
Initializing arrays, just as with variables, is essential. A recommended approach is to use the `forall` command to set all elements of the array to zero. This ensures each element is properly initialized before any further operations, helping to prevent unexpected behaviors due to uninitialized elements.

```
!******************************************************************************************
!                                        Giving Zero to Arrays
!******************************************************************************************
```

### 8. Displaying Input and Constant Values:
It is crucial for the user to verify the accuracy of all input values before executing the main segment of the program. To facilitate easy monitoring, it is recommended to display each variable name alongside its corresponding value. This allows the user to efficiently review and confirm the accuracy of the input before proceeding with the main operation of the program.

```
!******************************************************************************************
!                                        Printing Constants and Inputs
!******************************************************************************************
```

### 9. The Main Block of the Program:
In this section, we precisely establish the boundary conditions and define the problem-solving algorithm, ensuring that both are presented in a clear and straightforward manner.

```
!******************************************************************************************
!                                        Main Block of the Program
!******************************************************************************************
```

### 10. Displaying Results:
Although results are typically printed in a separate section, there are instances where it is necessary to store data directly into files within the main segment of the program. Carefully choosing the location for output retrieval within the program can help prevent `insufficient virtual memory` errors in most cases and, at times, significantly increase execution speed—occasionally by an order of magnitude.

```
!******************************************************************************************
!                                        Printing Results
!******************************************************************************************
```

### 11.  Closing the Opened Files and Terminating the Program:
At first glance, it may seem that leaving files opened at the start of a program, poses no issues. However, experience shows that to prevent potential complications, especially from interference in subsequent runs, it is best to close files upon completion of their tasks.

```
!******************************************************************************************
!                                        Closing Files and Ending the Program
!******************************************************************************************
```

### Conclusion
Based on years of experience, we have developed a template for writing programs that enables users—from beginners to professionals—to write optimal code and avoid many common errors. While modular programming also aims to achieve these objectives, we believe this format is simpler and more effective.
