# Installing IFORT (Intel FORTRAN) on Linux

## Step 1
Open the terminal and **copy/paste/run** the following two commands to download and install the first part. Then, follow the visual instructions:

`wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/7deeaac4-f605-4bcf-a81b-ea7531577c61/l_BaseKit_p_2023.1.0.46401_offline.sh`

`sudo sh ./l_BaseKit_p_2023.1.0.46401_offline.sh`


![Image 1](../Archive/images/1.%20ifort_installation%20images/image01.png)<p>&nbsp;</p>
![Image 2](../Archive/images/1.%20ifort_installation%20images/image02.png)<p>&nbsp;</p>
![Image 3](../Archive/images/1.%20ifort_installation%20images/image03.png)<p>&nbsp;</p>
![Image 4](../Archive/images/1.%20ifort_installation%20images/image04.png)<p>&nbsp;</p>
![Image 5](../Archive/images/1.%20ifort_installation%20images/image05.png)<p>&nbsp;</p>
![Image 6](../Archive/images/1.%20ifort_installation%20images/image06.png)<p>&nbsp;</p>
![Image 7](../Archive/images/1.%20ifort_installation%20images/image07.png)<p>&nbsp;</p>


## Step 2
Open the terminal and **copy/paste/run** the following two commands to download and install the second part. Then, follow the visual instructions:

`wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/1ff1b38a-8218-4c53-9956-f0b264de35a4/l_HPCKit_p_2023.1.0.46346_offline.sh`

`sudo sh ./l_HPCKit_p_2023.1.0.46346_offline.sh`

![Image 8](../Archive/images/1.%20ifort_installation%20images/image08.png)<p>&nbsp;</p>
![Image 9](../Archive/images/1.%20ifort_installation%20images/image09.png)<p>&nbsp;</p>
![Image 10](../Archive/images/1.%20ifort_installation%20images/image10.png)<p>&nbsp;</p>
![Image 11](../Archive/images/1.%20ifort_installation%20images/image11.png)<p>&nbsp;</p>
![Image 12](../Archive/images/1.%20ifort_installation%20images/image12.png)<p>&nbsp;</p>
![Image 13](../Archive/images/1.%20ifort_installation%20images/image13.png)<p>&nbsp;</p>

## Step 3
Install ZSH
- Open WSL terminal

`sudo apt install git zsh -y`

`zsh --version`

`sh -c "$(wget https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"`

- Close and reopen the WSL terminal

## Step 4

`nano .zshrc`

Add the bellow path to .zshrc

`source "/opt/intel/oneapi/compiler/2023.1.0/env/vars.sh" intel64`

Like below:


![Image 14](../Archive/images/1.%20ifort_installation%20images/image14.png)<p>&nbsp;</p>

Close and reopen the terminal and make sure that ifort is installed by this command:

`ifort --version`

## Sources
- Download ifort from Intel website [Link](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html?operatingsystem=linux&distributions=offline&version=2023.1)
- Installation tutorial from Fortran4duck YT Channel [Link](https://www.youtube.com/watch?v=ULR9Esf3sDc)
