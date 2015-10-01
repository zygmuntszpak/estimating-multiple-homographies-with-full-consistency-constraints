Sampson distance based joint estimation of multiple homographies with 
uncalibrated cameras Zygmunt L Szpak (c) 2014 

Version 0.1, 8 October 2014 

Please cite the following paper if you use this software in your own 
research. 

Z. L. Szpak, W. Chojnacki, A. Eriksson, and A. van den Hengel. Sampson 
distance based joint estimation of multiple homographies with 
uncalibrated cameras. Comput. Vis. Image Underst., 125:200-213, 2014. 

DISCLAIMER OF WARRANTY 

This source code is provided "as is" and without warranties as to 
performance or merchantability. The author and/or distributors of this 
source code may have made statements about this source code. Any such 
statements do not constitute warranties and shall not be relied on by 
the user in deciding whether to use this source code. 

This source code is provided without any express or implied warranties 
whatsoever. Because of the diversity of conditions and hardware under 
which this source code may be used, no warranty of fitness for a 
particular purpose is offered. The user is advised to test the source 
code thoroughly before relying on it. The user must assume the entire 
risk of using the source code. 

HOW TO USE 

Unzip the folder into a directory of your choice and add the folder and 
all its sub-folders to the matlab path. Then execute one of the runme_* 
scripts. 

This will generate several planar synthetic scenes (you can specify the 
number of scenes, data points, planes and noise levels) and estimate the 
corresponding homography matrices. The quality of the estimated 
homographies is evaluated on ground truth data and the program will 
produce an output of the root mean square error for each homography, as 
well as the overall mean root mean square error (averaged over all 
homographies and all scenes). Diagnostic information for each method 
(such as running time, iterations etc.) will also be computed. 

KNOWN ISSUES 

* None 

FUTURE RELEASE 

There are several functions that I have not included in the current 
release, because I am still in the process of cleaning the code. In 
particular, I have a standalone implementation of the 
Levenberg-Marquardt algorithm for optimising the Sampson distance, so 
that porting the code to a different language that doesn't have MATLAB's 
built in nonlinear least squares optimisation will not be a problem. 

Also, I have example scripts where outliers are added to the data and 
where RANSAC is used to estimate the initial homographies. I am also in 
the process of cleaning that code. 

In general I am in the process of compiling more example scripts 
(including tests on real data) and placing more comments throughout the 
code. 

There are several instances of duplicate code (in particular pertaining 
to normalising all data points so that they lie inside a unit box) that 
strictly speaking should be refractored and removed. 

If you require any assistance, have any suggestions or have questions, 
please feel free to contact me: zygmunt.szpak@adelaide.edu.au 



 


