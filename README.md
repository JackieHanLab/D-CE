# D-CE
D-CE (Developmental Coalescent Embedding): An effective landmark free and model free de novo 3D reconstruction method for single cell analysis.

# Introduction
De novo reconstruction of single cell 3D spatial tissue localization is hitherto landmark based so far, and de novo spatial reconstruction is a compelling computational open problem.Theoretically, cells adjacent in space have similar gene expression patterns, for example, the anterior and posterior determination during embryonic development.We rely merely on the spatial information encoded in the gene expression patterns, and we find that D-CE of cell-cell association DST-transcriptomic networks can reliably reconstructs the single cell samples 3D spatial tissue distribution. 

If you find D-CE useful, please give us a star at github and cite our paper.


# Folders description
- **visualization.m**
Matlab codes to visualize the D-CE reconstructed structure of E7.5 Geo-seq;  
- **expE7.5allsample.txt**
An example data of Geo-seq gene expression matrix, with genes on the columns and samples on the rows;  
- **Info.mat**
Annotation information of the example data, including layer, germ layer of E7.5 and color for visualization;  
- **D-CE_windows**
Windows version of the main function of spatial reconstruction, includes the codes for expression matrix normalization, pair-wise distance calculation and dimensionality reduction based on coalescent embedding and integrated into one function;  
&emsp;-**spatial_reconstruct.exe**    
&emsp;-**readme file**   

- **D-CE_linux**
Linux version of the main function of spatial reconstruction, includes the codes for expression matrix normalization, pair-wise distance calculation and dimensionality reduction based on coalescent embedding and integrated into one function;  
&emsp;-**spatial_reconstruct**  
&emsp;-**spatial_reconstruct.sh** (shell script for temporarily setting environment variables and executing the application)  
&emsp;-**readme file**  

# Requirements
**MATLAB(>=R2017b)**.

**MATLAB Runtime** please verify the MATLAB Runtime is installed and ensure you have installed version 9.2 (R2017a) for windows or version 9.3 (R2017b) for linux. If not, do as follow:  
&emsp;(1) enter 'mcrinstaller' at MATLAB prompt. The MCRINSTALLER command displays the location of the MATLAB Runtime installer  
&emsp;(2) run the MATLAB Runtime installer. Or download right version of the MATLAB Runtime from the MathWorks Web site  
&emsp;http://www.mathworks.com/products/compiler/mcr/index.html

**MatlabBGL library** The support functions of the MatlabBGL library can be downloaded at:
http://mathworks.com/matlabcentral/fileexchange/10922-matlabbgl

**PROPACK library** The support functions of the PROPACK library can be downloaded at:
https://github.com/mavenlin/PropackMatlab4Windows
https://github.com/epfl-lts2/unlocbox/tree/master/test_bench/private

# Usage
## For Windows User
To run the demo,first you need to paste the **expE7.5allsample.txt** into the same folder of **spatial_reconstruct.exe**. Next, just double click the **spatial_reconstruct.exe**, then type 'expE7.5allsample.txt', which claim the path of the expression data and 'Y' or 'N', which means you will use or don't use CSI matrix in the 3D reconstruction. Wait about 10 secounds for the output file '3Dcoordinates.txt',which is the reconstructed 3D coordinates of the example sample. Finally, run **visualization.m** in MATLAB for visualization.

## For Linux User
```
git clone https://github.com/JackieHanLab/D-CE
cd D-CE/D-CE_linux
chomd a+x spatial_reconstruct.sh
chomd a+x spatial_reconstruct
./spatial_reconstruct.sh <mcr_directory>
```
at Linux or Mac command prompt. <mcr_directory> is the directory where version 9.3 of the MATLAB Runtime is installed or the directory where MATLAB is installed on the machine.
for example, If you have MATLAB Runtime installed in */mathworks/home/application/v93*, run the shell script as:  
```  
./spatial_reconstruct.sh /mathworks/home/application/v93  
```  
If you have MATLAB installed in */mathworks/devel/application/matlab*, run the shell script as:  

```  
./spatial_reconstruct.sh /mathworks/devel/application/matlab
```  
Then, similar to Windows version, just type the path of matrix and 'Y' or 'N' to run the 3D reconstruction.  

# Reference
Muscoloni, A., Thomas, J. M., Ciucci, S., Bianconi, G. & Cannistraci, C. V. Machine learning meets complex networks via coalescent embedding in the hyperbolic space. Nat Commun 8, 1615, doi:10.1038/s41467-017-01825-5 (2017).  

Peng, G. et al. Molecular architecture of lineage allocation and tissue organization in early mouse embryo. Nature 572, 528-532, doi:10.1038/s41586-019-1469-8 (2019).  

Nitzan, M., Karaiskos, N., Friedman, N. & Rajewsky, N. Gene expression cartography. Nature 576, 132-137, doi:10.1038/s41586-019-1773-3 (2019).  


# Contact  
For any problems, please contact:  
Yuxuan Zhao: zhaoyuxuan2017@pku.edu.cn  
Jing-Dong J. Han: jackie.han@pku.edu.cn  
