
TIDD

TIDD is a universal post-processing tool which supports confident peptide identifications regardless of the search engine adopted. TIDD can work for any (including newly developed) search engines because it calculates universalfeatures that assess peptide-spectrum match quality while it allows additional featuresprovided by search engines (or users) as well.

Here, we support two types of TIDD version -- a simple GUI based TIDD and command version of TIDD

======================

Code download from Github: https://github.com/HanyangBISLab/TIDD.git

Download and test R&java code from R studio cloud: https://rstudio.cloud/spaces/178915/project/2994889

Test GUI from shinyapps.io: https://honglan-li.shinyapps.io/project/

<*This online GUI is just for testing, it may take very long time (or crash) to run real datasets becasuse of limited RAM/disk size, which supported by free R cloud studio. Therefore, we recommand you to run TIDD on local system.> 

TIDD user manual files: https://docs.google.com/document/d/168oGySS15xrobeqTY54_QLMgWQofnPIApx6ma5ntGa8/edit?usp=sharing

=====================



1. GUI 

It is based on R and java programming. 

So it required:

- install java jdk (17 or above) and R(4.1.2 or above)
- for R, the following packages are need to install.

   > install.packages("shiny")
   > install.packages("shinythemes")
   > install.packages("shinyFiles")
   > install.packages("shinydashboard")
   > install.packages("ROCR")
   > install.packages("e1071")
   
- run shinyApp 
 
  > PSM results:put psm result files to data/PSM/
  > mgf files: make a new dir. under data/MGF/ 
               put mgf files to new dir.
               
 Detail:
 TIDD: R and Java based program 

1. install R version 4.1.2 (or above)

   link: https://cran.r-project.org/
   
   E.g. download R-4.1.2* ver 

2. install R studio 

  link: https://www.rstudio.com/products/rstudio/download/

  E.g. download Rstudio-2021.09.1*

3.  install java 8 

    download: https://www.oracle.com/java/technologies/downloads/
    installation guide: https://docs.oracle.com/javase/10/install/installation-jdk-and-jre-microsoft-windows-platforms.htm#JSJIG-GUID-29333CFD-E7A6-498B-9317-97700C81D928

    E.g. jdk-17.01*

4. download TIDD code from git: https://github.com/HonglanLi/TIDD.git

5. unzip the TIDD*.zip 

6. open the R studio 
 
7. import TIDD project to R studio

   File --> Open project --> go to right position to select TIDD R project 

                                       for example: "~/Download/TIDD-main/TIDD/TIDD.rproj

8. install R packages

   required packages: shiny, shinythemes, shinydashboard, shinyFiles, e1071, ROCR

   E.g. install shiny package 

   console: install.packages("shiny")

   UI: click Packages --> install --> type "shiny" --> install 

9 Run shinyAPP

   open one of *R file in the TIDD rproject,  and click "Run 

2. commmand  

- Extract TIDD features 
 
   1. java -jar TIDD.jar <options> < attributes> 

   You can run "java -jar TIDD.jar -h" to check the options and their corresponding values

   (If you need to calculate XCorr, run "java -jar Xcorr_v2.jar " additionally)

   2. java -jar Xcorr_v2.jar <options> < attributes> 
   
   *detail: check usage (run java -jar Xcorr_v2.jar)
   
   3. run $Rscript TIDD_ModelFit.R
   
   * Change parameters inside the file: from line 30 to 41
   
   
