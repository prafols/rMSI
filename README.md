# rMSI
rMSI is an R package for mass spectrometry (MS) imaging data handling and visualization.
The package is a multi-platform tool that has been tested on Linux, MAC and Windows systems. 
It provides an optimized data model to allow loading large MS imaging datasets in low resource computers. MS data is stored in the hard disk drive (HDD) but rMSI is able to access data as if it where kept in computer’s memory using a virtual memory technology. The package also provides a graphical user interface (GUI) to facilitate MS imaging data exploration in R platform. The main rMSI GUI allows representing up to three MS ions spacial distribution, direct access to pixel spectrum and other usefull features. See and screenshot below.

![alt text](https://github.com/prafols/rMSI/blob/master/images/screenShotrMSI_RGB.png = 600x "rMSI Main GUI")

### Installation
##### Rtools (Windows only)
*rMSI* uses some compression methods that are not available by default in Windows operating system. To be able to run *rMSI*, Windows users must install Rtools available at: <https://cran.r-project.org/bin/windows/Rtools/> just download the *.exe installer and go trough the installation wizard.

##### RGtk2
rMSI provides a quite complex data model together with a graphical user interface (GUI), consequently rMSI depends on many other R packages that must be also installed. RGtk2 is one of these packages and is known to be problematic to install in some cases (specially on non-Linux systems). So, it is highly recommended installing and testing it before going through the process of installing *rMSI*. 
First install RGtk2 using R console:
```R
> install.packages("RGtk2")
```
Then test that it is working by loading the package and executing the demo in R:
```R
> library(RGtk2)
> demo(appWindow)
```
If it appears a Gtk Windows then continue with the *rMSI* installation process. Otherwise, please check out the RGtk2 website for solving issues related with Gtk installation: <http://http://www.ggobi.org/rgtk2/>

##### rMSI
The simplest way to install rMSI and keep it updated is using devtools package. Install devtools from CRAN into your R session:
```R
>  install.packages("devtools")
```
Then simply tell devtools to install rMSI from github latest release:
```R
> devtools::install_github("prafols/rMSI", ref = "0.3")
```
This will install rMSI package and all of its dependencies in your R environment. Then you can access to its functions by loading the rMSI package or through the `::` operator. For example, you can test the main rMSI GUI by executing:
```R
> foo <- rMSI::OpenMSI()
```
This will open a dialog that allows loading up to two MS images in imzML format or rMSI .tar format. Then the MS images will be displayed. A reference to each MS image will be returned and stored in *foo* variable as *foo$img1* and *foo$img2*. 
