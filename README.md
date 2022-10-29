This is the repo for VCF files plot

# Install plotVCF
Installing `plotVCF` is as simple as:
```
if (!require("devtools")) install.packages("devtools")
if (!require("BiocManager")) install.packages("BiocManager")
remotes::install_github(
    "cccnrc/plot-VCF",
    repos = BiocManager::repositories()
)
```
If you know a bit of R code (no worries, you don't really need to :wink:) you noticed that it only requires [devtools](https://devtools.r-lib.org/) and [BiocManager](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html).

The rest of dependencies you need will be installed directly with `plotVCF`: that's why installation will probably take a while (but you need to perform it only once :wink:)

# Run plotVCF
To use `plotVCF` you only need to point it to a VCF and a FASTA file. Then just call the `createVCFplot()` funtion on them and you will get your Manhattan VCF plot out!
```
library(plotVCF)

VCF <- <path-to-your-VCF-file>
FASTA <- <path-to-your-FASTA-file>

createVCFplot( VCF, FASTA )
```
There are multiple options you can specify, just have a look at them with:
```
?createVCFplot
```

# plotVCF usage
The default behavior of `plotVCF()` is to allow visualization of variants position. It will create random Y-values for variants just to allow their visualization:
```
plotVCF( VCF, FASTA )
```
![plotVCF() basic plot](plots/plotVCF.base.png)


# Save plotVCF
Once you created your plot, you can save it with any R graphic function ([png()](https://cran.r-project.org/web/packages/png/index.html),[pdf()](https://www.rdocumentation.org/packages/grDevices/versions/3.6.2/topics/pdf),[tiff()](https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/png.html), etc.).

Something I like:
```
VCF_PLOT <- createVCFplot( VCF, FASTA )

png( <path-to-your-PNG-output>,  width = 5000, height = 2500, res = 300 )
VCF_PLOT
whatever <- dev.off()
```
but fell free to use what you wish!
