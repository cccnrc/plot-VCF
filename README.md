This is the repo for VCF files plot

Main R script: [plot-VCF.R](plot-VCF.R)

# Install plotVCF
Installing `plotVCF` is as simple as:
```
if (!require("devtools")) install.packages("devtools")
remotes::install_github(
    "cccnrc/plot-VCF",
    repos = BiocManager::repositories()
)
```
If you know a bit of R code (no worries, you don't really need to :wink:) you noticed that it only requires [devtools](https://devtools.r-lib.org/).
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

# Save plotVCF
Once you created your plot, you can save it as you wish. Something I like:
```
VCF_PLOT <- createVCFplot( VCF, FASTA )

png( <path-to-your-PNG-output>,  width = 5000, height = 2500, res = 300 )
VCF_PLOT
whatever <- dev.off()
```
but fell free to use what you wish!
