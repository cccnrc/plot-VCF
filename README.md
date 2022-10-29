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
## FASTA file
To use `plotVCF()` you need a [FASTA file](https://en.wikipedia.org/wiki/FASTA_format) to let the software know chromosome boundaries etc.

You need to use the same FASTA format on which your VCF was aligned to. If you are not sure what are we talking about you should find this information in your [VCF header](https://samtools.github.io/hts-specs/VCFv4.2.pdf). In the examples below I used a [hg38](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/) aligned VCF. Thus, if you wish to reproduce the plots shown below [download a hg38 FASTA](https://www.gungorbudak.com/blog/2018/05/16/how-to-download-hg38-grch38-fasta-human-reference-genome/)!

<span style="color:red">*Important*</span>: make sure the chromosome names in your FASTA file exactly match chromosome names in your VCF file! (e.g. you can find VCF with `chr1` name format while FASTA may have `1` name format or reversal. Those <ins>VCF and FASTA chromosome names must match!</ins>)

# plotVCF usage
Let's have a look on how many different things you can achieve with this package!

All plots in examples below are generated from the same VCF file, that you can find (compressed) [here](inst/extdata/exampleVCF.vcf.gz).
<span style="color:red">*Important*</span>: you need to decompress it first!
```
gzip -d inst/extdata/exampleVCF.vcf.gz
```
## basic visualization plot
The default behavior of `plotVCF()` is to allow visualization of variants position. It will create random Y-values for variants just to allow their visualization:
```
VCF <- inst/extdata/exampleVCF.vcf.gz
FASTA <- <path-to-your-FASTA-file>

createVCFplot( VCF, FASTA )
```
![plotVCF() basic plot](plots/plotVCF.base.png)
## basic visualization plot - ordered
You can also choose to order your variants based on sample representation, this will allow you to graphically represent variant number differences across different chromosomes:
```
createVCFplot( VCF, FASTA, ORDERED=TRUE )
```
![plotVCF() ordered plot](plots/plotVCF.ordered.png)
## user-defined visualization plot
You can specify any numerical VCF variant flag (make sure you have that flag in your VCF) to use as Y-axis in your plot. This creates a sort of Manhattan VCF plot, based on that value.

In this example, I used `QUAL` flag of my variants as Y-axis coordinate:
```
createVCFplot( VCF, FASTA, VAR_FLAG="QUAL" )
```
![plotVCF() var-flag plot](plots/plotVCF.flag.png)
## user-defined visualization plot - threshold
You can also use a threshold that will be plotted on your Y-axis defined flag (such as you do with significance in a Manhattan plot):
```
createVCFplot( VCF, FASTA, VAR_FLAG="QUAL", THRESHOLD=22000 )
```
![plotVCF() var-flag threshold plot](plots/plotVCF.flag-threshold.png)
## focused plot
If you wish to graphically analyze just some chromosomes, you can use `CHR_NAMES` flag:
```
createVCFplot( VCF, FASTA, VAR_FLAG="QUAL", CHR_NAMES=c("chr3","chr4","chr5","chr13","chr19","chrX") )
```
![plotVCF() var-flag focused plot](plots/plotVCF.flag-focus.png)
## advanced plot
There are many more combinations you can create with this package, explore all possible values with:
```
?createVCFplot
```
and create your own plots!

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
