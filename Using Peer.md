**Statistical Considerations**

PEER is a powerful tool for finding hidden factors in gene expression data, but it has several limitations that must be understood if it is to be used properly.

* **Linear relationships** - PEER uses linear modeling for both the covariates you supply it with, as well as the hidden factors it finds.  If you have reason to believe that some of your covariates have non-linear relationships with gene expression, these effects would need to be removed before using PEER. If one of these covariates is one you wish to preserve for later analysis, PEER may not be the best choice for finding hidden factors.
* **Missing values** - Because PEER uses linear modeling, it does not support missing values in covariates under any circumstances.  Imputation methods may be used to estimate the missing values if the only a few are missing.  Otherwise samples missing these values must be dropped, or the covariate cannot be used.  If a covariate has many missing values, it may be possible to correlate it with the PEER factors to see if the effects of that covariate were captured in one of the factors.
* **Principal components** - The hidden factors are intially estimated using principal components analysis.  Because this is an abstract approximation, a hidden factor found by PEER may actually correspond to several underlying factors. This inherently less accurate than using known factors when they are available.  PEER factors should not be used as substitutes for known covariates; instead they should be used to complement these other covariates.  

**Data Preparation**

PEER should generally be used after more routine data processing steps have been completed.  This would include normalization of expression values, normalization between arrays, quality control selection of probes, removal of technical replicates and statistical outliers, and batch effect correction.  The gene expression data must be in the form of an <code>m x n</code> expression matrix with <code>n</code> genes and <code>m</code> samples.  The covariates should be in the form of an <code>m x c</code> matrix with <code>m</code> samples and <code>c</code> covariates.  Continuous covariates are supplied in a single column, while categorical covariates with <code>n</code> possible values must be encoded in <code>n - 1</code> dummy variables (like those used in regression models)

**Running PEER**

The code in this repository runs PEER inside a function.  This makes it easier to try different numbers of factors and covariates easily.  The full code for the function is provided below, but the individual steps inside the function will be explained first.

Instantiate the PEER object:

<code>model = PEER()</code>

Set the number of factors:

<code>PEER_setNk(model, num.factors)</code>

The next line sets the expression values.  Note that the matrix is transposed, so you should supply the function with the expression data as an <code>n x m</code> matrix.  The reason the matrix needs to be transposed is because other components of the pipeline such as limma and ComBat use an <code>n x m</code> matrix:

<code>PEER_setPhenoMean(model, as.matrix(t(intensities)))</code>

Enable the use of the mean expression values in the linear model:

<code>PEER_setAdd_mean(model, TRUE)</code>

If you have already corrected your expression data for all other covariates, you do not need to supply a known covariate matrix.  Otherwise, the covariates matrix is supplied:

```
if (use.covariates = TRUE)
{
    PEER_setCovariates(model, as.matrix(covariates))
}
```

The next line sets the maximum number of iterations.  Normally the model should converge within about 100 iterations, so a maximum of 1000 is far higher than should be needed:

<code>PEER_setNmax_iterations(model, 1000)</code>

Run the model with the given settings:

<code>PEER_update(model)</code>

If you wish to use the residuals which remain after fitting the known covariates and hidden factors, they can be extracted:

```
residuals.PEER = t(PEER_getResiduals(model))
rownames(residuals.PEER) = rownames(intensities)
colnames(residuals.PEER) = colnames(intensities)
```

There are several other outputs of PEER which are of importance.  The factor values from PEER can be directly supplied to limma as covariates for differential expression analysis.  The weights can be used to assess the relative importance of each factor.  The precision values for the factors can also be used to assess the importance of factors.  All of these outputs can be saved as tables:

```
write.csv(data.frame(residuals.PEER), file = paste("residuals_", num.factors, sep = "", ".csv"), row.names = FALSE)
write.csv(data.frame(PEER_getX(model)), file = paste("factor_", num.factors, sep = "", ".csv"), row.names = FALSE)
write.csv(data.frame(PEER_getW(model)), file = paste("weight_", num.factors, sep = "", ".csv"), row.names = FALSE)
write.csv(data.frame(PEER_getAlpha(model)), file = paste("precision_", num.factors, sep = "", ".csv"), row.names = FALSE)
```

Finally, several plots can be made of PEER output.  However, it may be necessary to remove the known covariates, which is not done in this code.  The <code>CairoPDF</code> driver from the <code>Cairo</code> package offers much better graphics output on most systems that the default <code>pdf</code> driver.  This code checks if the <code>Cairo</code> package is available before attempting to use it:

```
if (require(Cairo))
{
    CairoPDF(file = paste("model", num.factors, ".pdf", sep = ""), width = 10, height = 10)
    PEER_plotModel(model)
    dev.off()

    CairoPDF(file = paste("precision_", num.factors, ".pdf", sep = ""), width = 10, height = 10)
    plot(PEER_getAlpha(model), col = "red", lwd = 4, main = paste("precision", num.factors, "factor", sep = " "))
    dev.off()
}
else
{
    pdf(file = paste("model", num.factors, ".pdf", sep = ""), width = 10, height = 10)
    PEER_plotModel(model)
    dev.off()

    pdf(file = paste("precision_", num.factors, ".pdf", sep = ""), width = 10, height = 10)
    plot(PEER_getAlpha(model), col = "red", lwd = 4, main = paste("precision", num.factors, "factor", sep = " "))
    dev.off()
}    
```

The final function combines all of these steps together:

```
gen.peer <- function(num.factors, intensities, use.covariates, covariates)
{
    model = PEER() 
    PEER_setNk(model, num.factors) 
    PEER_setPhenoMean(model, as.matrix(t(intensities))) 
    PEER_setAdd_mean(model, TRUE) 
    if (use.covariates == TRUE)
    {
        PEER_setCovariates(model, as.matrix(covariates)) 
    }
    PEER_setNmax_iterations(model, 1000) 
    PEER_update(model) 
    residuals.PEER = t(PEER_getResiduals(model))
    rownames(residuals.PEER) = rownames(intensities)
    colnames(residuals.PEER) = colnames(intensities)

    write.csv(data.frame(residuals.PEER), file = paste("residuals_", num.factors, sep = "", ".csv"), row.names = FALSE)
    write.csv(data.frame(PEER_getX(model)), file = paste("factor_", num.factors, sep = "", ".csv"), row.names = FALSE)
    write.csv(data.frame(PEER_getW(model)), file = paste("weight_", num.factors, sep = "", ".csv"), row.names = FALSE)
    write.csv(data.frame(PEER_getAlpha(model)), file = paste("precision_", num.factors, sep = "", ".csv"), row.names = FALSE)

    if (require(Cairo))
    {
        CairoPDF(file = paste("model", num.factors, ".pdf", sep = ""), width = 10, height = 10)
        PEER_plotModel(model)
        dev.off()

        CairoPDF(file = paste("precision_", num.factors, ".pdf", sep = ""), width = 10, height = 10)
        plot(PEER_getAlpha(model), col = "red", lwd = 4, main = paste("precision", num.factors, "factor", sep = " "))
        dev.off()
    }
    else
    {
        pdf(file = paste("model", num.factors, ".pdf", sep = ""), width = 10, height = 10)
        PEER_plotModel(model)
        dev.off()

        pdf(file = paste("precision_", num.factors, ".pdf", sep = ""), width = 10, height = 10)
        plot(PEER_getAlpha(model), col = "red", lwd = 4, main = paste("precision", num.factors, "factor", sep = " "))
        dev.off()
    }    
}
```

**Using PEER in differential expression**

Once you have computed the desired number of factors, you can use these factors as additional covariates in a linear model fitted with limma. The PEER factors can simply be load from the .csv file and the known covariates removed.  Note that this uses the <code>%>%</code> from <code>magrittr</code> and the <code>select</code> function from <code>dplyr</code>:

<code>model.PEER_covariate <- read_csv("./factor_8.csv") %>% select(-(X1:X6))</code>

The factors can then be joined with the covariates matrix:

<code>model.full <- cbind(model.cov, model.PEER_covariate)</code>

This matrix can then be supplied to the lmFit function as the covariate matrix.

**Correcting gene expression data for effects of PEER factors**

Methods such as network analysis and time series analysis require that the effects of these factors be removed.  This is most easily done by fitting the factors to a linear model and getting the residuals:
