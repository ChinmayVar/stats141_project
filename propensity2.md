Propensity Scores Method for Treatment Pairings
================
Andre Hiwatig
2025-05-26

# Libraries

``` r
library(MatchIt)
library(tableone)
library(ggplot2)
library(cobalt)
```

    ##  cobalt (Version 4.6.0, Build Date: 2025-04-15)

    ## 
    ## Attaching package: 'cobalt'

    ## The following object is masked from 'package:MatchIt':
    ## 
    ##     lalonde

# Data Cleaning

Some changes made to the intial dataset were as follows:  
- changed response variable to binary (0/1) for matching function.  
- imputed gender by mode (didn’t end up using predictor).  
- imputed tissue resected (mm) by mean (didn’t end up using
predictor).  
- for patients with \> 1 surgery, remove the non-tarsectomy to avoid
being matched to themself.  
- changed character columns to factor.

``` r
df <- read.csv("new_updated.csv") # from github

# simplify colnames
colnames(df) <- gsub("[[:punct:]]", "", colnames(df))

# change resp var to binary for matching 
df['TarsusresectedYesNo'] <- sapply(df[['TarsusresectedYesNo']], function(x) {ifelse(x == "Y", 1, 0)})

# impute gender by mode
get_mode <- function(x) {
  ux <- na.omit(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}
mode_value <- get_mode(df[['GenderMale1Female2']])
df$GenderMale1Female2[is.na(df$GenderMale1Female2)] <- mode_value

# impute tissueresected by mean
df$TissueResectedmm[is.na(df$TissueResectedmm)] <- mean(df$TissueResectedmm, na.rm = TRUE)

# for patients with surgery in two eyes, remove the non-tarsectomy
df <- df[!(df$X %in% c(156, 157, 158)), ]

# datatypes to factor
df$GenderMale1Female2 <- as.factor(df$GenderMale1Female2)
df$ODOS <- as.factor(df$ODOS)

head(df)
```

    ##   X Dataset      Patientname  X                    SurgerytypeonSXdate
    ## 1 0     All  Aaron Greenberg  1 Ou upper lid ptosis posterior approach
    ## 2 1     All   Abby Hellwarth  3    OU Levator repair internal approach
    ## 3 2     All    Akram Hannani  6                    OU upper lid ptosis
    ## 4 3     All   Arleen Richman 10                    OU upper lid ptosis
    ## 5 4     All Barbara Einziger 18           OU upper lid ptosis internal
    ## 6 5     All       Barry Katz 21          OU upper lid ptosis posterior
    ##   MMCR0ELR1 GenderMale1Female2      Age EyeOD0OS1 Surgicaleyeyes1no0
    ## 1         0                  1 71.67397         0                  1
    ## 2         0                  2 76.71507         0                  1
    ## 3         0                  2 71.26027         1                  1
    ## 4         0                  2 65.44384         1                  1
    ## 5         0                  2 79.16438         1                  1
    ## 6         0                  1 76.34795         0                  1
    ##   TissueResectedmm TarsusresectedYesNo
    ## 1              8.0                   1
    ## 2              9.0                   0
    ## 3              5.5                   0
    ## 4              6.5                   0
    ## 5              7.5                   1
    ## 6              8.0                   1
    ##   Tarsusresected0fornotarsusformmoftarsusremoved Lengthoffollowupmo PreopMRD1
    ## 1                                            1.0           2.966667     0.000
    ## 2                                            0.0           2.600000     1.079
    ## 3                                            0.0          13.133333     0.690
    ## 4                                            0.0           3.000000     0.820
    ## 5                                            1.5           6.333333     0.000
    ## 6                                            1.0           4.233333     1.840
    ##   MostrecentPostOpMRD1 ChangeinMRDpostpre Unilateral X1eyetarsus1eyenormal
    ## 1                 2.92              2.920                              yes
    ## 2                 2.84              1.761                                 
    ## 3                 1.25              0.560                                 
    ## 4                 1.05              0.230                                 
    ## 5                 1.57              1.570                              yes
    ## 6                 3.18              1.340                                 
    ##   Botheyestarsus Othersurgeriesatthesametime ODOS SxDate
    ## 1             NA                          NA   OD       
    ## 2             NA                          NA   OD       
    ## 3             NA                          NA   OS       
    ## 4             NA                          NA   OS       
    ## 5             NA                          NA   OS       
    ## 6             NA                          NA   OD

# Propensity Matching (Ratio 1:1)

The following method calculates the probability (propensity) one would
receive a tarsectomy based on the predictors PreopMRD1, Age, and ODOS.
Then observations with similar propensity scores were paired with each
other.

``` r
# create a matched model with complete variables
match_model <- matchit(TarsusresectedYesNo ~ PreopMRD1 + Age + ODOS, 
                       data = df,
                       method = "optimal",     
                       ratio = 1)              
matched_data <- match.data(match_model)

# create a merged treatment control dataset for tests
treated <- matched_data[matched_data$TarsusresectedYesNo == 1, ]
control <- matched_data[matched_data$TarsusresectedYesNo == 0, ]
paired_data <- merge(treated, control, by = "subclass", suffixes = c("_treated", "_control"))
```

## Covariate Balance

To ensure that the treatment and control groups were even, we assessed
the distribution of covariates before and after matching. The solid line
represents the treatment group and the gray line is the control group.
As we can see, this matching method gives us roughly the same
distribution of covariates across the three predictors.

The plot of Covariate Balance shows Standardized Mean Differences before
and after matching. At a threshold of 0.2, all covariates are moderately
balanced. Pre-Operation MRD1 had an SMD of -0.170, but it’s difficult to
balance because we saw in EDA that tarsectomies tend to occur in
patients with lower Pre-Operation MRD1. Therefore, there aren’t that
many comparable patients in the control group for those with lower
Pre-Operation MRD1.

``` r
# covariate balance plots
plot(match_model, type = "density", interactive = FALSE,
     which.xs = ~ PreopMRD1 + Age + ODOS )
```

![](propensity2_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
love.plot(match_model, binary = "std", threshold = 0.2)
```

![](propensity2_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
summary(match_model)
```

    ## 
    ## Call:
    ## matchit(formula = TarsusresectedYesNo ~ PreopMRD1 + Age + ODOS, 
    ##     data = df, method = "optimal", ratio = 1)
    ## 
    ## Summary of Balance for All Data:
    ##           Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
    ## distance         0.3154        0.1369          0.8250     3.7507    0.2626
    ## PreopMRD1        0.4692        1.5722         -0.7223     2.0948    0.2139
    ## Age             69.6962       60.3601          0.5977     1.1607    0.2520
    ## ODOSOD           0.5385        0.4462          0.1852          .    0.0923
    ## ODOSOS           0.4615        0.5538         -0.1852          .    0.0923
    ##           eCDF Max
    ## distance    0.4538
    ## PreopMRD1   0.4462
    ## Age         0.4692
    ## ODOSOD      0.0923
    ## ODOSOS      0.0923
    ## 
    ## Summary of Balance for Matched Data:
    ##           Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
    ## distance         0.3154        0.2625          0.2445     2.0786    0.0218
    ## PreopMRD1        0.4692        0.7295         -0.1705     2.1850    0.0855
    ## Age             69.6962       69.0677          0.0402     1.4641    0.1141
    ## ODOSOD           0.5385        0.5000          0.0772          .    0.0385
    ## ODOSOS           0.4615        0.5000         -0.0772          .    0.0385
    ##           eCDF Max Std. Pair Dist.
    ## distance    0.2308          0.2465
    ## PreopMRD1   0.3077          0.5840
    ## Age         0.2692          0.6246
    ## ODOSOD      0.0385          1.4659
    ## ODOSOS      0.0385          1.4659
    ## 
    ## Sample Sizes:
    ##           Control Treated
    ## All           130      26
    ## Matched        26      26
    ## Unmatched     104       0
    ## Discarded       0       0

## Control/Treatment Comparisons

1)  Boxplot: Distributions of Postop MRD1 show a clear positive effect
    of tarsectomy on Post-Op MRD1. Range is similar between control and
    treatment groups.
2)  We only see 7/26 pairs where tarsectomy decreases postop mrd1.

``` r
# boxplot for postopmrd1
boxplot(MostrecentPostOpMRD1 ~ TarsusresectedYesNo, data = matched_data,
        main = "Post-op MRD1 by Tarsectomy with Propensity Score Matching",
        xlab = "Tarsectomy (0 = No, 1 = Yes)",
        ylab = "Post-op MRD1 (mm)",
        col = c("skyblue", "salmon"),
        border = "gray40")
```

![](propensity2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# plot of postop mrd differences (1)
ggplot(matched_data, aes(x = factor(subclass), y = MostrecentPostOpMRD1, color = factor(TarsusresectedYesNo), group = subclass)) +
  geom_point(size = 2) +
  geom_line() +
  labs(x = "Matched Pair", y = "Post-Op MRD1", color = "Group") +
  scale_color_manual(values = c("blue", "red"), labels = c("Control", "Treatment")) +
  theme_minimal() +
  theme(axis.text.x = element_blank())
```

![](propensity2_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
# plot of postop mrd differences (2)
paired_data$postopmrd1_change <- paired_data$MostrecentPostOpMRD1_treated - paired_data$MostrecentPostOpMRD1_control
mean_postopchange = mean(paired_data$postopmrd1_change)
paired_data$diff_color <- ifelse(paired_data$postopmrd1_change > 0, "Increased Post-Op MRD1", "Decreased Post-Op MRD1")
ggplot(paired_data, aes(x =factor(subclass), y = postopmrd1_change, color = diff_color), group = subclass) +
  geom_point(size = 2) +
  labs(x = "Matched Pairs", y = "Difference in Post-Op MRD1", color = "Group") +
  scale_color_manual(values = c("Increased Post-Op MRD1" = "green", "Decreased Post-Op MRD1" = "red")) +
  theme_minimal() 
```

![](propensity2_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

## Parametric Tests

Both paired and unpaired t-tests show that tarsectomy has a signficant
positive effect.

``` r
# Unpaired T-Test
t.test(paired_data$MostrecentPostOpMRD1_treated, paired_data$MostrecentPostOpMRD1_control)
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  paired_data$MostrecentPostOpMRD1_treated and paired_data$MostrecentPostOpMRD1_control
    ## t = 2.0175, df = 49.802, p-value = 0.04905
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  0.002726144 1.253043087
    ## sample estimates:
    ## mean of x mean of y 
    ##  2.880731  2.252846

``` r
# Paired T-Test
t.test(paired_data$MostrecentPostOpMRD1_treated, paired_data$MostrecentPostOpMRD1_control, paired = TRUE)
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  paired_data$MostrecentPostOpMRD1_treated and paired_data$MostrecentPostOpMRD1_control
    ## t = 2.1514, df = 25, p-value = 0.0413
    ## alternative hypothesis: true mean difference is not equal to 0
    ## 95 percent confidence interval:
    ##  0.02679895 1.22897028
    ## sample estimates:
    ## mean difference 
    ##       0.6278846

## Nonparametric Tests

Both paired and unpaired t-tests show that tarsectomy has a signficant
positive effect.

``` r
# Unpaired
wilcox.test(paired_data$MostrecentPostOpMRD1_treated, paired_data$MostrecentPostOpMRD1_control)
```

    ## Warning in wilcox.test.default(paired_data$MostrecentPostOpMRD1_treated, :
    ## cannot compute exact p-value with ties

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  paired_data$MostrecentPostOpMRD1_treated and paired_data$MostrecentPostOpMRD1_control
    ## W = 453.5, p-value = 0.03532
    ## alternative hypothesis: true location shift is not equal to 0

``` r
# Paired
wilcox.test(paired_data$MostrecentPostOpMRD1_treated, paired_data$MostrecentPostOpMRD1_control, paired = TRUE)
```

    ## Warning in wilcox.test.default(paired_data$MostrecentPostOpMRD1_treated, :
    ## cannot compute exact p-value with ties

    ## 
    ##  Wilcoxon signed rank test with continuity correction
    ## 
    ## data:  paired_data$MostrecentPostOpMRD1_treated and paired_data$MostrecentPostOpMRD1_control
    ## V = 264, p-value = 0.02541
    ## alternative hypothesis: true location shift is not equal to 0
