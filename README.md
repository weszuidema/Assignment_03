# Assignment 03

Instructions

1. [Fork this repository](https://help.github.com/articles/using-pull-requests/) to your GitHub account.
2. Write your solutions in R Markdown in a file named `solutions.Rmd`.
3. When you are ready to submit your assignment, [initiate a pull request](https://help.github.com/articles/using-pull-requests/#initiating-the-pull-request). Title your
pull request "Submission".

To update your fork from the upstream repository:

1. On your fork, e.g. `https://github.com/jrnold/Assignment_03` click on "New Pull request"
2. Set your fork `jrnold/Assignment_03` as the base fork on the left, and `UW-POLS503/Assignment_03` as the head fork on the right. In both cases the branch will be master. This means, compare any canes in the head fork that are not in the base fork. You will see differences between the `US-POLS503` repository and your fork. Click on "Create Pull Request", and if there are no issues, "Click Merge" A quick way is to use this link, but change the `jrnold` to your own username: `https://github.com/jrnold/Assignment_03/compare/master...UW-POLS503:master`.

We'll use these packages,

```r
library("foreign")
library("dplyr")
library("broom")
library("ggplot2")
library("DT")
```
Since we are going to do some simulation, we should set a seed, so the results are exactly replicable.

```r
set.seed(1234)
```
Since some of these computations will take time, we can cache the results so that knitr will
only run code that has changed.

```r
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE)
```

# Nunn and Wantchekon AER 2011 example

Let's run some regressions from 

> Nunn, Nathan and Leonard Wantchekon. 2011. "The Slave Trade and the Origins of Mistrust in Africa."	American Economic Review,
> 101(7):3221-52. [doi:10.1257/aer.101.7.3221](https://dx.doi.org/10.1257/aer.101.7.3221)

The replication data for the is available from its [AER site](https://dx.doi.org/10.1257/aer.101.7.3221), but the main
dataset is included in this repository.
Since the main dataset is a Stata `.dta` file load it using the `read.dta` function
and convert it to a **dplyr** `tbl` so the `print` function produces nicer output.

```r
nunn <- read.dta("Nunn_Wantchekon_AER_2011.dta") %>% tbl_df()
```

There are many variables in this data.
When `read.dta` converts a Stata data file the descriptions of the variables end up
in an R [attribute](https://stat.ethz.ch/R-manual/R-devel/library/base/html/attributes.html) `"var.labels"`.
Print out the variable labels to get the descriptions of the files[^datatable]

```r
data_frame(variable = names(nunn), description = attr(nunn, "var.labels")) %>%
  datatable(class = 'cell-border stripe')
```

<!--html_preserve--><div id="htmlwidget-1137" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-1137">{"x":{"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59"],["respno","ethnicity","murdock_name","isocode","region","district","townvill","location_id","trust_relatives","trust_neighbors","intra_group_trust","inter_group_trust","trust_local_council","ln_export_area","export_area","export_pop","ln_export_pop","age","age2","male","urban_dum","occupation","religion","living_conditions","education","near_dist","distsea","loc_murdock_name","loc_ln_export_area","local_council_performance","council_listen","corrupt_local_council","school_present","electricity_present","piped_water_present","sewage_present","health_clinic_present","district_ethnic_frac","frac_ethnicity_in_district","townvill_nonethnic_mean_exports","district_nonethnic_mean_exports","region_nonethnic_mean_exports","country_nonethnic_mean_exports","murdock_centr_dist_coast","centroid_lat","centroid_long","explorer_contact","railway_contact","dist_Saharan_node","dist_Saharan_line","malaria_ecology","v30","v33","fishing","exports","ln_exports","total_missions_area","ln_init_pop_density","cities_1400_dum"],["Respondent number in Afrobarometer dataset","Ethnicity name: from Afrobarometer q79","Ethnicity name: from Murdock","Country 3 digit iso code","Region: from Afrobarometer","District: from Afrobarometer","Town/village: from Afrobarometer","Unique location-of-repondent identifier - based on isocode region district townv","Trust of relatives: q84a","Trust of neighbors: q84b","Intra-group trust: q84d","Inter-group trust: q84c","Trust of local government council: q55d","Log [(total slave exports: Atlantic + Indian) / area (km^2)","(total slave exports: Atlantic + Indian) / area (km^2)","Exports divided by historic Murdock population","Ln (1+exports/Murdock historic population)","Age: q1","Age squared","Indicator for respondent being male: q101","Indicator for respondent living in urban area","Occupation categories: q95","Religion categories: q91","Living condition categories:q4b","Education categories: Afrobarometer q90","Current distance from coast 1000s kms","Historic distance from coast 1000s kms","Murdock identifier for the current location of the respondent","Slave exports measure based on current location of respondent","Perceived performance of local council: q68c","Does the local council listen: q62b","How much corruption in local council: q56c","Is there a school in the PSU: q116b","Is there electricity in the PSU: q116d","Is there piped water in the PSU: q116e","Is there sewage in the PSU: q116f","Is there a health clinic in the PSU: q116g","District-level ethnic fractionalization","Proportion of ethnic group in district","Avg slave exports of other ethnicities within town/village","Avg slave exports of other ethnicities within district","Avg slave exports of other ethnicities within region","Avg slave exports of other ethnicities within country","Historic distance of ethnicity's centroid from the coast (in kms)","Historic latitude of centroid of ethnic group","Historic longitude of centroid of ethnic group","Indicator for historic contact with European explorer","Indicator variable for historic integration into the colonial railway network","Historic distance of ethnicity's centroid from a centroid (town) in Saharan trad","Historic distance of ethnicity's centroid from a line (route) in Saharan trade (","Ethnic groups average malaria ecology measure","Pre-colonial settlement patterns of ethnicity: from Ethngraphic Atlas v30","Pre-colonial juris. hierarchy beyond the local community: Ethnographic Atlas v33","Pre-colonial reliance on fishing: Ethnographic Atlas v3","(Atlantic+Indian Exports)","ln(1+Atlantic+Indian Exports)","Total Catholic + Protestant mission per land area","Log population density during the colonial period - from Murdock","Indicator for existence of city among ethnic group in 1400"]],"container":"<table class=\"cell-border stripe\">\n  <thead>\n    <tr>\n      <th> \u003c/th>\n      <th>variable\u003c/th>\n      <th>description\u003c/th>\n    \u003c/tr>\n  \u003c/thead>\n\u003c/table>","options":{"order":[],"autoWidth":false,"orderClasses":false,"columnDefs":[{"orderable":false,"targets":0}]},"callback":null,"filter":"none"},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

[^datatable]: This uses the [DT](https://rstudio.github.io/DT/) package to produce pretty interactive tables in the HTML.

In Table 1, NW run several models with Trust in Neighbors as an outcome variable,
different measures of slave exports as the treatment variable, and the same set of controls variables.
Some of the relevant variables in the data are:

- `trust_neighbors`: Trust of neighbors
- `exports`: Slave exports in 1000s
- Individual controls: `age`, `age2`, `male`, `urban_dum`, `education`, `occupation`, `religion`, `living_conditions`
- District controls: `district_ethnic_frac`, `frac_ethnicity_in_district`
- Country-fixed effects: `isocode`

Note that NW use `education`, `occupation`, `religion`, and `living_conditions` as factor variables. Convert them accordingly,

```r
factor_vars <- c("education", "occupation", "religion", "living_conditions")
for (i in factor_vars) {
  nunn[[i]] <- factor(nunn[[i]])
}
```


## Bivariate regression

Run a regression of the Trust of Neighbors on Slave exports.
This is Table 1, Model 1, without any of the control variables.

```r
mod_1_0 <- lm(trust_neighbors ~ exports, data = nunn)
```

<div class="bs-callout bs-callout-info">
- Interpret the magnitude and statistical significance of the coefficient on `trust_neighbors`.
- Plot the fitted values and confidence interval of the fitted values of regression vs. `exports`.
- Plot the residuals of this regression against the fitted values of the regression. Do they appear to have
  constant variance? Are they approximately symmetric?
- What is the null hypothesis of the t-test reported by `summary()`? Explain the meaning of the p-value.
  Be precise. Is the p-value the probability that the null hypothesis is correct?
</div>

Example of using `augment` for fitted values:

```r
augment(mod_1_0, data_frame(exports = seq(min(nunn$exports, na.rm = TRUE), max(nunn$exports, na.rm = TRUE), length.out = 100)))
```
Example of plotting the confidence intervals

```r
ggplot() + geom_line(data = mod_1_0_fitted, mapping = aes(x = exports, y = .fitted)) + geom_ribbon(data = mod1_fitted, mapping = aes(x = exports, y = .fitted, ymin = .fitted - 2 * .se.fit, ymax = .fitted + 2 * .se.fit), alpha = 0.3) + geom_point(data = nunn, aes(x = exports, y = trust_neighbors))
```


## Probablities of Hypotheses

Frequentist statistics assigns no probabilities to hypotheses (parameter values).
They are either true or false, but they are unknown. Only samples are random variables, and have an associated probability.
But as scientists, we are generally interested in the probability that a hypothesis is correct.[^nature]
The probability that the research hypothesis ($H_0$) is correct can be calculated with Bayes law,
$$
p(H_0 | \text{data}) =
\frac{p(\text{data} | H_0) p(H_0)}{p(\text{data} | H_a) p(H_a) + p(\text{data} | H_0) p(H_0)} = \frac{p(\text{data} | H_0) p(H_0)}{p(\text{data})}
$$
Working somewhat informally, the p-value gives $p(\text{data} | H_0)$. An important missing piece of information is the baseline or prior probability that the null hypothesis is true, $p(H_0)$, which is the complement of the probability that the research hypothesis is true, $p(H_0) = 1 - p(H_a)$,[^h0] [^jeff]

<div class="bs-callout bs-callout-info">
- If more than the p-value is required to make sense of the research findings, what does the article do to increase your belief about the research hypothesis, $p(H_a)$?
- Suppose you believed that NW were p-value hacking (which I don't think they are!). What part of Bayes law is that 
  affecting? If you think that someone is p-value hacking, then you are saying that they will always produce significant p-values regardless of whether the null or alternative hypotheses are true.
</div>

[^nature]: See the discussion in (Scientific method: Statistical errors)[http://www.nature.com/news/scientific-method-statistical-errors-1.14700], *Nature*.

[^jeff]: This question is non-standard and idiosyncratic to the way I interpret research.

[^h0]: Assuming, for simplicity, that $H_0$ and $H_a$ are the only hypotheses so that $p(H_0) + p(H_a) = 1$.


## Multiple regression

In the models in Table 1, NW includes control variables to account for individual, district, and country-level 
variables that may explain differences. 

Run the model in Table 1, Model 1:

```r
mod_1_1 <- lm(trust_neighbors ~
              exports +
              # controls
              # individual level
              age + age2 + male + urban_dum + education +
              occupation + religion + living_conditions +
              # district-level 
              district_ethnic_frac + frac_ethnicity_in_district +
              # country-level
              isocode,
              data = nunn)
mod_1_1
```

```
## 
## Call:
## lm(formula = trust_neighbors ~ exports + age + age2 + male + 
##     urban_dum + education + occupation + religion + living_conditions + 
##     district_ethnic_frac + frac_ethnicity_in_district + isocode, 
##     data = nunn)
## 
## Coefficients:
##                (Intercept)                     exports  
##                  1.620e+00                  -6.791e-04  
##                        age                        age2  
##                  8.396e-03                  -5.473e-05  
##                       male                   urban_dum  
##                  4.550e-02                  -1.405e-01  
##                 education1                  education2  
##                  1.710e-02                  -5.225e-02  
##                 education3                  education4  
##                 -1.374e-01                  -1.890e-01  
##                 education5                  education6  
##                 -1.893e-01                  -2.401e-01  
##                 education7                  education8  
##                 -2.851e-01                  -1.232e-01  
##                 education9                 occupation1  
##                 -2.406e-01                   6.186e-02  
##                occupation2                 occupation3  
##                  7.392e-02                   3.356e-02  
##                occupation4                 occupation5  
##                  7.942e-03                   6.661e-02  
##                occupation6                 occupation7  
##                 -7.563e-02                   1.700e-02  
##                occupation8                 occupation9  
##                 -9.428e-02                  -9.981e-02  
##               occupation10                occupation11  
##                 -3.307e-02                  -2.300e-02  
##               occupation12                occupation13  
##                 -1.565e-01                  -1.441e-02  
##               occupation14                occupation15  
##                 -5.566e-02                  -2.344e-01  
##               occupation16                occupation18  
##                 -1.307e-02                  -1.730e-01  
##               occupation19                occupation20  
##                 -1.770e-01                  -2.458e-02  
##               occupation21                occupation22  
##                 -4.937e-02                  -1.069e-01  
##               occupation23                occupation24  
##                 -9.712e-02                   1.292e-02  
##               occupation25               occupation995  
##                  2.623e-02                  -1.195e-03  
##                  religion2                   religion3  
##                  5.396e-02                   7.888e-02  
##                  religion4                   religion5  
##                  4.749e-02                   4.318e-02  
##                  religion6                   religion7  
##                 -1.788e-02                  -3.617e-02  
##                 religion10                  religion11  
##                  6.015e-02                   2.238e-01  
##                 religion12                  religion13  
##                  2.627e-01                  -6.813e-02  
##                 religion14                  religion15  
##                  4.674e-02                   3.845e-01  
##                religion360                 religion361  
##                  3.657e-01                   3.416e-01  
##                religion362                 religion363  
##                  8.230e-01                   3.857e-01  
##                religion995          living_conditions2  
##                  4.161e-02                   4.396e-02  
##         living_conditions3          living_conditions4  
##                  8.627e-02                   1.197e-01  
##         living_conditions5        district_ethnic_frac  
##                  1.204e-01                  -1.554e-02  
## frac_ethnicity_in_district                  isocodeBWA  
##                  1.011e-01                  -4.259e-01  
##                 isocodeGHA                  isocodeKEN  
##                  1.135e-02                  -1.820e-01  
##                 isocodeLSO                  isocodeMDG  
##                 -5.511e-01                  -3.316e-01  
##                 isocodeMLI                  isocodeMOZ  
##                  7.528e-02                   8.224e-02  
##                 isocodeMWI                  isocodeNAM  
##                  3.062e-01                  -1.398e-01  
##                 isocodeNGA                  isocodeSEN  
##                 -2.382e-01                   3.867e-01  
##                 isocodeTZA                  isocodeUGA  
##                  2.079e-01                  -6.444e-02  
##                 isocodeZAF                  isocodeZMB  
##                 -2.179e-01                  -2.173e-01
```

<div class="bs-callout bs-callout-info">
- Interpret the coefficient on `exports`
- How much does the coefficient change with the addition of control variables? What does that suggest?
- Do the R^2 and number of observations match those reported in Table 1?
- Calculate the fitted values of the regression by multiplying the $\beta$ vector and the $\mat{X}$ matrix.
  Confirm that you get the same results as using `predict()`.
- How would you create a plot that shows the predicted values of `trust_neighbors` as the value of `exports` changes?
  What is different about the multiple regression case than the bivariate case?
</div>


### Understanding Multiple Regression

<div class="bs-callout bs-callout-info">
- Run the following regressions
    1. Regress regression of `trust_neighbors` on the controls.
        ```r
        lm(trust_neighbors ~ age + age2 + male + urban_dum +
           education + occupation + religion + living_conditions +
           district_ethnic_frac + frac_ethnicity_in_district +
           isocode, data = nunn)
        ```
        Save the residuals.
    2. Run the regression of `exports` on the controls. Save the residuals
        ```r
        lm(exports ~ age + age2 + male + urban_dum +
           education + occupation + religion + living_conditions +
           district_ethnic_frac + frac_ethnicity_in_district +
           isocode, data = nunn)
        ```
        Save the residuals.
    3. Regress the residuals from 1. on the residuals on 2.
- How does the coefficient from regression 3 compare the the coefficient on `exports` from the regression in Table 1, Model 1?
    What does that say about what multiple regression is doing?
- Are the steps 
</div>


## Validity of the standard errors

One of the assumptions necessary for OLS standard errors to be correct is homoskedasticity homoskedasticity (constant variance), and that the errors are uncorrelated.

<div class="bs-callout bs-callout-info">
- How might that assumption be violated?
- Plot the residuals of the regression by district. Do they appear to be uncorrelated? What does that say about the validity of the OLS standard errors?
- Do the standard errors match those reported in Table 1 of the article? What sort of standard errors does the article use?
</div>


## Regressions with log slave exports per capita

Run the regression in Table 1, model 6, which uses "log(1 + exports / pop)" as a measure of slave exports.


```r
mod_1_6 <- lm(trust_neighbors ~ ln_export_pop + 
              age + age2 + male + urban_dum + education +
              occupation + religion +
              living_conditions + district_ethnic_frac +
              frac_ethnicity_in_district + isocode,
              data = nunn)  
```


<div>
- Interpret the effect of `ln_export_pop` on `trust_neighbors`
- Why is "log(1 + exports / pop)" used as the measure instead of "log(exports / pop)"?
- Plot the fitted values of log(1 + exports / pop) and their confidence interval against "log(1 + exports / pop)" against the residuals of the controls only regression. Include the line, confidence intervals, and data points.
- Plot the fitted values of exports / pop against the residuals of the controls only regression. Include the line, confidence intervals, and data points. How does this relationship differ from the one which used the level of slave exports with out taking the logarithm or adjusting for population?
</div>


## Sampling distribution of OLS coefficients

Let's understand what the confidence intervals mean in terms of the sampling distribution.
Since we don't know the true parameter values for this, we will pretend that the OLS point estimates from the regression are the "true" population parameters.

The plan to generate a sampling distribution of $\beta$ is:

1. Draw a new sample $\tilde{y}_i \sim N(\hat{y}_i, \hat{\sigma}^2)$.
2. Estimate OLS estimates $\tilde{\vec{y}} = \tilde{vec{beta}} \mat{X} + \tilde{vec{varepsilon}} = \hat{y} + \tilde{\vec{\varepsilon}}$.
3. Repeat steps 1--2, `iter` times, store $\beta^*$ for each iteration, and
   return the estimates for all samples.
   
Then, the distribution of the $\beta^*$ is a sampling distribution of the parameters.

<div class="bs-callout bs-callout-info">
Why is only $\vec{y}$ being samples? Why is $\mat{X}$ fixed in these simulations? See Wooldridge Ch 2 and 3 discussion of the assumptions of OLS.
</div>

Let's take the results of the model on `ln_export_pop` and explore the sampling distribution of $\beta$ from that model.

First run the model,

```r
mod <- lm(trust_neighbors ~ ln_export_pop + 
          age + age2 + male + urban_dum + education +
          occupation + religion +
          living_conditions + district_ethnic_frac +
          frac_ethnicity_in_district + isocode,
          data = nunn, na.action = na.exclude)
```
The argument `na.action = na.exclude` ensures that when we calculate residuals, etc. they will be padded with missing
values so they are the same length as the original `nunn` data.
There are other ways to work around this, but this makes the code to run the simulations easier, especially the step
that draws `y_hat`.

Extract the values of the parameter estimates, $\hat{\beta}$, the model matrix, $X$,the regression standard error, $\hat{\sigma}$, and the number of observations, $N$.

```r
y_hat <- predict(mod, na.action = na.exclude)
sigma <- sqrt(sum(residuals(mod) ^ 2, na.rm = TRUE) / mod$df.residual)
n <- nrow(nunn)
```
Later we'll also need the original 
Choose a number of iterations to run. 
For this example use 1,024.

```r
iter <- 1024
```
Create a list to store the results

```r
results <- vector(mode = "list", length = iter)
```

For iterations `1 ... iter` we cant to 

- Draw the regression errors from i.i.d normal distributions, $\tilde\epsilon_i \sim N(0, \hat\sigma^2)$.
- Generate a new dependent variable, $\tilde{\vec{y}} = \mat{X} \hat{\vec{\beta}} + \tilde{\vec{epsilon}}$.
- Run an OLS regression to estimate the coefficients on the new data, $\tilde{\vec{y}} = \mat{X} \tilde{\vec{\beta}} + \tilde{\vec{\epsilon}}$
- Save the $\tilde{\beta}$


```r
p <- progress_estimated(iter, min_time = 2)
for (i in seq_len(iter)) {
  # draw errors
  errors <- rnorm(n, mean = 0, sd = sigma)
  # create new outcome variable from errors
  nunn[["trust_neighbors_new"]] <- y_hat + errors
  # Replace the dependent variable with the newly sampled y
  newmod <- lm(trust_neighbors_new ~ ln_export_pop + 
            age + age2 + male + urban_dum + education +
            occupation + religion +
            living_conditions + district_ethnic_frac +
            frac_ethnicity_in_district + isocode,
            data = nunn)
  # Formula objects are manipulable in R. So a more general
  # way to do the above is to alter the formula from the model
  # by replacing trust_neighbors ~ ... with trust_neighbors_new ~ ...
  # formula <- formula(mod$terms)
  # formula[[2]] <- "trust_neighbors_imputed"
  # # check that this worked
  # # print(formula)
  # # Also this should be put outside the loop, since it doesn't
  # # change
  # newmod <- lm(formula, data = nunn)
  # Save the coefficients as a data frame to the list
  results[[i]] <- tidy(newmod) %>% mutate(.iter = i)
  # Progress bar update
  p$tick()$print()
}
# clean up: remove the new variable
nunn[["trust_neighbors_new"]] <- NULL
```
Finally, since `results` is a list of data frames, stack the data frames in the list to form a single data frame that is easier to analyze:

```r
results <- bind_rows(results)
```
**Note:** this will take a few minutes.


<div class="bs-callout bs-callout-info">
Use the results  of this simulation to answer the following questions:

- Plot the distribution of the coefficients of `ln_export_pop`.
- What is the standard deviation of the sampling distribution of the coefficient of `ln_export_pop`? How does this compare to the standard error of this coefficient given by `lm()`?
- Calculate the correlation matrix of the coefficients. For the first step will need to create a data frame using `spread` in which the rows are iterations, the columns are coefficients, and the values are the esimates.
- Plot that correlation matrix using `geom_raster`. Are the coefficients of coefficients uncorrelated? In general, when would coefficients be more or less correlated?
- Why is this simulation not (directly) appropriate for calculating a $p$-value. What distribution would you have to simulate from to calculate a $p$-value?
</div>


## Non-parametric Bootstrap

The previous question was an example of parametric bootstrap.
It is a parametric bootstrap because you drew data from an assumed model (the OLS model that you estimated).

An alternative is a non-parametric bootstrap.
In a non-parametric bootstrap, instead of drawing samples from model, we are going to redraw samples from the sample.

An analogy is that the sample is to the population as the bootstrap is to the sample.
We are treating the sample distribution as an estimate of the population distribution and then drawing samples from
that estimated population distribution.

To do the bootstrapping we will use the `bootstrap` function in the **tidyr** package.
However, the [boot](https://cran.r-project.org/web/packages/boot/index.html) package supports many more advanced methods of bootstrapping.

Let's start by drawing a single bootstrap replication.
It is a sample of the same size as the original data, drawn from the data *with replacement*.

```r
nunn_bootstrapped <- bootstrap(nunn, 1)
```

So, in order to calculate bootstrap standard errors, we will need to draw a sample of 
To get bootstrap standard errors, we draw `B` replications, run an  regression, and save the estimates. 

```r
beta_bs <- 
  bootstrap(nunn, 1024) %>%
    do(tidy(lm(trust_neighbors ~ ln_export_pop, data = .)))
```

There are several ways to calculate standard errors from the bootstrap replications.
The following are two simple methods.

1. Calculate the standard error from these simulations by taking the standard deviation of the estimates.
   Suppose $\beta^{*b}_k$ is the estimated coefficient from replication $b \in 1:B$, and $\bar\beta^{*}_k = (\sum \beta^{*b}_k) / B$.
   Then the bootstrap standard error is,
   $$
   \se_k(\hat\beta_{k}) = \sqrt{\frac{1}{B - 1} \sum (\beta^{*b}_k - \bar\beta^{*b}_k)^2}
   $$
   The confidence interval is thus,
   $$
   \hat{\beta}_k \pm \se_{bs}(\hat\beta_k)
   $$
   Note that you use the estimate $\hat{\beta}_k$ from the original model, not the mean of the bootstrap estimates.
   This method works well if the sampling distribution of $\beta_k$ is symmetric.

2. The second method is to use the quantiles of the bootstrap estimates.
   E.g. a 95% confidence interval uses the 2.5% and 97.5% quantiles of the bootstrap estimates.
   This method allows for asymmetric confidence intervals. However, it takes more replications to get accurate 
   values of extreme quantiles than it does to calculate a standard deviation.

<div class="bs-callout bs-callout-info">
- Estimate the bootstrapped confidence intervals using those two methods.
- Compare the bootstrapped confidence intervals to the OLS confidence interval.
</div>

There are even more advanced methods such as the studentized bootstrap, and the adjusted bootstrap percentile (BCa) methods
included in `boot.ci`.

For bootstrapped standard errors to be valid, the samples from the data need to be taken in the same way as the sample
was taken from the population. 
For example, in a time series it would be inappropriate to sample observations without accounting for their order.

<div class="bs-callout bs-callout-info">
- What is the population in this paper?
- How was the sample drawn from this population?
- In the previous examples, did we draw the sample in the same way as it was drawn from the population? What would
  be a better way of drawing the bootstrapped samples?
  Try to implement it; see the `group_by` argument of `bootstrap`.
</div>

## F-test example

An $F$-test tests the null hypothesis that several coefficients in the regression are all 0 vs. the alternative that 
at least one of the coefficients is non-zero.
Suppose you want to test that the $q$ coefficients $\beta_j$ through $\beta_{j + q}$ are all 0,
$$
\begin{aligned}[t]
H_0: &\quad \beta_j = \dots = \beta_J = 0
H_a: &\quad \text{at least one $\beta_k \neq 0$}
\end{aligned}
$$

To run an F-test in R, use the `anova()` function to compare two models.
For example, to compare the regression of `trust_neighbors` on `exports` without controls to the regression with controls, use

```r
mod_1_0 <- lm(trust_neighbors ~ ln_export_pop, data = nunn)
mod_1_1 <- lm(trust_neighbors ~ ln_export_pop + age + age2 + male + urban_dum +
              education + occupation + religion + living_conditions +
              district_ethnic_frac + frac_ethnicity_in_district +
              isocode, data = nunn)
anova(mod_1_0, mod_1_1)
```

```
## Error in anova.lmlist(object, ...): models were not all fitted to the same size of dataset
```
We can't do it! At least not yet. The problem is that `lm()` drops all rows with at least one missing value.
So `mod_1_0` and `mod_1_1` run a regression on datasets with different numbers of observations,

```r
mod_1_0$df.residual + length(mod_1_0$coefficients)
```

```
## [1] 18112
```

```r
mod_1_1$df.residual + length(mod_1_1$coefficients)
```

```
## [1] 17644
```
The residual degrees of freedom is $N - K - 1$, and the length of the coefficient vector is $K + 1$, so
their sum is the number of observations in the regression, $(N - K - 1) + (K + 1) = N$.

To ensure that these models are run on the same set of data, create a subset of the `nunn` data that
has all the variables in the larger model, and drop an observations with missing values in any of those
variables, using `na.omit()`.

```r
nunn_nomiss <- nunn %>%
  select(trust_neighbors, ln_export_pop, age, age2, male, urban_dum,
         education, occupation, religion, living_conditions,
         district_ethnic_frac, frac_ethnicity_in_district,
         isocode) %>%
  na.omit()
```
Now, we can compare the models,

```r
mod_1_0 <- lm(trust_neighbors ~ ln_export_pop, data = nunn_nomiss)
mod_1_1 <- lm(trust_neighbors ~ ln_export_pop + age + age2 + male + urban_dum +
              education + occupation + religion + living_conditions +
              district_ethnic_frac + frac_ethnicity_in_district +
              isocode, data = nunn_nomiss)
anova(mod_1_0, mod_1_1)
```

```
## Analysis of Variance Table
## 
## Model 1: trust_neighbors ~ ln_export_pop
## Model 2: trust_neighbors ~ ln_export_pop + age + age2 + male + urban_dum + 
##     education + occupation + religion + living_conditions + district_ethnic_frac + 
##     frac_ethnicity_in_district + isocode
##   Res.Df   RSS Df Sum of Sq      F    Pr(>F)    
## 1  17642 17251                                  
## 2  17566 14777 76    2473.7 38.693 < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
The p-value for the F-test is approximately 0, so the test rejects the null 
hypothesis that all the control variables are equal to 0 at all commonly used
levels of significance.

<div class="bs-callout bs-callout-info">
- Run and interpret an F-test for a reasonable subset of coefficients.
- Why can't you use F-tests to compare different models in Table 1?
- Run an F-test comparing the model with only controls to the one in Model 1, Table 6.
    In other words, the null hypothesis is $\beta_{\mathtt{ln\_exports\_pop}} = 0$.
    ```r
    mod_controls <- lm(trust_neighbors ~ age + age2 + male + urban_dum +
              education + occupation + religion + living_conditions +
              district_ethnic_frac + frac_ethnicity_in_district +
              isocode, data = nunn)    
    mod_1_6 <- lm(trust_neighbors ~ ln_export_pop + age + age2 + male + urban_dum +
              education + occupation + religion + living_conditions +
              district_ethnic_frac + frac_ethnicity_in_district +
              isocode, data = nunn)
    ```
    How does the $p$-value of the F-test compare to the p-value of the regression
    coefficient on `ln_export_pop`? Square the t-statistic for the coefficient
    of `ln_export_pop`; how does it compare to the F-statistic? What is the
    relationship between a t-test and an F-test for a single parameter in a regression?
</div>
