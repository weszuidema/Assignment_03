# Assignment 03
$$
\DeclareMathOperator{\cor}{cor}
\DeclareMathOperator{\cov}{cov}
\DeclareMathOperator{\quantile}{quantile}
$$

## Instructions

1. [Fork this repository](https://help.github.com/articles/using-pull-requests/) to your GitHub account.
2. Write your solutions in R Markdown in a file named `solutions.Rmd`.
3. When you are ready to submit your assignment, [initiate a pull request](https://help.github.com/articles/using-pull-requests/#initiating-the-pull-request). Title your
pull request "Submission".

To update your fork from the upstream repository:

1. On your fork, e.g. `https://github.com/jrnold/Assignment_03` click on "New Pull reqest"
2. Set your fork `jrnold/Assignment_03` as the base fork on the left, and `UW-POLS503/Assignment_03` as the head fork on the right. In both cases the branch will be master. This means, compare any chanes in the head fork that are not in the base fork. You will see differences between the `US-POLS503` repo and your fork. Click on "Create Pull Request", and if there are no issues, "Click Merge" A quick way is to use this link, but change the `jrnold` to your own username: `https://github.com/jrnold/Assignment_03/compare/master...UW-POLS503:master`.

We'll use these packages,

```r
library("foreign")
library("dplyr")
library("broom")
library("ggplot2")
library("DT")
```
Since we are going to do some simulation, we shoudl set a seed, so the results are exactly replicable.

```r
set.seed(1234)
```
Since some of these computations will take time, we can cache the results so that knitr will
only run code that has changed.

```r
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE)
```

## Problem 1

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
Print out the variable lables to get the descriptions of the files[^datatable]

```r
data_frame(variable = names(nunn), description = attr(nunn, "var.labels")) %>%
  datatable(class = 'cell-border stripe')
```

<!--html_preserve--><div id="htmlwidget-1137" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-1137">{"x":{"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59"],["respno","ethnicity","murdock_name","isocode","region","district","townvill","location_id","trust_relatives","trust_neighbors","intra_group_trust","inter_group_trust","trust_local_council","ln_export_area","export_area","export_pop","ln_export_pop","age","age2","male","urban_dum","occupation","religion","living_conditions","education","near_dist","distsea","loc_murdock_name","loc_ln_export_area","local_council_performance","council_listen","corrupt_local_council","school_present","electricity_present","piped_water_present","sewage_present","health_clinic_present","district_ethnic_frac","frac_ethnicity_in_district","townvill_nonethnic_mean_exports","district_nonethnic_mean_exports","region_nonethnic_mean_exports","country_nonethnic_mean_exports","murdock_centr_dist_coast","centroid_lat","centroid_long","explorer_contact","railway_contact","dist_Saharan_node","dist_Saharan_line","malaria_ecology","v30","v33","fishing","exports","ln_exports","total_missions_area","ln_init_pop_density","cities_1400_dum"],["Respondent number in Afrobarometer dataset","Ethnicity name: from Afrobarometer q79","Ethnicity name: from Murdock","Country 3 digit iso code","Region: from Afrobarometer","District: from Afrobarometer","Town/village: from Afrobarometer","Unique location-of-repondent identifier - based on isocode region district townv","Trust of relatives: q84a","Trust of neighbors: q84b","Intra-group trust: q84d","Inter-group trust: q84c","Trust of local government council: q55d","Log [(total slave exports: Atlantic + Indian) / area (km^2)","(total slave exports: Atlantic + Indian) / area (km^2)","Exports divided by historic Murdock population","Ln (1+exports/Murdock historic population)","Age: q1","Age squared","Indicator for respondent being male: q101","Indicator for respondent living in urban area","Occupation categories: q95","Religion categories: q91","Living condition categories:q4b","Education categories: Afrobarometer q90","Current distance from coast 1000s kms","Historic distance from coast 1000s kms","Murdock identifier for the current location of the respondent","Slave exports measure based on current location of respondent","Perceived performance of local council: q68c","Does the local council listen: q62b","How much corruption in local council: q56c","Is there a school in the PSU: q116b","Is there electricity in the PSU: q116d","Is there piped water in the PSU: q116e","Is there sewage in the PSU: q116f","Is there a health clinic in the PSU: q116g","District-level ethnic fractionalization","Proportion of ethnic group in district","Avg slave exports of other ethnicities within town/village","Avg slave exports of other ethnicities within district","Avg slave exports of other ethnicities within region","Avg slave exports of other ethnicities within country","Historic distance of ethnicity's centroid from the coast (in kms)","Historic latitude of centroid of ethnic group","Historic longitude of centroid of ethnic group","Indicator for historic contact with European explorer","Indicator variable for historic integration into the colonial railway network","Historic distance of ethnicity's centroid from a centroid (town) in Saharan trad","Historic distance of ethnicity's centroid from a line (route) in Saharan trade (","Ethnic groups average malaria ecology measure","Pre-colonial settlement patterns of ethnicity: from Ethngraphic Atlas v30","Pre-colonial juris. hierarchy beyond the local community: Ethnographic Atlas v33","Pre-colonial reliance on fishing: Ethnographic Atlas v3","(Atlantic+Indian Exports)","ln(1+Atlantic+Indian Exports)","Total Catholic + Protestant mission per land area","Log population density during the colonial period - from Murdock","Indicator for existence of city among ethnic group in 1400"]],"container":"<table class=\"cell-border stripe\">\n  <thead>\n    <tr>\n      <th> \u003c/th>\n      <th>variable\u003c/th>\n      <th>description\u003c/th>\n    \u003c/tr>\n  \u003c/thead>\n\u003c/table>","options":{"order":[],"autoWidth":false,"orderClasses":false,"columnDefs":[{"orderable":false,"targets":0}]},"callback":null,"filter":"none"},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->


[^datatable]: This uses the [DT](https://rstudio.github.io/DT/) package to produce pretty interactive tables in the HTML.

In Table 1, NW run several models with Trust in Neighbors as an outcome variable,
different measures of slave exports as the treatment variable, and the same set of controls variables.
The relevant variables in the data are:

- `trust_neighbors`: Trust of neighbors
- `exports`: Slave exports in 1000s
- Individual controls: `age`, `age2`, `male`, `urban_dum`, `education`, `occupation`, `religion`, `living_conditions`
- District controls: `district_ethnic_frac`, `frac_ethnicity_in_district`
- Country-fixed effects: `isocode`
- 

Run a regression of the Trust of Neighbors on Slave exports.
This is Table 1, Model 1, witout any of the control variables.

```r
mod1 <- lm(trust_neighbors ~ exports, data = nunn)
```

<div class="bs-callout bs-callout-info">
- Interpret the magnitude and statistical significance of the coefficient on `trust_neighbors`.
- What is the null hypothesis of the t-test reported by `summary()`? Explain the meaning of the p-value.
  Be precise. Is the p-value the probability that the null hypothesis is correct?
</div>

Frequentist statistics assigns no probabilities to hypotheses (parameter values).
They are either true or false, but they are unknown. Only samples are random variables, and have an associated probability.
But as scientists, we are generally interested in the probability that a hypothesis is correct.[^nature]
This can be calculated using Bayes law as,
$$
p(H_0 | \text{data}) = \frac{p(\text{data | H_0} p(H_0))}{p(\text{data})}
$$
The p-value gives, teh The important missing piece of information is $p(! H_0)$, the probability that the null hypothesis is not correct, i.e. the probability that the null hypothesis is correct, which is the probability that the research hypothesis is incorrect.[^jeff]

<div class="bs-callout bs-callout-info">
- If more than the p-value is required to make sense of the research findings, how does the article increase your belief
  about $p(H_a)$?
- Suppose you believed that NW were p-value hacking (which I don't think they are!). What part of Bayes law is that 
  affecting?
</div>

[^nature]: See the discussion in (Scientific method: Statistical errors)[http://www.nature.com/news/scientific-method-statistical-errors-1.14700], *Nature*.

[^jeff]: This question is non-standard and idiosyncratic to the way I interpret research.

<div class="bs-callout bs-callout-info">
What are the control variables that NW include in the models in Table 1?
Run the model in Table 1, Regression 1.

Does the R^2 match that in Table 1?

Do the standard errors match those in Table 1? Any guesses why?

Run the regression in Table 1, model 6. Why is it "log(1 + exports / pop)" instead of "log(exports / pop)"?
</div>


```r
resampler_coef <- function(mod, .data, iter = 1) {
  # Remove missing values
  .data <- na.omit(.data)
  # Coefficients
  beta <- coef(mod)
  # mod$terms contains the formula used in the regression
  X <- model.matrix(mod$terms, data = .data)
  # estimate of std. dev. of errors
  sigma <- sqrt(sum(mod$residuals ^ 2) / mod$df.residual)
  # This produces the same result
  # sigma <- summary(mod1)$sigma  
  # Number of observations
  n <- nrow(X)
  # Name of dependent variable
  outcome_var_name <- all.vars(mod$terms)[1]
  # List to save results
  results <- vector(mode = "list", length = iter)
  for (i in seq_len(iter)) {
    # draw errors
    errors <- rnorm(n, mean = 0, sd = sigma)
    # create new outcome variable from errors
    y <- X %*% beta + errors
    # replace outcome variable
    .data[[outcome_var_name]] <- y
    # run regression
    newmod <- lm(mod$terms, data = .data)
    # Save coefficients as a data frame to the list
    results[[i]] <- tidy(newmod) %>% mutate(.iter = i)
  }
  # Convert the list of data frames to a single data frame by stacking the iterations
  bind_rows(results)
}
```

- Plot the distributions of the coefficients
- Calculate the correlation matrix of the coefficients. How similar is it to that from `vcov`?


### F-test example

- Run F-tests of the multiple regression model vs. the model with no controls.
- Run and interpet an F-test on some reasonable group of variables.

F-test simulations


```r
resampler_models <- function(mod, .data, iter = 1) {
  # Remove missing values
  .data <- na.omit(.data)
  # Coefficients
  beta <- coef(mod)
  # mod$terms contains the formula used in the regression
  X <- model.matrix(mod$terms, data = .data)
  # estimate of std. dev. of errors
  sigma <- sqrt(sum(mod$residuals ^ 2) / mod$df.residual)
  # This produces the same result
  # sigma <- summary(mod1)$sigma  
  # Number of observations
  n <- nrow(X)
  # Name of dependent variable
  outcome_var_name <- all.vars(mod$terms)[1]
  # List to save results
  results <- vector(mode = "list", length = iter)
  for (i in seq_len(iter)) {
    # draw errors
    errors <- rnorm(n, mean = 0, sd = sigma)
    # create new outcome variable from errors
    y <- X %*% beta + errors
    # replace outcome variable
    .data[[outcome_var_name]] <- y
    # run regression
    newmod <- lm(mod$terms, data = .data)
    # Save model stats as a data frame to the list
    results[[i]] <- glimpse(newmod) %>% mutate(.iter = i)
  }
  # Convert the list of data frames to a single data frame by stacking the iterations
  bind_rows(results)
}
```

### Bootstrap example

Example of a single bootstrap replication. We draw N observations *with replacement* 
from the original data.

```r
nunn_bootstrapped <- bootstrap(nunn, 1)
```

To get bootstrap standard errors, we draw `m` replications, run the regression, 
and save the estimates. 

```r
bootstrap(nunn, 1024) %>%
  do(tidy(lm(trust_neighbors ~ exports, data = nunn)))
```

There are several ways to calculate standard errors from bootstraped replications.

- Calculate the standard error from these simulations by taking the standard deviation of the estimates.
- Calculate the confidence interval using the 2.5% and 97.5% quantiles in the replications

However, in the bootstrap, we should draw the bootstrap samples the same way the sample
was drawn from the population. Why might this not be the case in what we just did? 

## Multiple comparisons and F-test


```r
noise <- data.frame(matrix(rnorm(2100), nrow = 100, ncol = 21))
summary(lm(noise))
```

- How many variables have t-tests that are significant?
- Is the F-test significant? 
- Explain the difference
