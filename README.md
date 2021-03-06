# Analysis of Covid-19 serology tests

This code implements the methods in the paper 
[Estimation of COVID-19 Prevalence from Serology Tests: A Partial Identification Approach](https://papers.ssrn.com/sol3/Papers.cfm?abstract_id=3587738) (Journal of Econometrics, 2020, forthcoming).
Compared to alternatives, the main method in this paper is valid in finite samples, and almost exact for all unknown parameters (see below).
As such, it does not require asymptotics or other approximations, like bootstrap and Bayesian methods do.

## Datasets

To load data we use `create_model`:

    m = create_model(1) 

Now `m` contains data from serology study 1 (= Santa Clara). Currently, there are 5 serology datasets in the script:
 * 1 = Santa Clara of [Bendavid et al (2020)](https://www.medrxiv.org/content/10.1101/2020.04.14.20062463v2).
 * 2 = LA county study [http://publichealth.lacounty.gov/phcommon/public/media/mediapubhpdetail.cfm?prid=2328]
 * 3 = Santa Clara + LA combined.
 * 4 = New York Study [https://www.nytimes.com/2020/04/23/nyregion/coronavirus-antibodies-test-ny.html]
 * 7 = SC + LA + NYC studies combined.

Each `model` contains data for `CalibrationStudy` and data for the `MainStudy`. The calibration studies contain results on known true positives and false positives.
The main study contains positive results of unknown nature (i.e., it is unknown how many true positives there are in all tests of the main study). The (unknown) proportion of true positives over the total number of tests in the main study is the (unknown) prevalence, which we want to estimate.

See paper (Section 2) for details.

## Testing a combination of parameters

Each study has three unknown parameters: false positive rate of antibody test (FPR), true positive rate (TPR) and prevalence.
Suppose we want to test whether in the Santa Clara study it is statistically plausible at the 5% level that 
FPR = 1.5%, TPR = 90% and prevalence = 0%. We can test this as follows:

    m = create_model(1)  # load Santa Clara data.
    num_infect = round((0/100) * m$MainStudy$total)  # how many true positives in main study correspond to specified prevalence. 
    theta0 = c(1.5/100, 90/100, num_infect)
    include_theta_Exact(theta0, m, alpha_level=0.05, vis=T)
    
This command will return a list with `incl` as an element for the test result (`TRUE` = the parameter values cannot be rejected).
With `vis=T` the results will be visualized (similar to Fig. 1 in the paper).

The 95% CIs in the paper (e.g., Figs. 3 and 4) can be obtained by running the above code for a grid of possible parameter values (could use a computing grid to do this in parallel).

## Other tests
To apply the likelihood ratio test of Section 4.3.2 in the paper we can run:

    LR_test(theta0, m)

This test is valid in finite samples as the main test described above, but it is much slower. The output is similar to `include_theta_Exact`.

To apply the MCMC-based partial identification methods (Procedure 1 and 3) of [Chen et al (2018)](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA14525) we run:

    Proc1_chen(m, num_mcmc=20000)
    Proc3_chen(m)


