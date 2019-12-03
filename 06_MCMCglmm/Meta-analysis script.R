
### Load packages

library("MCMCglmm") # for meta-analysis
library("dplyr") # for data manipulation

migrationdata <- read.csv(~"migration_metadata.csv", header = T) # import dataset
View(migrationdata) # have a look at the dataset. Check out the Predictor variable. There are two, time and temperature.

### Create dataset

migrationdata %>%
filter(Predictor == "year") -> migrationtime # this reduces the dataset to one predictor variable, time.

### Plot data

plot(migrationtime$Slope, I(1/migrationtime$SE)) # this makes the funnel plot.

### First model

randomtest <- MCMCglmm (Slope~1, random = ~Species+Location+Study, data = migrationtime)
summary(randomtest)

### Checking for significance

# Plot the posterior distribution as a histogram to check for significance and whether it's been well estimated or not
# Variance cannot be zero, and therefore if the mean value is pushed up against zero your effect is not significant
# The larger the spread of the histogram, the less well estimated the distribution is.

par(mfrow = c(1,3))

hist(mcmc(randomtest$VCV)[,"Study"])
hist(mcmc(randomtest$VCV)[,"Location"])
hist(mcmc(randomtest$VCV)[,"Species"])

### Assessing convergence

plot(randomtest$Sol) # Fixed effects
plot(randomtest$VCV) # Random effects

###### Priors ######

### Parameter expanded priors, model 1

a <- 1000
prior1<-list(R=list(V=diag(1),nu=0.002)
             , G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)
                      , G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)
                      , G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)))

randomprior <- MCMCglmm (Slope~1, random = ~Species+Location+Study, data = migrationtime, prior=prior1, nitt = 60000)
summary(randomprior)
plot(randomprior$Sol)
plot(randomprior$VCV)

#### Parameter expanded priors, variance of sampling error fixed at 1

prior2<-list(R=list(V=diag(1),nu=0.002)
             , G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)
                      , G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)
                      , G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)
                      , G1=list(V=diag(1),fix=1)))

randomerror2 <- MCMCglmm (Slope~1, random = ~Species+Location+Study+idh(SE):units, data = migrationtime, prior=prior2, nitt = 60000)
summary(randomerror)
plot(randomerror2$VCV)

# Model checks

xsim <- simulate(randomerror2) # reruns 100 new models, based around the same variance/covariance structures but with simulated data.

plot(migrationtime$Slope, I(1/migrationtime$SE))
points(xsim, I(1/migrationtime$SE), col = "red") # here you can plot the data from both your simulated and real datasets and compare them

#### Parameter expanded priors, estimating true variance in sampling error

prior3<-list(R=list(V=diag(1),nu=0.002)
             , G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)
                      , G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)
                      , G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)
                      , G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)))

randomerror3 <- MCMCglmm (Slope~1, random = ~Species+Location+Study+idh(SE):units, data = migrationtime, prior=prior3, nitt = 60000)

# Model checks and comparisons

xsim<-simulate(randomerror3)

plot(migrationtime$Slope, I(1/migrationtime$SE))
points(xsim, I(1/migrationtime$SE), col = "red") # here you can plot the data from both your simulated and real datasets and compare them

# A different way to check your model

xsim<-simulate(randomerror3, 1000) # 1000 represents the number of simulations, and for some reason needs to be higher than the default to work in this case
hist(apply(xsim, 2, max), breaks = 30) # plot your simulation data

abline(v = max(migration$Slope), col = "red") # check to see whether the max value of your real data falls within this histogram.


#### Fixed effects

fixedtest <- MCMCglmm (Slope~Migration_distance+Continent, random = ~Species+Location+Study+idh(SE):units, prior = prior3, data = migrationtime, nitt = 60000)

#### Saving posterior mode in $Sol

fixedtest <- MCMCglmm (Slope~Migration_distance+Continent, random = ~Species+Location+Study+idh(SE):units, data = migrationtime, prior = prior3, pr = TRUE, nitt = 60000)

# Extracting the information

posteriormode <- apply(fixedtest$Sol,2,mode)
names(posteriormode) # identify which posterior modes belong to Species
posteriormode[9:416] # there are a lot of species in this analysis!
sort(posteriormode[9:416]) # sort from smallest to largest, i.e. least responsive to most responsive species

#### Estimating Credible Intervals

HPDinterval(mcmc(fixedtest$Sol[,"(Intercept)"])) # this should look similar to the value in your summary

# particularly useful for combining effects

mean(mcmc(fixedtest$Sol[,"(Intercept)"])+fixedtest$Sol[,"Migration_distanceshort"]+fixedtest$Sol[,"ContinentEurope"])
HPDinterval(mcmc(fixedtest$Sol[,"(Intercept)"])+fixedtest$Sol[,"Migration_distanceshort"]+fixedtest$Sol[,"ContinentEurope"])

#### (Co)variance structures

levels(migrationtime$Response_variable) # check how many levels there are in this variable to account for the heterogeneity of variance in responses

prior4<-list(R=list(V=diag(3),nu=0.002) # changing the matrix here to 3x3 means you can have a separate variance for each response variable
             , G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)
                      , G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)
					  , G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)
                      , G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)))

fixedtest <- MCMCglmm (Slope~Migration_distance+Continent, random = ~Species+Location+Study+idh(SE):units, data = migrationtime, rcov = ~idh(Response_variable):units, prior = prior4, pr = TRUE, nitt = 60000)

####### Now it's your turn #######

#### Now it’s your turn!

# 1. Filter the data by rows which have temperature as the predictor
# 2. Plot the data using a funnel plot
# 3. Run a basic random effects model. Save the posterior mode.
# 4. Plot VCV (random) and Sol (fixed) and check for autocorrelation
# 5. Increase the number of iterations and burn in, check your priors
# 6. Do model checks
# 7. Interpret your model!
# 8. You might want to include some fixed effects, or use different variance structures for your residual as well.

