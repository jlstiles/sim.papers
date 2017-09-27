# getting a random sample from a randomly drawn dgp 
dgp = get.dgp(n = 1000, d = 4)

# population blip variance
dgp$BV0
# population average treatment effect
dgp$ATE0
# sample proportion with treatment
mean(dgp$DF$A)
# sample proportion who died
mean(dgp$DF$Y)
# min sample pscore
min(dgp$PGn)
# max sample pscore
max(dgp$PGn)
# sample dataframe
head(dgp$DF)
# histogram of blips
hist(dgp$blip_n)

