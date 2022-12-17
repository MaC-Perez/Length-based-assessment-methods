library(TropFishR)
graphics.off()
rm(list=ls(all=TRUE))

# 1

data("alba")
alba

# 2 
#length frequency data set for the mollusc Abra alba
#Determine the number of bins to include in your moving average
#TropFishR the default value is 5 but we can specify different values with the argunet MA 
#MA should approximate the number of bins spanning the width covered by the smallest
#(i.e. youngest) cohort
# if we check in the data we could use 7 bins 

# explore different values for MA 
alba <- lfqRestructure(alba, MA=7)
#alba <- lfqRestructure(alba, MA=5)
#alba <- lfqRestructure(alba, MA=3)
plot(alba)


#Identify appropriate limits for vonB parameters to constrain the search
#space in ELEFAN_GA. the lower and upper bounds of Linf and K as an argument when calling the function.
# start with Linf .... use a catch curve and the data 
PW <- powell_wetherall(alba, catch_columns = 1:7)

#$Linf_est
#[1] 11.3413

#$se_Linf
#[1] 0.9659902

#$confidenceInt_Linf
#[1]  8.977609 13.704995

#$ZK
#[1] 0.8569509

#$se_ZK
#[1] 0.03007039


# Compare the Linf and K values you generated with Response Surface Analysis (RSA)
# but also we can use response surface analysis 

#Linf between 9 and 12 and K between 2 and 3

alba2 <- ELEFAN(
  lfq = alba, MA = 7,
  Linf_range = seq(7, 20, length.out = 30),
  K_range = exp(seq(log(0.1),log(4), length.out = 30)),
  method = "cross",
  cross.date = alba$dates[3],
  cross.midLength = alba$midLengths[5],
  contour = TRUE, add.values = FALSE,
  hide.progressbar = TRUE # change to 'TRUE' to follow algorithm's progression
)
points(alba2$par["Linf"], alba2$par["K"], pch="*", cex=2, col=2)
unlist(alba2$par)
alba2$Rn_max
plot(alba2)

# Linf        K t_anchor        C       ts     phiL 
#9.689655 2.400000 0.370000 0.000000 0.000000 2.352828 

#Run ELEFAN_GA using your values for constraining Linf and K from Step 4. Do this by changing the Linf
#and K values in the low_par and up_par arguments.
#Generate growth curve parameter value estimates and plot those curves on your restructured data.

#$confidenceInt_Linf
#[1]  8.977609 13.704995
Linf=PW$confidenceInt_Linf[1]
Linf=PW$confidenceInt_Linf[2]

set.seed(1)
alba4 <- ELEFAN_GA(
  lfq = alba,
  seasonalised = FALSE,
  low_par = list(Linf=PW$confidenceInt_Linf[1], K=1, t_anchor=0, ts=0, C=0),
  up_par = list(Linf=PW$confidenceInt_Linf[2], K=4, t_anchor=1, ts=1, C=1),
  popSize = 60,
  pmutation = 0.2,
  maxiter = 100,
  run = 20,
  MA = 7,
  plot.score = TRUE,
  monitor = FALSE,
  parallel = FALSE
)
unlist(alba4$par)
#Linf          K   t_anchor       phiL 
#10.1301838  1.9618729  0.2643963  2.3039055 
alba4$Rn_max
plot(alba4)


