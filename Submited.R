## Code for Dream Challenge mCRPC ###
require('randomForestSRC')
require('survival')
require('caret')
require('MASS')
require('ROCR')

######## Data Process ########
set.seed(0)
trainCRPC <- read.csv('CoreTable_training.csv', na.string = c(' ', '.',''))
testCRPC1 <- read.csv('CoreTable_validation.csv', na.string = c(' ', '.',''))
testCRPC2 <-  read.csv('CoreTable_leaderboard.csv', na.string = c(' ', '.',''))
na.fill <- function(mydata, na.col = 1:ncol(mydata), fill = "NO"){
  # Fill na with certain rules.
  # 
  for (col.idx in na.col)
  {
    mydata[, col.idx] <- factor(mydata[, col.idx], 
                                levels = c(levels(mydata[, col.idx]), fill))
    mydata[is.na(mydata[, col.idx]), col.idx] <- fill   
  }
  mydata
}
trainCRPC <- na.fill(trainCRPC, na.col = c(5,55:131), fill="NO")
death <- rep(1,1600)
death[trainCRPC$DEATH == 'NO'] = 0
trainCRPC$DEATH = death


# apply coverts data.frame into character matrix, use lapply or sapply
one.level <- which(sapply(X = trainCRPC,  function(x) length(levels(x))) == 1)
CRPC.Q1 <- trainCRPC[, c(-1,-2,-3,-6,-7,-8,-9,-10,-11,-18,-20, -one.level,
                         -which(colSums(is.na(trainCRPC)) / 1600 > 0.5) )]

CRPC.Q2 <- trainCRPC[, c(-1,-2,-3,-4,-5,-7,-9,-10,-11,-12, -18,-20,-one.level,
                         -which(colSums(is.na(trainCRPC)) / 1600 > 0.5) )]
CRPC.Q2 <- CRPC.Q2[!is.na(CRPC.Q2$DISCONT), ]   # remove patient without responsible varaible 
CRPC.test1 <- testCRPC1[, colnames(CRPC.Q1)] 
CRPC.test2 <- testCRPC2[, colnames(CRPC.Q1)]
CRPC.test1 <- na.fill(CRPC.test1, na.col = c(5,30:102), fill="NO")
CRPC.test2 <- na.fill(CRPC.test2, na.col = c(5,30:102), fill="NO")
# remove column 1 and 2
CRPC.all <- rbind(CRPC.Q1)
CRPC.test1 <- CRPC.test1[, c(-1,-2)]
CRPC.test2 <- CRPC.test2[, c(-1,-2)]
# impute data
CRPC.Q1.rfimpute <- impute.rfsrc(Surv(LKADT_P, DEATH) ~ ., data = CRPC.Q1,
                                 ntree = 500, 
                                 na.action = "na.impute",
                                 nimpute = 5)
CRPC.Q2.rfimpute <- impute.rfsrc(Surv(ENTRT_PC, DISCONT) ~ ., data = CRPC.Q2,
                                 ntree = 500, 
                                 na.action = "na.impute",
                                 nimpute = 5)
CRPC.test <- rbind(CRPC.test1, CRPC.test2)
CRPC.test.rfimpute <- impute.rfsrc(data = CRPC.test, ntree = 500,
                                    na.action = 'na.impute',
                                   nimpute = 5)
CRPC.validation.rfimpute <- CRPC.test.rfimpute[1:313,]

######## Cross Validation Tools and Help Function ########
cv.rf <- function(k = 3, ntree = 100, nodesize = 3, mtry = NULL, 
                  data = CRPC.Q1.rfimpute, xvar.wt = NULL, myformula = formula(Surv(LKADT_P, DEATH) ~ .))
{
  # k-fold cross validation
  # use default parameters of RSF model.
  
  if(is.null(mtry)) mtry = floor(sqrt(ncol(data)))
  result <- data.frame(idx = 0, ntree = 0, nodesize, mtry = mtry,
                       train.err = 0, test.err= 0)
  
  folds <- createFolds(1:nrow(data), k)
  for (fold.idx in 1:k)
  {
    data.test <- data[folds[[fold.idx]], ]
    data.train <- data[-folds[[fold.idx]], ]
    rsf <- rfsrc(myformula, data = data.train,
                 ntree = ntree,
                 mtry = mtry, 
                 nodesize = nodesize,
                 importance = 'none', xvar.wt = xvar.wt)
    prediction.rsf <- predict(rsf, data.test) 
    result <- rbind(result, c(fold.idx, ntree, nodesize, mtry,
                              rsf$err.rate[ntree], 
                              prediction.rsf$err.rate[ntree]))
  }
  result[-1, ]
}

tune.cv <- function(k = 3, ntree.range = 100, nodesize.range = 3, mtry.range = 36, 
                    myformula = formula(Surv(LKADT_P, DEATH) ~ .), mydata = CRPC.Q1.rfimpute)
{
  # Tuning prameters use cross-validation.
  # Considered parameters are ntree, mtry, nodesize 
  result.lst <- list()
  for (ntree in ntree.range)
  {
    for (mtry in mtry.range)
    {
      for(nodesize in nodesize.range)
      {
        kfold.df <- cv.rf(k= k, ntree = ntree, nodesize = nodesize, data = mydata, myformula = myformula)
        result.lst <- c(result.lst, list(kfold.df))
      }
    }
  }
  result.lst
}

mut.kfold <- function(k = 3, n = 20, ntree.range = 100, mtry.range = 36, nodesize.range = 3, 
                      myformula = formula(Surv(LKADT_P, DEATH) ~ .), data = CRPC.Q1.rfimpute)
{
  # Do n tiral of k-fold cross validation. 
  result.lst <- c()
  for (i in 1:n)
  {
    ntree.kfold <- tune.cv(k = k, ntree.range = ntree.range, nodesize.range = 3, mtry.range = mtry.range, 
                           myformula = myformula, mydata = data)
    result.lst <- c(result.lst, list(ntree.kfold))
  }
  result.lst
}


######## Q1 a ########
### Use Cox model to select variable and decide weights. ####
#### Help Function ####
convert.to.num <- function(mydata){
  # covert a data frame to 0,1 .
  # 'Y' or 'Yes' to 1
  res <- matrix(FALSE, nrow = nrow(mydata), ncol = ncol(mydata))
  res <- (mydata == 'Y') + (mydata == 'Yes') + res
  as.data.frame(res)
}




cox.one <- function(data = CRPC.Q1.rfimpute, alpha = .1, myformula = formula(Surv(LKADT_P, DEATH) ~ .))
{
  # Build cox model fol every variable and return the value is smaller than .05 
  res <- data.frame(variable = NULL, pvalue = NULL, concordance = NULL)
  for (var.name in colnames(data)[c(-1,-2)])
  {
    newdata <- cbind(data[, c(1,2)], data[, var.name])
    colnames(newdata)[3] <- var.name
    cox1 <- coxph(formula = myformula, data = newdata, x = TRUE)    
    sum.cox <- summary(cox1)
    if (sum.cox$logtest[3] < alpha)
    {
      res <- rbind(res,data.frame(name = var.name, pvalue = sum.cox$logtest[3], 
                              concordance = sum.cox$concordance[1]))
    }
  }
  res
}
newdata <- CRPC.Q1.rfimpute[, c(1,2,4)]

#### Subset data and create new features ##### 
visceral <- c('RECTAL','KIDNEYS','LUNGS','LIVER','PLEURA','PROSTATE','ADRENAL','BLADDER','PERITONEUM','COLON','STOMACH','PANCREAS','THYROID')
visceral <- visceral[visceral %in% colnames(CRPC.lesion)]
CRPC.dmg <- CRPC.Q1.rfimpute[, 3:10]
CRPC.LDH.upper <- CRPC.Q1.rfimpute$LDH > 200
CRPC.lab <- CRPC.Q1.rfimpute[, 11:29]
CRPC.lesion <- CRPC.Q1.rfimpute[,30:46]
CRPC.surgery <- CRPC.Q1.rfimpute[, 47:65]
CRPC.medical <- CRPC.Q1.rfimpute[, 66:102]
# CRPC.lesion.visceral <- rowSums(convert.to.num(CRPC.lesion[, visceral]))
CRPC.lymph.only <-  (convert.to.num(CRPC.lesion)$LYM == 1 & 
                     rowSums(convert.to.num(CRPC.lesion)) <= 1)
CRPC.bone.lym <-  (convert.to.num(CRPC.lesion)$LYM == 1 & 
                     convert.to.num(CRPC.lesion)$BONE == 1 )
CRPC.lesion.other <- CRPC.lesion[, setdiff(colnames(CRPC.lesion), visceral)]


######    Step Wise  build model ##############
# Use one variable cox model to select variables. 
# model should have more than 10 variables.
# 5 fold cross validation to avoid overfitting.

step.wise <- function(fulldata, k = 5, ntree = 500, weight_fun = NULL, var.set = NULL)
{
  # this function is specific for Q1.
  # Backward stepwise.
  step.path <- data.frame(idx = NULL, train.err.mean = NULL, test.err.mean = NULL, train.sd = NULL, test.sd = NULL)
  if(is.null(var.set)) var.set <- colnames(full.data)[c(-1,-2)]
  mydata <- fulldata
  idx = 1
  var.set.list <- list()
  while( length(var.set) > 9)
  {
    var.set.list <- c(var.set.list, list(var.set)) 
    model.res <- cv.rf(k = k, data = mydata, ntree = ntree)
    step.path <- rbind(step.path, data.frame(idx = idx, train.err.mean = mean(model.res[, 5]), 
                                             test.err.mean = mean(model.res[, 6]),
	                                           train.err.sd = sd(model.res[, 5]), 
                                             test.err.sd = sd(model.res[, 6])))
	  var.set <- setdiff(var.set, best.to.delete(fulldata, k = k, ntree = ntree, var.set)[[1]] )
    mydata <- cbind(fulldata[, c(1,2)], fulldata[, var.set])
    idx <- idx + 1
  }
  list(path = step.path, var.list = var.set.list)
}

best.to.delete <- function(fulldata, k = k, ntree = ntree, var.set, myformula = formula(Surv(LKADT_P, DEATH) ~ .)) 
{
  res.df <- data.frame(var = NULL, test.err = NULL)
  for (var.candidate in var.set)
  {   
	  mydata <- cbind(fulldata[, c(1,2)], fulldata[, setdiff(var.set, var.candidate)])
	  res <- cv.rf(k = k, data = mydata, ntree = ntree,myformula = myformula)
	  res.df <- rbind(res.df, data.frame(var = var.candidate, test.err = mean(res[, 6])))
  }
  list(as.character(res.df$var[which.min(res.df$test.err)]), res.df)
}
## Select features using backward step-wise elimination
# var.set <- cox.one(data = data, alpha = 0.2)$name
# res <- step.wise(data, k = 5, ntree =500, var.set = var.set)
# var.chosen.Q1 <- res[[2]][[8]]  # No new feature selected
# To reduce running time the process of step-wise elimination has been commented.
var.chosen.Q1 <- c('AGEGRP2','BMI','WGTBLCAT','REGION_C','ECOG_C','ALP','ALT',
                   'AST','CA','CREAT','HB','LDH','NEU','PLT','WBC','ALB','TARGET',
                   'LYMPH_NODES','PLEURA','ADRENAL','PROSTATECTOMY','TURP','ANALGESICS',
                   'GONADOTROPIN','ESTROGENS','CHF','MI','MHCARD','MHGEN','MHNEOPLA',
                   'MHPSYCH','MHRESP','MHVASC') # selected from step-wise process
Q1a.data <- cbind(CRPC.Q1.rfimpute[, c(1, 2)], CRPC.Q1.rfimpute[, var.chosen.Q1])
Q1a.data <- rbind(Q1a.data, cbind(LKADT_P = NA, DEATH = NA, CRPC.validation.rfimpute[, var.chosen.Q1]))
Q1a.rf <- rfsrc(Surv(LKADT_P, DEATH) ~ ., data = Q1a.data[1:1600, ], ntree = 500)
Q1a.predict <- predict(Q1a.rf, Q1a.data[1601:1913, c(-1,-2)])

# save Q1a resutl
risk <- Q1a.predict$predicted / max(Q1a.predict$predicted)
Q1a.submission <- data.frame(RPT = testCRPC1$RPT,
                             riskScoreGlobal = risk,
                             riskScore12 = risk,
                             riskScore18 = risk,
                             riskScore24 = risk)
write.csv(Q1a.submission, 'finalq1a.csv', row.names = F)

#######   Q1b Exact Death Time ######
Q1b.data = CRPC.Q1.rfimpute[CRPC.Q1.rfimpute$DEATH == 1, ]
Q1b.data <- Q1b.data[, -2]

#cv.rf(k = 5, ntree = 500, data =Q1b.data, myformula = formula(LKADT_P ~ .))
lm.model <- glm(formula = LKADT_P ~ MHPSYCH, data = Q1b.data)
summary(lm.model)$coefficients[2,4]
linear_rg <- function(data = Q1b.data, alpha = 0.2, myformula = formula(LKADT_P ~ .)){
  # one variable regression to select variables
  p.vector <- rep(0, ncol(data))
  for (idx in 2:ncol(data))
  {
    mydata <- data[, c(1,idx)]
    mymodel <- glm(myformula, data = mydata)
    p.vector[idx] <- summary(mymodel)$coefficients[2,4]
  }
  p.vector < alpha
}

Q1b.all <- rbind(Q1b.data, cbind(LKADT_P = NA, CRPC.test[1:313,])) 
Q1b.all <- Q1b.all[, linear_rg(alpha = 0.2)]
Q1b.rf <- rfsrc(LKADT_P ~., Q1b.all[1:663, ])
Q1b.predict <- predict(Q1b.rf, Q1b.all[664:976, -1], na.action = 'na.impute')
Q1b.submission <- data.frame(RPT = as.character(testCRPC1$RPT),
                            TIMETOEVENT = Q1b.predict$predicted)
# write csv file for Q1b
write.csv(Q1b.submission, 'finalq1b.csv', row.names = F)

#### Q2    ########
# var.set.q2 <- cox.one(myformula = formula(Surv(ENTRT_PC, DISCONT) ~ .), data = CRPC.Q2.rfimpute, alpha = 0.4)$name
# res <- step.wise(data, k = 5, ntree =500, var.set = var.set, myformula = formula(Surv(ENTRT_PC, DISCONT) ~ .))
# var.set.q2 <- res[[2]][[12]]
var.set.q2 <- c('AGEGRP2','REGION_C','ECOG_C','ALP','AST','CREAT','HB','LDH',
                'PSA','NA.','PHOS','ALB','CCRC','GLU','NON_TARGET','TARGET',
                'RECTAL','LYMPH_NODES','LUNGS','PLEURA','PROSTATE','ADRENAL',
                'BLADDER','ORCHIDECTOMY','PROSTATECTOMY','LYMPHADENECTOMY',
                'SPINAL_CORD_SURGERY','BILATERAL_ORCHIDECTOMY','PRIOR_RADIOTHERAPY',
                'ANALGESICS','ANTI_ANDROGENS','CORTICOSTEROID','BETA_BLOCKING',
                'ESTROGENS','ARTTHROM','PULMEMB','PATHFRAC','SPINCOMP','COPD',
                'MHBLOOD','MHEAR','MHENDO','MHEYE','MHGASTRO','MHGEN','MHHEPATO',
                'MHINVEST','MHMETAB','MHMUSCLE','MHNEOPLA','MHPSYCH','MHRENAL',
                'MHSOCIAL','MHSURG','MHVASC')
# choose the threshold to maxmize the F1 score. 
# res <- c()
# for (k in 1:10)
# {
# Q2.data <- cbind(CRPC.Q2.rfimpute[, c(1,2)], CRPC.Q2.rfimpute[ , var.set.q2])
# sample.idx <- sample(1:1489, 300)
# Q2.train <- Q2.data[-sample.idx, ]
# Q2.test <- Q2.data[sample.idx, ]
# Q2.rf <- rfsrc(Surv(ENTRT_PC, DISCONT) ~ ., data = Q2.train, importance = 'none')
# Q2.predict <- predict(Q2.rf, Q2.test)
# rf.prediction <- prediction(predictions = Q2.predict$predicted, labels = Q2.test$DISCONT)
# rf.pf <- performance(rf.prediction, measure="f")
# quant <- sum(slot(rf.pf, 'x.values')[[1]] < (slot(rf.pf, 'x.values')[[1]])[which.max(slot(rf.pf, 'y.values')[[1]])], na.omit = T) / 300
# res <- c(res, quant)
# }
# mean(res)
threshold.q <- 0.563
# Build Q2 model
Q2.data <- cbind(CRPC.Q2.rfimpute[, c(1,2)], CRPC.Q2.rfimpute[ , var.set.q2])
Q2.all <- rbind(Q2.data, cbind(DISCONT = NA, ENTRT_PC = NA, CRPC.test[,colnames(Q2.data)[c(-1,-2)]]))
Q2.rf <- rfsrc(Surv(ENTRT_PC, DISCONT) ~ ., data = Q2.all[1:1489, ])
Q2.predict <- predict(Q2.rf, Q2.all[1490:1959, c(-1,-2)], na.action = 'na.impute')
threshold <- quantile(Q2.predict$predicted, prob = threshold.q)
Q2.submission <- data.frame(RPT = c(as.character(testCRPC1$RPT), as.character(testCRPC2$RPT)), RISK = Q2.predict$predicted, 
                            DISCONT = as.numeric(Q2.predict$predicted < threshold))
# write the csv file
write.csv(Q2.submission, 'finalq2.csv', row.names = F)
