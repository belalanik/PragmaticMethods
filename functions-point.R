function.point <- function(data = data, Z = Z, A = A, vars = vars, method = "ITT"){
  
  # ITT
  if(method == "ITT"){
    fit <- glm(Y ~ Z, data = data, family = binomial("logit"))
  }
  
  # Naive PP
  if(method == "Naive-PP"){
    fit <- glm(Y ~ A, data = data, family = binomial("logit"), subset=(Z==A))
  }
  
  # Naive AT
  if(method == "Naive-AT"){
    fit <- glm(Y ~ A, data = data, family = binomial("logit"))
  }
  
  # Baseline-adjusted ITT
  if(method == "B-adj-ITT"){
    fit <- glm(formula(paste("Y ~ Z +", paste(vars, collapse=" + "))), 
                data = data, family = binomial("logit"))
  }
  
  # Baseline-adj PP
  if(method == "B-adj-PP"){
    fit <- glm(formula(paste("Y ~ A +", paste(vars, collapse=" + "))), 
                data = data, family = binomial("logit"), subset=(Z==A))
  }
  
  # Baseline-adj AT
  if(method == "B-adj-AT"){
    fit <- glm(formula(paste("Y ~ A +", paste(vars, collapse=" + "))), 
               data = data, family = binomial("logit"))
  }
  
  # IPW-PP UnStabilized 
  if(method == "IPW-PP-Un"){
    dat1 <- data[data$Z==data$A,]
    ps.model <- glm(formula(paste("A ~ ", paste(vars, collapse=" + "))), 
                    data = dat1, family = "binomial")
    dat1$ps <- predict(ps.model, type = "response")
    wgt.un <- with(dat1, ifelse(A == 1, 1/ps, 1/(1-ps)))
    fit <- glm(Y ~ A, weights = wgt.un, data = dat1, family = quasibinomial("logit"))
  }
  
  # IPW-PP Stabilized 
  if(method == "IPW-PP-S"){
    dat1 <- data[data$Z==data$A,]
    ps.model <- glm(formula(paste("A ~ ", paste(vars, collapse=" + "))), 
                     data = dat1, family = "binomial")
    dat1$ps <- predict(ps.model, type = "response")
    wgt <- with(dat1, ifelse(A == 1, mean(A)/ps, (1-mean(A))/(1-ps)))
    fit <- glm(Y ~ A, weights = wgt, data = dat1, family = quasibinomial("logit"))
  }
  
  # IPW-AT UnStabilized 
  if(method == "IPW-AT-Un"){
    dat1 <- data
    ps.model <-  glm(formula(paste("A ~ ", paste(vars, collapse=" + "))), 
                     data = dat1, family = "binomial")
    dat1$ps <- predict(ps.model, type = "response")
    wgt.un <- with(dat1, ifelse(A == 1, 1/ps, 1/(1-ps)))
    fit <- glm(Y ~ A, weights = wgt.un, data = dat1, family = quasibinomial("logit"))
  }
  
  # IPW-AT Stabilized 
  if(method == "IPW-AT-S"){
    dat1 <- data
    ps.model <-  glm(formula(paste("A ~ ", paste(vars, collapse=" + "))), 
                     data = dat1, family = "binomial")
    dat1$ps <- predict(ps.model, type = "response")
    wgt <- with(dat1, ifelse(A == 1, mean(A)/ps, (1-mean(A))/(1-ps)))
    fit <- glm(Y ~ A, weights = wgt, data = dat1, family = quasibinomial("logit"))
  }
  
  # 2SLS
  if(method == "2SLS"){
    s1 <- glm(formula(paste("A ~ Z + ", paste(vars, collapse=" + "))), 
              data = data, family = binomial)
    data$pr.A <- predict(s1, type = "response")
    fit <- glm(formula(paste("Y ~ pr.A + ", paste(vars, collapse=" + "))), 
               data = data, family = binomial("logit"))
  }
  
  # 2SRI 
  if(method == "2SRI"){
    s1 <- glm(formula(paste("A ~ Z + ", paste(vars, collapse=" + "))), 
              data = data, family = binomial)
    data$pr.A <- predict(s1, type = "response")
    data$e <- data$A - data$pr.A
    fit <- glm(formula(paste("Y ~ A + e + ", paste(vars, collapse=" + "))), 
               data = data, family = binomial("logit"))
  }
  
  # NPCB
  if(method == "NPCB"){
    NPCB <- function(rawData){
      xt <- xtabs(~ X + Y + Z, data = rawData)
      p <- prop.table(xt, margin = 3)
      bpres <- bpbounds(p)
      sbp <- summary(bpres)
      
      #If the inequality constraint is not satisfied, bounds will not be calculated
      ace.LB = NA
      ace.UB = NA
      if (sbp$inequality){
        ace.LB = round(sbp$bounds$`Lower bound`[4], 4) 
        ace.UB = round(sbp$bounds$`Upper bound`[4], 4)
      }
      res <- c(LB = ace.LB, UB = ace.UB)
      return(res)
    }
    dat.npcb <- with(data, data.frame(X = A, Y = Y, Z = Z))
    fit <- NPCB(rawData = dat.npcb)
  }
  return(fit)
}