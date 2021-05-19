setwd('/Users/godongyoung/Dropbox/MyFiles/Research/quantile_regression/Rcode')
source('fn_wo_ME.R')
source('fn_w_ME.R')

load_data  = function(){
  data = read.csv('../Data/income_short.csv',skip = 1,header = F)
  colnames(data) = c('year','index','asset','salary_income','property_income') #조사연도,가구고유번호,자산,근로소득(보완),재산소득(보완)

  y = data$asset
  W1 = data$property_income
  W2 = data$salary_income
  
  # truncate_data--------------------------------------------------------
  inval_idx = (W2==0)|(W1==0)
  y = y[!inval_idx]
  W1 = W1[!inval_idx]
  W2 = W2[!inval_idx]
  
  inval_W2_idx = (W2 < quantile(W2,probs = c(0.001))) | (W2 > quantile(W2,probs = c(0.999))) # or 0.999 might be better solution
  inval_W1_idx = (W1 < quantile(W1,probs = c(0.001))) | (W1 > quantile(W1,probs = c(0.999))) 
  inval_y_idx = (y < quantile(y,probs = c(0.001))) | (y > quantile(y,probs = c(0.999))) 
  sum(inval_W2_idx)
  y = y[(!inval_W2_idx)&(!inval_W1_idx)&(!inval_y_idx)]
  W1 = W1[(!inval_W2_idx)&(!inval_W1_idx)&(!inval_y_idx)]
  W2 = W2[(!inval_W2_idx)&(!inval_W1_idx)&(!inval_y_idx)]
  
  # plot(W2,y)  
  
  # preprocess_data--------------------------------------------------------
  W1.const = 5000
  W2.const = 5000
  y.const = 100000
  
  log_W1 = log(W1+W1.const)
  # Knots.transform = seq(min(log_W1),max(log_W1),length.out = 30)
  # hist(W1,nclass=100)
  # abline(v= exp(Knots.transform)-W1.const,col=2,lwd=3)
  
  log_W2 = log(W2+W2.const)
  # Knots.transform = seq(min(log_W2),max(log_W2),length.out = 30)
  # hist(W2,nclass=100)
  # abline(v= exp(Knots.transform)-W2.const,col=2,lwd=3)
  
  log_y = log(y+y.const)
  # Knots.transform = seq(min(log_y),max(log_y),length.out = 30)
  # hist(y,nclass=100)
  # abline(v= exp(Knots.transform)-y.const,col=2,lwd=3)
  
  res_list = list()
  res_list[['W1']] = W1
  res_list[['W2']] = W2
  res_list[['y']] = y
  
  res_list[['log_W1']] = log_W1
  res_list[['log_W2']] = log_W2
  res_list[['log_y']] = log_y
  
  res_list[['log_y']] = log_y
  res_list[['log_y']] = log_y
  res_list[['log_y']] = log_y
  
  res_list[['W1.const']] = W1.const
  res_list[['W2.const']] = W2.const
  res_list[['y.const']] = y.const
  
  Knots.driect = exp(seq(min(log_W2),max(log_W2),length.out = 30))-W2.const
  
  res_list[['Knots.driect']] = Knots.driect
  return(res_list)
}


