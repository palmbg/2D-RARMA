simu.2drarma = function (alpha = NA, phi = NA, theta = NA, 
                         p = NA, q = NA, n = NA, k = NA)
{

  source("functions.r")
  ar = p
  ma = q
    
    ###### ARMA model
    if(any(is.na(phi)==F) && any(is.na(theta)==F))
    {
      m = max(p,q)
      
      ynew = matrix(1,ncol=(k+m),nrow=(n+m)) * alpha 
      
      ynew1 = c()
      
      y = matrix(0, ncol=(k+m),nrow=(n+m))
      eta = error = mu = y 
      
      for (i in (m+1):(n+m))
      {
        for (j in (m+1):(k+m))
        {
          xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
          y_new1 = as.matrix(xx[(length(xx)-1):1])
          
          yy = as.vector(t(error[(i-ma):i,(j-ma):j]))
          error_new = as.matrix(yy[(length(yy)-1):1])
          
          eta[i,j]  = alpha + phi%*%y_new1 + theta%*%error_new
          mu[i,j] = exp(eta[i,j])
          y[i,j] = rr(mu[i,j])
          ynew[i,j] = log(y[i,j])
          
          error[i,j] = ynew[i,j]-log(mu[i,j])
        }
        
      }
      
      y = y[(m+1):(n+m),(m+1):(k+m)]
      
      return(y)
    } # end ARMA model
    
    
    ###### AR model
    if(any(is.na(phi)==F) && any(is.na(theta)==T))
    {
      
      m = p
      
      ynew  =  matrix(1,ncol=(k+m),nrow=(n+m)) * alpha 
      
      ynew1 = c()
      
      y = matrix(0, ncol=(k+m),nrow=(n+m))
      eta = mu = y 
      
      for (i in (m+1):(n+m))
      {
        for (j in (m+1):(k+m))
        {
          xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
          y_new1 = as.matrix(xx[(length(xx)-1):1])
          
          eta[i,j]  = alpha + phi%*%y_new1 
          mu[i,j] = exp(eta[i,j])
          y[i,j] = rr(mu[i,j])
          ynew[i,j] = log(y[i,j])
          
        }
        
      }
      
      y = y[(m+1):(n+m),(m+1):(k+m)]
      
      return(y)
      
    } # end AR model
  
    ###### MA model
    if(any(is.na(phi)==T) && any(is.na(theta)==F))
    {
      m = q
      
      ynew  = matrix(1,ncol=(k+m),nrow=(n+m)) * alpha 
      
      y = matrix(0, ncol=(k+m),nrow=(n+m))
      eta = error = mu = y 
      
      for (i in (m+1):(n+m))
      {
        for (j in (m+1):(k+m))
        {
          
          yy = as.vector(t(error[(i-ma):i,(j-ma):j]))
          error_new = as.matrix(yy[(length(yy)-1):1])

          eta[i,j]  = alpha + theta%*%error_new
          mu[i,j] = exp(eta[i,j])
          y[i,j] = rr(mu[i,j])
          ynew[i,j] = log(y[i,j])
          
          error[i,j] = ynew[i,j]-log(mu[i,j])
        }
        
      }
      
      y = y[(m+1):(n+m),(m+1):(k+m)]
      
      return(y)
      
    
    } # end MA model  

    
  
  
  
}

