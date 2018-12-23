TestModel <- function(){
  
  Nsim =200
  for (i in (1:Nsim))
  {  tryCatch(
  
      ## This is what I want to do:
    {
      source("NB.exch.sim.r")
      
      
      source("JGEE1.R")
      mydata = NB.exch.sim()
      print("hello")
      res1 = JGEE1(mydata)
      
      write.table(t(as.matrix(res1)),"res_0501_rho09.csv",row.names=FALSE, 
                  sep=',', col.names=FALSE,append=TRUE)
      
    print(res1)
    
    # source("jointGEE_2.R")
    # res2 = jointGEE_2(mydata)
    # 
    # write.table(t(as.matrix(res2)),"res2.csv",row.names=FALSE, 
    #             sep=',', col.names=FALSE,append=TRUE)
    
    cat("This is simulation number :", i)
    
    
    },
    ## But if an error occurs, do the following: 
    error=function(error_message) {
      message("Yet another error message.")
      message("Here is the actual R error message:")
      message(error_message)
      return(NA)
    }
  )}
  return(res1)
}
