pm <- function(nde, nie){
  result = (exp(nde)*(exp(nie)-1))/(exp(nde)*exp(nie)-1)
  return(result)
}
pm(0.3525641, -0.01972528)
pm(0.311737991, -0.005869267)
