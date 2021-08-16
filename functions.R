# Functions for drawing threshold lines on a volcano plot

# Function to compute the smooth threshold curve
# Four arguments:
# 1 - 
# x = All the values of x from xmin to xmax

# 2 - 
# s0, the fudge factor, according to Tusher, V., Tibshirani, R., and Chu, G. (2001)
# ff0 = fudge2(FC, sdd)$s.zero
# or ff = 0.5

# 3 - 
# ta = the quantile function of the t-distribution at a specific confidence level and degrees of freedom.
# i.e. the critical value for alpha as per qt(confidence_level, df = df)

# 4 - 
# df = (n1 + n2) - 2, degrees of freedom for t test
# Where, e.g.
# n1 = 3
# n2 = 3

smooth.threshold = function(x, ta, s0, df) {
  xp = x[x > (ta*s0)]
  xn = x[x<(-ta*s0)]
  dp = xp/ta-s0
  dn = xn/(-ta)-s0
  dp = s0/dp
  dp = ta*(1+dp)
  dn = s0/dn
  dn = ta*(1+dn)
  fp = pt(dp, df = df)
  fn = pt(dn, df = df)
  yp = -log10(2*(1-fp))
  yn = -log10(2*(1-fn))
  
  return(cbind(
    c(xn,xp),
    c(yn,yp)
  ))
}

# As a tibble, for drawing with ggplot
ggthreshold <- function(x, ta, s0, df) {
  xp = x[x > (ta*s0)]
  xn = x[x<(-ta*s0)]
  dp = xp/ta-s0
  dn = xn/(-ta)-s0
  dp = s0/dp
  dp = ta*(1+dp)
  dn = s0/dn
  dn = ta*(1+dn)
  fp = pt(dp, df = df)
  fn = pt(dn, df = df)
  yp = -log10(2*(1-fp))
  yn = -log10(2*(1-fn))
  
  # return(cbind(
  #   c(xn,xp),
  #   c(yn,yp)
  # ))
  
  return(tibble(x = c(xn,xp),
                y = c(yn,yp)))
}


# As a vector of y coordinates
# As a tibble, for drawing with ggplot
ggthres_func <- function(x, ta, s0, df) {
  xp = x[x > (ta*s0)] # x positive
  # xn = x[x<(-ta*s0)]  # x negative
  
  dp = xp/ta-s0
  dp = s0/dp
  dp = ta*(1+dp)
  fp = pt(dp, df = df)
  yp = -log10(2*(1-fp))

  # dn = xn/(-ta)-s0
  # dn = s0/dn
  # dn = ta*(1+dn)
  # fn = pt(dn, df = df)
  # yn = -log10(2*(1-fn))
  
  
  # return(cbind(
  #   c(xn,xp),
  #   c(yn,yp)
  # ))
  
  return(yp)
}
