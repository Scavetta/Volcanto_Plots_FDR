# Get data ----
tabl_2 <- LFQRatio2

# Columns A.R1, A.R2 and A.R3 correspond to the (numeric) abundance values of proteins in the three replicates of condition A.
# Columns B.R1, B.R2 and B.R3 correspond to the (numeric) abundance values of proteins in the three replicates of condition B.
# Column Welch.test.pval contains the p-values of the Welch t-test between condition A and condition B computed with the Perseus software.

# Calculate the Fold-Change of the averages of each group ----
tabl_2 %>% 
  as_tibble() %>% 
  mutate(ID = seq_len(nrow(.))) %>% 
  select(ID, matches("^(A|B)")) %>% 
  pivot_longer(-ID,
    names_to = "condition",
    values_to = "abundance") %>% 
  mutate(condition = str_remove(condition, "\\..*$")) %>%
  group_by(ID, condition) %>%
  # na.omit() %>% 
  summarize(n = n(),                                                                                # Ideally, this should be consistent among all tests
            avg = mean(abundance),                                                                  # The mean of the 3 biological replicates in each condition 
            stdev = sd(abundance)) %>%                                                              # The standard deviation of the 3 biological replicates in each condition
  ungroup() %>% 
  group_by(ID) %>%                                                                                  # For each ID:
  summarise(FC = avg[condition == "A"] - avg[condition == "B"] ,                                    # the signal
            sdd = sqrt(((stdev[condition == "A"]^2) + (stdev[condition == "B"]^2) )/ sum(n)),       # the noise
            df = sum(n) - 2,                                                                        # the degrees of freedom
            t.stat = FC/sdd,                                                                        # the test statistic
            p.value = 2 * (1 - pt(abs(t.stat), df = df))                                            # the unadjusted p.value
            ) -> tabl_2_summary
            
# The fudge factor
ff0 = fudge2(tabl_2_summary$FC, tabl_2_summary$sdd)$s.zero

# How to not calculte the x in the function and let ggplot2 define the x?
test2 <- ggthres_func(x = seq(-4, 4, 0.05),     # x value
             ta = qt(0.975, df = 4),   # The critical value at alpha = 0.05
             s0 = ff0,                 # The fudge factor
             df = 4                    # The degrees of freedom (n1 + n2 - 2) 
             )

ggplot() +
  xlim(-5, 5) + 
  geom_function(fun = dnorm, args = list(mean = 2, sd = .5))

ggplot() +
  xlim(-4, 4) +
  geom_function(fun = ggthres_func, args = list(ta = qt(0.975, df = 4),   # The critical value at alpha = 0.05
                                                s0 = ff0,                 # The fudge factor
                                                df = 4                    # The degrees of freedom (n1 + n2 - 2) 
  ))