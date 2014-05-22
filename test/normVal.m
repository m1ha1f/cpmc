function new_val = normVal(val, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA)
      a = 0.5 / (SIGMA * SIGMA);
      val = single(val);
      new_val = (CONTRAST_SENSITIVE_WEIGHT*exp(-(val)*a)) + POTTS_WEIGHT+ 0.007;           
      duh = hist(new_val);
      if(duh(2)>duh(1))
          a = a + 0.07;
          new_val = (CONTRAST_SENSITIVE_WEIGHT*exp(-(val)*a)) + POTTS_WEIGHT+ 0.007;  
      end
    end