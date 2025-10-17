rTaylor <- function(rNbar, betapi, betay, inflation, outputgap, meaninflation, meanoutputgap) {
  # ------------------------------------------------------------
  # rTaylor:
  # Computes the Taylor rule implied nominal interest rate
  #
  # Inputs:
  #   rNbar          - natural nominal interest rate
  #   betapi         - inflation response coefficient
  #   betay          - output gap response coefficient
  #   inflation      - current inflation
  #   outputgap      - current output gap
  #   meaninflation  - mean inflation over sample
  #   meanoutputgap  - mean output gap over sample
  #
  # Output:
  #   rval           - model-implied interest rate (non-negative)
  # ------------------------------------------------------------
  
  rval <- max(rNbar + betapi * (inflation - meaninflation) +
                betay * (outputgap - meanoutputgap), 0)
  
  return(rval)
}
