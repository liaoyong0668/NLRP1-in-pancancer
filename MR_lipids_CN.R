
library(TwoSampleMR)


getwd()


library(vroom)
exp_dat_clumped<-vroom(file.choose("exposure_clumped.gz"))

View(exp_dat_clumped)


outcome_dat<-vroom("outcome.gz",comment = "##")

dat<-harmonise_data(exposure_dat = exp_dat_clumped,
                    outcome_dat = outcome_dat)

dat$EAF2 <- (1 - dat$eaf.exposure)
dat$MAF <- pmin(dat$eaf.exposure, dat$EAF2)
PVEfx <- function(BETA, MAF, SE, N){
  pve <- (2*(BETA^2)*MAF*(1 - MAF))/((2*(BETA^2)*MAF*(1 - MAF)) + ((SE^2)*2*N*MAF*(1 - MAF)))
  return(pve) 
}
dat$PVE <- mapply(PVEfx, dat$beta.exposure, dat$MAF, dat$se.exposure, N = dat$samplesize.exposure)
dat$FSTAT <- ((dat$samplesize.exposure - 1 - 1)/1)*(dat$PVE/(1 - dat$PVE))


mr(dat)

generate_odds_ratios(mr_res = mr(dat))


mr(dat,method_list = "mr_ivw")

mr_heterogeneity(dat)

mr_pleiotropy_test(dat)

steiger_data(dat)




