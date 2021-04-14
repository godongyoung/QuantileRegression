echo `date`
Rscript NQR_w_MME_wGamma.R 1 50 &
Rscript NQR_w_MME_wGamma.R 51 100 &
Rscript NQR_w_MME_wGamma.R 101 150 &
Rscript NQR_w_MME_wGamma.R 151 200 &
Rscript NQR_w_MME_wGamma.R 201 250 &
Rscript NQR_w_MME_wGamma.R 251 300 &
Rscript NQR_w_MME_wGamma.R 301 350 &
Rscript NQR_w_MME_wGamma.R 351 400 &
Rscript NQR_w_MME_wGamma.R 401 450 &
Rscript NQR_w_MME_wGamma.R 451 500 &
wait # do not return before background tasks are complete
echo `date`
