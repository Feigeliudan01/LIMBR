FILES += $(shell for i in 1 2 3; do echo output/simdata/simulated_data_baseline_$$i__jtkout_GammaP.txt; done)
FILES += $(shell for i in 1 2 3; do echo output/simdata/simulated_data_with_noise_$$i__jtkout_GammaP.txt; done)
FILES += $(shell for i in 1 2 3; do for j in 1 2 3 4 5 6 7 8 9 10; do echo output/simdata/denoised_circ_lowess_$$i\_$$j.txt; done; done)

all: ${FILES}
	@echo ${FILES}
output/simdata/simulated_data_baseline_%__jtkout_GammaP.txt: output/simdata/simulated_data_baseline_%.txt
	@echo Running eJTK on $<
	@sed -i -e 's/_[[:digit:]]//g' $<
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/simulated_data_with_noise_%__jtkout_GammaP.txt: output/simdata/simulated_data_with_noise_%.txt
	@echo Running eJTK on $<
	@sed -i -e 's/_[[:digit:]]//g' $<
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_%_1.txt output/simdata/denoised_circ_lowess_%_2.txt output/simdata/denoised_circ_lowess_%_3.txt output/simdata/denoised_circ_lowess_%_4.txt output/simdata/denoised_circ_lowess_%_5.txt output/simdata/denoised_circ_lowess_%_6.txt output/simdata/denoised_circ_lowess_%_7.txt output/simdata/denoised_circ_lowess_%_8.txt output/simdata/denoised_circ_lowess_%_1.txt output/simdata/denoised_circ_lowess_%_9.txt output/simdata/denoised_circ_lowess_%_10.txt :  output/simdata/simulated_data_with_noise_for_sva_%.txt
	@echo --- SVA Normalizing Simulated Data with Noise—
	@sed -i -e 's/_[[:digit:]]//g' $<
	python src/sva_normalize.py -i $< -o $@ -s 25 -p 1000 -a .05 -d c

clean:
	@echo Cleaning up files..
	rm output/simdata/*__jtknull1000.txt
	rm output/simdata/*__jtkout.txt