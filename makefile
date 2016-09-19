all: output/simdata/simulated_data_baseline_1__jtkout_GammaP.txt output/simdata/simulated_data_baseline_2__jtkout_GammaP.txt output/simdata/simulated_data_baseline_3__jtkout_GammaP.txt output/simdata/simulated_data_baseline_4__jtkout_GammaP.txt output/simdata/simulated_data_baseline_5__jtkout_GammaP.txt output/simdata/simulated_data_baseline_6__jtkout_GammaP.txt output/simdata/simulated_data_baseline_7__jtkout_GammaP.txt output/simdata/simulated_data_baseline_8__jtkout_GammaP.txt output/simdata/simulated_data_baseline_9__jtkout_GammaP.txt output/simdata/simulated_data_baseline_10__jtkout_GammaP.txt output/simdata/simulated_data_with_noise_1__jtkout_GammaP.txt output/simdata/simulated_data_with_noise_2__jtkout_GammaP.txt output/simdata/simulated_data_with_noise_3__jtkout_GammaP.txt output/simdata/simulated_data_with_noise_4__jtkout_GammaP.txt output/simdata/simulated_data_with_noise_5__jtkout_GammaP.txt output/simdata/simulated_data_with_noise_6__jtkout_GammaP.txt output/simdata/simulated_data_with_noise_7__jtkout_GammaP.txt output/simdata/simulated_data_with_noise_8__jtkout_GammaP.txt output/simdata/simulated_data_with_noise_9__jtkout_GammaP.txt output/simdata/simulated_data_with_noise_10__jtkout_GammaP.txt  output/simdata/denoised_block_1.txt output/simdata/denoised_block_2.txt output/simdata/denoised_block_3.txt output/simdata/denoised_block_4.txt output/simdata/denoised_block_5.txt output/simdata/denoised_block_6.txt output/simdata/denoised_block_7.txt output/simdata/denoised_block_8.txt output/simdata/denoised_block_9.txt output/simdata/denoised_block_10.txt output/simdata/denoised_circ_lowess_1.txt output/simdata/denoised_circ_lowess_2.txt output/simdata/denoised_circ_lowess_3.txt output/simdata/denoised_circ_lowess_4.txt output/simdata/denoised_circ_lowess_5.txt output/simdata/denoised_circ_lowess_6.txt output/simdata/denoised_circ_lowess_7.txt output/simdata/denoised_circ_lowess_8.txt output/simdata/denoised_circ_lowess_9.txt output/simdata/denoised_circ_lowess_10.txt output/simdata/denoised_circ_lowess_1__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_2__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_3__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_4__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_5__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_6__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_7__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_8__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_9__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_10__jtkout_GammaP.txt output/simdata/denoised_block_1__jtkout_GammaP.txt output/simdata/denoised_block_2__jtkout_GammaP.txt output/simdata/denoised_block_3__jtkout_GammaP.txt output/simdata/denoised_block_4__jtkout_GammaP.txt output/simdata/denoised_block_5__jtkout_GammaP.txt output/simdata/denoised_block_6__jtkout_GammaP.txt output/simdata/denoised_block_7__jtkout_GammaP.txt output/simdata/denoised_block_8__jtkout_GammaP.txt output/simdata/denoised_block_9__jtkout_GammaP.txt output/simdata/denoised_block_10__jtkout_GammaP.txt output/simdata/block_pvals.txt output/simdata/block_classes.txt output/simdata/block_classifications.txt output/simdata/circ_lowess_pvals.txt output/simdata/circ_lowess_classes.txt output/simdata/circ_lowess_classifications.txt

output/simdata/simulated_data_baseline_1__jtkout_GammaP.txt : output/simdata/simulated_data_baseline_1.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/simulated_data_baseline_2__jtkout_GammaP.txt : output/simdata/simulated_data_baseline_2.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/simulated_data_baseline_3__jtkout_GammaP.txt : output/simdata/simulated_data_baseline_3.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/simulated_data_baseline_4__jtkout_GammaP.txt : output/simdata/simulated_data_baseline_4.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/simulated_data_baseline_5__jtkout_GammaP.txt : output/simdata/simulated_data_baseline_5.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/simulated_data_baseline_6__jtkout_GammaP.txt : output/simdata/simulated_data_baseline_6.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/simulated_data_baseline_7__jtkout_GammaP.txt : output/simdata/simulated_data_baseline_7.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/simulated_data_baseline_8__jtkout_GammaP.txt : output/simdata/simulated_data_baseline_8.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/simulated_data_baseline_9__jtkout_GammaP.txt : output/simdata/simulated_data_baseline_9.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/simulated_data_baseline_10__jtkout_GammaP.txt : output/simdata/simulated_data_baseline_10.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/simulated_data_with_noise_1__jtkout_GammaP.txt : output/simdata/simulated_data_with_noise_1.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_1.txt :  output/simdata/simulated_data_with_noise_for_sva_1.txt
	@echo --- SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_circ_lowess_1.txt -s 25 -p 1000 -a .05 -d c

output/simdata/denoised_block_1.txt :  output/simdata/simulated_data_with_noise_for_sva_1.txt
	@echo --- Block SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_block_1.txt -s 25 -p 1000 -a .05 -d b -b output/block_design.p

output/simdata/simulated_data_with_noise_2__jtkout_GammaP.txt : output/simdata/simulated_data_with_noise_2.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_2.txt : output/simdata/simulated_data_with_noise_for_sva_2.txt
	@echo --- SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_circ_lowess_2.txt -s 25 -p 1000 -a .05 -d c

output/simdata/denoised_block_2.txt :  output/simdata/simulated_data_with_noise_for_sva_2.txt
	@echo --- Block SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_block_2.txt -s 25 -p 1000 -a .05 -d b -b output/block_design.p


output/simdata/simulated_data_with_noise_3__jtkout_GammaP.txt : output/simdata/simulated_data_with_noise_3.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_3.txt : output/simdata/simulated_data_with_noise_for_sva_3.txt
	@echo --- SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_circ_lowess_3.txt -s 25 -p 1000 -a .05 -d c

output/simdata/denoised_block_3.txt :  output/simdata/simulated_data_with_noise_for_sva_3.txt
	@echo --- Block SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_block_3.txt -s 25 -p 1000 -a .05 -d b -b output/block_design.p


output/simdata/simulated_data_with_noise_4__jtkout_GammaP.txt : output/simdata/simulated_data_with_noise_4.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_4.txt : output/simdata/simulated_data_with_noise_for_sva_4.txt
	@echo --- SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_circ_lowess_4.txt -s 25 -p 1000 -a .05 -d c

output/simdata/denoised_block_4.txt :  output/simdata/simulated_data_with_noise_for_sva_4.txt
	@echo --- Block SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_block_4.txt -s 25 -p 1000 -a .05 -d b -b output/block_design.p


output/simdata/simulated_data_with_noise_5__jtkout_GammaP.txt : output/simdata/simulated_data_with_noise_5.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_5.txt : output/simdata/simulated_data_with_noise_for_sva_5.txt
	@echo --- SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_circ_lowess_5.txt -s 25 -p 1000 -a .05 -d c

output/simdata/denoised_block_5.txt :  output/simdata/simulated_data_with_noise_for_sva_5.txt
	@echo --- Block SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_block_5.txt -s 25 -p 1000 -a .05 -d b -b output/block_design.p


output/simdata/simulated_data_with_noise_6__jtkout_GammaP.txt : output/simdata/simulated_data_with_noise_6.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_6.txt : output/simdata/simulated_data_with_noise_for_sva_6.txt
	@echo --- SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_circ_lowess_6.txt -s 25 -p 1000 -a .05 -d c

output/simdata/denoised_block_6.txt :  output/simdata/simulated_data_with_noise_for_sva_6.txt
	@echo --- Block SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_block_6.txt -s 25 -p 1000 -a .05 -d b -b output/block_design.p


output/simdata/simulated_data_with_noise_7__jtkout_GammaP.txt : output/simdata/simulated_data_with_noise_7.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_7.txt : output/simdata/simulated_data_with_noise_for_sva_7.txt
	@echo --- SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_circ_lowess_7.txt -s 25 -p 1000 -a .05 -d c

output/simdata/denoised_block_7.txt :  output/simdata/simulated_data_with_noise_for_sva_7.txt
	@echo --- Block SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_block_7.txt -s 25 -p 1000 -a .05 -d b -b output/block_design.p


output/simdata/simulated_data_with_noise_8__jtkout_GammaP.txt : output/simdata/simulated_data_with_noise_8.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_8.txt : output/simdata/simulated_data_with_noise_for_sva_8.txt
	@echo --- SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_circ_lowess_8.txt -s 25 -p 1000 -a .05 -d c

output/simdata/denoised_block_8.txt :  output/simdata/simulated_data_with_noise_for_sva_8.txt
	@echo --- Block SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_block_8.txt -s 25 -p 1000 -a .05 -d b -b output/block_design.p


output/simdata/simulated_data_with_noise_9__jtkout_GammaP.txt : output/simdata/simulated_data_with_noise_9.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_9.txt : output/simdata/simulated_data_with_noise_for_sva_9.txt
	@echo --- SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_circ_lowess_9.txt -s 25 -p 1000 -a .05 -d c

output/simdata/denoised_block_9.txt :  output/simdata/simulated_data_with_noise_for_sva_9.txt
	@echo --- Block SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_block_9.txt -s 25 -p 1000 -a .05 -d b -b output/block_design.p


output/simdata/simulated_data_with_noise_10__jtkout_GammaP.txt : output/simdata/simulated_data_with_noise_10.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_10.txt : output/simdata/simulated_data_with_noise_for_sva_10.txt
	@echo --- SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_circ_lowess_10.txt -s 25 -p 1000 -a .05 -d c

output/simdata/denoised_block_10.txt :  output/simdata/simulated_data_with_noise_for_sva_10.txt
	@echo --- Block SVA Normalizing Simulated Data with Noise—
	python src/sva_normalize.py -i $< -o output/simdata/denoised_block_10.txt -s 25 -p 1000 -a .05 -d b -b output/block_design.p

output/simdata/denoised_circ_lowess_1__jtkout_GammaP.txt : output/simdata/denoised_circ_lowess_1.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_2__jtkout_GammaP.txt : output/simdata/denoised_circ_lowess_2.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_3__jtkout_GammaP.txt : output/simdata/denoised_circ_lowess_3.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_4__jtkout_GammaP.txt : output/simdata/denoised_circ_lowess_4.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_5__jtkout_GammaP.txt : output/simdata/denoised_circ_lowess_5.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_6__jtkout_GammaP.txt : output/simdata/denoised_circ_lowess_6.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_7__jtkout_GammaP.txt : output/simdata/denoised_circ_lowess_7.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_8__jtkout_GammaP.txt : output/simdata/denoised_circ_lowess_8.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_9__jtkout_GammaP.txt : output/simdata/denoised_circ_lowess_9.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_circ_lowess_10__jtkout_GammaP.txt : output/simdata/denoised_circ_lowess_10.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_block_1__jtkout_GammaP.txt : output/simdata/denoised_block_1.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_block_2__jtkout_GammaP.txt : output/simdata/denoised_block_2.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_block_3__jtkout_GammaP.txt : output/simdata/denoised_block_3.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_block_4__jtkout_GammaP.txt : output/simdata/denoised_block_4.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_block_5__jtkout_GammaP.txt : output/simdata/denoised_block_5.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_block_6__jtkout_GammaP.txt : output/simdata/denoised_block_6.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_block_7__jtkout_GammaP.txt : output/simdata/denoised_block_7.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_block_8__jtkout_GammaP.txt : output/simdata/denoised_block_8.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_block_9__jtkout_GammaP.txt : output/simdata/denoised_block_9.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/denoised_block_10__jtkout_GammaP.txt : output/simdata/denoised_block_10.txt
	@echo --- Running eJTK on Simulated Data ---
	eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

output/simdata/block_pvals.txt : output/simdata/denoised_block_1__jtkout_GammaP.txt output/simdata/denoised_block_2__jtkout_GammaP.txt output/simdata/denoised_block_3__jtkout_GammaP.txt output/simdata/denoised_block_4__jtkout_GammaP.txt output/simdata/denoised_block_5__jtkout_GammaP.txt output/simdata/denoised_block_6__jtkout_GammaP.txt output/simdata/denoised_block_7__jtkout_GammaP.txt output/simdata/denoised_block_8__jtkout_GammaP.txt output/simdata/denoised_block_9__jtkout_GammaP.txt output/simdata/denoised_block_10__jtkout_GammaP.txt
	tail -q -n +2 $^ |cut -f 19 > output/simdata/block_pvals.txt

output/simdata/block_classes.txt : output/simdata/simulated_data_key_1.txt output/simdata/simulated_data_key_2.txt output/simdata/simulated_data_key_3.txt output/simdata/simulated_data_key_4.txt output/simdata/simulated_data_key_5.txt output/simdata/simulated_data_key_6.txt output/simdata/simulated_data_key_7.txt output/simdata/simulated_data_key_8.txt output/simdata/simulated_data_key_9.txt output/simdata/simulated_data_key_10.txt
	tail -q -n +2 $^ |cut -f 2 > output/simdata/block_classes.txt

output/simdata/block_classifications.txt : output/simdata/block_pvals.txt output/simdata/block_classes.txt
	paste $^ > output/simdata/block_classifications.txt

output/simdata/circ_lowess_pvals.txt : output/simdata/denoised_circ_lowess_1__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_2__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_3__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_4__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_5__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_6__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_7__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_8__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_9__jtkout_GammaP.txt output/simdata/denoised_circ_lowess_10__jtkout_GammaP.txt
	 tail -q -n +2 $^ |cut -f 19 > output/simdata/circ_lowess_pvals.txt

output/simdata/circ_lowess_classes.txt : output/simdata/simulated_data_key_1.txt output/simdata/simulated_data_key_2.txt output/simdata/simulated_data_key_3.txt output/simdata/simulated_data_key_4.txt output/simdata/simulated_data_key_5.txt output/simdata/simulated_data_key_6.txt output/simdata/simulated_data_key_7.txt output/simdata/simulated_data_key_8.txt output/simdata/simulated_data_key_9.txt output/simdata/simulated_data_key_10.txt
	tail -q -n +2 $^ |cut -f 2 > output/simdata/circ_lowess_classes.txt

output/simdata/circ_lowess_classifications.txt : output/simdata/circ_lowess_pvals.txt output/simdata/circ_lowess_classes.txt
	paste $^ > output/simdata/circ_lowess_classifications.txt


clean:
	rm output/simdata/*__jtknull1000.txt
	rm output/simdata/*__jtkout.txt