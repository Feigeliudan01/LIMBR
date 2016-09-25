BLINE = $(shell for i in 1 2 3; do echo output/simdata/simulated_data_baseline_$$i.txt; done)
BLINE += $(shell for i in 1 2 3; do echo output/simdata/simulated_data_with_noise_$$i.txt; done)
BLINE += $(shell for i in 1 2 3; do echo output/simdata/simulated_data_with_noise_for_sva_$$i.txt; done)
FILES += $(BLINE)
FILES += $(shell for i in 1 2 3; do echo output/simdata/simulated_data_baseline_$$i\__jtkout_GammaP.txt; done)
FILES += $(shell for i in 1 2 3; do echo output/simdata/simulated_data_with_noise_$$i\__jtkout_GammaP.txt; done)
CSVA := $(shell for i in 1 2 3 4 5 6 7 8 9 10; do echo output/simdata/denoised_circ_lowess_%_$$i.txt; done)
BSVA := $(shell for i in 1 2 3 4 5 6 7 8 9 10; do echo output/simdata/denoised_block_%_$$i.txt; done)
FILES += $(shell for i in 1 2 3; do for j in 1 2 3 4 5 6 7 8 9 10; do echo output/simdata/denoised_circ_lowess_$$i\_$$j\.txt; done; done)
FILES += $(shell for i in 1 2 3; do for j in 1 2 3 4 5 6 7 8 9 10; do echo output/simdata/denoised_block_$$i\_$$j\.txt; done; done)
FILES += $(shell for i in 1 2 3; do for j in 1 2 3 4 5 6 7 8 9 10; do echo output/simdata/denoised_circ_lowess_$$i\_$$j\__jtkout_GammaP.txt; done; done)
FILES += $(shell for i in 1 2 3; do for j in 1 2 3 4 5 6 7 8 9 10; do echo output/simdata/denoised_block_$$i\_$$j\__jtkout_GammaP.txt; done; done)
FILES += $(shell for i in 1 2 3; do echo output/simdata/block_$$i\_pvals.txt; done)
FILES += $(shell for i in 1 2 3; do echo output/simdata/block_$$i\_classes.txt; done)
FILES += $(shell for i in 1 2 3; do echo output/simdata/circ_lowess_$$i\_pvals.txt; done)
FILES += $(shell for i in 1 2 3; do echo output/simdata/circ_lowess_$$i\_classes.txt; done)
FILES += $(shell for i in 1 2 3; do echo output/simdata/circ_lowess_$$i\_classifications.txt; done)
FILES += $(shell for i in 1 2 3; do echo output/simdata/block_$$i\_classifications.txt; done)

all: ${FILES}

$(BLINE) :
	@echo Generating Simulated Data
	@python src/generate_simulated_data.py


output/simdata/simulated_data%__jtkout_GammaP.txt: output/simdata/simulated_data%.txt
	@echo Running eJTK on $<
	@sed -i '' -e 's/_[[:digit:]]//g' $<
	@eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
	@rm output/simdata/*__jtknull1000.txt
	@rm output/simdata/*__jtkout.txt

$(CSVA) :  output/simdata/simulated_data_with_noise_for_sva_%.txt
	@echo Circadian Lowess SVA Normalizing $<
	@for i in 1 2 3 4 5 6 7 8 9 10; do python src/sva_normalize.py -i $< -o $(subst _1.txt,_$$i.txt,$@) -s 25 -p 1000 -a .05 -d c; done

$(BSVA) :  output/simdata/simulated_data_with_noise_for_sva_%.txt
	@echo Circadian Block SVA Normalizing $<
	@for i in 1 2 3 4 5 6 7 8 9 10; do python src/sva_normalize.py -i $< -o $(subst _1.txt,_$$i.txt,$@) -s 25 -p 1000 -a .05 -d b -b output/block_design.p; done

output/simdata/denoised_circ_lowess_%__jtkout_GammaP.txt output/simdata/denoised_block_%__jtkout_GammaP.txt : output/simdata/denoised_circ_lowess_%.txt output/simdata/denoised_block_%.txt
	@echo Running eJTK on $<
	@sed -i '' -e 's/_[[:digit:]]//g' $<
	@sed -i '' -e 's/_[[:digit:]]//g' $(word 2, $^)
	@eJTK-CalcP.py -f $< -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
	@echo Running eJTK on $(word 2, $^)
	@eJTK-CalcP.py -f $(word 2, $^) -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
	@rm output/simdata/*__jtknull1000.txt
	@rm output/simdata/*__jtkout.txt

output/simdata/block_%_pvals.txt : $(shell for i in 1 2 3 4 5 6 7 8 9 10; do echo output/simdata/denoised_block_%_$$i\__jtkout_GammaP.txt; done)
	@tail -q -n +2 $^ |cut -f 19 > $@

output/simdata/block_%_classes.txt : output/simdata/simulated_data_key_%.txt
	@$(shell for i in 1 2 3 4 5 6 7 8 9 10; do tail -q -n +2 $< |cut -f 2 > $@; done)

output/simdata/circ_lowess_%_pvals.txt : $(shell for i in 1 2 3 4 5 6 7 8 9 10; do echo output/simdata/denoised_circ_lowess_%_$$i\__jtkout_GammaP.txt; done)
	@tail -q -n +2 $^ |cut -f 19 > $@

output/simdata/circ_lowess_%_classes.txt : output/simdata/simulated_data_key_%.txt
	@$(shell for i in 1 2 3 4 5 6 7 8 9 10; do tail -q -n +2 $< |cut -f 2 > $@; done)

output/simdata/%_classifications.txt : output/simdata/%_pvals.txt output/simdata/%_classes.txt
	@paste $^ > $@


clean:
	@echo Cleaning up files..
	@rm output/simdata/*.txt
	

test: 
	@echo $(PVALS)

