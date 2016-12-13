FILES += output/wt_lowess_normed.txt output/csp_lowess_normed.txt output/rna_lowess_normed.txt

all: ${FILES}

output/rna_qnormed.txt output/wt_pool_qnormed.txt output/csp_pool_qnormed.txt :

output/rna_lowess_normed.txt : output/rna_qnormed.txt
	@echo Lowess normalizing $<
	@python src/sva_normalize.py -i $< -s 25 -p 10000 -a .05 -d c

output/%_lowess_normed.txt : output/%_pool_qnormed.txt
	@echo Lowess normalizing $<
	@python src/sva_normalize.py -i $< -s 25 -p 10000 -a .05 -d c