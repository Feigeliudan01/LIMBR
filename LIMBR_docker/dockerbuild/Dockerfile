FROM ubuntu:latest AS base
MAINTAINER Alec Crowell <alexander.m.crowell@gmail.com>

# Install all the software needed to run the pipeline
ENV DEBIAN_FRONTEND noninteractive
ENV LANG C.UTF-8
RUN apt-get -qq update
RUN apt-get update && apt-get install -y  software-properties-common curl wget git python python3 nano python-tk python3-pip python-pip python3-tk
RUN pip2 install --upgrade pip
RUN pip3 install --upgrade pip
RUN pip3 install -U numpy cython matplotlib seaborn astroML jinja2 limbr
RUN pip2 install -U numpy cython matplotlib limbr
RUN add-apt-repository ppa:git-core/ppa
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt-get install git-lfs
RUN git lfs install

# clone the most recent version of the pipeline
ARG expire_after=3-28-18
RUN git lfs clone https://github.com/aleccrowell/LIMBR.git
RUN mv LIMBR temp
RUN mv temp/LIMBR_docker LIMBR
RUN rm -rf temp
WORKDIR /LIMBR/
RUN mkdir ./output/
RUN mkdir ./output/actual/
RUN mkdir ./output/simdata/
RUN mkdir ./output/figs/
WORKDIR /LIMBR/src/bin/
RUN python3 setup.py build_ext --inplace
RUN python setup.py build_ext --inplace
WORKDIR /LIMBR/
COPY ./data/Formatted_input.txt /LIMBR/data/Formatted_input.txt
COPY ./data/Jen_rnaseq_formatted_raw_counts.txt /LIMBR/data/Jen_rnaseq_formatted_raw_counts.txt
COPY ./data/wt_classes.csv /LIMBR/data/wt_classes.csv
COPY ./data/csp_classes.csv /LIMBR/data/csp_classes.csv

FROM base AS actual
#make actual data files
RUN python3 src/knn_impute.py -i data/Formatted_input.txt -o output/actual/imputed_peptide.txt -m 0.3 -p 2
RUN python3 src/separate.py
RUN python3 src/sva_normalize.py -i output/actual/wt_for_sva.txt -o output/actual/wt_lowess_normed.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p -j output/actual/wt_pool_design.p
RUN python3 src/sva_normalize.py -i output/actual/csp_for_sva.txt -o output/actual/csp_lowess_normed.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p -j output/actual/csp_pool_design.p
RUN python3 src/sva_normalize.py -i output/actual/rna_for_sva.txt -o output/actual/rna_lowess_normed.txt -s 25 -n 12 -p 1000 -a .05 -d c -e r
RUN python3 src/old_normalize.py -i output/actual/wt_for_sva.txt -o output/actual/wt_old_normed.txt -e p -j output/actual/wt_pool_design.p
RUN python3 src/old_normalize.py -i output/actual/csp_for_sva.txt -o output/actual/csp_old_normed.txt -e p -j output/actual/csp_pool_design.p
RUN python3 src/old_normalize.py -i output/actual/rna_for_sva.txt -o output/actual/rna_old_normed.txt -e r
RUN sed -e 's/_[[:digit:]]//g' output/actual/wt_lowess_normed.txt > temp.txt
RUN cut -f 1,5- temp.txt > output/actual/wt_lowess_normed.txt
RUN python src/eJTK-CalcP.py -f output/actual/wt_lowess_normed.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN sed -e 's/_[[:digit:]]//g' output/actual/csp_lowess_normed.txt > temp.txt
RUN cut -f 1,5- temp.txt > output/actual/csp_lowess_normed.txt
RUN python src/eJTK-CalcP.py -f output/actual/csp_lowess_normed.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN sed -e 's/_[[:digit:]]//g' output/actual/rna_lowess_normed.txt > temp.txt
RUN mv temp.txt output/actual/rna_lowess_normed.txt
RUN python src/eJTK-CalcP.py -f output/actual/rna_lowess_normed.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN sed -e 's/_[[:digit:]]//g' output/actual/wt_old_normed.txt > temp.txt
RUN cut -f 1,5- temp.txt > output/actual/wt_old_normed.txt
RUN python src/eJTK-CalcP.py -f output/actual/wt_old_normed.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN sed -e 's/_[[:digit:]]//g' output/actual/csp_old_normed.txt > temp.txt
RUN cut -f 1,5- temp.txt > output/actual/csp_old_normed.txt
RUN python src/eJTK-CalcP.py -f output/actual/csp_old_normed.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN sed -e 's/_[[:digit:]]//g' output/actual/rna_old_normed.txt > temp.txt
RUN mv temp.txt output/actual/rna_old_normed.txt
RUN python src/eJTK-CalcP.py -f output/actual/rna_old_normed.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN rm output/actual/*__jtknull1000.txt
RUN rm output/actual/*__jtkout.txt

From base AS sim_mb_lowess
COPY --from=actual /LIMBR/output /LIMBR/output
#Simulation testing
RUN echo Generating Simulated Data
RUN python3 src/generate_simulated_data_mb.py

RUN echo Running eJTK
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do sed -e 's/_[[:digit:]]//g' output/simdata/mb_simulated_data_baseline_$i.txt > temp_baseline_$i.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do mv temp_baseline_$i.txt output/simdata/mb_simulated_data_baseline_$i.txt; done
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_1.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_2.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_3.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_4.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_5.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_6.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_7.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_8.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_9.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_10.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_11.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_12.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_13.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_14.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_15.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_16.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_17.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_18.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_19.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_baseline_20.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do sed -e 's/_[[:digit:]]//g' output/simdata/mb_simulated_data_with_noise_$i.txt > temp_with_noise_$i.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do mv temp_with_noise_$i.txt output/simdata/mb_simulated_data_with_noise_$i.txt; done
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_1.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_2.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_3.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_4.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_5.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_6.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_7.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_8.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_9.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_10.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_11.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_12.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_13.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_14.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_15.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_16.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_17.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_18.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_19.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt
RUN python src/eJTK-CalcP.py -f output/simdata/mb_simulated_data_with_noise_20.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt

RUN echo Circadian Lowess SVA Normalizing
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_1.txt -o output/simdata/mb_denoised_circ_lowess_1.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_2.txt -o output/simdata/mb_denoised_circ_lowess_2.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_3.txt -o output/simdata/mb_denoised_circ_lowess_3.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_4.txt -o output/simdata/mb_denoised_circ_lowess_4.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_5.txt -o output/simdata/mb_denoised_circ_lowess_5.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_6.txt -o output/simdata/mb_denoised_circ_lowess_6.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_7.txt -o output/simdata/mb_denoised_circ_lowess_7.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_8.txt -o output/simdata/mb_denoised_circ_lowess_8.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_9.txt -o output/simdata/mb_denoised_circ_lowess_9.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_10.txt -o output/simdata/mb_denoised_circ_lowess_10.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_11.txt -o output/simdata/mb_denoised_circ_lowess_11.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_12.txt -o output/simdata/mb_denoised_circ_lowess_12.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_13.txt -o output/simdata/mb_denoised_circ_lowess_13.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_14.txt -o output/simdata/mb_denoised_circ_lowess_14.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_15.txt -o output/simdata/mb_denoised_circ_lowess_15.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_16.txt -o output/simdata/mb_denoised_circ_lowess_16.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_17.txt -o output/simdata/mb_denoised_circ_lowess_17.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_18.txt -o output/simdata/mb_denoised_circ_lowess_18.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_19.txt -o output/simdata/mb_denoised_circ_lowess_19.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_20.txt -o output/simdata/mb_denoised_circ_lowess_20.txt -s 25 -n 12 -p 1000 -a .05 -d c -e p

From base AS sim_mb_sse
COPY --from=sim_mb_lowess /LIMBR/output /LIMBR/output

RUN echo Circadian Lowess SSE Normalizing
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_1.txt -o output/simdata/mb_denoised_circ_lowess_sse_1.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_2.txt -o output/simdata/mb_denoised_circ_lowess_sse_2.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_3.txt -o output/simdata/mb_denoised_circ_lowess_sse_3.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_4.txt -o output/simdata/mb_denoised_circ_lowess_sse_4.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_5.txt -o output/simdata/mb_denoised_circ_lowess_sse_5.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_6.txt -o output/simdata/mb_denoised_circ_lowess_sse_6.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_7.txt -o output/simdata/mb_denoised_circ_lowess_sse_7.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_8.txt -o output/simdata/mb_denoised_circ_lowess_sse_8.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_9.txt -o output/simdata/mb_denoised_circ_lowess_sse_9.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_10.txt -o output/simdata/mb_denoised_circ_lowess_sse_10.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_11.txt -o output/simdata/mb_denoised_circ_lowess_sse_11.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_12.txt -o output/simdata/mb_denoised_circ_lowess_sse_12.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_13.txt -o output/simdata/mb_denoised_circ_lowess_sse_13.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_14.txt -o output/simdata/mb_denoised_circ_lowess_sse_14.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_15.txt -o output/simdata/mb_denoised_circ_lowess_sse_15.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_16.txt -o output/simdata/mb_denoised_circ_lowess_sse_16.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_17.txt -o output/simdata/mb_denoised_circ_lowess_sse_17.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_18.txt -o output/simdata/mb_denoised_circ_lowess_sse_18.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_19.txt -o output/simdata/mb_denoised_circ_lowess_sse_19.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_20.txt -o output/simdata/mb_denoised_circ_lowess_sse_20.txt -s 25 -n 12 -p 1000 -a .05 -d l -e p

From base AS sim_mb_block
COPY --from=sim_mb_sse /LIMBR/output /LIMBR/output

RUN echo Circadian Block SVA Normalizing
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_1.txt -o output/simdata/mb_denoised_block_1.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_2.txt -o output/simdata/mb_denoised_block_2.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_3.txt -o output/simdata/mb_denoised_block_3.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_4.txt -o output/simdata/mb_denoised_block_4.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_5.txt -o output/simdata/mb_denoised_block_5.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_6.txt -o output/simdata/mb_denoised_block_6.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_7.txt -o output/simdata/mb_denoised_block_7.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_8.txt -o output/simdata/mb_denoised_block_8.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_9.txt -o output/simdata/mb_denoised_block_9.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_10.txt -o output/simdata/mb_denoised_block_10.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_11.txt -o output/simdata/mb_denoised_block_11.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_12.txt -o output/simdata/mb_denoised_block_12.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_13.txt -o output/simdata/mb_denoised_block_13.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_14.txt -o output/simdata/mb_denoised_block_14.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_15.txt -o output/simdata/mb_denoised_block_15.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_16.txt -o output/simdata/mb_denoised_block_16.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_17.txt -o output/simdata/mb_denoised_block_17.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_18.txt -o output/simdata/mb_denoised_block_18.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_19.txt -o output/simdata/mb_denoised_block_19.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p
RUN python3 src/sva_normalize.py -i output/simdata/mb_simulated_data_with_noise_for_sva_20.txt -o output/simdata/mb_denoised_block_20.txt -s 25 -n 12 -p 1000 -a .05 -d b -b data/block_design.p -e p


From base AS ejtk
COPY --from=sim_mb_block /LIMBR/output /LIMBR/output
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do sed -e 's/_[[:digit:]]//g' output/simdata/mb_denoised_block_$i.txt > temp_block_$i.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do sed -e 's/_[[:digit:]]//g' output/simdata/mb_denoised_circ_lowess_sse_$i.txt > temp_circ_lowess_sse_$i.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do sed -e 's/_[[:digit:]]//g' output/simdata/mb_denoised_circ_lowess_$i.txt > temp_circ_lowess_$i.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do mv temp_block_$i.txt output/simdata/mb_denoised_block_$i.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do mv temp_circ_lowess_sse_$i.txt output/simdata/mb_denoised_circ_lowess_sse_$i.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do mv temp_circ_lowess_$i.txt output/simdata/mb_denoised_circ_lowess_$i.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do python src/eJTK-CalcP.py -f output/simdata/mb_denoised_block_$i.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do python src/eJTK-CalcP.py -f output/simdata/mb_denoised_circ_lowess_sse_$i.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do python src/eJTK-CalcP.py -f output/simdata/mb_denoised_circ_lowess_$i.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done

From base AS class
COPY --from=ejtk /LIMBR/output /LIMBR/output
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do tail -q -n +2 output/simdata/mb_simulated_data_with_noise_${i}__jtkout_GammaP.txt |cut -f 20 > output/simdata/mb_simulated_data_with_noise_${i}_pvals.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do tail -q -n +2 output/simdata/mb_simulated_data_baseline_${i}__jtkout_GammaP.txt |cut -f 20 > output/simdata/mb_simulated_data_baseline_${i}_pvals.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do tail -q -n +2 output/simdata/mb_denoised_block_${i}__jtkout_GammaP.txt |cut -f 20 > output/simdata/mb_denoised_block_${i}_pvals.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do tail -q -n +2 output/simdata/mb_denoised_circ_lowess_sse_${i}__jtkout_GammaP.txt |cut -f 20 > output/simdata/mb_denoised_circ_lowess_sse_${i}_pvals.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do tail -q -n +2 output/simdata/mb_denoised_circ_lowess_${i}__jtkout_GammaP.txt |cut -f 20 > output/simdata/mb_denoised_circ_lowess_${i}_pvals.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do python3 src/get_classes.py -l output/simdata/mb_simulated_data_key_${i}.txt output/simdata/mb_simulated_data_baseline_${i}.txt output/simdata/mb_simulated_data_baseline_${i}_classes.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do python3 src/get_classes.py -l output/simdata/mb_simulated_data_key_${i}.txt output/simdata/mb_simulated_data_with_noise_${i}.txt output/simdata/mb_simulated_data_with_noise_${i}_classes.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do python3 src/get_classes.py -l output/simdata/mb_simulated_data_key_${i}.txt output/simdata/mb_denoised_block_${i}.txt output/simdata/mb_denoised_block_${i}_classes.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do python3 src/get_classes.py -l output/simdata/mb_simulated_data_key_${i}.txt output/simdata/mb_denoised_circ_lowess_sse_${i}.txt output/simdata/mb_denoised_circ_lowess_sse_${i}_classes.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do python3 src/get_classes.py -l output/simdata/mb_simulated_data_key_${i}.txt output/simdata/mb_denoised_circ_lowess_${i}.txt output/simdata/mb_denoised_circ_lowess_${i}_classes.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do paste output/simdata/mb_simulated_data_baseline_${i}_pvals.txt output/simdata/mb_simulated_data_baseline_${i}_classes.txt > output/simdata/mb_simulated_data_baseline_${i}_classifications.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do paste output/simdata/mb_simulated_data_with_noise_${i}_pvals.txt output/simdata/mb_simulated_data_with_noise_${i}_classes.txt > output/simdata/mb_simulated_data_with_noise_${i}_classifications.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do paste output/simdata/mb_denoised_block_${i}_pvals.txt output/simdata/mb_denoised_block_${i}_classes.txt > output/simdata/mb_denoised_block_${i}_classifications.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do paste output/simdata/mb_denoised_circ_lowess_sse_${i}_pvals.txt output/simdata/mb_denoised_circ_lowess_sse_${i}_classes.txt > output/simdata/mb_denoised_circ_lowess_sse_${i}_classifications.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do paste output/simdata/mb_denoised_circ_lowess_${i}_pvals.txt output/simdata/mb_denoised_circ_lowess_${i}_classes.txt > output/simdata/mb_denoised_circ_lowess_${i}_classifications.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do paste output/simdata/mb_simulated_data_baseline_${i}_pvals.txt output/simdata/mb_simulated_data_baseline_${i}_classes.txt > output/simdata/mb_simulated_data_baseline_${i}_classifications.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do paste output/simdata/mb_simulated_data_with_noise_${i}_pvals.txt output/simdata/mb_simulated_data_with_noise_${i}_classes.txt > output/simdata/mb_simulated_data_with_noise_${i}_classifications.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do python3 src/analysis/ROC_curves_mb.py -l output/simdata/mb_simulated_data_with_noise_${i}_classifications.txt output/simdata/mb_simulated_data_baseline_${i}_classifications.txt output/simdata/mb_denoised_block_${i}_classifications.txt output/simdata/mb_denoised_circ_lowess_${i}_classifications.txt output/simdata/mb_denoised_circ_lowess_sse_${i}_classifications.txt output/simdata/mb_ROC_curves_${i}.pdf; done

From base AS figs
COPY --from=class /LIMBR/output /LIMBR/output
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-get -y build-dep r-base
RUN mkdir /usr/local/src/R
WORKDIR /usr/local/src/R
RUN curl -O http://cran.utstat.utoronto.ca/src/base/R-3/R-3.4.1.tar.gz
RUN tar xzvf R-3.4.1.tar.gz
WORKDIR R-3.4.1
RUN ./configure --enable-R-shlib  --with-blas --with-lapack
RUN make
RUN make install
WORKDIR /usr/local/src/R
RUN rm -f R Rscript
RUN ln -s R-3.4.1/bin/R R
RUN ln -s R-3.4.1/bin/Rscript Rscript
WORKDIR /
RUN su - -c "R -e \"install.packages('circular', repos = 'http://cran.rstudio.com/')\""
WORKDIR /LIMBR/
RUN pip3 install -U rpy2
COPY ./src/Karpeivitch_comparison.r /LIMBR/src/Karpeivitch_comparison.r
COPY ./src/analysis/gen_figs.py /LIMBR/src/analysis/gen_figs.py
RUN Rscript src/Karpeivitch_comparison.r
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do sed -e 's/_[[:digit:]]//g' output/simdata/mb_simdata_eigenMS_$i.txt > temp_eigenMS_$i.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do mv temp_eigenMS_$i.txt output/simdata/mb_denoised_eigenMS_$i.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do python src/eJTK-CalcP.py -f output/simdata/mb_denoised_eigenMS_$i.txt -w src/ref_files/waveform_cosine.txt -a src/ref_files/asymmetries_02-22_by2.txt -s src/ref_files/phases_00-22_by2.txt -p src/ref_files/period24.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do tail -q -n +2 output/simdata/mb_denoised_eigenMS_${i}__jtkout_GammaP.txt |cut -f 20 > output/simdata/mb_denoised_eigenMS_${i}_pvals.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do python3 src/get_classes.py -l output/simdata/mb_simulated_data_key_${i}.txt output/simdata/mb_denoised_eigenMS_${i}.txt output/simdata/mb_denoised_eigenMS_${i}_classes.txt; done
RUN for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do paste output/simdata/mb_denoised_eigenMS_${i}_pvals.txt output/simdata/mb_denoised_eigenMS_${i}_classes.txt > output/simdata/mb_denoised_eigenMS_${i}_classifications.txt; done
COPY ./src/analysis/gen_figs.py /LIMBR/src/analysis/gen_figs.py
RUN python3 src/analysis/gen_figs.py 
