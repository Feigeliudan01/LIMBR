FILES = make_checkpoints/install
FILES += make_checkpoints/simulate1 make_checkpoints/simulate2 make_checkpoints/simulate3 make_checkpoints/simulate4
FILES += make_checkpoints/simulate_mb1 make_checkpoints/simulate_mb2 make_checkpoints/simulate_mb3 make_checkpoints/simulate_mb4
FILES += make_checkpoints/actual make_checkpoints/figures make_checkpoints/figures_mb

all: ${FILES}

#make_checkpoints/install : dockerbuild/Dockerfile_install
make_checkpoints/install :
	@echo installing base image and downloading data
	@docker pull ubuntu
	@docker build --force-rm --no-cache --squash -f dockerbuild/Dockerfile_install -t acrowell/limbr .
	@docker tag acrowell/limbr acrowell/limbr:1
	@touch $@

make_checkpoints/simulate1 : dockerbuild/Dockerfile_sim1 make_checkpoints/install
	@echo running single bias trend simulations
	@docker build --force-rm --no-cache --squash -f dockerbuild/Dockerfile_sim1 -t acrowell/limbr .
	@docker tag acrowell/limbr acrowell/limbr:2
	@touch $@

make_checkpoints/simulate2 : dockerbuild/Dockerfile_sim2 make_checkpoints/simulate1
	@echo running single bias trend simulations 2
	@docker build --force-rm --no-cache --squash -f dockerbuild/Dockerfile_sim2 -t acrowell/limbr .
	@docker tag acrowell/limbr acrowell/limbr:3
	@touch $@

make_checkpoints/simulate3 : dockerbuild/Dockerfile_sim3 make_checkpoints/simulate2
	@echo running single bias trend simulations 3
	@docker build --force-rm --no-cache --squash -f dockerbuild/Dockerfile_sim3 -t acrowell/limbr .
	@docker tag acrowell/limbr acrowell/limbr:4
	@touch $@

make_checkpoints/simulate4 : dockerbuild/Dockerfile_sim4 make_checkpoints/simulate3
	@echo running single bias trend simulations 4
	@docker build --force-rm --no-cache --squash -f dockerbuild/Dockerfile_sim4 -t acrowell/limbr .
	@docker tag acrowell/limbr acrowell/limbr:5
	@touch $@

make_checkpoints/simulate_mb1 : dockerbuild/Dockerfile_sim_mb1 make_checkpoints/simulate4
	@echo running multiple bias trend simulations
	@docker build --force-rm --no-cache --squash -f dockerbuild/Dockerfile_sim_mb1 -t acrowell/limbr .
	@docker tag acrowell/limbr acrowell/limbr:6
	@touch $@

make_checkpoints/simulate_mb2 : dockerbuild/Dockerfile_sim_mb2 make_checkpoints/simulate_mb1
	@echo running multiple bias trend simulations 2
	@docker build --force-rm --no-cache --squash -f dockerbuild/Dockerfile_sim_mb2 -t acrowell/limbr .
	@docker tag acrowell/limbr acrowell/limbr:7
	@touch $@

make_checkpoints/simulate_mb3 : dockerbuild/Dockerfile_sim_mb3 make_checkpoints/simulate_mb2
	@echo running multiple bias trend simulations 3
	@docker build --force-rm --no-cache --squash -f dockerbuild/Dockerfile_sim_mb3 -t acrowell/limbr .
	@docker tag acrowell/limbr acrowell/limbr:8
	@touch $@

make_checkpoints/simulate_mb4 : dockerbuild/Dockerfile_sim_mb4 make_checkpoints/simulate_mb3
	@echo running multiple bias trend simulations 4
	@docker build --force-rm --no-cache --squash -f dockerbuild/Dockerfile_sim_mb4 -t acrowell/limbr .
	@docker tag acrowell/limbr acrowell/limbr:9
	@touch $@

make_checkpoints/actual : dockerbuild/Dockerfile_actual make_checkpoints/simulate_mb4
	@echo running actual data
	@docker build --force-rm --no-cache --squash -f dockerbuild/Dockerfile_actual -t acrowell/limbr .
	@docker tag acrowell/limbr acrowell/limbr:10
	@touch $@

make_checkpoints/figures : dockerbuild/Dockerfile_figs make_checkpoints/actual
	@echo making figures
	@docker build --force-rm --no-cache --squash -f dockerbuild/Dockerfile_figs -t acrowell/limbr .
	@docker tag acrowell/limbr acrowell/limbr:11
	@touch $@

make_checkpoints/figures_mb : dockerbuild/Dockerfile_figs_mb make_checkpoints/figures
	@echo making figures mb
	@docker build --force-rm --no-cache --squash -f dockerbuild/Dockerfile_figs_mb -t acrowell/limbr .
	@docker tag acrowell/limbr acrowell/limbr:12
	@touch $@
