from LIMBR import simulations, imputation, batch_fx

for i in range(1, 21):
    simulation = simulations.simulate()
    simulation.generate_pool_map()
    simulation.write_output('standard_'+str(i))

for i in range(1, 21):
    simulation = simulations.simulate(effect_size=4)
    simulation.generate_pool_map()
    simulation.write_output('double_noise_'+str(i))

for i in range(1, 21):
    simulation = simulations.simulate(effect_size=1)
    simulation.generate_pool_map()
    simulation.write_output('half_noise_'+str(i))

for i in range(1, 21):
    simulation = simulations.simulate(effect_size=.2)
    simulation.generate_pool_map()
    simulation.write_output('tenth_noise_'+str(i))