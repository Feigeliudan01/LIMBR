from LIMBR import simulations, imputation, batch_fx

for i in range(1, 21):
    simulation = simulations.simulate(effect_size=0.25, p_miss=.4, lam_miss=15, phase_noise=.5, amp_noise=1.5, rseed=i)
    simulation.generate_pool_map()
    simulation.write_output('standard_'+str(i))

for i in range(1, 21):
    simulation = simulations.simulate(effect_size=0.5, p_miss=.4, lam_miss=15, phase_noise=.5, amp_noise=1.5, rseed=i)
    simulation.generate_pool_map()
    simulation.write_output('double_noise_'+str(i))

for i in range(1, 21):
    simulation = simulations.simulate(effect_size=0.125, p_miss=.4, lam_miss=15, phase_noise=.5, amp_noise=1.5, rseed=i)
    simulation.generate_pool_map()
    simulation.write_output('half_noise_'+str(i))

for i in range(1, 21):
    simulation = simulations.simulate(effect_size=.0625, p_miss=.4, lam_miss=15, phase_noise=.5, amp_noise=1.5, rseed=i)
    simulation.generate_pool_map()
    simulation.write_output('tenth_noise_'+str(i))
