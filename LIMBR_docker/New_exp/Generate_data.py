from LIMBR import simulations, imputation, batch_fx

for i in range(1, 21):
    simulation = simulations.simulate(p_miss=.4, lam_miss=15, phase_noise=.5, amp_noise=1.5, rseed=i)
    simulation.generate_pool_map()
    simulation.write_output('standard_'+str(i))

for i in range(1, 21):
    simulation = simulations.simulate(effect_size=4,p_miss=.4, lam_miss=15, phase_noise=.5, amp_noise=1.5, rseed=i)
    simulation.generate_pool_map()
    simulation.write_output('double_noise_'+str(i))

for i in range(1, 21):
    simulation = simulations.simulate(effect_size=1,p_miss=.4, lam_miss=15, phase_noise=.5, amp_noise=1.5, rseed=i)
    simulation.generate_pool_map()
    simulation.write_output('half_noise_'+str(i))

for i in range(1, 21):
    simulation = simulations.simulate(effect_size=.2,p_miss=.4, lam_miss=15, phase_noise=.5, amp_noise=1.5, rseed=i)
    simulation.generate_pool_map()
    simulation.write_output('tenth_noise_'+str(i))
