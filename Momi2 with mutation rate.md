```python
import momi		## momi2 analysis
import logging		## create log file

logging.basicConfig(level=logging.INFO,
                    filename="new_muts_momi_log.txt")
```


```sh
%%sh
python -m momi.extract_sfs new_final_sfs.gz 100 /scratch/fm2371/momi_new/counts_new/*.snpAlleleCounts.gz
```


```python
## Load the sfs data
sfs = momi.Sfs.load("new_final_sfs.gz")
```


```bash
%%bash
gunzip -c new_final_sfs.gz | head
```

    {
    	"sampled_pops": ["US"],
    	"folded": true,
    	"length": 2252044789.0,
    	"configs": [
    		[[186, 2]],
    		[[185, 3]],
    		[[184, 4]],
    		[[183, 5]],
    		[[182, 6]],



```python
#load the sfs data,avg number of hets per ind per pop (assuming Hardy-Weinberg)
#print("Avg pairwise heterozygosity", sfs.avg_pairwise_hets)
print("populations", sfs.populations)
print("percent missing data per population", sfs.p_missing)
```

    populations ('US',)
    percent missing data per population [0.02991884]



```python
# Common parameters
Ne = 2.2e7  # Effective population size, calculated based on He
mu = 2.8e-9  # Mutation rate
gen_time =1 #Generation Time

# Time of introduction into Korea and the US, reported by people, the actual invasive time could be ealier
t_korea = 2022 - 2004  # 18 years ago
t_us = 2022 - 2014  # 8 years ago

```


```python
# Model 1: Constant population size

# Initialize a demographic model with a constant population size
model_constant = momi.DemographicModel(N_e=Ne, gen_time=gen_time,muts_per_gen=mu)

# Add a size parameter for the contemporary population
model_constant.add_size_param(name="n_contemp", lower=1e1, upper=1e8)

# Add the SLF population as a leaf in the demographic model
model_constant.add_leaf("US")

model_constant.add_time_param("t_inital", lower=0, upper=100) 

# Set the data for the model (replace 'sfs' with your actual SFS data)
model_constant.set_data(sfs)

# Optimize the model to fit the data
result = model_constant.optimize(method="L-BFGS-B",options={"maxiter": 500}) #Truncated Newton Conjugate-gradient not get converged.This is very accurate but can be slow for large models; for large models, method="L-BFGS-B" is a good choice.

# Display the results
print(result)
```

                fun: 1796.693253099363
           hess_inv: <2x2 LbfgsInvHessProduct with dtype=float64>
                jac: array([0., 0.])
      kl_divergence: 1796.693253099363
     log_likelihood: -538950785.7185186
            message: b'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
               nfev: 1
                nit: 0
               njev: 1
         parameters: ParamsDict({'n_contemp': 7908971.1354903625, 't_inital': 66.41102486116002})
             status: 0
            success: True
                  x: array([15.88350826,  0.68166519])



```python
# Model 2: Expansion population size

model_expansion = momi.DemographicModel(N_e=Ne, gen_time=gen_time,muts_per_gen=mu)#add mutation rate makes very recent time period converged

# Add size parameters
model_expansion.add_size_param("n_init", lower=1e1, upper=1e4)  # Initial population size before expansion, the range has to be lower than current size
model_expansion.add_size_param("n_contemp", lower=1e4, upper=1e8)  # Contemporary population size

# Add time parameters
model_expansion.add_time_param("t_initial", lower=0, upper=100)  # Time of expansion, even could happen before SLF were there

# Add the leaf population (SLF in the US)
model_expansion.add_leaf("US", N="n_init")

# Add the demographic events
model_expansion.set_size("US", t="t_initial", N="n_contemp")

# Set the data
model_expansion.set_data(sfs)

# Print parameters
#print(model_expansion.get_params())
# Run the optimization
result_model_expansion=model_expansion.optimize(method="TNC")
print(result_model_expansion)

```

                fun: 0.7765045295758033
                jac: array([3.94560153e-05, 6.20471412e-01, 1.37118862e-08])
      kl_divergence: 0.7765045295758033
     log_likelihood: -1317964.1996729425
            message: 'Local minimum reached (|pg| ~= 0)'
               nfev: 32
                nit: 13
         parameters: ParamsDict({'n_init': 10.000000000000002, 'n_contemp': 10000.00000000001, 't_initial': 0.1774388543541894})
             status: 0
            success: True
                  x: array([ 2.30258509,  9.21034037, -6.33252343])



```python
# Model 3: Decline population size

model_decline = momi.DemographicModel(N_e=Ne, gen_time=gen_time, muts_per_gen=mu)

# Add size parameters
model_decline.add_size_param("n_init", lower=1e4, upper=1e8)  # Initial population size before decline
model_decline.add_size_param("n_contemp", lower=1e1, upper=1e4)  # Contemporary population size after decline, the same with inital size for Expansion model.

# Add time parameters
model_decline.add_time_param("t_initial", lower=0, upper=100)  # keep time of decline the same as in expansion model

# Add the leaf population (SLF in the US)
model_decline.add_leaf("US", N="n_init")

# Add the demographic events
#model_decline.set_size("US", t="t_init", N="n_init")
model_decline.set_size("US", t="t_initial", N="n_contemp")

# Set the data
model_decline.set_data(sfs)

# Print parameters (Optional)
print(model_decline.get_params())

# Run the optimization
result_model_decline = model_decline.optimize(method="TNC", options={"maxiter": 500})

# Print the result
print(result_model_decline)

```

    ParamsDict({'n_init': 15873754.571860114, 'n_contemp': 5668.009283561372, 't_initial': 58.02568233928568})
                fun: 0.4679806211029051
                jac: array([2.08397511e-10, 6.13176710e-09, 1.72392567e-06])
      kl_divergence: 0.4679806211029051
     log_likelihood: -1225603.2483368618
            message: 'Converged (|f_n-f_(n-1)| ~= 0)'
               nfev: 28
                nit: 11
         parameters: ParamsDict({'n_init': 19194206.536149774, 'n_contemp': 2319.9630149611294, 't_initial': 0.0004899486454302199})
             status: 1
            success: True
                  x: array([ 16.77011905,   7.74930652, -12.22637526])



```python
# Model 4: Decline then Expansion Model(arrived the US 10 years ago when size decline and recently expand)

model_decline_expand = momi.DemographicModel(N_e=Ne, gen_time=gen_time, muts_per_gen=mu)

# Add size parameters
model_decline_expand.add_size_param("n_init", lower=1e4, upper=1e8)
model_decline_expand.add_size_param("n_us_decline", lower=1e1, upper=1e4)
model_decline_expand.add_size_param("n_us_expand", lower=1e4, upper=1e8)

# Add time parameters
model_decline_expand.add_time_param("t_initial", lower=12, upper=100)
# Time of decline: assuming the decline might have started up to 5 years ago and lasted until the expansion phase
model_decline_expand.add_time_param("t_us_decline", lower=6, upper=11)
# Time of expansion: assuming the expansion started sometime in the last 10 years
model_decline_expand.add_time_param("t_us_expand", lower=0, upper=5)

# Add the leaf population (SLF in the US)
model_decline_expand.add_leaf("US", N="n_init")

# Add the demographic events
# First decline
#model_decline_expand.set_size("US", t="t_initial", N="n_init")
model_decline_expand.set_size("US", t="t_us_decline", N="n_us_decline")
# Then expansion
model_decline_expand.set_size("US", t="t_us_expand", N="n_us_expand")

# Set the data
model_decline_expand.set_data(sfs)

# Print parameters
print(model_decline_expand.get_params())

# Run the optimization
result_model_decline_expand = model_decline_expand.optimize(method="TNC", options={"maxiter": 500})

# Print the result
print(result_model_decline_expand)

```

    ParamsDict({'n_init': 98262553.56485932, 'n_us_decline': 2555.0362282148985, 'n_us_expand': 9540334.050488358, 't_initial': 75.56962793576072, 't_us_decline': 7.350189457985172, 't_us_expand': 0.9192705923136802})
                fun: 0.48219874903490745
                jac: array([4.35826433e-12, 7.30480477e-06, 5.37817353e-03, 0.00000000e+00,
           4.14090617e-09, 2.59520357e-10])
      kl_divergence: 0.48219874903490745
     log_likelihood: -1229859.6439870978
            message: 'Converged (|f_n-f_(n-1)| ~= 0)'
               nfev: 34
                nit: 8
         parameters: ParamsDict({'n_init': 605470.462180968, 'n_us_decline': 2536.528197155776, 'n_us_expand': 10000.00000000001, 't_initial': 75.56962793576072, 't_us_decline': 6.000001987639386, 't_us_expand': 2.412783957699757e-07})
             status: 1
            success: True
                  x: array([ 13.31376106,   7.83855157,   9.21034037,   0.95630869,
           -14.73800038, -16.84675227])



```python
# Model 5: Complex Migration and Size Change Model(Assume SLF expericed size decline when they reach to Korea, which is around 20 years ago, and then expand in Korea, and then they were introduced to the US around 10 years ago, which leads size decline in their size and recently expand again)
# As we don’t have data from Korea and China, basically we are testing a model that has size changed four times, which is decline and then expasion and then decline again and fianlly expansion.

model_complex_us = momi.DemographicModel(N_e=Ne, gen_time=gen_time, muts_per_gen=mu)

# Add size parameters
# Add size parameters
model_complex_us.add_size_param("n_initial", lower=1e4, upper=1e8)   # Initial population size
model_complex_us.add_size_param("n_korea_decline1", lower=1e1, upper=1e4)  # Size after decline in Korea
model_complex_us.add_size_param("n_korea_expand1", lower=1e4, upper=1e8)   # Size after expansion in Korea
model_complex_us.add_size_param("n_us_decline2", lower=1e1, upper=1e4)  # Size after decline in the US
model_complex_us.add_size_param("n_us_expand2", lower=1e4, upper=1e8)   # Size after expansion in the US

# Add time parameters for each size change
# First decline in Korea, they were reported in Korean since 2004, but it could be ealier
model_complex_us.add_time_param("t_initial", lower=23, upper=100) #  before in Korea
model_complex_us.add_time_param("t_korea_decline1", lower=17, upper=22) #  2000 to 2006
# First expansion in Korea, must be after the first decline
model_complex_us.add_time_param("t_korea_expand1", lower=12, upper=16) #  2006 to 2010
# Second decline in the US
model_complex_us.add_time_param("t_us_decline2", lower=6, upper=11) # upper limit in the US: could be ealier than 2014, (2011 to 2017)
# Recent expansion in the US
model_complex_us.add_time_param("t_us_expand2", lower=0, upper=5) #0 means the year of 2022, (2017 to 2022)

# Add the leaf population (SLF in the US)
model_complex_us.add_leaf("US", N="n_initial")

# Add the demographic events
#model_complex_us.set_size("US", t="t_initial", N="n_initial")
model_complex_us.set_size("US", t="t_korea_decline1", N="n_korea_decline1")
model_complex_us.set_size("US", t="t_korea_expand1", N="n_korea_expand1")
model_complex_us.set_size("US", t="t_us_decline2", N="n_us_decline2")
model_complex_us.set_size("US", t="t_us_expand2", N="n_us_expand2")

# Set the data
model_complex_us.set_data(sfs)

# Run the optimization
result_model_complex_us = model_complex_us.optimize(method="TNC", options={"maxiter": 500})

# Print the result
print(result_model_complex_us)

```

                fun: 0.45449585092505757
                jac: array([ 1.58650882e-12,  2.62094940e-11, -2.39564311e-08, -5.53164596e-10,
            1.72821720e-03,  0.00000000e+00, -2.11532155e-08,  1.65867718e-08,
            2.18113930e-11,  7.29355001e-10])
      kl_divergence: 0.45449585092505757
     log_likelihood: -1221566.3935973407
            message: 'Converged (|f_n-f_(n-1)| ~= 0)'
               nfev: 81
                nit: 18
         parameters: ParamsDict({'n_initial': 4609542.905089804, 'n_korea_decline1': 2594.983302081868, 'n_korea_expand1': 100000000.00000018, 'n_us_decline2': 307.75347207479257, 'n_us_expand2': 10000.00000000001, 't_initial': 34.029693804484324, 't_korea_decline1': 21.999798828336836, 't_korea_expand1': 12.000024748703398, 't_us_decline2': 6.000000009485941, 't_us_expand2': 1.4234980233590303e-06})
             status: 1
            success: True
                  x: array([ 15.34363926,   7.86133536,  18.42068074,   5.72929905,
             9.21034037,  -1.78861366,  10.12074965, -11.99302563,
           -20.08289285, -15.07183095])



```python
# Model 6: Complex Migration and Size Change Model through Shanghai
#In this model, we assume the SLF in shanghai gone through population size change beacuse urbanization or small population went to shanghai urban areas, which has population size decline. 
# Then some of them were introduced to Korea and had population size change again.
#Finally, some of them were introduced to the US and had population size change again.

model_complex_global = momi.DemographicModel(N_e=Ne, gen_time=gen_time)
                                             
                                             #, muts_per_gen=mu)

# Add size parameters
model_complex_global.add_size_param("n_initial", lower=1e4, upper=1e8)   # Initial population size
model_complex_global.add_size_param("n_shanghai_decline1", lower=1e1, upper=1e4)   # Size after decline in Shanghai
model_complex_global.add_size_param("n_shanghai_expand1", lower=1e4, upper=1e8)   # Size after expansion in Shanghai
model_complex_global.add_size_param("n_korea_decline2", lower=1e1, upper=1e4)  # Size after decline in Korea
model_complex_global.add_size_param("n_korea_expand2", lower=1e4, upper=1e8)   # Size after expansion in Korea
model_complex_global.add_size_param("n_us_decline3", lower=1e1, upper=1e4)  # Size after decline in the US
model_complex_global.add_size_param("n_us_expand3", lower=1e4, upper=1e8)   # Size after expansion in the US

# Add time parameters for each event
model_complex_global.add_time_param("t_initial", lower=61, upper=100)  # Time of initial population
model_complex_global.add_time_param("t_shanghai_decline1", lower=36, upper=60)  # Time of first decline in Shanghai
model_complex_global.add_time_param("t_shanghai_expand1", lower=23, upper=35)   # Time of first expansion in Shanghai
model_complex_global.add_time_param("t_korea_decline2", lower=17, upper=22)  # Time of second decline in Korea
model_complex_global.add_time_param("t_korea_expand2", lower=12, upper=16)   # Time of second expansion in Korea
model_complex_global.add_time_param("t_us_decline3", lower=6, upper=11)  # Time of third decline in the US
model_complex_global.add_time_param("t_us_expand3", lower=0, upper=5)   # Time of third expansion in the US

# Add the leaf population (SLF in the US)
model_complex_global.add_leaf("US", N="n_initial")

# Add the demographic events
#model_complex_global.set_size("US", t="t_initial", N="n_initial")
model_complex_global.set_size("US", t="t_us_expand3", N="n_us_expand3")
model_complex_global.set_size("US", t="t_us_decline3", N="n_us_decline3")
model_complex_global.set_size("US", t="t_korea_expand2", N="n_korea_expand2")
model_complex_global.set_size("US", t="t_korea_decline2", N="n_korea_decline2")
model_complex_global.set_size("US", t="t_shanghai_expand1", N="n_shanghai_expand1")
model_complex_global.set_size("US", t="t_shanghai_decline1", N="n_shanghai_decline1")

# Set the data
model_complex_global.set_data(sfs)

# Run the optimization
result_model_complex_global = model_complex_global.optimize(method="TNC", options={"maxiter": 500})

# Print the result
print(result_model_complex_global)

```

                fun: 0.43688243168785423
                jac: array([ 1.64194377e-08, -5.55022829e-03, -5.22867795e-08, -9.26232007e-06,
           -1.98238367e-08,  1.36774601e-08,  7.02238941e-08,  0.00000000e+00,
           -4.85439705e-12,  1.12804294e-08, -2.59826203e-08,  5.51787568e-15,
            1.52780374e-18,  8.67610316e-09])
      kl_divergence: 0.43688243168785423
     log_likelihood: -1216287.1776836417
            message: 'Converged (|f_n-f_(n-1)| ~= 0)'
               nfev: 53
                nit: 12
         parameters: ParamsDict({'n_initial': 88204171.97201204, 'n_shanghai_decline1': 10000.00000000001, 'n_shanghai_expand1': 100000000.00000018, 'n_korea_decline2': 10000.00000000001, 'n_korea_expand2': 39133760.33574585, 'n_us_decline3': 314.18069667530483, 'n_us_expand3': 50848723.67364753, 't_initial': 82.24808044796693, 't_shanghai_decline1': 59.9999997432368, 't_shanghai_expand1': 23.00120595379491, 't_korea_decline2': 21.997140015179433, 't_korea_expand2': 12.000000000026972, 't_us_decline3': 6.000000000000002, 't_us_expand3': 1.079899133008733})
             status: 1
            success: True
                  x: array([ 18.29516482,   9.21034037,  18.42068074,   9.21034037,
            17.48249609,   5.74996829,  17.74436558,   0.17977291,
            18.3531654 ,  -9.20529064,   7.46580471, -25.72253649,
           -35.37331523,  -1.28924974])






```python
#Try to plot for the best model, 

import matplotlib.pyplot as plt
import numpy as np

# Parameters for the population sizes and time points
n_initial = 3689520
n_shanghai_decline1 = 2587
n_shanghai_expand1 = 13799597
n_korea_decline2 = 10000.00
n_korea_expand2 = 21436465
n_us_decline3 = 231
n_us_expand3 = 6748614
t_initial = 81
t_shanghai_decline1 = 60
t_shanghai_expand1 = 25
t_korea_decline2 = 18
t_korea_expand2 = 12
t_us_decline3 = 6
t_us_expand3 = 1.7

# Time points for the plot
time_points = np.linspace(-100, 0, 1000)  # Covering 100 years to present

# Population size at each time point
population_size = np.zeros_like(time_points)
for i, t in enumerate(time_points):
    if t <= -t_shanghai_decline1:
        population_size[i] = n_initial
    elif -t_shanghai_decline1 < t <= -t_shanghai_expand1:
        population_size[i] = n_shanghai_decline1
    elif -t_shanghai_expand1 < t <= -t_korea_decline2:
        population_size[i] = n_shanghai_expand1
    elif -t_korea_decline2 < t <= -t_korea_expand2:
        population_size[i] = n_korea_decline2
    elif -t_korea_expand2 < t <= -t_us_decline3:
        population_size[i] = n_korea_expand2
    elif -t_us_decline3 < t <= -t_us_expand3:
        population_size[i] = n_us_decline3
    else:
        population_size[i] = n_us_expand3

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(time_points, population_size, color='Orange', linewidth=3)
plt.xlabel('Years before present in USA', labelpad=10)
plt.ylabel('Population Size')
plt.title('SLF Population Dynamics Over Time (model_complex_global)')
plt.grid(False)
plt.yscale('log')  # Log scale for better visualization

# Adding labels for specific events
plt.axvline(x=-t_us_decline3, color='blue', linestyle='--', linewidth=1)
plt.text(-t_us_decline3 -1, n_us_decline3, 'Arrived USA', horizontalalignment='right', verticalalignment='bottom')

plt.axvline(x=-t_korea_decline2, color='blue', linestyle='--', linewidth=1)
plt.text(-t_korea_decline2 - 2 , n_korea_decline2, 'Arrived South Korea', horizontalalignment='right', verticalalignment='bottom')

plt.axvline(x=-t_shanghai_decline1, color='blue', linestyle='--', linewidth=1)
plt.text(-t_shanghai_decline1 - 2, n_shanghai_decline1, 'Urbanization of Shanghai', horizontalalignment='right', verticalalignment='bottom')

# Saving the custom plot as a PDF
plt.savefig("SLF_Complex_Global_Model_Custom_Plot.pdf", format='pdf')

plt.show()


```


    
![png](output_13_0.png)
    



```python
### Model selection with AIC

import numpy as np
AICs = []
model_list = [model_constant,model_expansion,model_decline,model_decline_expand,model_complex_us,model_complex_global]
model_names = ["Constant", "Expansion", "Decline", "Decline and Expand", "Complex US", "Complex Global"]

for model in model_list:
    lik = model.log_likelihood()
    nparams = len(model.get_params())
    aic = 2 * nparams - 2 * lik
    print(f"AIC for {model.__class__.__name__}: {aic}")
    AICs.append(aic)

minv = np.min(AICs)
delta_aic = np.array(AICs) - minv
print("Delta AIC per model: ", delta_aic)
print("AIC weight per model: ", np.exp(-0.5 * delta_aic))
```

    AIC for DemographicModel: 1077901575.4370372
    AIC for DemographicModel: 2635934.399345885
    AIC for DemographicModel: 2451212.4966737237
    AIC for DemographicModel: 2459731.2879741956
    AIC for DemographicModel: 2443152.7871946814
    AIC for DemographicModel: 2441627.182442166
    Delta AIC per model:  [1.07545995e+09 1.94307217e+05 9.58531423e+03 1.81041055e+04
     1.52560475e+03 0.00000000e+00]
    AIC weight per model:  [0. 0. 0. 0. 0. 1.]



```python
### Based on the AIC result, Delta AIC, AIC weight, the the "model_complex_migration_revised_time" has the best fit to the data among the models.
### Its AIC value is the lowest, and its AIC weight is 1, which suggests a very strong likelihood that this model is the best.
```


```python
# run 1000 replicates for each model and choose the best model

import numpy as np

np.random.seed(318)

model_list = [model_constant, model_expansion, model_decline, model_decline_expand, model_complex_us, model_complex_global]
model_names = ["Constant", "Expansion", "Decline", "Decline and Expand", "Complex US", "Complex Global"]

num_replicates = 1100

with open('318aic_and_likelihoods.txt', 'w') as aic_file, open('318best_models_based_on_weights.txt', 'w') as best_models_file, open('318best_model_params.txt', 'w') as params_file:

    for i in range(num_replicates):
        print(f"Running replicate {i+1} out of {num_replicates}")

        try:
            resampled_sfs = sfs.resample()

            aics = []
            log_likelihoods = []

            for model in model_list:
                model.set_data(resampled_sfs)
                model.set_params(randomize=True)
                model.optimize(method="TNC")
                log_likelihood = model.log_likelihood()
                nparams = len(model.get_params())
                aic = 2 * nparams - 2 * log_likelihood
                aics.append(aic)
                log_likelihoods.append(log_likelihood)

            # Delta AIC calculation
            min_aic = np.min(aics)
            delta_aics = np.array(aics) - min_aic

            # AIC weights calculation
            aic_weights = np.exp(-0.5 * delta_aics) / np.sum(np.exp(-0.5 * delta_aics))

            # Choose the best model based on the highest AIC weight
            best_model_index = np.argmax(aic_weights)
            best_model = model_list[best_model_index]
            params = best_model.get_params()

            # Write detailed info for each replicate to '555aic_and_likelihoods.txt'
            aic_file.write(f'Replicate {i+1}:\n')
            for model_name, aic, ll, delta_aic, aic_weight in zip(model_names, aics, log_likelihoods, delta_aics, aic_weights):
                aic_file.write(f'{model_name} - AIC: {aic}, Log Likelihood: {ll}, Delta AIC: {delta_aic}, AIC Weight: {aic_weight}\n')
            aic_file.write('\n')
            aic_file.flush()

            # Record the index of the best model based on AIC weight for this replicate
            best_models_file.write(f'{best_model_index}\n')
            best_models_file.flush()

            # Save the parameters of the best model for this replicate
            params_file.write(f'Replicate {i+1} Best Model: {model_names[best_model_index]}\n')
            for key, value in params.items():
                params_file.write(f'{key}: {value}\n')
            params_file.write('\n')
            params_file.flush()

        except Exception as e:
            print(f"An error occurred in replicate {i+1}: {str(e)}")
            print("Skipping this replicate and moving on to the next one.")

```

    Running replicate 1 out of 1100


    /ext3/miniconda3/envs/momi-py36/lib/python3.6/site-packages/autograd/tracer.py:14: UserWarning: Output seems independent of input.
      warnings.warn("Output seems independent of input.")


    Running replicate 2 out of 1100
    Running replicate 3 out of 1100
    Running replicate 4 out of 1100
    Running replicate 5 out of 1100
    Running replicate 6 out of 1100
    Running replicate 7 out of 1100
    Running replicate 8 out of 1100
    Running replicate 9 out of 1100
    Running replicate 10 out of 1100
    Running replicate 11 out of 1100
    Running replicate 12 out of 1100
    Running replicate 13 out of 1100
    Running replicate 14 out of 1100
    Running replicate 15 out of 1100
    Running replicate 16 out of 1100
    Running replicate 17 out of 1100
    Running replicate 18 out of 1100
    Running replicate 19 out of 1100
    Running replicate 20 out of 1100
    Running replicate 21 out of 1100
    Running replicate 22 out of 1100
    Running replicate 23 out of 1100
    Running replicate 24 out of 1100
    Running replicate 25 out of 1100
    Running replicate 26 out of 1100
    Running replicate 27 out of 1100
    Running replicate 28 out of 1100
    Running replicate 29 out of 1100
    Running replicate 30 out of 1100
    Running replicate 31 out of 1100
    Running replicate 32 out of 1100
    Running replicate 33 out of 1100


    /ext3/miniconda3/envs/momi-py36/lib/python3.6/site-packages/autograd/numpy/numpy_vjps.py:53: RuntimeWarning: overflow encountered in double_scalars
      lambda ans, x, y : unbroadcast_f(y, lambda g: - g * x / y**2))


    Running replicate 34 out of 1100
    Running replicate 35 out of 1100
    Running replicate 36 out of 1100
    Running replicate 37 out of 1100
    Running replicate 38 out of 1100
    Running replicate 39 out of 1100
    Running replicate 40 out of 1100
    Running replicate 41 out of 1100
    Running replicate 42 out of 1100
    Running replicate 43 out of 1100
    Running replicate 44 out of 1100
    Running replicate 45 out of 1100
    Running replicate 46 out of 1100
    Running replicate 47 out of 1100
    Running replicate 48 out of 1100
    Running replicate 49 out of 1100
    Running replicate 50 out of 1100
    Running replicate 51 out of 1100
    Running replicate 52 out of 1100
    Running replicate 53 out of 1100
    Running replicate 54 out of 1100
    Running replicate 55 out of 1100
    Running replicate 56 out of 1100
    Running replicate 57 out of 1100
    Running replicate 58 out of 1100
    Running replicate 59 out of 1100
    Running replicate 60 out of 1100
    Running replicate 61 out of 1100
    Running replicate 62 out of 1100
    Running replicate 63 out of 1100
    Running replicate 64 out of 1100
    Running replicate 65 out of 1100
    Running replicate 66 out of 1100
    Running replicate 67 out of 1100
    Running replicate 68 out of 1100
    Running replicate 69 out of 1100
    Running replicate 70 out of 1100
    Running replicate 71 out of 1100
    Running replicate 72 out of 1100
    Running replicate 73 out of 1100
    Running replicate 74 out of 1100
    Running replicate 75 out of 1100
    Running replicate 76 out of 1100
    Running replicate 77 out of 1100
    Running replicate 78 out of 1100
    Running replicate 79 out of 1100
    Running replicate 80 out of 1100
    Running replicate 81 out of 1100
    Running replicate 82 out of 1100
    Running replicate 83 out of 1100
    Running replicate 84 out of 1100
    Running replicate 85 out of 1100
    Running replicate 86 out of 1100
    Running replicate 87 out of 1100
    Running replicate 88 out of 1100
    Running replicate 89 out of 1100
    Running replicate 90 out of 1100
    Running replicate 91 out of 1100
    Running replicate 92 out of 1100
    Running replicate 93 out of 1100
    Running replicate 94 out of 1100
    Running replicate 95 out of 1100
    Running replicate 96 out of 1100
    Running replicate 97 out of 1100
    Running replicate 98 out of 1100
    Running replicate 99 out of 1100
    Running replicate 100 out of 1100
    Running replicate 101 out of 1100
    Running replicate 102 out of 1100
    Running replicate 103 out of 1100
    Running replicate 104 out of 1100
    Running replicate 105 out of 1100
    Running replicate 106 out of 1100
    Running replicate 107 out of 1100
    Running replicate 108 out of 1100
    Running replicate 109 out of 1100
    Running replicate 110 out of 1100
    Running replicate 111 out of 1100
    Running replicate 112 out of 1100
    Running replicate 113 out of 1100
    Running replicate 114 out of 1100
    Running replicate 115 out of 1100
    Running replicate 116 out of 1100
    Running replicate 117 out of 1100
    Running replicate 118 out of 1100
    Running replicate 119 out of 1100
    Running replicate 120 out of 1100
    Running replicate 121 out of 1100
    Running replicate 122 out of 1100


    /ext3/miniconda3/envs/momi-py36/lib/python3.6/site-packages/autograd/tracer.py:48: RuntimeWarning: overflow encountered in exp
      return f_raw(*args, **kwargs)
    /ext3/miniconda3/envs/momi-py36/lib/python3.6/site-packages/autograd/numpy/numpy_vjps.py:75: RuntimeWarning: invalid value encountered in double_scalars
      defvjp(anp.exp,    lambda ans, x : lambda g: ans * g)


    An error occurred in replicate 122: Exception encountered at parameter values {'n_initial': 'Autograd ArrayBox with value nan', 'n_shanghai_decline1': 'Autograd ArrayBox with value nan', 'n_shanghai_expand1': 'Autograd ArrayBox with value nan', 'n_korea_decline2': 'Autograd ArrayBox with value nan', 'n_korea_expand2': 'Autograd ArrayBox with value nan', 'n_us_decline3': 'Autograd ArrayBox with value nan', 'n_us_expand3': 'Autograd ArrayBox with value nan', 't_initial': 'Autograd ArrayBox with value nan', 't_shanghai_decline1': 'Autograd ArrayBox with value nan', 't_shanghai_expand1': 'Autograd ArrayBox with value nan', 't_korea_decline2': 'Autograd ArrayBox with value nan', 't_korea_expand2': 'Autograd ArrayBox with value nan', 't_us_decline3': 'Autograd ArrayBox with value nan', 't_us_expand3': 'Autograd ArrayBox with value nan'} (internal scaling x = {'n_initial': 'Autograd ArrayBox with value nan', 'n_shanghai_decline1': 'Autograd ArrayBox with value nan', 'n_shanghai_expand1': 'Autograd ArrayBox with value nan', 'n_korea_decline2': 'Autograd ArrayBox with value nan', 'n_korea_expand2': 'Autograd ArrayBox with value nan', 'n_us_decline3': 'Autograd ArrayBox with value nan', 'n_us_expand3': 'Autograd ArrayBox with value nan', 't_initial': 'Autograd ArrayBox with value nan', 't_shanghai_decline1': 'Autograd ArrayBox with value nan', 't_shanghai_expand1': 'Autograd ArrayBox with value nan', 't_korea_decline2': 'Autograd ArrayBox with value nan', 't_korea_expand2': 'Autograd ArrayBox with value nan', 't_us_decline3': 'Autograd ArrayBox with value nan', 't_us_expand3': 'Autograd ArrayBox with value nan'})
    Skipping this replicate and moving on to the next one.
    Running replicate 123 out of 1100
    Running replicate 124 out of 1100
    Running replicate 125 out of 1100
    Running replicate 126 out of 1100
    Running replicate 127 out of 1100
    Running replicate 128 out of 1100
    Running replicate 129 out of 1100
    Running replicate 130 out of 1100
    Running replicate 131 out of 1100
    Running replicate 132 out of 1100
    Running replicate 133 out of 1100
    Running replicate 134 out of 1100
    Running replicate 135 out of 1100
    Running replicate 136 out of 1100
    Running replicate 137 out of 1100
    Running replicate 138 out of 1100
    Running replicate 139 out of 1100
    Running replicate 140 out of 1100
    Running replicate 141 out of 1100
    Running replicate 142 out of 1100
    Running replicate 143 out of 1100
    Running replicate 144 out of 1100
    Running replicate 145 out of 1100
    Running replicate 146 out of 1100
    Running replicate 147 out of 1100
    Running replicate 148 out of 1100
    Running replicate 149 out of 1100
    Running replicate 150 out of 1100
    Running replicate 151 out of 1100
    Running replicate 152 out of 1100
    Running replicate 153 out of 1100
    Running replicate 154 out of 1100
    Running replicate 155 out of 1100
    Running replicate 156 out of 1100
    Running replicate 157 out of 1100
    Running replicate 158 out of 1100
    Running replicate 159 out of 1100
    Running replicate 160 out of 1100
    Running replicate 161 out of 1100
    Running replicate 162 out of 1100
    Running replicate 163 out of 1100
    Running replicate 164 out of 1100
    Running replicate 165 out of 1100
    Running replicate 166 out of 1100
    Running replicate 167 out of 1100
    Running replicate 168 out of 1100
    Running replicate 169 out of 1100
    Running replicate 170 out of 1100
    Running replicate 171 out of 1100
    Running replicate 172 out of 1100
    Running replicate 173 out of 1100
    Running replicate 174 out of 1100
    Running replicate 175 out of 1100
    Running replicate 176 out of 1100
    Running replicate 177 out of 1100
    Running replicate 178 out of 1100
    Running replicate 179 out of 1100
    Running replicate 180 out of 1100
    Running replicate 181 out of 1100
    Running replicate 182 out of 1100
    Running replicate 183 out of 1100
    Running replicate 184 out of 1100
    Running replicate 185 out of 1100
    Running replicate 186 out of 1100
    Running replicate 187 out of 1100
    Running replicate 188 out of 1100
    Running replicate 189 out of 1100
    Running replicate 190 out of 1100
    Running replicate 191 out of 1100
    Running replicate 192 out of 1100
    Running replicate 193 out of 1100
    Running replicate 194 out of 1100
    Running replicate 195 out of 1100
    Running replicate 196 out of 1100
    Running replicate 197 out of 1100
    Running replicate 198 out of 1100
    Running replicate 199 out of 1100
    Running replicate 200 out of 1100
    Running replicate 201 out of 1100
    Running replicate 202 out of 1100
    Running replicate 203 out of 1100
    Running replicate 204 out of 1100
    Running replicate 205 out of 1100
    Running replicate 206 out of 1100
    Running replicate 207 out of 1100
    Running replicate 208 out of 1100
    Running replicate 209 out of 1100
    Running replicate 210 out of 1100
    Running replicate 211 out of 1100
    Running replicate 212 out of 1100
    Running replicate 213 out of 1100
    Running replicate 214 out of 1100
    Running replicate 215 out of 1100
    Running replicate 216 out of 1100
    Running replicate 217 out of 1100
    Running replicate 218 out of 1100
    Running replicate 219 out of 1100
    Running replicate 220 out of 1100
    Running replicate 221 out of 1100
    Running replicate 222 out of 1100
    Running replicate 223 out of 1100
    Running replicate 224 out of 1100
    Running replicate 225 out of 1100
    Running replicate 226 out of 1100
    Running replicate 227 out of 1100
    Running replicate 228 out of 1100
    Running replicate 229 out of 1100
    Running replicate 230 out of 1100
    Running replicate 231 out of 1100
    Running replicate 232 out of 1100
    Running replicate 233 out of 1100
    Running replicate 234 out of 1100
    Running replicate 235 out of 1100
    Running replicate 236 out of 1100
    Running replicate 237 out of 1100
    Running replicate 238 out of 1100
    Running replicate 239 out of 1100
    Running replicate 240 out of 1100
    Running replicate 241 out of 1100
    Running replicate 242 out of 1100
    Running replicate 243 out of 1100
    Running replicate 244 out of 1100
    Running replicate 245 out of 1100
    Running replicate 246 out of 1100
    Running replicate 247 out of 1100
    Running replicate 248 out of 1100
    Running replicate 249 out of 1100
    Running replicate 250 out of 1100
    Running replicate 251 out of 1100
    Running replicate 252 out of 1100
    Running replicate 253 out of 1100
    Running replicate 254 out of 1100
    Running replicate 255 out of 1100
    Running replicate 256 out of 1100
    Running replicate 257 out of 1100
    Running replicate 258 out of 1100
    Running replicate 259 out of 1100
    Running replicate 260 out of 1100
    Running replicate 261 out of 1100
    Running replicate 262 out of 1100
    Running replicate 263 out of 1100
    Running replicate 264 out of 1100
    Running replicate 265 out of 1100
    Running replicate 266 out of 1100
    Running replicate 267 out of 1100
    Running replicate 268 out of 1100
    Running replicate 269 out of 1100
    Running replicate 270 out of 1100
    Running replicate 271 out of 1100
    Running replicate 272 out of 1100
    Running replicate 273 out of 1100
    Running replicate 274 out of 1100
    Running replicate 275 out of 1100
    Running replicate 276 out of 1100
    Running replicate 277 out of 1100
    Running replicate 278 out of 1100
    Running replicate 279 out of 1100
    Running replicate 280 out of 1100
    Running replicate 281 out of 1100
    Running replicate 282 out of 1100
    Running replicate 283 out of 1100
    Running replicate 284 out of 1100
    Running replicate 285 out of 1100
    Running replicate 286 out of 1100
    Running replicate 287 out of 1100
    Running replicate 288 out of 1100
    Running replicate 289 out of 1100
    Running replicate 290 out of 1100
    Running replicate 291 out of 1100
    Running replicate 292 out of 1100
    Running replicate 293 out of 1100
    Running replicate 294 out of 1100
    Running replicate 295 out of 1100
    Running replicate 296 out of 1100
    Running replicate 297 out of 1100
    Running replicate 298 out of 1100
    Running replicate 299 out of 1100
    Running replicate 300 out of 1100
    Running replicate 301 out of 1100
    Running replicate 302 out of 1100
    Running replicate 303 out of 1100
    Running replicate 304 out of 1100
    Running replicate 305 out of 1100
    Running replicate 306 out of 1100
    Running replicate 307 out of 1100
    Running replicate 308 out of 1100
    Running replicate 309 out of 1100
    Running replicate 310 out of 1100
    Running replicate 311 out of 1100
    Running replicate 312 out of 1100
    Running replicate 313 out of 1100
    Running replicate 314 out of 1100
    Running replicate 315 out of 1100
    Running replicate 316 out of 1100
    Running replicate 317 out of 1100
    Running replicate 318 out of 1100
    Running replicate 319 out of 1100
    Running replicate 320 out of 1100
    Running replicate 321 out of 1100
    Running replicate 322 out of 1100
    Running replicate 323 out of 1100
    Running replicate 324 out of 1100
    Running replicate 325 out of 1100
    Running replicate 326 out of 1100
    Running replicate 327 out of 1100
    Running replicate 328 out of 1100
    Running replicate 329 out of 1100
    Running replicate 330 out of 1100
    Running replicate 331 out of 1100
    Running replicate 332 out of 1100
    Running replicate 333 out of 1100
    Running replicate 334 out of 1100
    Running replicate 335 out of 1100
    Running replicate 336 out of 1100
    Running replicate 337 out of 1100
    Running replicate 338 out of 1100
    Running replicate 339 out of 1100
    Running replicate 340 out of 1100
    Running replicate 341 out of 1100
    Running replicate 342 out of 1100
    Running replicate 343 out of 1100
    Running replicate 344 out of 1100
    Running replicate 345 out of 1100
    Running replicate 346 out of 1100
    Running replicate 347 out of 1100
    Running replicate 348 out of 1100
    Running replicate 349 out of 1100
    Running replicate 350 out of 1100
    Running replicate 351 out of 1100
    Running replicate 352 out of 1100
    Running replicate 353 out of 1100
    Running replicate 354 out of 1100
    Running replicate 355 out of 1100
    Running replicate 356 out of 1100
    Running replicate 357 out of 1100
    Running replicate 358 out of 1100
    Running replicate 359 out of 1100
    Running replicate 360 out of 1100
    Running replicate 361 out of 1100
    Running replicate 362 out of 1100
    Running replicate 363 out of 1100
    Running replicate 364 out of 1100
    Running replicate 365 out of 1100
    Running replicate 366 out of 1100
    Running replicate 367 out of 1100
    Running replicate 368 out of 1100
    Running replicate 369 out of 1100
    Running replicate 370 out of 1100
    Running replicate 371 out of 1100
    Running replicate 372 out of 1100
    Running replicate 373 out of 1100
    Running replicate 374 out of 1100
    Running replicate 375 out of 1100
    Running replicate 376 out of 1100
    Running replicate 377 out of 1100
    Running replicate 378 out of 1100
    Running replicate 379 out of 1100
    An error occurred in replicate 379: Exception encountered at parameter values {'n_initial': 'Autograd ArrayBox with value nan', 'n_shanghai_decline1': 'Autograd ArrayBox with value nan', 'n_shanghai_expand1': 'Autograd ArrayBox with value nan', 'n_korea_decline2': 'Autograd ArrayBox with value nan', 'n_korea_expand2': 'Autograd ArrayBox with value nan', 'n_us_decline3': 'Autograd ArrayBox with value nan', 'n_us_expand3': 'Autograd ArrayBox with value nan', 't_initial': 'Autograd ArrayBox with value nan', 't_shanghai_decline1': 'Autograd ArrayBox with value nan', 't_shanghai_expand1': 'Autograd ArrayBox with value nan', 't_korea_decline2': 'Autograd ArrayBox with value nan', 't_korea_expand2': 'Autograd ArrayBox with value nan', 't_us_decline3': 'Autograd ArrayBox with value nan', 't_us_expand3': 'Autograd ArrayBox with value nan'} (internal scaling x = {'n_initial': 'Autograd ArrayBox with value nan', 'n_shanghai_decline1': 'Autograd ArrayBox with value nan', 'n_shanghai_expand1': 'Autograd ArrayBox with value nan', 'n_korea_decline2': 'Autograd ArrayBox with value nan', 'n_korea_expand2': 'Autograd ArrayBox with value nan', 'n_us_decline3': 'Autograd ArrayBox with value nan', 'n_us_expand3': 'Autograd ArrayBox with value nan', 't_initial': 'Autograd ArrayBox with value nan', 't_shanghai_decline1': 'Autograd ArrayBox with value nan', 't_shanghai_expand1': 'Autograd ArrayBox with value nan', 't_korea_decline2': 'Autograd ArrayBox with value nan', 't_korea_expand2': 'Autograd ArrayBox with value nan', 't_us_decline3': 'Autograd ArrayBox with value nan', 't_us_expand3': 'Autograd ArrayBox with value nan'})
    Skipping this replicate and moving on to the next one.
    Running replicate 380 out of 1100
    Running replicate 381 out of 1100
    Running replicate 382 out of 1100
    Running replicate 383 out of 1100
    Running replicate 384 out of 1100
    Running replicate 385 out of 1100
    Running replicate 386 out of 1100
    Running replicate 387 out of 1100
    Running replicate 388 out of 1100
    Running replicate 389 out of 1100
    Running replicate 390 out of 1100
    Running replicate 391 out of 1100
    Running replicate 392 out of 1100
    Running replicate 393 out of 1100
    Running replicate 394 out of 1100
    Running replicate 395 out of 1100
    Running replicate 396 out of 1100
    Running replicate 397 out of 1100
    Running replicate 398 out of 1100
    Running replicate 399 out of 1100
    Running replicate 400 out of 1100
    Running replicate 401 out of 1100
    Running replicate 402 out of 1100
    Running replicate 403 out of 1100
    Running replicate 404 out of 1100
    Running replicate 405 out of 1100
    Running replicate 406 out of 1100
    Running replicate 407 out of 1100
    Running replicate 408 out of 1100
    Running replicate 409 out of 1100
    Running replicate 410 out of 1100
    Running replicate 411 out of 1100
    Running replicate 412 out of 1100
    Running replicate 413 out of 1100
    Running replicate 414 out of 1100
    Running replicate 415 out of 1100
    Running replicate 416 out of 1100
    Running replicate 417 out of 1100
    Running replicate 418 out of 1100
    Running replicate 419 out of 1100
    Running replicate 420 out of 1100
    Running replicate 421 out of 1100
    Running replicate 422 out of 1100
    Running replicate 423 out of 1100
    Running replicate 424 out of 1100
    Running replicate 425 out of 1100
    Running replicate 426 out of 1100
    Running replicate 427 out of 1100
    Running replicate 428 out of 1100
    Running replicate 429 out of 1100
    Running replicate 430 out of 1100
    Running replicate 431 out of 1100
    Running replicate 432 out of 1100
    Running replicate 433 out of 1100
    Running replicate 434 out of 1100
    Running replicate 435 out of 1100
    Running replicate 436 out of 1100
    Running replicate 437 out of 1100
    Running replicate 438 out of 1100
    Running replicate 439 out of 1100
    Running replicate 440 out of 1100
    Running replicate 441 out of 1100
    Running replicate 442 out of 1100
    Running replicate 443 out of 1100
    Running replicate 444 out of 1100
    Running replicate 445 out of 1100
    Running replicate 446 out of 1100
    Running replicate 447 out of 1100
    Running replicate 448 out of 1100
    Running replicate 449 out of 1100
    Running replicate 450 out of 1100
    Running replicate 451 out of 1100
    Running replicate 452 out of 1100
    Running replicate 453 out of 1100
    Running replicate 454 out of 1100
    Running replicate 455 out of 1100
    Running replicate 456 out of 1100
    Running replicate 457 out of 1100
    Running replicate 458 out of 1100
    Running replicate 459 out of 1100
    Running replicate 460 out of 1100
    Running replicate 461 out of 1100
    Running replicate 462 out of 1100
    Running replicate 463 out of 1100
    Running replicate 464 out of 1100
    Running replicate 465 out of 1100
    Running replicate 466 out of 1100
    Running replicate 467 out of 1100
    Running replicate 468 out of 1100
    Running replicate 469 out of 1100
    Running replicate 470 out of 1100
    Running replicate 471 out of 1100
    Running replicate 472 out of 1100
    Running replicate 473 out of 1100
    Running replicate 474 out of 1100
    Running replicate 475 out of 1100
    Running replicate 476 out of 1100
    Running replicate 477 out of 1100
    Running replicate 478 out of 1100
    Running replicate 479 out of 1100
    Running replicate 480 out of 1100
    Running replicate 481 out of 1100
    Running replicate 482 out of 1100
    Running replicate 483 out of 1100
    Running replicate 484 out of 1100
    Running replicate 485 out of 1100
    Running replicate 486 out of 1100
    Running replicate 487 out of 1100
    Running replicate 488 out of 1100
    Running replicate 489 out of 1100
    Running replicate 490 out of 1100
    Running replicate 491 out of 1100
    Running replicate 492 out of 1100
    Running replicate 493 out of 1100
    Running replicate 494 out of 1100
    Running replicate 495 out of 1100
    Running replicate 496 out of 1100
    Running replicate 497 out of 1100
    Running replicate 498 out of 1100
    Running replicate 499 out of 1100
    Running replicate 500 out of 1100
    Running replicate 501 out of 1100
    Running replicate 502 out of 1100
    Running replicate 503 out of 1100
    Running replicate 504 out of 1100
    Running replicate 505 out of 1100
    Running replicate 506 out of 1100
    Running replicate 507 out of 1100
    Running replicate 508 out of 1100
    Running replicate 509 out of 1100
    Running replicate 510 out of 1100
    Running replicate 511 out of 1100
    Running replicate 512 out of 1100
    Running replicate 513 out of 1100
    Running replicate 514 out of 1100
    Running replicate 515 out of 1100
    Running replicate 516 out of 1100
    Running replicate 517 out of 1100
    Running replicate 518 out of 1100
    Running replicate 519 out of 1100
    Running replicate 520 out of 1100
    Running replicate 521 out of 1100
    Running replicate 522 out of 1100
    Running replicate 523 out of 1100
    Running replicate 524 out of 1100
    Running replicate 525 out of 1100
    Running replicate 526 out of 1100
    Running replicate 527 out of 1100
    Running replicate 528 out of 1100
    Running replicate 529 out of 1100
    Running replicate 530 out of 1100
    Running replicate 531 out of 1100
    Running replicate 532 out of 1100
    Running replicate 533 out of 1100
    Running replicate 534 out of 1100
    Running replicate 535 out of 1100
    Running replicate 536 out of 1100
    Running replicate 537 out of 1100
    Running replicate 538 out of 1100
    Running replicate 539 out of 1100
    Running replicate 540 out of 1100
    Running replicate 541 out of 1100
    Running replicate 542 out of 1100
    Running replicate 543 out of 1100
    Running replicate 544 out of 1100
    Running replicate 545 out of 1100
    Running replicate 546 out of 1100
    Running replicate 547 out of 1100
    Running replicate 548 out of 1100
    Running replicate 549 out of 1100
    Running replicate 550 out of 1100
    Running replicate 551 out of 1100
    Running replicate 552 out of 1100
    Running replicate 553 out of 1100
    Running replicate 554 out of 1100
    Running replicate 555 out of 1100
    Running replicate 556 out of 1100
    Running replicate 557 out of 1100
    Running replicate 558 out of 1100
    Running replicate 559 out of 1100
    Running replicate 560 out of 1100
    Running replicate 561 out of 1100
    Running replicate 562 out of 1100
    Running replicate 563 out of 1100
    Running replicate 564 out of 1100
    Running replicate 565 out of 1100
    Running replicate 566 out of 1100
    Running replicate 567 out of 1100
    Running replicate 568 out of 1100
    Running replicate 569 out of 1100
    Running replicate 570 out of 1100
    Running replicate 571 out of 1100
    Running replicate 572 out of 1100
    Running replicate 573 out of 1100
    Running replicate 574 out of 1100
    Running replicate 575 out of 1100
    Running replicate 576 out of 1100
    Running replicate 577 out of 1100
    Running replicate 578 out of 1100
    Running replicate 579 out of 1100
    Running replicate 580 out of 1100
    Running replicate 581 out of 1100
    Running replicate 582 out of 1100
    Running replicate 583 out of 1100
    Running replicate 584 out of 1100
    Running replicate 585 out of 1100
    Running replicate 586 out of 1100
    Running replicate 587 out of 1100
    Running replicate 588 out of 1100
    Running replicate 589 out of 1100
    Running replicate 590 out of 1100
    Running replicate 591 out of 1100
    Running replicate 592 out of 1100
    Running replicate 593 out of 1100
    Running replicate 594 out of 1100
    Running replicate 595 out of 1100
    Running replicate 596 out of 1100
    Running replicate 597 out of 1100
    Running replicate 598 out of 1100
    Running replicate 599 out of 1100
    Running replicate 600 out of 1100
    Running replicate 601 out of 1100
    Running replicate 602 out of 1100
    Running replicate 603 out of 1100
    Running replicate 604 out of 1100
    Running replicate 605 out of 1100
    Running replicate 606 out of 1100
    Running replicate 607 out of 1100
    Running replicate 608 out of 1100
    Running replicate 609 out of 1100
    Running replicate 610 out of 1100
    Running replicate 611 out of 1100
    Running replicate 612 out of 1100
    Running replicate 613 out of 1100
    Running replicate 614 out of 1100
    Running replicate 615 out of 1100
    Running replicate 616 out of 1100
    Running replicate 617 out of 1100
    Running replicate 618 out of 1100
    Running replicate 619 out of 1100
    Running replicate 620 out of 1100
    Running replicate 621 out of 1100
    Running replicate 622 out of 1100
    Running replicate 623 out of 1100
    Running replicate 624 out of 1100
    Running replicate 625 out of 1100
    Running replicate 626 out of 1100
    Running replicate 627 out of 1100
    Running replicate 628 out of 1100
    Running replicate 629 out of 1100
    Running replicate 630 out of 1100
    Running replicate 631 out of 1100
    Running replicate 632 out of 1100
    Running replicate 633 out of 1100
    Running replicate 634 out of 1100
    Running replicate 635 out of 1100
    Running replicate 636 out of 1100
    Running replicate 637 out of 1100
    Running replicate 638 out of 1100
    Running replicate 639 out of 1100
    Running replicate 640 out of 1100
    Running replicate 641 out of 1100
    Running replicate 642 out of 1100
    Running replicate 643 out of 1100
    Running replicate 644 out of 1100
    Running replicate 645 out of 1100
    Running replicate 646 out of 1100
    Running replicate 647 out of 1100
    Running replicate 648 out of 1100
    Running replicate 649 out of 1100
    Running replicate 650 out of 1100
    Running replicate 651 out of 1100
    Running replicate 652 out of 1100
    Running replicate 653 out of 1100
    Running replicate 654 out of 1100
    Running replicate 655 out of 1100
    Running replicate 656 out of 1100
    Running replicate 657 out of 1100
    Running replicate 658 out of 1100
    Running replicate 659 out of 1100
    Running replicate 660 out of 1100
    Running replicate 661 out of 1100
    Running replicate 662 out of 1100
    Running replicate 663 out of 1100
    Running replicate 664 out of 1100
    Running replicate 665 out of 1100
    Running replicate 666 out of 1100
    Running replicate 667 out of 1100
    Running replicate 668 out of 1100
    Running replicate 669 out of 1100
    Running replicate 670 out of 1100
    Running replicate 671 out of 1100
    Running replicate 672 out of 1100
    Running replicate 673 out of 1100
    Running replicate 674 out of 1100
    Running replicate 675 out of 1100
    Running replicate 676 out of 1100
    Running replicate 677 out of 1100
    Running replicate 678 out of 1100
    Running replicate 679 out of 1100
    Running replicate 680 out of 1100
    Running replicate 681 out of 1100
    Running replicate 682 out of 1100
    Running replicate 683 out of 1100
    Running replicate 684 out of 1100
    Running replicate 685 out of 1100
    Running replicate 686 out of 1100
    Running replicate 687 out of 1100
    Running replicate 688 out of 1100
    Running replicate 689 out of 1100
    Running replicate 690 out of 1100
    Running replicate 691 out of 1100
    Running replicate 692 out of 1100
    Running replicate 693 out of 1100
    Running replicate 694 out of 1100
    Running replicate 695 out of 1100
    Running replicate 696 out of 1100
    Running replicate 697 out of 1100
    Running replicate 698 out of 1100
    Running replicate 699 out of 1100
    Running replicate 700 out of 1100
    Running replicate 701 out of 1100
    Running replicate 702 out of 1100
    Running replicate 703 out of 1100
    Running replicate 704 out of 1100
    Running replicate 705 out of 1100
    Running replicate 706 out of 1100
    Running replicate 707 out of 1100
    Running replicate 708 out of 1100
    Running replicate 709 out of 1100
    Running replicate 710 out of 1100
    Running replicate 711 out of 1100
    Running replicate 712 out of 1100
    Running replicate 713 out of 1100
    Running replicate 714 out of 1100
    Running replicate 715 out of 1100
    Running replicate 716 out of 1100
    Running replicate 717 out of 1100
    Running replicate 718 out of 1100
    Running replicate 719 out of 1100
    Running replicate 720 out of 1100
    Running replicate 721 out of 1100
    Running replicate 722 out of 1100
    Running replicate 723 out of 1100
    Running replicate 724 out of 1100
    Running replicate 725 out of 1100
    Running replicate 726 out of 1100
    Running replicate 727 out of 1100
    Running replicate 728 out of 1100
    Running replicate 729 out of 1100
    Running replicate 730 out of 1100
    Running replicate 731 out of 1100
    Running replicate 732 out of 1100
    Running replicate 733 out of 1100
    Running replicate 734 out of 1100
    Running replicate 735 out of 1100
    Running replicate 736 out of 1100
    Running replicate 737 out of 1100
    Running replicate 738 out of 1100
    Running replicate 739 out of 1100
    Running replicate 740 out of 1100
    Running replicate 741 out of 1100
    Running replicate 742 out of 1100
    Running replicate 743 out of 1100
    Running replicate 744 out of 1100
    Running replicate 745 out of 1100
    Running replicate 746 out of 1100
    Running replicate 747 out of 1100
    Running replicate 748 out of 1100
    Running replicate 749 out of 1100
    Running replicate 750 out of 1100
    Running replicate 751 out of 1100
    Running replicate 752 out of 1100
    Running replicate 753 out of 1100
    Running replicate 754 out of 1100
    Running replicate 755 out of 1100
    Running replicate 756 out of 1100
    Running replicate 757 out of 1100
    Running replicate 758 out of 1100
    Running replicate 759 out of 1100
    Running replicate 760 out of 1100
    Running replicate 761 out of 1100
    Running replicate 762 out of 1100
    Running replicate 763 out of 1100
    Running replicate 764 out of 1100
    Running replicate 765 out of 1100
    Running replicate 766 out of 1100
    Running replicate 767 out of 1100
    Running replicate 768 out of 1100
    Running replicate 769 out of 1100
    Running replicate 770 out of 1100
    Running replicate 771 out of 1100
    Running replicate 772 out of 1100
    Running replicate 773 out of 1100
    Running replicate 774 out of 1100
    Running replicate 775 out of 1100
    Running replicate 776 out of 1100
    Running replicate 777 out of 1100
    Running replicate 778 out of 1100
    Running replicate 779 out of 1100
    Running replicate 780 out of 1100
    Running replicate 781 out of 1100
    Running replicate 782 out of 1100
    Running replicate 783 out of 1100
    Running replicate 784 out of 1100
    Running replicate 785 out of 1100
    Running replicate 786 out of 1100
    Running replicate 787 out of 1100
    Running replicate 788 out of 1100
    Running replicate 789 out of 1100
    Running replicate 790 out of 1100
    Running replicate 791 out of 1100
    Running replicate 792 out of 1100
    Running replicate 793 out of 1100
    Running replicate 794 out of 1100
    Running replicate 795 out of 1100
    Running replicate 796 out of 1100
    Running replicate 797 out of 1100
    Running replicate 798 out of 1100
    Running replicate 799 out of 1100
    Running replicate 800 out of 1100
    Running replicate 801 out of 1100
    Running replicate 802 out of 1100
    Running replicate 803 out of 1100
    Running replicate 804 out of 1100
    Running replicate 805 out of 1100
    Running replicate 806 out of 1100
    Running replicate 807 out of 1100
    Running replicate 808 out of 1100
    Running replicate 809 out of 1100
    Running replicate 810 out of 1100
    Running replicate 811 out of 1100
    Running replicate 812 out of 1100
    Running replicate 813 out of 1100
    Running replicate 814 out of 1100
    Running replicate 815 out of 1100
    Running replicate 816 out of 1100
    Running replicate 817 out of 1100
    Running replicate 818 out of 1100
    Running replicate 819 out of 1100
    Running replicate 820 out of 1100
    Running replicate 821 out of 1100
    Running replicate 822 out of 1100
    Running replicate 823 out of 1100
    Running replicate 824 out of 1100
    Running replicate 825 out of 1100
    Running replicate 826 out of 1100
    Running replicate 827 out of 1100
    Running replicate 828 out of 1100
    Running replicate 829 out of 1100
    Running replicate 830 out of 1100
    Running replicate 831 out of 1100
    Running replicate 832 out of 1100
    Running replicate 833 out of 1100
    Running replicate 834 out of 1100
    Running replicate 835 out of 1100
    Running replicate 836 out of 1100
    Running replicate 837 out of 1100
    Running replicate 838 out of 1100
    Running replicate 839 out of 1100
    Running replicate 840 out of 1100
    Running replicate 841 out of 1100
    Running replicate 842 out of 1100
    Running replicate 843 out of 1100
    Running replicate 844 out of 1100
    Running replicate 845 out of 1100
    Running replicate 846 out of 1100
    Running replicate 847 out of 1100
    Running replicate 848 out of 1100
    Running replicate 849 out of 1100
    Running replicate 850 out of 1100
    Running replicate 851 out of 1100
    Running replicate 852 out of 1100
    Running replicate 853 out of 1100
    Running replicate 854 out of 1100
    Running replicate 855 out of 1100
    Running replicate 856 out of 1100
    Running replicate 857 out of 1100
    Running replicate 858 out of 1100
    Running replicate 859 out of 1100
    Running replicate 860 out of 1100
    Running replicate 861 out of 1100
    Running replicate 862 out of 1100
    Running replicate 863 out of 1100
    Running replicate 864 out of 1100
    Running replicate 865 out of 1100
    Running replicate 866 out of 1100
    Running replicate 867 out of 1100
    Running replicate 868 out of 1100
    Running replicate 869 out of 1100
    Running replicate 870 out of 1100
    Running replicate 871 out of 1100
    Running replicate 872 out of 1100
    Running replicate 873 out of 1100
    Running replicate 874 out of 1100
    Running replicate 875 out of 1100
    Running replicate 876 out of 1100
    Running replicate 877 out of 1100
    Running replicate 878 out of 1100
    Running replicate 879 out of 1100
    Running replicate 880 out of 1100
    Running replicate 881 out of 1100
    Running replicate 882 out of 1100
    Running replicate 883 out of 1100
    Running replicate 884 out of 1100
    Running replicate 885 out of 1100
    Running replicate 886 out of 1100
    Running replicate 887 out of 1100
    Running replicate 888 out of 1100
    Running replicate 889 out of 1100
    Running replicate 890 out of 1100
    Running replicate 891 out of 1100
    Running replicate 892 out of 1100
    Running replicate 893 out of 1100
    Running replicate 894 out of 1100
    Running replicate 895 out of 1100
    Running replicate 896 out of 1100
    Running replicate 897 out of 1100
    Running replicate 898 out of 1100
    Running replicate 899 out of 1100
    Running replicate 900 out of 1100
    Running replicate 901 out of 1100
    Running replicate 902 out of 1100
    Running replicate 903 out of 1100
    Running replicate 904 out of 1100
    Running replicate 905 out of 1100
    Running replicate 906 out of 1100
    Running replicate 907 out of 1100
    Running replicate 908 out of 1100
    Running replicate 909 out of 1100
    Running replicate 910 out of 1100
    Running replicate 911 out of 1100
    Running replicate 912 out of 1100
    Running replicate 913 out of 1100
    Running replicate 914 out of 1100
    Running replicate 915 out of 1100
    Running replicate 916 out of 1100
    Running replicate 917 out of 1100
    Running replicate 918 out of 1100
    Running replicate 919 out of 1100
    Running replicate 920 out of 1100
    Running replicate 921 out of 1100
    Running replicate 922 out of 1100
    Running replicate 923 out of 1100
    Running replicate 924 out of 1100
    Running replicate 925 out of 1100
    Running replicate 926 out of 1100
    Running replicate 927 out of 1100
    Running replicate 928 out of 1100
    Running replicate 929 out of 1100
    Running replicate 930 out of 1100
    Running replicate 931 out of 1100
    Running replicate 932 out of 1100
    Running replicate 933 out of 1100
    Running replicate 934 out of 1100
    Running replicate 935 out of 1100
    Running replicate 936 out of 1100
    Running replicate 937 out of 1100
    Running replicate 938 out of 1100
    Running replicate 939 out of 1100
    Running replicate 940 out of 1100
    Running replicate 941 out of 1100
    Running replicate 942 out of 1100
    Running replicate 943 out of 1100
    Running replicate 944 out of 1100
    Running replicate 945 out of 1100
    Running replicate 946 out of 1100
    Running replicate 947 out of 1100
    Running replicate 948 out of 1100
    An error occurred in replicate 948: Exception encountered at parameter values {'n_initial': 'Autograd ArrayBox with value nan', 'n_korea_decline1': 'Autograd ArrayBox with value nan', 'n_korea_expand1': 'Autograd ArrayBox with value nan', 'n_us_decline2': 'Autograd ArrayBox with value nan', 'n_us_expand2': 'Autograd ArrayBox with value nan', 't_korea_decline1': 'Autograd ArrayBox with value nan', 't_korea_expand1': 'Autograd ArrayBox with value nan', 't_us_decline2': 'Autograd ArrayBox with value nan', 't_us_expand2': 'Autograd ArrayBox with value nan'} (internal scaling x = {'n_initial': 'Autograd ArrayBox with value nan', 'n_korea_decline1': 'Autograd ArrayBox with value nan', 'n_korea_expand1': 'Autograd ArrayBox with value nan', 'n_us_decline2': 'Autograd ArrayBox with value nan', 'n_us_expand2': 'Autograd ArrayBox with value nan', 't_korea_decline1': 'Autograd ArrayBox with value nan', 't_korea_expand1': 'Autograd ArrayBox with value nan', 't_us_decline2': 'Autograd ArrayBox with value nan', 't_us_expand2': 'Autograd ArrayBox with value nan'})
    Skipping this replicate and moving on to the next one.
    Running replicate 949 out of 1100
    Running replicate 950 out of 1100
    Running replicate 951 out of 1100
    Running replicate 952 out of 1100
    Running replicate 953 out of 1100
    Running replicate 954 out of 1100
    Running replicate 955 out of 1100
    Running replicate 956 out of 1100
    Running replicate 957 out of 1100
    Running replicate 958 out of 1100
    Running replicate 959 out of 1100
    Running replicate 960 out of 1100
    Running replicate 961 out of 1100
    Running replicate 962 out of 1100
    Running replicate 963 out of 1100
    Running replicate 964 out of 1100
    Running replicate 965 out of 1100
    Running replicate 966 out of 1100
    Running replicate 967 out of 1100
    Running replicate 968 out of 1100
    Running replicate 969 out of 1100
    Running replicate 970 out of 1100
    Running replicate 971 out of 1100
    Running replicate 972 out of 1100
    Running replicate 973 out of 1100
    Running replicate 974 out of 1100
    Running replicate 975 out of 1100
    Running replicate 976 out of 1100
    Running replicate 977 out of 1100
    Running replicate 978 out of 1100
    Running replicate 979 out of 1100
    Running replicate 980 out of 1100
    Running replicate 981 out of 1100
    Running replicate 982 out of 1100
    Running replicate 983 out of 1100
    Running replicate 984 out of 1100
    Running replicate 985 out of 1100
    Running replicate 986 out of 1100
    Running replicate 987 out of 1100
    Running replicate 988 out of 1100
    Running replicate 989 out of 1100
    Running replicate 990 out of 1100
    Running replicate 991 out of 1100
    Running replicate 992 out of 1100
    Running replicate 993 out of 1100
    Running replicate 994 out of 1100
    Running replicate 995 out of 1100
    Running replicate 996 out of 1100
    Running replicate 997 out of 1100
    Running replicate 998 out of 1100
    Running replicate 999 out of 1100
    Running replicate 1000 out of 1100
    Running replicate 1001 out of 1100
    Running replicate 1002 out of 1100
    Running replicate 1003 out of 1100
    Running replicate 1004 out of 1100
    Running replicate 1005 out of 1100
    Running replicate 1006 out of 1100
    Running replicate 1007 out of 1100
    Running replicate 1008 out of 1100
    Running replicate 1009 out of 1100
    Running replicate 1010 out of 1100
    Running replicate 1011 out of 1100
    Running replicate 1012 out of 1100
    Running replicate 1013 out of 1100
    Running replicate 1014 out of 1100
    Running replicate 1015 out of 1100
    Running replicate 1016 out of 1100
    Running replicate 1017 out of 1100
    Running replicate 1018 out of 1100
    Running replicate 1019 out of 1100
    Running replicate 1020 out of 1100
    Running replicate 1021 out of 1100
    Running replicate 1022 out of 1100
    Running replicate 1023 out of 1100
    Running replicate 1024 out of 1100
    Running replicate 1025 out of 1100
    Running replicate 1026 out of 1100
    Running replicate 1027 out of 1100
    Running replicate 1028 out of 1100
    Running replicate 1029 out of 1100
    Running replicate 1030 out of 1100
    Running replicate 1031 out of 1100
    Running replicate 1032 out of 1100
    Running replicate 1033 out of 1100
    Running replicate 1034 out of 1100
    Running replicate 1035 out of 1100
    Running replicate 1036 out of 1100
    Running replicate 1037 out of 1100
    Running replicate 1038 out of 1100
    Running replicate 1039 out of 1100
    Running replicate 1040 out of 1100
    Running replicate 1041 out of 1100
    Running replicate 1042 out of 1100
    Running replicate 1043 out of 1100
    Running replicate 1044 out of 1100
    Running replicate 1045 out of 1100
    Running replicate 1046 out of 1100
    Running replicate 1047 out of 1100
    Running replicate 1048 out of 1100
    Running replicate 1049 out of 1100
    Running replicate 1050 out of 1100
    Running replicate 1051 out of 1100
    Running replicate 1052 out of 1100
    Running replicate 1053 out of 1100
    Running replicate 1054 out of 1100
    Running replicate 1055 out of 1100
    Running replicate 1056 out of 1100
    Running replicate 1057 out of 1100
    Running replicate 1058 out of 1100
    Running replicate 1059 out of 1100
    Running replicate 1060 out of 1100
    Running replicate 1061 out of 1100
    Running replicate 1062 out of 1100
    Running replicate 1063 out of 1100
    Running replicate 1064 out of 1100
    Running replicate 1065 out of 1100
    Running replicate 1066 out of 1100
    Running replicate 1067 out of 1100
    Running replicate 1068 out of 1100
    Running replicate 1069 out of 1100
    Running replicate 1070 out of 1100
    Running replicate 1071 out of 1100
    Running replicate 1072 out of 1100
    Running replicate 1073 out of 1100
    Running replicate 1074 out of 1100
    Running replicate 1075 out of 1100
    Running replicate 1076 out of 1100
    Running replicate 1077 out of 1100
    Running replicate 1078 out of 1100
    Running replicate 1079 out of 1100
    Running replicate 1080 out of 1100
    Running replicate 1081 out of 1100
    Running replicate 1082 out of 1100
    Running replicate 1083 out of 1100
    Running replicate 1084 out of 1100
    Running replicate 1085 out of 1100
    Running replicate 1086 out of 1100
    Running replicate 1087 out of 1100
    Running replicate 1088 out of 1100
    Running replicate 1089 out of 1100
    Running replicate 1090 out of 1100
    Running replicate 1091 out of 1100
    Running replicate 1092 out of 1100
    Running replicate 1093 out of 1100
    Running replicate 1094 out of 1100
    Running replicate 1095 out of 1100
    Running replicate 1096 out of 1100
    Running replicate 1097 out of 1100
    Running replicate 1098 out of 1100
    Running replicate 1099 out of 1100
    Running replicate 1100 out of 1100



```python
model_complex_global.fit_within_pop_diversity() #Estimating mutation rate,To evaluate the total number of mutations in the data, e.g. to fit the mutation rate, use the method DemographicModel.fit_within_pop_diversity, which computes the within-population nucleotide diversity, i.e. the heterozygosity of a random individual in that population assuming Hardy-Weinberg Equilibrium:
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Pop</th>
      <th>EstMutRate</th>
      <th>JackknifeSD</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>US</td>
      <td>6.538848e-10</td>
      <td>6.857533e-12</td>
    </tr>
  </tbody>
</table>
</div>




```python
import numpy as np
from momi import SfsModelFitStats

fit_stats = momi.SfsModelFitStats(model_complex_global)

num_samples = sfs.sampled_n[0]  # Number of samples for the "US" population'
print("num_samples",num_samples)

weights = np.linspace(1, 0, num_samples + 1)  # Generate weights for n+1 categories

# Compute the tensor product with the defined weights
tensor_stat = fit_stats.tensor_prod(derived_weights_dict={"US": weights})

# Extract standard deviation and Z-score
print("Standard Deviation:", tensor_stat.sd)
print("Z-Score:", tensor_stat.z_score)

```

    num_samples 196
    Standard Deviation: 0.01732579974556093
    Z-Score: -0.5127106187943128



```python
#Standard Deviation (SD): This value (0.017) represents the variability or dispersion of the jackknifed values of the tensor product statistic from their mean. A lower standard deviation indicates that the jackknifed values are closer to their mean, suggesting more certainty in the estimate. In demographic inference, this can imply a more consistent pattern across different parts of the genome, thus more confidence in the demographic model's fit to the data.
#A Z-score near zero indicates that the observed data aligns closely with the model's expectations. 
```


```python
#Plot CI interval for 1100 bootstrap for the best model(three bottlenecks)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Read the statistics data
stats_df = pd.read_csv('stats.momi2.wide.csv')

# Set the 'parameter' column as the index
stats_pivot = stats_df.set_index('parameter')

# Time points for the plot
time_points = np.linspace(-80, 0, 1000)  # Covering 100 years to 2022

# Initialize population size arrays for mean, lower CI, and upper CI
population_size_mean = np.zeros_like(time_points)
population_size_lower = np.zeros_like(time_points)
population_size_upper = np.zeros_like(time_points)

# Define a function to determine the population size at each time point
def get_population_size(t):
    if t <= -stats_pivot.at['t_shanghai_decline1', 'Mean']:
        return ('n_initial',)
    elif -stats_pivot.at['t_shanghai_decline1', 'Mean'] < t <= -stats_pivot.at['t_shanghai_expand1', 'Mean']:
        return ('n_shanghai_decline1',)
    elif -stats_pivot.at['t_shanghai_expand1', 'Mean'] < t <= -stats_pivot.at['t_korea_decline2', 'Mean']:
        return ('n_shanghai_expand1',)
    elif -stats_pivot.at['t_korea_decline2', 'Mean'] < t <= -stats_pivot.at['t_korea_expand2', 'Mean']:
        return ('n_korea_decline2',)
    elif -stats_pivot.at['t_korea_expand2', 'Mean'] < t <= -stats_pivot.at['t_us_decline3', 'Mean']:
        return ('n_korea_expand2',)
    elif -stats_pivot.at['t_us_decline3', 'Mean'] < t <= -stats_pivot.at['t_us_expand3', 'Mean']:
        return ('n_us_decline3',)
    else:
        return ('n_us_expand3',)

# Calculate the mean, lower CI, and upper CI for each time point
for i, t in enumerate(time_points):
    param = get_population_size(t)
    population_size_lower[i] = stats_pivot.at[param[0], '95% CI Lower']
    population_size_upper[i] = stats_pivot.at[param[0], '95% CI Upper']
    population_size_mean[i] = stats_pivot.at[param[0], 'Mean']

# Plotting
plt.figure(figsize=(10, 6))
#plt.plot(time_points, population_size_mean, color='indianred', linewidth=2, label='Mean Values')
plt.fill_between(time_points, population_size_lower, population_size_upper, color='cornflowerblue', alpha=0.5, label='95% CI')
plt.plot(time_points, population_size_mean, color='indianred', linewidth=2, label='Mean Values')
plt.legend(loc='lower left')

# Plot CI for time parameters (excluding the initial time point)
for param in ['t_shanghai_decline1', 't_shanghai_expand1', 't_korea_decline2', 't_korea_expand2', 't_us_decline3', 't_us_expand3']:
    mean = -stats_pivot.at[param, 'Mean']
    lower = -stats_pivot.at[param, '95% CI Lower']
    upper = -stats_pivot.at[param, '95% CI Upper']
    plt.axvspan(lower, upper, alpha=0.5, color='cornflowerblue')

plt.xlabel('Years before 2022', labelpad=15)
plt.ylabel('Population Size')
plt.title('SLF Population Dynamics Over Time (1100 bootstraps)')
plt.grid(False)
plt.yscale('log')  # Log scale for better visualization

plt.legend(loc='lower left')

# Saving the custom plot as a PDF
plt.savefig("SLF_Final_Three_Bottlenecks_CI_Plot.png", format='png')
plt.show()
```


    
![png](output_20_0.png)
    



```python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Read the statistics data
stats_df = pd.read_csv('stats.momi2.txt')

# Convert the DataFrame into a pivot table with parameters as index, statistics as columns
stats_pivot = stats_df.pivot(index='parameter', columns='statistic', values='value')

# Time points for the plot
time_points = np.linspace(-100, 0, 1000)  # Covering 100 years to present

# Initialize population size arrays for mean, lower CI, and upper CI
population_size_mean = np.zeros_like(time_points)
population_size_lower = np.zeros_like(time_points)
population_size_upper = np.zeros_like(time_points)

# Define a function to determine the population size at each time point
def get_population_size(t):
    if t <= -stats_pivot.at['t_shanghai_decline1', 'mean']:
        return ('n_initial',)
    elif -stats_pivot.at['t_shanghai_decline1', 'mean'] < t <= -stats_pivot.at['t_shanghai_expand1', 'mean']:
        return ('n_shanghai_decline1',)
    elif -stats_pivot.at['t_shanghai_expand1', 'mean'] < t <= -stats_pivot.at['t_korea_decline2', 'mean']:
        return ('n_shanghai_expand1',)
    elif -stats_pivot.at['t_korea_decline2', 'mean'] < t <= -stats_pivot.at['t_korea_expand2', 'mean']:
        return ('n_korea_decline2',)
    elif -stats_pivot.at['t_korea_expand2', 'mean'] < t <= -stats_pivot.at['t_us_decline3', 'mean']:
        return ('n_korea_expand2',)
    elif -stats_pivot.at['t_us_decline3', 'mean'] < t <= -stats_pivot.at['t_us_expand3', 'mean']:
        return ('n_us_decline3',)
    else:
        return ('n_us_expand3',)

# Calculate the mean, lower CI, and upper CI for each time point
for i, t in enumerate(time_points):
    param = get_population_size(t)
    population_size_mean[i] = stats_pivot.at[param[0], 'mean']
    population_size_lower[i] = stats_pivot.at[param[0], 'lower']
    population_size_upper[i] = stats_pivot.at[param[0], 'upper']

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(time_points, population_size_mean, color='red', linewidth=3, label='Mean Values')
plt.fill_between(time_points, population_size_lower, population_size_upper, color='blue', alpha=0.5, label='95% CI')
plt.xlabel('Years before present in USA', labelpad=10)
plt.ylabel('Population Size')
plt.title('SLF Population Dynamics Over Time (1300 bootstraps)')
plt.grid(False)
plt.yscale('log')  # Log scale for better visualization

# Add the vertical lines and text annotations for specified events
params = ['t_us_decline3', 't_korea_decline2', 't_shanghai_decline1']
labels = ['Arrived USA', 'Arrived South Korea', 'Urbanization of Shanghai']
n_params = ['n_us_decline3', 'n_korea_decline2', 'n_shanghai_decline1']

for param, label, n_param in zip(params, labels, n_params):
    mean = -stats_pivot.at[param, 'mean']
    pop_mean = stats_pivot.at[n_param, 'mean']
    plt.axvline(x=mean, color='lightgreen', linestyle='--', linewidth=1)
    plt.text(mean, pop_mean, label, horizontalalignment='left', verticalalignment='bottom')

# Saving the custom plot as a PDF
plt.savefig("SLF_Complex_Global_Model_CI_Plot.pdf", format='pdf')
plt.show()
```


    
![png](output_21_0.png)
    



```python

```
