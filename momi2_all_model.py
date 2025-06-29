#!/usr/bin/env python3
"""
Bootstrap Analysis for 8 Spotted Lanternfly Invasion Models:
1. Urban origin with Korea (Bridgehead Model): SH_rural -> SH_urban -> Korea -> USA
2. Urban origin direct (Direct Invasion Model): SH_rural -> SH_urban -> USA
3. Rural origin with Korea: SH_rural -> Korea -> USA and SH_rural -> SH_urban
4. Rural origin with direct paths: SH_rural -> USA and SH_rural -> SH_urban
5. Independent Urban with Korea: SH_rural and SH_urban (parallel), SH_urban -> Korea -> USA
6. Independent Urban Direct: SH_rural and SH_urban (parallel), SH_urban -> USA
7. Independent Rural with Korea: SH_rural and SH_urban (parallel), SH_rural -> Korea -> USA
8. Independent Rural Direct: SH_rural and SH_urban (parallel), SH_rural -> USA
"""

import momi
import numpy as np
import logging
import sys
import os
import time
from datetime import datetime

# Setup logging
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
output_dir = f"bootstrap_results_{timestamp}"
os.makedirs(output_dir, exist_ok=True)

log_file = f"{output_dir}/slf_models_bootstrap_{timestamp}.log"

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout)
    ]
)

# Set random seed for reproducibility
np.random.seed(3113)

# Constants
GEN_TIME = 1  # One generation per year for spotted lanternfly
MU = 2.8e-9   # Mutation rate per site per generation
NE_INIT = 1e4  # Initial effective population size

def load_sfs_data(sfs_path):
    """Load Site Frequency Spectrum data"""
    try:
        logging.info(f"Loading SFS data from {sfs_path}...")
        sfs_data = momi.Sfs.load(sfs_path)
        logging.info(f"SFS loaded successfully with populations: {sfs_data.sampled_pops}")
        logging.info(f"Number of SNPs: {sfs_data.n_snps()}")
        return sfs_data
    except Exception as e:
        logging.error(f"Failed to load SFS: {e}")
        raise

# Model 1: Full Bridgehead Model (Shanghai Rural → Shanghai Urban → South Korea (ghost) → USA)
def create_bridgehead_model(sfs_data, randomize=False):
    """Create the full bridgehead model with South Korea as a ghost population"""
    model = momi.DemographicModel(N_e=NE_INIT, gen_time=GEN_TIME, muts_per_gen=MU)
    
    # Add time parameters with chronological constraints
    model.add_time_param("t_usa_split", t0=10, lower=8, upper=20)       # USA invasion ~2014
    model.add_time_param("t_korea_split", t0=16, lower_constraints=["t_usa_split"], upper=30)  # Korea invasion ~2006
    model.add_time_param("t_urban_split", t0=50, lower_constraints=["t_korea_split"], upper=1e6)  # Urbanization
    
    # Add contemporary population size parameters
    model.add_size_param("N_rural", lower=1e3, upper=1e7)     # Rural Shanghai current size
    model.add_size_param("N_urban", lower=1e3, upper=1e7)     # Urban Shanghai current size
    model.add_size_param("N_usa", lower=1e3, upper=1e7)       # USA current size
    
    # Add bottleneck size parameters
    model.add_size_param("N_urban_founder", lower=10, upper=1e3)  # Bottleneck during urbanization
    model.add_size_param("N_korea_founder", lower=10, upper=1e3)  # Bottleneck during Korea invasion
    model.add_size_param("N_usa_founder", lower=10, upper=1e3)    # Bottleneck during USA invasion
    
    # Define contemporary populations 
    model.add_leaf("SH_rural", N="N_rural")
    model.add_leaf("SH_urban", N="N_urban")
    model.add_leaf("USA", N="N_usa")
    
    # Set the data
    model.set_data(sfs_data)
    
    # Define population splits - working BACKWARDS in time
    # First, USA split from Korea with a bottleneck
    model.set_size("USA", t="t_usa_split", N="N_usa_founder")
    model.move_lineages("USA", "Korea", t="t_usa_split")  # Korea implicitly created as ghost
    
    # Set properties of the ghost Korea population
    model.set_size("Korea", t=0, N="N_korea_founder")  # Small founder population for Korea
    
    # Then, Korea split from SH_urban
    model.move_lineages("Korea", "SH_urban", t="t_korea_split")
    
    # Finally, SH_urban split from SH_rural with a bottleneck
    model.set_size("SH_urban", t="t_urban_split", N="N_urban_founder")
    model.move_lineages("SH_urban", "SH_rural", t="t_urban_split")
    
    # Randomize parameters if requested (for bootstrapping)
    if randomize:
        model.set_params(randomize=True)
    
    return model

# Model 2: Direct Invasion Model (Shanghai Rural → Shanghai Urban → USA)
def create_direct_model(sfs_data, randomize=False):
    """Create the direct invasion model without South Korea"""
    model = momi.DemographicModel(N_e=NE_INIT, gen_time=GEN_TIME, muts_per_gen=MU)
    
    # Add time parameters with chronological constraints
    model.add_time_param("t_usa_split", t0=10, lower=8, upper=20)       # USA invasion ~2014
    model.add_time_param("t_urban_split", t0=50, lower_constraints=["t_usa_split"], upper=1e6)  # Urbanization
    
    # Add contemporary population size parameters
    model.add_size_param("N_rural", lower=1e3, upper=1e7)     # Rural Shanghai current size
    model.add_size_param("N_urban", lower=1e3, upper=1e7)     # Urban Shanghai current size
    model.add_size_param("N_usa", lower=1e3, upper=1e7)       # USA current size
    
    # Add bottleneck size parameters
    model.add_size_param("N_urban_founder", lower=10, upper=1e3)  # Bottleneck during urbanization
    model.add_size_param("N_usa_founder", lower=10, upper=1e3)    # Bottleneck during USA invasion
    
    # Define contemporary populations 
    model.add_leaf("SH_rural", N="N_rural")
    model.add_leaf("SH_urban", N="N_urban")
    model.add_leaf("USA", N="N_usa")
    
    # Set the data
    model.set_data(sfs_data)
    
    # Define population splits - working BACKWARDS in time
    # First, USA split directly from SH_urban with a bottleneck
    model.set_size("USA", t="t_usa_split", N="N_usa_founder")
    model.move_lineages("USA", "SH_urban", t="t_usa_split")
    
    # Then, SH_urban split from SH_rural with a bottleneck
    model.set_size("SH_urban", t="t_urban_split", N="N_urban_founder")
    model.move_lineages("SH_urban", "SH_rural", t="t_urban_split")
    
    # Randomize parameters if requested (for bootstrapping)
    if randomize:
        model.set_params(randomize=True)
    
    return model

# Model 3: Rural Origin with Korea (Shanghai Rural → Korea → USA, Shanghai Rural → Shanghai Urban)
def create_rural_origin_with_korea_model(sfs_data, randomize=False):
    """Create model with rural origin, Korea as intermediary to USA, and direct path to urban"""
    model = momi.DemographicModel(N_e=NE_INIT, gen_time=GEN_TIME, muts_per_gen=MU)
    
    # Add time parameters
    model.add_time_param("t_usa_split", t0=10, lower=8, upper=20)       # USA invasion ~2014
    model.add_time_param("t_korea_split", t0=16, lower_constraints=["t_usa_split"], upper=30)  # Korea invasion ~2006
    model.add_time_param("t_urban_split", t0=50, lower_constraints=["t_korea_split"], upper=1e6)    # Urbanization
    
    # Add contemporary population size parameters
    model.add_size_param("N_rural", lower=1e3, upper=1e7)     # Rural Shanghai current size
    model.add_size_param("N_urban", lower=1e3, upper=1e7)     # Urban Shanghai current size
    model.add_size_param("N_usa", lower=1e3, upper=1e7)       # USA current size
    
    # Add bottleneck size parameters
    model.add_size_param("N_urban_founder", lower=10, upper=1e3)  # Bottleneck during urbanization
    model.add_size_param("N_korea_founder", lower=10, upper=1e3)  # Bottleneck during Korea invasion
    model.add_size_param("N_usa_founder", lower=10, upper=1e3)    # Bottleneck during USA invasion
    
    # Define contemporary populations 
    model.add_leaf("SH_rural", N="N_rural")
    model.add_leaf("SH_urban", N="N_urban")
    model.add_leaf("USA", N="N_usa")
    
    # Set the data
    model.set_data(sfs_data)
    
    # Define population splits - working BACKWARDS in time
    # USA split from Korea with a bottleneck
    model.set_size("USA", t="t_usa_split", N="N_usa_founder")
    model.move_lineages("USA", "Korea", t="t_usa_split")  # Korea implicitly created as ghost
    
    # Set properties of the ghost Korea population
    model.set_size("Korea", t=0, N="N_korea_founder")
    
    # Korea split from SH_rural
    model.move_lineages("Korea", "SH_rural", t="t_korea_split")
    
    # SH_urban split from SH_rural (independently of Korea/USA)
    model.set_size("SH_urban", t="t_urban_split", N="N_urban_founder")
    model.move_lineages("SH_urban", "SH_rural", t="t_urban_split")
    
    # Randomize parameters if requested
    if randomize:
        model.set_params(randomize=True)
    
    return model

# Model 4: Rural Origin with Direct Paths (Shanghai Rural → USA, Shanghai Rural → Shanghai Urban)
def create_rural_origin_model(sfs_data, randomize=False):
    """Create model with rural origin and direct paths to USA and urban"""
    model = momi.DemographicModel(N_e=NE_INIT, gen_time=GEN_TIME, muts_per_gen=MU)
    
    # Add time parameters
    model.add_time_param("t_usa_split", t0=10, lower=8, upper=20)       # USA invasion ~2014
    model.add_time_param("t_urban_split", t0=50, lower_constraints=["t_usa_split"], upper=1e6)    # Urbanization
    
    # Add contemporary population size parameters
    model.add_size_param("N_rural", lower=1e3, upper=1e7)     # Rural Shanghai current size
    model.add_size_param("N_urban", lower=1e3, upper=1e7)     # Urban Shanghai current size
    model.add_size_param("N_usa", lower=1e3, upper=1e7)       # USA current size
    
    # Add bottleneck size parameters
    model.add_size_param("N_urban_founder", lower=10, upper=1e3)  # Bottleneck during urbanization
    model.add_size_param("N_usa_founder", lower=10, upper=1e3)    # Bottleneck during USA invasion
    
    # Define contemporary populations 
    model.add_leaf("SH_rural", N="N_rural")
    model.add_leaf("SH_urban", N="N_urban")
    model.add_leaf("USA", N="N_usa")
    
    # Set the data
    model.set_data(sfs_data)
    
    # Define population splits - working BACKWARDS in time
    # USA split directly from SH_rural with a bottleneck
    model.set_size("USA", t="t_usa_split", N="N_usa_founder")
    model.move_lineages("USA", "SH_rural", t="t_usa_split")
    
    # SH_urban split from SH_rural with a bottleneck (independently of USA split)
    model.set_size("SH_urban", t="t_urban_split", N="N_urban_founder")
    model.move_lineages("SH_urban", "SH_rural", t="t_urban_split")
    
    # Randomize parameters if requested
    if randomize:
        model.set_params(randomize=True)
    
    return model

# Model 5: Independent Urban with Korea (SH_rural and SH_urban parallel, SH_urban -> Korea -> USA)
def create_independent_urban_with_korea_model(sfs_data, randomize=False):
    """Create model with independent rural and urban populations, with Korea as a bridgehead"""
    model = momi.DemographicModel(N_e=NE_INIT, gen_time=GEN_TIME, muts_per_gen=MU)
    
    # Add time parameters with chronological constraints
    model.add_time_param("t_usa_split", t0=10, lower=8, upper=20)       # USA invasion ~2014
    model.add_time_param("t_korea_split", t0=16, lower_constraints=["t_usa_split"], upper=30)  # Korea invasion ~2006
    model.add_time_param("t_rural_urban_split", t0=50, lower_constraints=["t_korea_split"], upper=1e6)  # Ancient split between rural/urban
    
    # Add contemporary population size parameters
    model.add_size_param("N_rural", lower=1e3, upper=1e7)     # Rural Shanghai current size
    model.add_size_param("N_urban", lower=1e3, upper=1e7)     # Urban Shanghai current size
    model.add_size_param("N_usa", lower=1e3, upper=1e7)       # USA current size
    model.add_size_param("N_ancestral", lower=1e3, upper=1e7)  # Ancestral population size
    
    # Add bottleneck size parameters
    model.add_size_param("N_korea_founder", lower=10, upper=1e3)  # Bottleneck during Korea invasion
    model.add_size_param("N_usa_founder", lower=10, upper=1e3)    # Bottleneck during USA invasion
    
    # Define contemporary populations 
    model.add_leaf("SH_rural", N="N_rural")
    model.add_leaf("SH_urban", N="N_urban")
    model.add_leaf("USA", N="N_usa")
    
    # Set the data
    model.set_data(sfs_data)
    
    # Define population splits - working BACKWARDS in time
    # First, USA split from Korea with a bottleneck
    model.set_size("USA", t="t_usa_split", N="N_usa_founder")
    model.move_lineages("USA", "Korea", t="t_usa_split")  # Korea implicitly created as ghost
    
    # Set properties of the ghost Korea population
    model.set_size("Korea", t=0, N="N_korea_founder")  # Small founder population for Korea
    
    # Then, Korea split from SH_urban
    model.move_lineages("Korea", "SH_urban", t="t_korea_split")
    
    # Finally, SH_urban and SH_rural merge in ancestral population
    model.move_lineages("SH_urban", "Ancestral", t="t_rural_urban_split")  # Create ancestral ghost population
    model.move_lineages("SH_rural", "Ancestral", t="t_rural_urban_split")
    model.set_size("Ancestral", t=0, N="N_ancestral")
    
    # Randomize parameters if requested (for bootstrapping)
    if randomize:
        model.set_params(randomize=True)
    
    return model

# Model 6: Independent Urban Direct (SH_rural and SH_urban parallel, SH_urban -> USA)
def create_independent_urban_direct_model(sfs_data, randomize=False):
    """Create model with independent rural and urban populations, direct invasion from urban"""
    model = momi.DemographicModel(N_e=NE_INIT, gen_time=GEN_TIME, muts_per_gen=MU)
    
    # Add time parameters with chronological constraints
    model.add_time_param("t_usa_split", t0=10, lower=8, upper=20)       # USA invasion ~2014
    model.add_time_param("t_rural_urban_split", t0=50, lower_constraints=["t_usa_split"], upper=1e6)  # Ancient split between rural/urban
    
    # Add contemporary population size parameters
    model.add_size_param("N_rural", lower=1e3, upper=1e7)     # Rural Shanghai current size
    model.add_size_param("N_urban", lower=1e3, upper=1e7)     # Urban Shanghai current size
    model.add_size_param("N_usa", lower=1e3, upper=1e7)       # USA current size
    model.add_size_param("N_ancestral", lower=1e3, upper=1e7)  # Ancestral population size
    
    # Add bottleneck size parameters
    model.add_size_param("N_usa_founder", lower=10, upper=1e3)    # Bottleneck during USA invasion
    
    # Define contemporary populations 
    model.add_leaf("SH_rural", N="N_rural")
    model.add_leaf("SH_urban", N="N_urban")
    model.add_leaf("USA", N="N_usa")
    
    # Set the data
    model.set_data(sfs_data)
    
    # Define population splits - working BACKWARDS in time
    # First, USA split directly from SH_urban with a bottleneck
    model.set_size("USA", t="t_usa_split", N="N_usa_founder")
    model.move_lineages("USA", "SH_urban", t="t_usa_split")
    
    # Finally, SH_urban and SH_rural merge in ancestral population
    model.move_lineages("SH_urban", "Ancestral", t="t_rural_urban_split")  # Create ancestral ghost population
    model.move_lineages("SH_rural", "Ancestral", t="t_rural_urban_split")
    model.set_size("Ancestral", t=0, N="N_ancestral")
    
    # Randomize parameters if requested (for bootstrapping)
    if randomize:
        model.set_params(randomize=True)
    
    return model

# Model 7: Independent Rural with Korea (SH_rural and SH_urban parallel, SH_rural -> Korea -> USA)
def create_independent_rural_with_korea_model(sfs_data, randomize=False):
    """Create model with independent rural and urban populations, Korea as intermediary from rural"""
    model = momi.DemographicModel(N_e=NE_INIT, gen_time=GEN_TIME, muts_per_gen=MU)
    
    # Add time parameters with chronological constraints
    model.add_time_param("t_usa_split", t0=10, lower=8, upper=20)       # USA invasion ~2014
    model.add_time_param("t_korea_split", t0=16, lower_constraints=["t_usa_split"], upper=30)  # Korea invasion ~2006
    model.add_time_param("t_rural_urban_split", t0=50, lower_constraints=["t_korea_split"], upper=1e6)  # Ancient split between rural/urban
    
    # Add contemporary population size parameters
    model.add_size_param("N_rural", lower=1e3, upper=1e7)     # Rural Shanghai current size
    model.add_size_param("N_urban", lower=1e3, upper=1e7)     # Urban Shanghai current size
    model.add_size_param("N_usa", lower=1e3, upper=1e7)       # USA current size
    model.add_size_param("N_ancestral", lower=1e3, upper=1e7)  # Ancestral population size
    
    # Add bottleneck size parameters
    model.add_size_param("N_korea_founder", lower=10, upper=1e3)  # Bottleneck during Korea invasion
    model.add_size_param("N_usa_founder", lower=10, upper=1e3)    # Bottleneck during USA invasion
    
    # Define contemporary populations 
    model.add_leaf("SH_rural", N="N_rural")
    model.add_leaf("SH_urban", N="N_urban")
    model.add_leaf("USA", N="N_usa")
    
    # Set the data
    model.set_data(sfs_data)
    
    # Define population splits - working BACKWARDS in time
    # First, USA split from Korea with a bottleneck
    model.set_size("USA", t="t_usa_split", N="N_usa_founder")
    model.move_lineages("USA", "Korea", t="t_usa_split")  # Korea implicitly created as ghost
    
    # Set properties of the ghost Korea population
    model.set_size("Korea", t=0, N="N_korea_founder")  # Small founder population for Korea
    
    # Korea split from SH_rural (instead of urban in model 5)
    model.move_lineages("Korea", "SH_rural", t="t_korea_split")
    
    # Finally, SH_urban and SH_rural merge in ancestral population
    model.move_lineages("SH_urban", "Ancestral", t="t_rural_urban_split")  # Create ancestral ghost population
    model.move_lineages("SH_rural", "Ancestral", t="t_rural_urban_split")
    model.set_size("Ancestral", t=0, N="N_ancestral")
    
    # Randomize parameters if requested (for bootstrapping)
    if randomize:
        model.set_params(randomize=True)
    
    return model

# Model 8: Independent Rural Direct (SH_rural and SH_urban parallel, SH_rural -> USA)
def create_independent_rural_direct_model(sfs_data, randomize=False):
    """Create model with independent rural and urban populations, direct invasion from rural"""
    model = momi.DemographicModel(N_e=NE_INIT, gen_time=GEN_TIME, muts_per_gen=MU)
    
    # Add time parameters with chronological constraints
    model.add_time_param("t_usa_split", t0=10, lower=8, upper=20)       # USA invasion ~2014
    model.add_time_param("t_rural_urban_split", t0=50, lower_constraints=["t_usa_split"], upper=1e6)  # Ancient split between rural/urban
    
    # Add contemporary population size parameters
    model.add_size_param("N_rural", lower=1e3, upper=1e7)     # Rural Shanghai current size
    model.add_size_param("N_urban", lower=1e3, upper=1e7)     # Urban Shanghai current size
    model.add_size_param("N_usa", lower=1e3, upper=1e7)       # USA current size
    model.add_size_param("N_ancestral", lower=1e3, upper=1e7)  # Ancestral population size
    
    # Add bottleneck size parameters
    model.add_size_param("N_usa_founder", lower=10, upper=1e3)    # Bottleneck during USA invasion
    
    # Define contemporary populations 
    model.add_leaf("SH_rural", N="N_rural")
    model.add_leaf("SH_urban", N="N_urban")
    model.add_leaf("USA", N="N_usa")
    
    # Set the data
    model.set_data(sfs_data)
    
    # Define population splits - working BACKWARDS in time
    # First, USA split directly from SH_rural with a bottleneck (instead of urban in model 6)
    model.set_size("USA", t="t_usa_split", N="N_usa_founder")
    model.move_lineages("USA", "SH_rural", t="t_usa_split")
    
    # Finally, SH_urban and SH_rural merge in ancestral population
    model.move_lineages("SH_urban", "Ancestral", t="t_rural_urban_split")  # Create ancestral ghost population
    model.move_lineages("SH_rural", "Ancestral", t="t_rural_urban_split")
    model.set_size("Ancestral", t=0, N="N_ancestral")
    
    # Randomize parameters if requested (for bootstrapping)
    if randomize:
        model.set_params(randomize=True)
    
    return model

# Main execution
def main():
    logging.info("Starting bootstrap analysis for Spotted Lanternfly demographic models")
    
    try:
        # Load SFS data
        sfs_path = "/scratch/fm2371/Momi2/momi_new/allele_counts/All_pop_sfs.gz"
        sfs_data = load_sfs_data(sfs_path)
        
        # Define all model creation functions with their names
        model_creators = [
        (create_bridgehead_model, "Model1_Urban_Korea"),
        (create_direct_model, "Model2_Urban_Direct"),
    (create_rural_origin_with_korea_model, "Model3_Rural_Korea"),
    (create_rural_origin_model, "Model4_Rural_Direct"),
    (create_independent_urban_with_korea_model, "Model5_Independent_Urban_Korea"),
    (create_independent_urban_direct_model, "Model6_Independent_Urban_Direct"),
    (create_independent_rural_with_korea_model, "Model7_Independent_Rural_Korea"),
    (create_independent_rural_direct_model, "Model8_Independent_Rural_Direct")
]
        
        # Files to store results
        aic_file_path = os.path.join(output_dir, 'model_aic_and_likelihoods.txt')
        best_models_file_path = os.path.join(output_dir, 'best_models_based_on_weights.txt')
        params_file_path = os.path.join(output_dir, 'best_model_params.txt')
        
        # Number of bootstrap replicates
        num_replicates = 100
        
        with open(aic_file_path, 'w') as aic_file, \
             open(best_models_file_path, 'w') as best_models_file, \
             open(params_file_path, 'w') as params_file:
            
            # Write headers
            aic_file.write("Replicate,Model,AIC,LogLikelihood,DeltaAIC,AICWeight\n")
            best_models_file.write("Replicate,BestModelIndex,BestModelName\n")
            params_file.write("Replicate,BestModel,Parameter,Value\n")
            
            successful_replicates = 0
            
            # For storing all best model indexes
            best_model_counts = {i: 0 for i in range(len(model_creators))}
            
            for i in range(num_replicates):
                logging.info(f"Running bootstrap replicate {i+1} out of {num_replicates}")
                
                try:
                    # Resample the SFS
                    resampled_sfs = sfs_data.resample()
                    
                    aics = []
                    log_likelihoods = []
                    model_instances = []
                    
                    # Run optimization for each model with the resampled SFS
                    for j, (creator_func, model_name) in enumerate(model_creators):
                        logging.info(f"Optimizing model: {model_name}")
                        
                        try:
                            # Create model instance with resampled data
                            model = creator_func(resampled_sfs, randomize=True)
                            model_instances.append(model)
                            
                            start_time = time.time()
                            # Increased maxiter from 100 to 200 for better optimization
                            result = model.optimize(
                                method="L-BFGS-B", 
                                options={
                                    "maxiter": 200,  # Increased from 100
                                    "ftol": 1e-5,    # More stringent tolerance
                                    "gtol": 1e-5     # More stringent tolerance
                                }
                            )
                            end_time = time.time()
                            
                            # Calculate model fit statistics
                            log_likelihood = model.log_likelihood()
                            nparams = len(model.get_params())
                            aic = 2 * nparams - 2 * log_likelihood
                            
                            aics.append(aic)
                            log_likelihoods.append(log_likelihood)
                            
                            logging.info(f"Model {model_name} optimization completed in {end_time - start_time:.2f} seconds")
                            logging.info(f"AIC: {aic}, Log-likelihood: {log_likelihood}")
                        
                        except Exception as e:
                            logging.error(f"Error optimizing model {model_name}: {e}")
                            # If a model fails, assign it a very high AIC
                            aics.append(float('inf'))
                            log_likelihoods.append(float('-inf'))
                            model_instances.append(None)
                    
                    # Calculate delta AICs and AIC weights
                    if not all(np.isinf(aic) for aic in aics):
                        valid_aics = [aic for aic in aics if not np.isinf(aic)]
                        if valid_aics:
                            min_aic = min(valid_aics)
                            delta_aics = [aic - min_aic if not np.isinf(aic) else float('inf') for aic in aics]
                            
                            # Calculate AIC weights, handling infinity
                            exp_delta = [np.exp(-0.5 * d) if not np.isinf(d) else 0.0 for d in delta_aics]
                            sum_exp_delta = sum(exp_delta)
                            
                            if sum_exp_delta > 0:
                                aic_weights = [exp / sum_exp_delta for exp in exp_delta]
                            else:
                                aic_weights = [0.0] * len(aics)
                                
                            # Identify the best model
                            best_model_index = aics.index(min(aics))
                            best_model = model_instances[best_model_index]
                            best_model_name = model_creators[best_model_index][1]
                            best_model_counts[best_model_index] += 1
                            
                            # Get parameters of the best model
                            if best_model is not None:
                                params = best_model.get_params()
                            else:
                                params = {}
                            
                            # Write detailed info for this replicate
                            for j, (creator_func, model_name) in enumerate(model_creators):
                                aic_file.write(f"{i+1},{model_name},{aics[j]},{log_likelihoods[j]},{delta_aics[j]},{aic_weights[j]}\n")
                            
                            # Record the best model for this replicate
                            best_models_file.write(f"{i+1},{best_model_index},{best_model_name}\n")
                            
                            # Save parameters of the best model
                            for param_name, param_value in params.items():
                                params_file.write(f"{i+1},{best_model_name},{param_name},{param_value}\n")
                            
                            successful_replicates += 1
                            logging.info(f"Best model for replicate {i+1}: {best_model_name}")
                        else:
                            logging.warning(f"No valid AICs found for replicate {i+1}, skipping")
                    else:
                        logging.warning(f"All models failed for replicate {i+1}, skipping")
                
                except Exception as e:
                    logging.error(f"Error in bootstrap replicate {i+1}: {e}")
                    continue
                
                # Flush files periodically to ensure data is written
                if (i + 1) % 5 == 0:
                    aic_file.flush()
                    best_models_file.flush()
                    params_file.flush()
            
            # Write summary at the end
            aic_file.write("\n\nSummary:\n")
            aic_file.write(f"Successful replicates: {successful_replicates} out of {num_replicates}\n")
            
            for model_idx, count in best_model_counts.items():
                model_name = model_creators[model_idx][1]
                if successful_replicates > 0:
                    percentage = (count / successful_replicates) * 100
                else:
                    percentage = 0
                aic_file.write(f"Model {model_name}: {count} replicates ({percentage:.2f}%)\n")
            
            logging.info(f"Bootstrap analysis completed. Successful replicates: {successful_replicates}/{num_replicates}")
            
            for model_idx, count in best_model_counts.items():
                model_name = model_creators[model_idx][1]
                if successful_replicates > 0:
                    percentage = (count / successful_replicates) * 100
                else:
                    percentage = 0
                logging.info(f"Model {model_name}: {count} replicates ({percentage:.2f}%)")
    
    except Exception as e:
        logging.error(f"An error occurred in the main execution: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    start_time = time.time()
    logging.info("Script started")
    
    try:
        main()
    except Exception as e:
        logging.error(f"Unhandled exception: {e}")
        sys.exit(1)
    
    end_time = time.time()
    logging.info(f"Script completed in {(end_time - start_time)/60:.2f} minutes")