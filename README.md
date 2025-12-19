A motif detection framework using SeqNMF with Split Bregman Iteration (SBI) for neural activity.

## Setup

**Download data first:** Before running any scripts, create a `Data/` directory and download the required data files to the `Data/` directory.

## Repository Structure

### Core Components

- **FlexMF/** - Main motif detection scripts using SeqNMF with Split Bregman Iteration (SBI)
- **TFOCS/** - Supporting optimization toolbox for FlexMF
- **Utils/** - Additional supporting functions

### Demo Scripts

- **FlexMF/demo_synthetic.m** - Demonstration using synthetic data
- **find_motifs_demo.m** - Demonstration for detecting motifs on real data

### Analysis Scripts

- **analyze_motif_encodings.m** - Analyze motif encoding properties for a given session
- **analyze_single_neuron_encoding.m** - Analyze single neuron encoding properties for a given session
- **analyze_motifs_trials.m** - Visualize detected motifs relative to dF/F traces for a given session

## Quick Start

1. Download and place data files in the `Data/` directory
2. Run `FlexMF/demo_synthetic.m` to test the method on synthetic data
3. Run `find_motifs_demo.m` to test the method on real data
4. Use the analysis scripts to examine motif properties and visualizations

## Requirements

- MATLAB (version 2021 or higher)
- TFOCS (https://github.com/cvxr/TFOCS/tree/master)

## Citation

If you use this code, please cite: 

[1] Xiang, Liang, et al. "Precise Extraction of Neural Motifs Reveals Multiscale, Parallel Encoding Schemes in Auditory Cortex." bioRxiv (2025).

[2] Becker, Stephen R., Emmanuel J. Candès, and Michael C. Grant. "Templates for convex cone problems with applications to sparse signal recovery." Mathematical programming computation 3.3 (2011): 165-218.
