# HeavyHex

This is a Julia package is for running the MLD algorithm in [`Matching and maximum likelihood decoding of a multi-round subsystem quantum error correction experiment`](https://arxiv.org/pdf/2203.07205.pdf) on the MLD experiments data set `HHCData` 

Please note, the `results.csv` file contains the results for all 44 jobs, including each job's runtime on a large, shared computing cluster, with up to 80 threads.

## Installation / Usage

### Some dependencies in a registry

The packages [`DictTools.jl`](https://github.com/jlapeyre/DictTools.jl) and [`BitsX.jl`](https://github.com/jlapeyre/BitsX.jl) have been submitted
for inclusion in the the Julia [General Registry](https://github.com/JuliaRegistries/General) soon. They will probably be available there around May 1, 2023.
In the meantime, they are available in a small registry which can be accessed easily like this:
```julia
> julia # Hit the ']' key to go to pkg mode
(@v1.9) pkg> registry add https://github.com/jlapeyre/LapeyreRegistry
```

### Installing HeavyHex.jl
It might be best to use this package in a development mode. So do the following.

* Clone this repo and cd into the top level
* Start Julia in terminal and do the following

```julia
julia>
 # Now hit the ']' key to go to pkg mode

(@v1.7) pkg> activate .
 # Now hit backspace to return to Julia mode

julia> using HeavyHex
```

### Using HeavyHex.jl 

# Notes 
*Remember, when running Julia, to start Julia with `julia -t auto` to allow for multithreading*

This package can run any job from [HHCData](https://github.com/qiskit-community/HHCData). *As a prerequisite for running any job, you must first download a local copy of the HHCData set by doing a `git clone` for this repo:  [HHCData](https://github.com/qiskit-community/HHCData).*  

*Some jobs can get quite large so we recommend using jobs with a small round number, such as job 656. All jobs and their round numbers can be found in the `results.csv` file.*

# Code Example

```
import HeavyHex:JobRunner
import HeavyHex:DecoderIO

### Basic Job Example
# We begin by giving `HeavyHex.jl` the path to our  local HHCData repo. This sets the path globally.
path_to_hhcdata = "<local-path>/HHCData" 
DecoderIO.setup_hhc_directory_path(path_to_hhcdata)

# We call `JobRunner.run_job` and pass in the `job_id`, which in our example is 656. This runs our job.
println(JobRunner.run_job(656)) # don't use logging


### Setting Up Logging Example

import HeavyHex:LoggingUtils # setup logging, globally

log_folder_path = "<local-path-to-where-you-want-logs-stored>"
LoggingUtils.setup_logging(true, false, false, log_folder_path) # We only have to specify the logging information once (unless we want to change it).
# We toggle `write_res` to `true`, `write_prefix_store` to `false`, `write_log` to `false`. 
# write_res tells HeavyHex to write a `results.txt`  file containing the results of each job's run. 
# write_prefix_store tells the logger, for each job, to write the prefix_store to a seperate file after the decoder finishes running
# write_log tells HeavyHex to create a seperate log file and write updates as the decoder runs

# Now, because we have logs, when we run `JobRunner.run_job`, we must also include our own custom name for this run, `test_job`. This custom name is used create a unique sub-directory in which to store the logs for this run. 
println(JobRunner.run_job(656, "test_job")) # run job 

### Multiple Jobs Example

# Running multiple jobs but storing the logging/results in the same subdirectory "my_jobs"
for job_id in [656, 660]
JobRunner.run_job(job_id, "my_jobs")
end 

```
# HeavyHex.jl
