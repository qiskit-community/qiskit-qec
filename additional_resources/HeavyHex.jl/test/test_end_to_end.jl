import HeavyHex: JobRunner
import HeavyHex: Experiment
using Test

@testset "HeavyHex.jl" begin
    hhc_data_folder_path = joinpath("resources", "HHCData")
    DecoderIO.setup_hhc_directory_path(hhc_data_folder_path)
    job_656_correct_res = (
        weight_good = 356567,
        count_good = 4,
        num_obs = 364203,
        num_obs_strings = 8,
        error_ratio = 0.02096632921749686,
    )
    job_656_test_res = JobRunner.run_job(656, "")
    @test job_656_correct_res == job_656_test_res
end
