# importing more than we need here

using HeavyHex.Bits
import HeavyHex: Experiment
import HeavyHex: Decoders
import HeavyHex: DecoderIO
using Dictionaries: Dictionary

using Dates
using Logging



function run_job(job_id, logio)

    println("running... $(job_id)")
    @info("running... $job_id")
    job_info = Experiment.generate_job_decoding_info(job_id)
    preprocessed_data = DecoderIO.read_preprocessed(job_id)

    _hyperedges = preprocessed_data.numerical_hyper_edges_dict
    hyperedges = Dictionary(_hyperedges.bits, _hyperedges.prob)

    prep_info = preprocessed_data.preprocessed_info
    observed_counts = prep_info.counts


    observations_all = reverse.(prep_info.bits)

    function bit_from_start(observations_all)
        observations = [x[2:end] for x in observations_all]
        expected_bits = [x[1:1] for x in observations_all]
        return (observations, expected_bits)
    end

    (observations, expected_bits) = bit_from_start(observations_all)

    println("no obser: $(length(observations))")

    function run_observations(decoder, observations, expected_bits, observed_counts, logio)
        cur_j = 1
        weight_good = 0
        count_good = 0
        for (expected, observation, bcount) in zip(expected_bits, observations, observed_counts)
            cur_j += 1
            result = Decoders.decode(decoder, observation)
            predicted = first(result)
            if expected == predicted
                weight_good += bcount
                count_good += 1
            end
            println("yellow")
            if cur_j % 1 == 0
                @info("currently running... $cur_j at $(Dates.format(now(), "DD:HH:MM:SS")) aka $(DateTime(now())) \n\n current:  bcount: $bcount ;; weight_good: $weight_good ;; count_good: $count_good \n\ncurrently:  predict: $predicted ;; expected: $expected")
                flush(logio)
                println("currently running... $cur_j at $(Dates.format(now(), "DD:HH:MM:SS")) aka $(DateTime(now()))")
            end
        end
        num_obs = sum(observed_counts)
        num_obs_strings = length(observations)
        error_ratio = (1 - weight_good / num_obs)
        return (; weight_good, count_good, num_obs, num_obs_strings, error_ratio)
    end

    decoder = Decoders.OnlineMLDecoder(job_info; hyperedges)

    # To run the decoding
    results = run_observations(decoder, observations, expected_bits, observed_counts, logio)
    println(results)
    @info("$job_id results: $(results)")
    open(joinpath("exp_results", "recovery_res", "$(job_id)_results.txt"), "w") do io
        write(io, "$(results)")
    end
end

function run_basis(basis, log_file)
    logio = open(log_file, "w+")
    logger = SimpleLogger(logio, Logging.Debug)
    global_logger(logger)
    for job in basis
        run_job(job, logio)
    end
end


Z0 = [656,
    661,
    665,
    669,
    673,
    677,
    681,
    685,
    689,
    693,
    697]

Z1 = [
    657,
    662,
    666,
    670,
    674,
    678,
    682,
    686,
    690,
    694,
    698
]

Xplus = [
    659,
    663,
    667,
    671,
    675,
    679,
    683,
    687,
    691,
    695,
    699
]

Xminus = [
    660,
    664,
    668,
    672,
    676,
    680,
    684,
    688,
    692,
    696,
    700
]

missing_rounds = [
    # 700,
    # # 699,
    656, # test 
    695, #x+9
    696, #x-9
    698, # z 1 10 
    697, # z 0 10
]


log_file_path = joinpath("exp_results", "recovery_res")
