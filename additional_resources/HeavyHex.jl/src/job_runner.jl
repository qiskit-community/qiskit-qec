"""
    module JobRunner

    This is the only module the user needs to interact with in order to run the decoder.
"""
module JobRunner
using HeavyHex.Bits
import HeavyHex: LoggingUtils
import HeavyHex: Experiment
import HeavyHex: OmlDecoding
import HeavyHex: DecoderIO
using Dates



function load_job(job_id)
    job_info = Experiment.generate_job_decoding_info(job_id)
    hyperedges, observed_counts, observations, expected_bits =
        DecoderIO.load_job_run_from_id(job_id)
    return job_info, hyperedges, observed_counts, observations, expected_bits
end


"""
    run_job(job_id, unique_info)
    `job_id` is the JobID labelling of the jobs in HHCData, clearly laid out in `results.csv`
    `unique_info` is a custom name used create a unique sub-directory in which to store the logs/results/prefix_store for this run. 
"""

function run_job(job_id, unique_info = "")
    function run_observations(
        decoder::OmlDecoding.OnlineMLDecoder,
        observations,
        expected_bits,
        observed_counts,
    )

        safe_lock = ReentrantLock() # writing to files should probably be threadsafe
        weight_good = 0
        count_good = 0

        cur_j = 1
        lbits = length(observations)
        num_obs = sum(observed_counts)
        Threads.@threads for i = 1:lbits # using length over eachindex is discourage but I'm not sure how eachindex works w threads

            expected = expected_bits[i]
            observation = observations[i]
            bcount = observed_counts[i]
            result = OmlDecoding.decode(decoder, observation) # since all attempts at updating key k will be with the same value v, maybe locking isn't necessary
            predicted = first(result)

            lock(safe_lock)
            cur_j += 1

            if expected == predicted
                weight_good += bcount
                count_good += 1
            end
            if cur_j % 100 == 0 || cur_j == lbits
                # there's definitely a better way to do this https://instance-factory.com/?p=181
                info_writer.write_logs(
                    "currently running... $cur_j at $(Dates.format(Dates.now(), "DD:HH:MM:SS")) aka $(DateTime(Dates.now())) \n\n current:  bcount: $bcount ;; weight_good: $weight_good ;; count_good: $count_good \n\ncurrently:  predict: $predicted ;; expected: $expected) current thread is: $(Threads.threadid()) ",
                )
                info_writer.flush_logs()
            end

            unlock(safe_lock)

        end
        info_writer.write_logs(
            "final cur_j: $(cur_j);  current thread is: $(Threads.threadid()) ",
        )
        info_writer.flush_logs()
        num_obs_strings = length(observations)
        error_ratio = (1 - weight_good / num_obs)

        return (; weight_good, count_good, num_obs, num_obs_strings, error_ratio)
    end

    info_writer = LoggingUtils.setup_info_writer(job_id, unique_info)
    info_writer.write_logs("running... $job_id")
    start_time = Dates.now()

    job_info, hyperedges, observed_counts, observations, expected_bits = load_job(job_id)
    decoder = OmlDecoding.OnlineMLDecoder(job_info; hyperedges)
    results = run_observations(decoder, observations, expected_bits, observed_counts)
    finish_time = Dates.now()
    runtime = finish_time - start_time

    info_writer.write_logs("Finished with runtime: $(runtime)")
    info_writer.write_logs("$job_id results: $(results)")
    info_writer.flush_logs()
    info_writer.write_results(results, runtime)
    info_writer.write_prefix_store(decoder.prefix_store)
    info_writer.close()

    return results
end
end
