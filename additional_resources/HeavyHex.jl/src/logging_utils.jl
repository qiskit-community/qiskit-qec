"""
    module LoggingUtils

    Tools for more cleanly logging the decoder. 
"""
module LoggingUtils
using HeavyHex.Bits
export convert_prefix_store_to_cross_platform_dict,
    convert_prefix_store_to_non_concurrent_dict, setup_info_writer
using Logging
using Dates
using Serialization


global WRITE_RES = false
global WRITE_PREFIX_STORE = false
global WRITE_DYNAMIC_LOG = false
global WRITTEN_DIR = ""

# Set up what logs get written 
""" 
    setup_logging(write_res=false, write_prefix_store=false, write_log=false, written_dir=")
    `write_res` if true will write the results of each experiment to a  <written_dir>/<job_id>/mld_results/results.txt file 
    `write_prefix_store` if true will write the prefix_store (in 3 formats) to   <written_dir>/<job_id>/mld_dicts/  
        This is used *primarily* for debugging and is unlikely to be useful to a user.
        The 3 formats are concurrent_dictionary, the Julia native dictionary, and the Julia String dictionary. 
        Because `prefix_store` has custom bit packaging, not all dictionary types are compatable across platforms,
        hence the 3 types 
    `write_log` if true will write logs to <written_dir>/<job_id>/mld_logs/
        while the job is being decoded. Periodically checking the log outputs is
         an easy way to make sure a job is still running. However, they are unlikely to be useful for small jobs.
    `written_dir` is the full path to the directory that will store any results, prefix_stores, or logs. 

"""
function setup_logging(write_res, write_prefix_store, write_log, written_dir)
    global WRITE_RES = write_res # after job finishes, write results to a separate file 
    global WRITE_PREFIX_STORE = write_prefix_store # after job finishes, store the prefix_store dictionary
    global WRITE_DYNAMIC_LOG = write_log # write logs as job is running 
    global WRITTEN_DIR = written_dir # directory to contain logs, results, and prefix_store
end


####
#### INFO WRITER
####

function setup_info_writer(job_id, unique_info)

    ### CREATE FOLDER PATHS (won't be used if all logging vars are set to false)
    unique_info_folder_path = joinpath(WRITTEN_DIR, "$(unique_info)")
    folder_path = joinpath(unique_info_folder_path, "$(job_id)")

    if (WRITE_RES || WRITE_DYNAMIC_LOG || WRITE_PREFIX_STORE)
        if !(isdir(unique_info_folder_path))
            mkdir(unique_info_folder_path)
        end
        mkdir(folder_path) # make a few folder with unique_info where everything will be stored
    end

    if WRITE_RES
        results_folder_path = joinpath(folder_path, "mld_results")
        results_file_path = joinpath(results_folder_path, "$(job_id)_results.txt")
        mkdir(results_folder_path)
        write_results = write_results_wrapper(results_file_path)
    else
        write_results = write_results_wrapper(nothing)
    end

    if WRITE_DYNAMIC_LOG
        log_folder_path = joinpath(folder_path, "mld_logs")
        log_file_path = joinpath(log_folder_path, "$(unique_info)_$(job_id).log")
        mkdir(log_folder_path)
        logio, logger = setup_logio(log_file_path)
        write_log = write_log_wrapper(logger)
        flush_logs = flush_log_wrapper(logio)
        close_logs = close_wrapper(logio)
    else
        write_log = write_log_wrapper(nothing)
        flush_logs = flush_log_wrapper(nothing)
        close_logs = close_wrapper(nothing)
    end

    if WRITE_PREFIX_STORE
        dicts_folder_path = joinpath(folder_path, "mld_dicts")
        con_file_path =
            joinpath(dicts_folder_path, "$(unique_info)_$(job_id)_concurrent_dictionary")
        julia_dict_file_path =
            joinpath(dicts_folder_path, "$(unique_info)_$(job_id)_julia_dictionary")
        julia_string_dict_filepath =
            joinpath(dicts_folder_path, "$(unique_info)_$(job_id)_julia_string_dictionary")
        mkdir(dicts_folder_path)
        write_prefix_store = write_prefix_store_wrapper(
            con_file_path,
            julia_dict_file_path,
            julia_string_dict_filepath,
        )
    else
        write_prefix_store = write_prefix_store_wrapper(nothing, nothing, nothing)
    end

    info_writer =
        InfoWriter(write_prefix_store, write_results, write_log, flush_logs, close_logs)
    return info_writer
end

struct InfoWriter
    write_prefix_store::Base.Callable
    write_results::Base.Callable
    write_logs::Base.Callable
    flush_logs::Base.Callable
    close::Base.Callable
end

#### WRAPPERS are needed because Julia won't let you define the same local function in different blocks -> https://github.com/JuliaLang/julia/issues/15602

function write_results_wrapper(res_file_store_path)
    if !(res_file_store_path === nothing)
        function write_results(results, runtime)
            open(res_file_store_path, "w") do io
                write(io, "$(results)")
                write(io, "Runtime: $(runtime)")
            end
        end
        return write_results
    else
        function write_no_results(results, runtime) end
        return write_no_results
    end
end

function write_prefix_store_wrapper(
    con_file_path,
    julia_dict_file_path,
    julia_string_dict_filepath,
)
    if !(con_file_path === nothing)
        function save_concurrent_dictionary(con_dictionary)

            Serialization.serialize(con_file_path, con_dictionary)

            non_con_dict =
                LoggingUtils.convert_prefix_store_to_non_concurrent_dict(con_dictionary)
            Serialization.serialize(julia_dict_file_path, non_con_dict)

            non_con_str_dict =
                LoggingUtils.convert_prefix_store_to_cross_platform_dict(con_dictionary)
            Serialization.serialize(julia_string_dict_filepath, non_con_str_dict)
        end
        return save_concurrent_dictionary
    else
        function save_nothing(con_dictionary) end
        return save_nothing
    end
end

function write_log_wrapper(logger)
    if !(logger === nothing)
        function write_log(data)
            with_logger(logger) do
                @info(data)
            end
        end
        return write_log
    else
        function write_no_log(data) end
        return write_no_log
    end
end

function flush_log_wrapper(logio)
    if !(logio === nothing)
        function flush_logs()
            flush(logio)
        end
        return flush_logs
    else
        function flush_no_logs() end
        return flush_no_logs
    end
end

function close_wrapper(logio)
    if !(logio === nothing)
        function close_log()
            close(logio)
        end
        return close_log
    else
        function do_nothing() end
        return do_nothing
    end
end

### LOGGING 
function setup_logio(log_file)
    logio = open(log_file, "w+")
    logger = SimpleLogger(logio, Logging.Debug)
    return logio, logger
end


####
#### Dictionary Conversions 
####

function convert_prefix_store_to_cross_platform_dict(con_dictionary)
    storage_dict = Dict{String,Any}()
    for key in keys(con_dictionary)
        string_key = Bits.bitvector1mask_2_string(key)
        bv1mtup, bv1mdict = con_dictionary[key]
        edges_tuple = Tuple(Bits.bitvector1mask_2_string(bv1m) for bv1m in bv1mtup)
        string_dict = Dict{}(
            Bits.bitvector1mask_2_string(key) => val for (key, val) in pairs(bv1mdict)
        )
        value_tuple = (edges_tuple, string_dict)
        storage_dict[string_key] = value_tuple
    end
    return storage_dict
end


function convert_prefix_store_to_non_concurrent_dict(con_dictionary)
    storage_dict = Dict{BitVector1Mask{UInt128},Any}()
    for key in keys(con_dictionary)
        storage_dict[key] = con_dictionary[key]
    end
    return storage_dict
end
end
