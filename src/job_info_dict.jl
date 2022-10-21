## THIS FILE IS OBSOLETE !

function bit_order_func end
function Hd3zExperiment end
function Hd3zExperimentProcessing end

# Copied from python output of gjl_run_experiment.py
# Setting up the job creates this dict

const _test_job_info = Dict(
    :name => "post_selection_job661",
    :filename => "hyperedges_PSall_and_last_init_both_fav2_false_3_3_3_counts_661.txt",
    :round_no => 2,
    :experiment => Hd3zExperiment,
    :experiment_processing => Hd3zExperimentProcessing,
    :deflag => true,
    :hyperedge_params => [(:p1, 0.00028),
                          (:pw, 0.001),
                          (:p2, 0.01),
                          (:pm, 0.0005),
                          (:pr, 1e-05),
                          (:pi, 0)],
    :bit_order_func => bit_order_func,
    :pre_processed_file => "../examples/test_data/220107_Skywalker_data/post_select_data/ps_counts_GMM_v2/hyperedges_PSall_and_last_init_both_fav2_false_3_3_3_counts_661_preprocessed_online.txt",
    :prefix_len => 3,
    :abs_path => "../examples/test_data/220107_Skywalker_data/post_select_data/ps_counts_GMM_v2/hyperedges_PSall_and_last_init_both_fav2_false_3_3_3_counts_661_preprocessed_online.txt"
)

const test_job_info = (; _test_job_info...)

# The keys won't be ordered the same every time
"""
    JobInfo

`NamedTuple` with info about the job. You can get a pretty display
with `pairs(jobinfo::JobInfo)`.
"""
const JobInfo = NamedTuple{keys(test_job_info), Ttup} where Ttup

#     ks = (keys(_test_job_info)...,)
#     vs = (values(_test_job_info)...,)
#     nt = NamedTuple{ks, typeof(vs)}(vs)
#     JobInfo(nt)
# end

# const test_job_info = let
#     ks = (keys(_test_job_info)...,)
#     vs = (values(_test_job_info)...,)
#     nt = NamedTuple{ks, typeof(vs)}(vs)
#     JobInfo(nt)
# end

# Some original python entries
#    :filename =>  "PSall_and_last_init_both_fav2_false_3_3_3_counts_661.pkl",
#<function yoderdecoder.preprocessor_helper.Hd3zExperiment(rounds, deflag=False, gateSpecific=False, fineOne=False, fineTwo=False, logicalSupport=(-9, -6, -3)) -> Dict[Tuple, sympy.core.add.Add]>,
# function yoderdecoder.experiment.generate_job_decoding_info.<locals>.<lambda>(x)>,
#  :experiment_processing" =>  <function yoderdecoder.preprocessor_helper.Hd3zExperimentProcessing(bstring, deflag=False, flatten=False, isString=True, logicalSupport=(-9, -6, -3))>,
#  :abs_path" =>  "../tests/test_data/220107_Skywalker_data/post_select_data/ps_counts_GMM_v2/PSall_and_last_init_both_fav2_false_3_3_3_counts_661.pkl"
