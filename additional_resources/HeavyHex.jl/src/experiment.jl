"""
module Experiment

This module if for getting metadata/circuit information about each job that
    is  needed by the decoder in order to decode that job correctly.
"""
module Experiment
using ..ExperimentHelper

const POST_SELECTION_HYPEREDGE_PARAMS =
    [(:p1, 0.00028), (:pw, 0.001), (:p2, 0.01), (:pm, 0.0005), (:pr, 0.00001), (:pi, 0.0)]

const JOB_INFO_DICT = Dict(
    656 => ("Z(|0>)", 0),
    661 => ("Z(|0>)", 1),
    665 => ("Z(|0>)", 2),
    669 => ("Z(|0>)", 3),
    673 => ("Z(|0>)", 4),
    677 => ("Z(|0>)", 5),
    681 => ("Z(|0>)", 6),
    685 => ("Z(|0>)", 7),
    689 => ("Z(|0>)", 8),
    693 => ("Z(|0>)", 9),
    697 => ("Z(|0>)", 10),
    657 => ("Z(|1>)", 0),
    662 => ("Z(|1>)", 1),
    666 => ("Z(|1>)", 2),
    670 => ("Z(|1>)", 3),
    674 => ("Z(|1>)", 4),
    678 => ("Z(|1>)", 5),
    682 => ("Z(|1>)", 6),
    686 => ("Z(|1>)", 7),
    690 => ("Z(|1>)", 8),
    694 => ("Z(|1>)", 9),
    698 => ("Z(|1>)", 10),
    659 => ("X(|+>)", 0),
    663 => ("X(|+>)", 1),
    667 => ("X(|+>)", 2),
    671 => ("X(|+>)", 3),
    675 => ("X(|+>)", 4),
    679 => ("X(|+>)", 5),
    683 => ("X(|+>)", 6),
    687 => ("X(|+>)", 7),
    691 => ("X(|+>)", 8),
    695 => ("X(|+>)", 9),
    699 => ("X(|+>)", 10),
    660 => ("X(|->)", 0),
    664 => ("X(|->)", 1),
    668 => ("X(|->)", 2),
    672 => ("X(|->)", 3),
    676 => ("X(|->)", 4),
    680 => ("X(|->)", 5),
    684 => ("X(|->)", 6),
    688 => ("X(|->)", 7),
    692 => ("X(|->)", 8),
    696 => ("X(|->)", 9),
    700 => ("X(|->)", 10),
)

function generate_job_decoding_info(job_id::Integer)
    job_info = JOB_INFO_DICT[job_id]
    str_experiment_type = first(job_info)
    round_no::Int = job_info[2] + 1  # Ted

    if any(in(str_experiment_type), ('0', '1'))
        bof = bit_order_function(round_no) #  @eval (x -> bit_order_function($round_no, x))
        prefix_len = round_no < 2 ? 10 : 4 * (round_no - 1) + 2 * round_no - 5
    end

    if any(in(str_experiment_type), ('+', '-'))
        bof = bit_order_function_x(round_no)
        prefix_len = round_no < 2 ? 10 : 4 * round_no + 2 * (round_no - 1) - 5
    end

    # Note that the preprocessed data is in a spc-separated text file.
    decoding_info_tuple = (
        name = "post_selection_job$(job_id)",
        round_no = round_no,
        deflag = true,
        hyperedge_params = POST_SELECTION_HYPEREDGE_PARAMS,
        bit_order_func = bof,
        prefix_len = prefix_len,
    )
    return decoding_info_tuple
end


const Z0 = [656, 661, 665, 669, 673, 677, 681, 685, 689, 693, 697]

const Z1 = [657, 662, 666, 670, 674, 678, 682, 686, 690, 694, 698]

const XPLUS = [659, 663, 667, 671, 675, 679, 683, 687, 691, 695, 699]

const XMINUS = [660, 664, 668, 672, 676, 680, 684, 688, 692, 696, 700]

const ALL_JOBS = [Z0, Z1, XPLUS, XMINUS]

end # module Experiment
