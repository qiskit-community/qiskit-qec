using HeavyHex: Experiment, JobRunner, DecoderIO, OmlDecoding
using HeavyHex.Bits: bits
using Dictionaries: Dictionary

@testset "Experiment" begin
    hhc_data_folder_path = joinpath("resources", "HHCData")
    DecoderIO.setup_hhc_directory_path(hhc_data_folder_path)

    job_id = 660
    job_info = Experiment.generate_job_decoding_info(job_id)
    @test job_info isa NamedTuple
    @test fieldnames(typeof(job_info)) ==
          (:name, :round_no, :deflag, :hyperedge_params, :bit_order_func, :prefix_len)
    @test job_info.name == "post_selection_job660"
    @test job_info.round_no == 1
    @test job_info.prefix_len == 10

    preproc_data = DecoderIO.read_preprocessed(job_id)
    @test preproc_data isa NamedTuple
    @test fieldnames(typeof(preproc_data)) ==
          (:preprocessed_info, :numerical_hyperedges_dict)
    @test length(preproc_data.preprocessed_info) == 32
    @test length(preproc_data.numerical_hyperedges_dict) == 8
    ppinfo1 = first(preproc_data.preprocessed_info)
    hyperedge1 = first(preproc_data.numerical_hyperedges_dict)
    @test ppinfo1 == (bits = bits("00000"), counts = 1224)
    @test hyperedge1 == (bits = bits("00000"), prob = 0.19625333333333334)
    @test ppinfo1.bits == HeavyHex.Bits.bits("00000")
    @test hyperedge1.bits == HeavyHex.Bits.bits("00000")

    _hyperedges = preproc_data.numerical_hyperedges_dict
    hyperedges = Dictionary(_hyperedges.bits, _hyperedges.prob)

    preproc_info = preproc_data.preprocessed_info
    observed_counts = preproc_info.counts
    observations_all = reverse.(preproc_info.bits)
    function bit_from_start(observations_all)
        observations = [x[2:end] for x in observations_all]
        expected_bits = [x[1:1] for x in observations_all]
        return (observations, expected_bits)
    end
    (observations, expected_bits) = bit_from_start(observations_all)
    obs = observations[1]
    decoder = OmlDecoding.OnlineMLDecoder(job_info; hyperedges)
    @test decoder isa OmlDecoding.OnlineMLDecoder
    result = OmlDecoding.decode(decoder, obs)
    @test result isa Tuple
    @test length(result) == 2
    @test first(result) == bits("0")
end
