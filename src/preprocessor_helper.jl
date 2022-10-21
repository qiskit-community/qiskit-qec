module PreprocessorHelper

# Make dummy functions
for func in (:Hd3zExperiment, :Hd3zExperimentProcessing,
             :Hd3pExperimentProcessing, :Hd3pExperiment,
             :Hd3zCircExperiment, :Hd3pCircExperiment)
    @eval function $func end
    @eval export $func
end

end # module PreprocessorHelper
