from cosap import (
    MDUP,
    BamReader,
    Mapper,
    Pipeline,
    PipelineRunner,
    Recalibrator,
    Trimmer,
    VariantCaller,
)

tumor_bam = BamReader("/workdir/unprocessed_tumor_bowtie.bam")
normal_bam = BamReader("/workdir/unprocessed_normal_bowtie.bam")

# Creating preprocessors
mdup_normal = MDUP(input_step=normal_bam)
mdup_tumor = MDUP(input_step=tumor_bam)

# Creating variant callers
caller = VariantCaller(library="strelka", germline=mdup_normal, tumor=mdup_tumor)

# Creating pipeline and adding steps to it
pipeline = (
    Pipeline()
    .add(mdup_normal)
    .add(mdup_tumor)
    .add(caller)
)

# Creating the config that contains all information about the pipeline
config = pipeline.build()

# Create a pipeline runner and run the config file
pipeline_runner = PipelineRunner()
pipeline_runner.run_pipeline(config)
