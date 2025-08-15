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

tumor_bam = BamReader("/home/rengin/cosap_data/unprocessed_tumor_bwa.bam")
normal_bam = BamReader("/home/rengin/cosap_data/unprocessed_normal_bwa.bam")

# Creating preprocessors
mdup_normal = MDUP(input_step=normal_bam)
mdup_tumor = MDUP(input_step=tumor_bam)

basecal_normal = Recalibrator(input_step=mdup_normal)
basecal_tumor = Recalibrator(input_step=mdup_tumor)

# Creating variant callers
caller = VariantCaller(library="somaticsniper", normal_sample=basecal_normal, tumor_sample=basecal_tumor)

# Creating pipeline and adding steps to it
pipeline = (
    Pipeline()
    .add(mdup_normal)
    .add(mdup_tumor)
    .add(basecal_normal)
    .add(basecal_tumor)
    .add(caller)
)

# Creating the config that contains all information about the pipeline
config = pipeline.build()

# Create a pipeline runner and run the config file
pipeline_runner = PipelineRunner()
pipeline_runner.run_pipeline(config)
