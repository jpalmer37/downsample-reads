# downsample-reads

## Usage

```
nextflow run BCCDC-PHL/downsample-reads \
  -profile conda \
  --cache ~/.conda/envs \
  --coverage 30 \
  --fastq_input </path/to/fastqs> \
  --outdir </path/to/output_dir>
```

```
nextflow run BCCDC-PHL/downsample-reads \
  -profile conda \
  --cache ~/.conda/envs \
  --fastq_input </path/to/fastqs> \
  --outdir </path/to/output_dir>
```

The default genome size is 5 megabases (`5m`). To specify another genome size, use the `--genome_size` flag:

```
nextflow run BCCDC-PHL/downsample-reads \
  -profile conda \
  --cache ~/.conda/envs \
  --genome_size '30k' \
  --coverage 100 \
  --fastq_input </path/to/fastqs> \
  --outdir </path/to/output_dir>
```

### SampleSheet Input

If sample-specific downsampling parameters are needed, they can be provided via samplesheet input.

Prepare a `samplesheet.csv` file with the following fields:

```
ID
R1
R2
GENOME_SIZE
COVERAGE
```

...for example:

```csv
ID,R1,R2,GENOME_SIZE,COVERAGE
sample-01,/path/to/sample-01_R1.fastq.gz,/path/to/sample-01_R2.fastq.gz,5.0m,100
sample-02,/path/to/sample-02_R1.fastq.gz,/path/to/sample-02_R2.fastq.gz,3.0m,100
sample-03,/path/to/sample-03_R1.fastq.gz,/path/to/sample-03_R2.fastq.gz,3.0m,50
```

...then run the pipeline using the `--samplesheet_input` flag as follows:

```
nextflow run BCCDC-PHL/downsample-reads \
  -profile conda \
  --cache ~/.conda/envs \
  --samplesheet_input samplesheet.csv \
  --fastq_input </path/to/fastqs> \
  --outdir </path/to/output_dir>
```

Note that you can include multiple entries for each sample:

```csv
ID,R1,R2,GENOME_SIZE,COVERAGE
sample-01,/path/to/sample-01_R1.fastq.gz,/path/to/sample-01_R2.fastq.gz,5.0m,25
sample-01,/path/to/sample-01_R1.fastq.gz,/path/to/sample-01_R2.fastq.gz,5.0m,50
sample-01,/path/to/sample-01_R1.fastq.gz,/path/to/sample-01_R2.fastq.gz,5.0m,100
sample-02,/path/to/sample-02_R1.fastq.gz,/path/to/sample-02_R2.fastq.gz,3.0m,10
sample-02,/path/to/sample-02_R1.fastq.gz,/path/to/sample-02_R2.fastq.gz,3.0m,100
sample-03,/path/to/sample-03_R1.fastq.gz,/path/to/sample-03_R2.fastq.gz,3.0m,50
sample-03,/path/to/sample-03_R1.fastq.gz,/path/to/sample-03_R2.fastq.gz,3.0m,100
sample-03,/path/to/sample-03_R1.fastq.gz,/path/to/sample-03_R2.fastq.gz,3.0m,200
```


## Output

A pair of fastq.gz files will be produced for each target coverage, for each sample.
Filenames are appended with `-downsample-Nx`, where `N` is the target coverage for that file pair.

```
outdir
|-- sample-01-downsample-10x_R1.fastq.gz
|-- sample-01-downsample-10x_R2.fastq.gz
|-- sample-01-downsample-20x_R1.fastq.gz
|-- sample-01-downsample-20x_R2.fastq.gz
|-- sample-01-downsample-2x_R1.fastq.gz
|-- sample-01-downsample-2x_R2.fastq.gz
|-- sample-01-downsample-30x_R1.fastq.gz
|-- sample-01-downsample-30x_R2.fastq.gz
|-- sample-01-downsample-40x_R1.fastq.gz
|-- sample-01-downsample-40x_R2.fastq.gz
|-- sample-01-downsample-5x_R1.fastq.gz
|-- sample-01-downsample-5x_R2.fastq.gz
`-- fastp.csv
```

The `fastp.csv` file will include a summary of the number of reads and bases, for each sample for each target coverage, plus the `original` input files.
The depth of coverage of the output files is estimated based on the `total_bases` divided by the `genome_size`.

```csv
sample_id  total_reads  total_bases  q30_rate  genome_size  target_coverage  estimated_coverage
sample-01  2549698      383269283    0.884585  5.0m         original         76.654
sample-02  2500548      375831165    0.877859  5.0m         original         75.166
sample-03  3432324      515552128    0.887493  5.5m         original         93.737
sample-01  166352       25000128     0.884852  5.0m         5                5.0
sample-02  332658       50000175     0.877714  5.0m         10               10.0
sample-03  366112       55000224     0.887422  5.5m         10               10.0
sample-01  1663158      250000208    0.884677  5.0m         50               50.0
sample-02  1830796      275000177    0.887517  5.5m         50               50.0
sample-03  2500548      375831165    0.877859  5.0m         100              75.166
```