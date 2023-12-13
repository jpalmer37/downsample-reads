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

Alternatively, samples can be downsampled to multiple depths using a `coverages.csv` file:

```
nextflow run BCCDC-PHL/downsample-reads \
  -profile conda \
  --cache ~/.conda/envs \
  --coverages </path/to/coverages.csv> \
  --fastq_input </path/to/fastqs> \
  --outdir </path/to/output_dir>
```

The `coverages.csv` file should be a single-column list of target depths of coverage, with no header. For example:

```
30
20
10
5
```

The default genome size is 5 megabases (`5m`). To specify another genome size, use the `--genome_size` flag:

```
nextflow run BCCDC-PHL/downsample-reads \
  -profile conda \
  --cache ~/.conda/envs \
  --genome_size '30k' \
  --coverages </path/to/coverages.csv> \
  --fastq_input </path/to/fastqs> \
  --outdir </path/to/output_dir>
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

```
sample_id  total_reads  total_bases  q30_rate  target_coverage
sample-01  42158        10000353     0.881622  2
sample-01  105510       25000351     0.883446  5
sample-01  210932       50000154     0.883396  10
sample-01  421748       100000287    0.882951  20
sample-01  632628       150000289    0.88305   30
sample-01  843538       200000061    0.883072  40
sample-01  1172548      278027194    0.882882  original
```