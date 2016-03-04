## Running the compiled software ##
There are 2 ways to run the compiled software: 1) executable or 2) shell script. These options are offered by the Matlab as a result of using the MCR compiler. If you have MCR already installed and added to your path (specifically the LD_LIBRARY_PATH environment variable) then you can use the executable; otherwise, use the shell script as it allows you to manually specify the MCR install path.  

```shell
### In Linux ###
cd <$install_dir>/APOLLOH_0.1.1/bin/

### Option 1) Run APOLLOH: Executable ###
export LD_LIBRARY_PATH=<$install_dir>/MATLAB_Component_Runtime/v77/runtime/glnxa64/:<$install_dir>/MATLAB_Component_Runtime/v77/sys/os/glnxa64/:<$install_dir>/MATLAB_Component_Runtime/v77/bin/glnxa64/
./apolloh <$input_allelic_ratio_datafile> <$input_copy_number_datafile> ../parameters/K18/apolloh_K18_params_Illumina_stromalRatio_Hyper10k.mat <$output_converged_parameter_file> <$output_results_file>


### Option 2) Run APOLLOH: Shell Script ###
./run_apolloh.sh <$install_dir>/MATLAB_Component_Runtime/v77/ <$input_allelic_ratio_datafile> <$input_copy_number_datafile> ../parameters/K18/apolloh_K18_params_Illumina_stromalRatio_Hyper10k.mat <$output_converged_parameter_file> <$output_results_file>


### Sort data by chromosome and position ###
sort -k1,1n -k2,2n <$output_results_file> tmp.txt
mv tmp.txt <$output_results_file>

### Generate segment file ###
../scripts/createSingleSegFileFromAPOLLOH.pl -i <$output_results_file> -o <$output_segment_results_file> -calls 1;
```

## Example: Running the compiled software on test data ##
We have provided some data and a script to test that APOLLOH is properly installed and able to run. These files can also serve as a guide to preparing and formatting your input data files.

```shell
> cd test/
> head -5 cnfile.txt
ID chrom start end num.mark seg.mean state.num state.name
test 1 1 340000 340 0.128103 3 NEUT
test 1 340001 354000 14 0.672409 5 AMP
test 1 354001 1509000 1155 0.135829 3 NEUT
test 1 1509001 1542000 33 0.833131 6 HLAMP
test 1 1542001 5288000 3746 0.137101 3 NEUT
> head -5 infile.txt
1 1 N 27 N 3
1 698 N 16 N 3
1 4803 N 2 N 9
1 7652 N 37 N 6
1 7997 N 1 N 8
> ./run_APOLLOH_test.sh
```

## Running the software in Matlab ##
If you have Matlab installed and wish to run within the Matlab environment, then you can start up Matlab and add the source files before executing the main function.

~~~~~
> ### In Linux ###
> cd <$install_dir>/APOLLOH_0.1.0/bin/
> <$matlab_dir>/bin/matlab
~~~~~

```Matlab
>> % In Matlab
>> addpath(genpath("<$install_dir>/APOLLOH_0.1.0/hmm"))
>> addpath(genpath("<$install_dir>/APOLLOH_0.1.0/util"))
>> % assuming your have all the arguments to the main function assigned...
>> apolloh(infile,cnfile,paramset,outparam,outfile)

```

## Inputs and Outputs ##
When using the compiled executable or running within the Matlab environment, the same input/output files and parameters are required:
```matlab
Function apolloh(infile,cnfile,paramset,outparam,outfile)

INPUTS:  
infile         Tab-delimited input file containing allelic counts
                 from the tumour at positions determined as heterozygous
                 from the normal genome. 
                 No header line is assumed.
                6 columns:
                    1) chr (integer; 'X' and 'Y' strings can be used)
                    2) position
                    3) reference base (can be arbitrary; not used)
                    4) referenc count
                    5) non-reference base (can be arbitrary; not used)
                    6) non-reference count 

cnfile         Tab-delimited input copy number segment prior file.
                  The accepted format is the output from HMMcopy,
                  a read-depth for analyzing copy number in tumour-
                  normal sequenced genomes.  
                  However, copy number segments from any source can be used.                
                 8-columns:
                    1) id (can be arbitrary, not used)
                    2) chr
                    3) start  
                    4) stop
                    5) Number of 1kb intervals (can be arbitrary; not used)
                    6) median log2 ratio (normal and tumour) for segment
                    7) HMM state: 1=HOMD (0 copies), 2=HEMD (1 copy),
                        3=NEUT (2 copies), 4=GAIN (3 copies),
                        5=AMP (4 copies), 6=HLAMP (5+ copies)
                        Note that for AMP and HLAMP, relative numbers of copies
                        can be used (i.e. GAIN is 3-4 copies, AMP is 5-6 copies,
                        HLAMP is 7+ copies)
                    8 ) CN state (can be arbitrary; not used)
                  If cnfile='0' is used, then copy number of 2 (diploid) is used
                    for all positions.  

paramset       Parameter intialization file is a matlab binary (.mat) file.
                      This file contains model and setting paramters necessary
                      to run the program.
                      See examples in "<$install_dir>/APOLLOH_0.1.0/parameters/".

outfile        Tab-delimited output file for position-level results. 
                  9-columns:
                     1) chr ('X' and 'Y' will be output as 23 and 24)
                     2) position
                     3) reference count
                     4) non-reference count
                     5) total depth
                     6) allelic ratio
                     7) copy number (from input)
                     8 ) APOLLOH genotype state
                     9) Zygosity state.
                  N additional columns:
                     posterior marginal probabilities (responsibilities) for
                     each APOLLOH genotype state. 
                 Zygosity states are:
                    DLOH=deletion-LOH (state 1)
                    NLOH=copy-neutral-LOH (states 2,4)
                    ALOH=amplified-LOH (states 5,8,9,13,14,19)
                    HET=heterozygous (states 3,6,7)
                    ASCNA=allele-specific-amplification (states 10,12,15,18)
                    BCNA=balanced-amplification (states 11,16,17)

                  Segment boundaries are determined as consecutive
                   marginal states of DLOH, NLOH, ALOH, HET, BCNA,
                   ASCNA; this implementation does not output this
                   information. An external Perl script handles this: "
                   <$install_dir>/APOLLOH_0.1.0/scripts/createSingleSegFileFromAPOLLOH.pl"

outparam       Tab-delimited output file storing converged parameters
                       after model training using Expectation Maximization (EM)
                       algorithm.
                     1) Number of iterations 
                     2) Global normal contamination parameter 
                     3) Binomial parameters for each HMM class/state.
```