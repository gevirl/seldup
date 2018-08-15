seldup
======

seldup reads short read alignments from a BAM file,
estimates the PCR duplicate read rate, and
selects apparent duplicate reads for removal
by the rmdupsam program.

## Notes

seldup
- appears to estimate PCR duplicate rates adequately when run on paired end reads but estimates poorly on single end reads.
- is written for unstranded libraries. It will not process correctly strand-specific libraries.
- does not consider read discrepancies in its analysis except to identify read ends in genomic coordinates.
- was written for unstranded, paired end read sets from C. elegans sequencing. It is untested on other data sets. Please test it thoroughly before using it for any purpose.
- please do not use this program to remove duplicates unless you are certain that it will work on your data sets.

## Building seldup

seldup.c requires a C compiler in order to make
an executable file. There is a Makefile to simplify
compiling and linking.

## Running seldup

Usage: seldup <input_bam_filename> <output_root_filename> <check_sl_flag> <single_end/paired_end_flag>
- <check_sl_flag> 0=do not check for spurious SL matches; 1=check for spurious SL matches
- <single_end/paired_end_flag> 1=single end only; 2=paired end only; 3=mixed single and paired end

See [seldup.c](seldup.c) for additional information.

## License

See the [license](LICENSE.md) document.
