/** seldup.c **/

/*
*|   File: seldup.c
*|                                                                            |*
*|   Copyright (c) 2016-2018 University of Washington All rights reserved.    |*
*|                                                                            |*
*|   Redistribution and use in source and binary forms, with or without       |*
*|   modification, are permitted provided that the following conditions are   |*
*|   met:                                                                     |*
*|                                                                            |*
*|   Redistributions of source code must retain the above copyright notice,   |*
*|   this list of conditions and the following disclaimer.                    |*
*|                                                                            |*
*|   Redistributions in binary form must reproduce the above copyright        |*
*|   notice, this list of conditions and the following disclaimer in the      |*
*|   documentation and/or other materials provided with the distribution.     |*
*|                                                                            |*
*|   Neither the name of the University of Washington nor the names of its    |*
*|   contributors may be used to endorse or promote products derived from     |*
*|   this software without specific prior written permission.                 |*
*|                                                                            |*
*|   This software is provided by the university of washington and            |*
*|   contributors "as is" and any express or implied warranties, including,   |*
*|   but not limited to, the implied warranties of merchantability and        |*
*|   fitness for a particular purpose are disclaimed. In no event shall the   |*
*|   University of Washington or contributors be liable for any direct,       |*
*|   indirect, incidental, special, exemplary, or consequential damages       |*
*|   (including, but not limited to, procurement of substitute goods or       |*
*|   services; loss of use, data, or profits; or business interruption)       |*
*|   however caused and on any theory of liability, whether in contract,      |*
*|   strict liability, or tort (including negligence or otherwise) arising    |*
*|   in any way out of the use of this software, even if advised of the       |*
*|   possibility of such damage.                                              |*
*/

/*
** PCR duplicate removal description
** 
** steps
** 1. read alignments
**     •o skip alignments flagged as chimeric or alternate
**     •o adjust read ends for indels and clipping. This 
**        gives a read start and end in the chromosome 
**        coordinates.
** 
** 2. find regions with high density of read starts (high coverage regions)
**     •o count and record read starts at each chromosome base
**     •o consider a chromosome base with one or more read 
**        starts to be occupied
**     •o detect high coverage regions using the d-segment dynamic 
**        programming algorithm to scan each chromosome for 
**        regions with relatively high densities of occupied bases. 
**        In short, the algorithm finds maximal scoring segments by 
**        calculating the sum of positive scores for occupied bases 
**        and negative scores for unoccupied bases. For the purpose 
**        of evaluating the high coverage regions, we find four 
**        sets of high coverage regions where the sets differ by 
**        the target minimum occupied base density. The target 
**        densities and corresponding d-segment parameters are 
**        0.01, 0.1, 1.0 and 10.0 occupied bases per chromosome 
**        base, the occupied base scores are 500, 50, 5, and 2, 
**        respectively, the drop-off, D, values are -500, -50,
**        -5, and -2, respectively, and the minimum ratios of 
**        occupied to total chr bases are .01, .1, 1.0, and 
**        10.0, respectively. The unoccupied base score is -1 and 
**        the minimum segment length is 10. Regions identifed by 
**        the d-segment algorithm as having on the order of >= 0.1 
**        occupied bases per chromosome base are flagged as 
**        high coverage regions for these paired end read data. 
**        Regions not flagged as high coverage are considered 
**        to be low coverage regions.
** 
** 3. estimate PCR duplicate rate
**     •o assuming that finding by chance two or more reads 
**        with the same start and end in the low coverage 
**        region is near zero, we consider such reads to be 
**        PCR duplicates
**     •o these data consist of paired end reads so we use the 
**        starts and ends of each read pair to identify observed 
**        duplicate read pairs. Incidentally, we found that 
**        estimating PCR duplicate rates from single-end read 
**        data using this method is unreliable because of 
**        substantially increased likelihood of chance 
**        duplicates.
**     •o we count the number of reads and duplicate reads in 
**        low coverage regions with length >= 200 bases and 
**        record the estimated PCR duplicate rate as the ratio 
**        of the two counts
**   
** 4. remove read pairs that are likely to be PCR duplicates
**     •o in each high coverage region individually, count the 
**        numbers of read starts 
**     •o calculate the number of duplicates to remove in each 
**        high coverage region using the estimated PCR duplicate 
**        rate and the number of read starts
**     •o choose randomly for removal read pairs in each high 
**        coverage region from the list of its observed duplicates
**     •o remove all read pair duplicates in the low coverage 
**        regions
*/

/*
** Notes:
**   o  filter out reads that are not primary but keep multi-mapped reads.
**      (This avoids multiply counting reads but preserves counts of
**      genuine reads. The remaining problem is assigning the read to
**      the correct gene... How often does this occur?)
**
**   o  identify duplicate reads to remove
**        o  tabulate high coverage regions (not low) by region: store
**           read start rate and duplicate rate for region
**        o  sort high coverage regions by chromosome and region start
**        o  sort duplicates by region and read start
**        o  move through the high coverage regions and mark duplicate
**           in high coverage regions and store high coverage region
**           index
**
**   o  this program knows nothing about strandedness so it cannot
**      process properly strand-specific data sets
**
**   o  this program treats implicitly reads that start and end
**      on the same bases of different diploid chromosomes to be
**      chance duplicates. That is, it makes no provision for
**      diploid genomes but this should not affect the results
**      greatly. The PCR duplicate rate estimate should be unaffected
**      by diploid genomes because it looks for regions with very
**      low numbers of read starts.
**
**   o  I want to remove very short spurious SL matches...I do not care
**      so much about very short spurious non-SL matches because they
**      will tend to be random, and so will not cluster at SL-like
**      locations.
**
**   o  things to think about
**        o  I use type 'int' and 'int32_t' in places where, someday,
**           wrap-around may happen. Test for wrap-around or change
**           type to int64_t. An example, in the number BAM file
**           alignments. The chromosome lengths may exceed 2Gbases,
**           which will affect BAM files too.
**         
**   o  d-segment scoring is somewhat ad-hoc and tends to be conservative
**      in order to avoid merging nearby higher coverage regions. I was
**      advised to think about the d-segment algorithm as a two-state
**      HMM. I look at the HMM model and find that for the chromosome
**      regions with coverage below 0.01, the occupancy rate is 4.046892e-05
**      in chromosome I of 20120411_EMB-0. The occupied scores for
**      the cutoffs are
**
**        cutoff    occupied score     my score
**         .01       553.4             500
**         .1         77.3              50
**        1.0          9.1               5
**
**      The calculated transition frequencies based on WormBase annotated
**      protein coding exons are
**        a11: 0.995359  a10: 0.004641
**        a00: 0.996489  a01: 0.003511
**      so |D|=4.82. In all cases -D > S so set -D to S. My scores and
**      dropoffs look reasonable so I leave them in place.
**
**   o  d-segment finds currently low-read occupancy base regions in the
**      regions where the base occupancy is less than 1. That is regions
**      where there are few bases with any read starts. The issue here
**      is not coverage because the coverage varies with PCR duplication
**      rate. However, for higher coverage rate regions, one may use
**      coverage, possibly corrected for the estimated PCR duplicate
**      rate. 
**
**   o  add command line options
**        o  histograms of terminal match lengths by occupancy region
**        o  add usage warning at startup that points to a command line
**           option with extended PCR duplication rate information
**
**   o  thoughts on PCR duplicate removal in general
**        o  reasons for removing PCR duplicates where the problem
**           may affect the calculated values
**             o  an analysis depends on absolute read counts
**             o  PCR amplification has context biases that may
**                affect adversely the analysis and duplicates are
**                likely to worsen the problem
**             o  duplicates are distributed non-uniformly creating
**                a patchiness or clumping of reads
**        o  reasons for removing PCR duplicates where the problem
**           may affect the significance of statistical tests
**             o  PCR duplicates are not independent data so significance
**                tests will over-state the significance of a result
**                (Phil Green points this out and adds that the
**                significance test may be modified to adjust for
**                duplicate reads but it's better probably to remove
**                the duplicates.)
**        o  in a recent conversation, Phil expressed a clear preference
**           for removing duplicates because
**             o  it's difficult to assess whether and how duplicates
**                affect an analysis (e.g, MISO) so it makes sense to
**                try to remove the duplicates
**             o  in cases where some analyses are duplicate sensitive
**                and others are insensitive such that one might be
**                tempted to run some analyses without removing duplicates
**                and others with (selected) duplicates removed, one might
**                invite criticism for using more than one data set.
**             o  his intuition suggests that the merging of blocks
**                with different expression levels is not a serious
**                problem (I compared gene, transcript, and exon dcpm
**                value differences for with and without duplicate
**                alignment sets and found the, for dcpm values >= .1,
**                the (relative) differences were small for projects with
**                lower estimated PCR duplicate rates and increased
**                modestly as the estimated PCR duplicate rates increased
**                but remained reasonable (with the possible exception
**                of a project with about 67% estimated PCR duplicate
**                rate). These data sets are C. elegans RNA-seq,
**                not strand-specific, and with between about 15M and
**                300M paired end reads.)
**        o  reasons for not removing all duplicate reads
**             o  highly expressed transcripts have numerous chance
**                duplicates that are independent data so removing
**                all duplicates erroneously depresses read counts
**        o  reasons for not removing selected duplicates
**             o  many analyses give results that are normalized by
**                the total number of reads in which case there may
**                be no need to remove duplicates
**             o  there is no way to distinguish between PCR and chance
**                duplicates (without adding unique identifiers to the
**                read sequence) so any method for attempting to remove
**                PCR duplicates while leaving chance duplicates will
**                be a compromise (of some sort) that creates artifacts.
**
**   o  thoughts on PCR duplicate removal by seldup
**        o  seldup distinguishes between low coverage and high
**           coverage chromosome bases where low coverage bases
**           have a low likelihood of having chance duplicate
**           reads. It uses a d-segment algorithm to find identify
**           the high coverage regions and the remaining chromosome
**           bases are low coverage.
**        o  the d-segment algorithm is a modified maximal scoring
**           segment detector. In this implementation, chromosome
**           bases with at least one read start (occupied) contribute
**           a positive score and chromosome bases with no read start
**           contribute a negative score to the score sum.
**        o  a candidate high coverage region begins with an occupied
**           chromosome base and ends when the score sum falls below
**           a threshold value.
**        o  a candidate high coverage region is accepted if its length
**           is at least a minimum value and the average coverage is at
**           least a minimum value
**        o  regarding duplicate removal, an important concern about
**           high coverage regions is that a region can contain
**           'blocks' of chromosome bases with differing expression
**           levels. These blocks will be merged when there are too few
**           unoccupied chromosome bases between the blocks, for
**           example when genes (nearly) overlap and when genes have
**           alternatively spliced transcripts. If the blocks have
**           different numbers of chance duplicates, the uniform
**           random selection of duplicates for removal will be more
**           likely to remove duplicates from the blocks with larger
**           numbers of chance duplicates, which results in block
**           coverage inaccuracies.
**        o  the threshold for high coverage regions is set low
**           in order to minimize the number of chance duplicates
**           falling in the low coverage regions. It's important to
**           keep chance duplicates out of low coverage regions because
**           they elevate the PCR duplicate rate estimate, and because
**           seldup removes all duplicates in the low coverage regions.
**           However, using a low threshold for detecting high coverage
**           regions increases the likelihood of merging chromosome
**           blocks with different coverage.
**        o  one can reduce the effect of merged chromosome blocks by
**           modifying this program to divide high coverage regions
**           into two or more sets of regions that differ by their
**           'coverage' (or read start rate), and removing duplicates
**           from these regions independently. Additionally, one can
**           use a less sensitive search for higher coverage regions
**           for the purpose of removing duplicates while keeping the
**           more sensitive search for estimating the PCR duplicate
**           rate.
**        o  it's possible to add diagnostic monitors to the program,
**           for example,
**             o  check the PCR duplicate estimate value by finding
**                the empirical distribution duplicate counts at each
**                chromosome base and comparing it to the corresponding
**                Poisson distribution.
**             o  report the high coverage region length distribution
**                (I added this although it can be extended and have
**                details added)
**
**   o  some additional thoughts on PCR duplicate removal
**        o  genes below .1 dcpm are more affected by duplicate
**           removal than those above. One explanation for this is
**           that genes with dcpm < .1 are most likely to be in
**           very low coverage regions where all duplicates are
**           removed. If a region has a single read with a duplicate,
**           then the numerator is halved by the duplicate removal.
**           It's unlikely that the denominator is halved as well
**           so the removal is likely to create a substantial
**           difference.
**        o  projects with higher PCR duplicate rates are more
**           affected by duplicate removal than those with low
**           duplicate rates. I think that a uniform chance of
**           removing each duplicate cannot compensate for the
**           increasing numbers of PCR duplicate read starts at
**           various chromosome bases (with a Poisson distribution).
**           This affects genes with small dcpm values more
**           strongly.
**        o  the relative dcpm differences can be used to
**             o  identify projects with unacceptably large numbers/values
**                of relative dcpm differences
**             o  verify analysis results for specific genes, transcripts,
**                and exons if there is a concern that the results are
**                affected by duplicate removal artifacts
**             o  evaluate modifications to the duplicate removal program
**                during development (should we decide to improve it)
**
**   o  possible improvements
**        o  divide current high occupancy region into two or more
**           occupancy ranges and remove duplicates in regions+level
**           Note: I believe that higher level regions are short and
**           relatively sparse so they're less likely to merge with
**           some exceptions such as certain alternatively spliced
**           exons. Check this supposition. It was suggested that
**           the program scan the initial high coverage regions for
**           block merging.
**        o  for high coverage regions with occupancy 1.0 (and high
**           averge coverage such as 10+), I can use a d-segment
**           algorithm that uses scores based on coverage rather than
**           occupancy. I can adjust the coverage rate by the estimated
**           PCR duplicate rate, if necessary (think about this).
*/

#define SE_VERY_LOW_REG
#define PE_VERY_LOW_REG


// #define VERIFY_LONG
#define FASTER_PROC_REG

#define VERSION 	"20180628"


#define MXLENREAD	1024
#define MXCOUNTER	128
#define MXHIST		8192

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <limits.h>

#include "util.h"
#include "bamutil.h"
#include "seldup.h"


#define BF1_B00  ( (uint8_t)( 1 << 0 ) )
#define BF1_B01  ( (uint8_t)( 1 << 1 ) )
#define BF1_B02  ( (uint8_t)( 1 << 2 ) )
#define BF1_B03  ( (uint8_t)( 1 << 3 ) )
#define BF1_B04  ( (uint8_t)( 1 << 4 ) )
#define BF1_B05  ( (uint8_t)( 1 << 5 ) )
#define BF1_B06  ( (uint8_t)( 1 << 6 ) )
#define BF1_B07  ( (uint8_t)( 1 << 7 ) )


#define BF2_B00	( (uint16_t)( 1 <<  0 ) )
#define BF2_B01	( (uint16_t)( 1 <<  1 ) )
#define BF2_B02	( (uint16_t)( 1 <<  2 ) )
#define BF2_B03	( (uint16_t)( 1 <<  3 ) )
#define BF2_B04	( (uint16_t)( 1 <<  4 ) )
#define BF2_B05	( (uint16_t)( 1 <<  5 ) )
#define BF2_B06	( (uint16_t)( 1 <<  6 ) )
#define BF2_B07	( (uint16_t)( 1 <<  7 ) )
#define BF2_B08	( (uint16_t)( 1 <<  8 ) )
#define BF2_B09	( (uint16_t)( 1 <<  9 ) )
#define BF2_B10	( (uint16_t)( 1 << 10 ) )
#define BF2_B11	( (uint16_t)( 1 << 11 ) )
#define BF2_B12	( (uint16_t)( 1 << 12 ) )
#define BF2_B13	( (uint16_t)( 1 << 13 ) )
#define BF2_B14	( (uint16_t)( 1 << 14 ) )
#define BF2_B15	( (uint16_t)( 1 << 15 ) )


#ifdef SE_VERY_LOW_REG
#define SE_BIT_TEST	( BF2_B08 | BF2_B09 | BF2_B10 | BF2_B11 )
#else
#define SE_BIT_TEST	( BF2_B09 | BF2_B10 | BF2_B11 )
#endif

#ifdef PE_VERY_LOW_REG
#define PE_BIT_TEST	( BF2_B09 | BF2_B10 | BF2_B11 )
#else
#define PE_BIT_TEST	( BF2_B10 | BF2_B11 )
#endif


/*
** xflag
**    bit    description
**      0    0=not proper paired read; 1=proper paired read (proper pair by test in xprocBam)
**      1    0=not detected duplicate read; 1=detected duplicate read
**      2    0=not in high coverage region; 1=in high coverage region
**      3    0=not removed duplicate; 1=removed duplicate
*/
typedef struct
{
  int32_t  refid;
  char    *read_name;
  int32_t  pos;
  int32_t  mapq;
  uint32_t flag;   /* sam flag */
  int32_t  pos5p;
  int32_t  tlen; 
  int32_t  index;  /* alignment index */
  int32_t  sbasq;
  uint8_t  xflag;
  int32_t  mateindex;  /* AlData index of proper pair mate alignment */
  int32_t  hcrindex; /* high coverage region index, if relevant */
} AlData;


/*
** flag: 
**   bit     description
**     2     detected duplicate read
*/
typedef struct
{
  int32_t  index[2];  /* alignment indices */
  char    *read_name;
  int32_t  refid[2];
  int32_t  pos[2];
  int32_t  pos5p[2];
  int32_t  mapq[2];
  int32_t  sbasq[2];
  int32_t  tlen;
  uint8_t  xflag;
  uint8_t  nalign;
} RedPair;


/*
** flag: same as RedPair.
*/
typedef struct
{
  int32_t  index;    /* alignment index */
  char    *read_name;
  int32_t  refid;
  int32_t  pos;
  int32_t  pos5p;
  int32_t  mapq;
  int32_t  sbasq;
  int32_t  tlen;
  uint8_t  xflag;
} RedOne;


typedef struct
{
  int refid;
  int beg;
  int end;
  int nstart;
  int ndup;
} HiCovReg;


typedef struct
{
  char       *nameBamFile;
  char       *starttime;
  int         numHiCovReg;
  HiCovReg   *hiCovReg;
  int         checkSLFlag; 
} Bundle;


typedef struct
{
  int sfalse;      /* no read start score */
  int strue;       /* yes read start score */
  int dscore;      /* D-segment score difference */
  int lmin;        /* minimum region length for reporting */
  double rmin;     /* minimum read start rate for reporting */
  uint16_t bflag;  /* bit flag value for chromosome bases */
} RegPar;


typedef struct
{
  int index;
  int flag;
  int mateidup;
} DupSet;


static int xallocGenome( BamHeader *bamHeader, Genome *genome, int *fstatus );
static int xinitGenome( Genome *genome, int *fstatus );
static int xprocBam( BGZF *fp, Bundle *bundle, BamHeader *bamHeader, AlData **falData, int *fnumAlData, int64_t counts[MXCOUNTER], int *fstatus );
static int xfindStartReadSE( AlData *alData, int numAlData, Genome *genome, int *fnumAlignSE, int64_t counts[MXCOUNTER], int *fstatus );
static int xfindStartReadPE( AlData *alData, int numAlData, Genome *genome, int *fnumAlignPE, int64_t counts[MXCOUNTER], int *fstatus );
static int xcheckSLMatch( BamAlign *bamAlign, BamDsc *bamDsc, int numBamDsc, int clip5[2], int clip3[2], int *fslFlag, int64_t *fiadj, int64_t counts[MXCOUNTER], int *fstatus );
static int xcheckShortMatch( BamAlign *bamAlign, BamDsc *bamDsc, int numBamDsc, int clip5[2], int clip3[2], int *fshortFlag, int64_t *fiadj, int64_t counts[MXCOUNTER], int *fstatus );
static int xtestUniqMap( BamAlign *align, int *funiqMapFlag, int *fstatus );
static int xbamCalcEndsRead( BamHeader *bamHeader, BamAlign *bamAlign, BamDsc *bamDsc, int numBamDsc, int clip5[2], int clip3[2], int checkSLFlag, int64_t *fbpos, int64_t *fepos, int64_t counts[MXCOUNTER], int *fstatus );
static int xmarkCovReg( RegPar *regPar, Genome *genome, int rsFlag, int64_t counts[MXCOUNTER], int *fstatus );
static int xmarkCovRegSet( Genome *genome, int rsFlag, int64_t counts[MXCOUNTER], int *fstatus );
static int xcountCovRegion( Genome *genome, uint16_t bit_true, uint16_t bit_false, int minSpc, int minLen, int *fnstart, int *fndup, int64_t counts[MXCOUNTER], int *fstatus );
static int xcountLoCovRegion( Genome *genome, uint16_t bit_false, int minSpc, int minLen, int *fnstart, int *fndup, int32_t rsFlag, int64_t hist[MXHIST], int64_t counts[MXCOUNTER], int *fstatus );
static int xestimatePCRDupRate( Bundle *bundle, BamHeader *bamHeader, Genome *genome, int rsFlag, double *fpcrDupRate, int64_t hist[MXHIST], int64_t counts[MXCOUNTER], int *fstatus );
static int xprocRegAlign( Genome *genome, HiCovReg *hiCovReg, AlData *alData, int ialign0, int ialign1, int nstart, int ndup, double pcrDupRate, int *fnrem, int *fstatus );
static int xreportHiCovRegDist( HiCovReg *hiCovReg, int nreg, int *fstatus );
static int xremoveDuplicateReads( Bundle *bundle, BamHeader *bamHeader, AlData *alData, int numAlData, Genome *genome, int rsFlag, double pcrDupRate, int64_t counts[MXCOUNTER], int *fstatus );
static int xreportCounts( FILE *ofp, Bundle *bundle, double pcrDupRateSE, double pcrDupRatePE, int64_t hist[2][MXHIST], int64_t counts[3][MXCOUNTER], int *fstatus );
static int xdumpStartRead( char *nameRoot, Genome *genome, int *fstatus );
static int xdumpCovReg( Genome *genome, int *fstatus );
static int xdumpCovNonReg( Genome *genome, int *fstatus );
static int xcheckAlign( Bundle *bundle, BamHeader *bamHeader, AlData *alData, int numAlData, Genome *genome, int64_t counts[MXCOUNTER], int *fstatus );
static int xcheckRegAlign( Genome *genome, AlData *alData, int ialign0, int ialign1, HiCovReg *hiCovReg, int nstart, int ndup, int *fstatus );
static int xcheckDupSet( AlData *alData, DupSet *dupSet, int ndup, HiCovReg *hiCovReg, int ireg, int *fstatus );
static int xwriteReads( FILE *ofp, Bundle *bundle, BamHeader *bamHeader, AlData *alData, int numAlData, Genome *genome, int *fstatus );
static int xinitXorShift1024( void );
uint64_t xorShift1024(void);


/*
** xorshift1024 internal variable.
*/
extern uint64_t xrs[16];


int main( int argc, char **argv )
{
  int i;
  int numAlData;
  int numAlignSE;
  int numAlignPE;
  int sepeFlag;
  int status;

  char string[8192];
  char *nameBamFile;
  char *nameRootOutFile;
  char starttime[8192];

  double pcrDupRateSE;
  double pcrDupRatePE;

  int64_t counts[3][MXCOUNTER];
  int64_t hist[2][MXHIST];

  Bundle    bundle;
  BGZF     *fp;
  BamHeader bamHeader;
  Genome    genome;
  AlData   *alData;

  FILE     *lfp;
  FILE     *rfp;

  fprintf( stderr, "seldup reads a BAM file and estimates the PCR duplicate\n" );
  fprintf( stderr, "read rate and then selects duplicates for removal. It is\n" );
  fprintf( stderr, "impossible to distinguish between PCR and chance\n" );
  fprintf( stderr, "duplicate reads without barcode labeled reads so the\n" );
  fprintf( stderr, "estimation of a PCR duplicate rate and the selection of\n" );
  fprintf( stderr, "duplicates for removal involve methods that depend on\n" );
  fprintf( stderr, "compromises and chance. Factors that may affect adversely\n" );
  fprintf( stderr, "either or both the PCR duplicate rate estimate and the\n" );
  fprintf( stderr, "selection of duplicates for removal include but are likely\n" );
  fprintf( stderr, "not limited to\n" );
  fprintf( stderr, "  o  the distributions of expressed regions in the genome,\n" );
  fprintf( stderr, "     for example, the lengths of low and high expression\n" );
  fprintf( stderr, "     level regions and how they are interspersed and\n" );
  fprintf( stderr, "     overlap\n" );
  fprintf( stderr, "  o  insert length distribution for paired end reads\n" );
  fprintf( stderr, "  o  the number of reads in the data set\n" );
  fprintf( stderr, "  o  the mis-alignment rate and characteristics of the\n" );
  fprintf( stderr, "     alignment program\n" );
  fprintf( stderr, "  o  characteristics of the genome that might contribute\n" );
  fprintf( stderr, "     to misalignments such as trans-splicing\n" );
  fprintf( stderr, "Due to the much greater rate of chance duplicates in single\n" );
  fprintf( stderr, "end read data sets in comparison to paired end read sets,\n" );
  fprintf( stderr, "this program is more likely to give less accurate estimates\n" );
  fprintf( stderr, "of PCR duplicate rates and selection of duplicates for removal\n" );
  fprintf( stderr, "when run on single end read data sets. This program does not\n" );
  fprintf( stderr, "account for strand-specific reads or heterozygosity in its\n" );
  fprintf( stderr, "current form.\n" );
  fprintf( stderr, "This program is unsuitable for strand-specific reads in its\n" );
  fprintf( stderr, "form.\n" );
  fprintf( stderr, "Please see the seldup.c source code file header for additional\n" );
  fprintf( stderr, "information and evaluate carefully the performance of this\n" );
  fprintf( stderr, "program on your data before using it for important analyses.\n" );

  if( argc == 5 )
  {
    nameBamFile        = mstrcpy( argv[1] );
    nameRootOutFile    = mstrcpy( argv[2] );
    bundle.checkSLFlag = atoi( argv[3] );
    sepeFlag           = atoi( argv[4] );
  }
  else
  {
    printf( "\n" );

    printf( "enter BAM filename: " );
    gets( string );
    nameBamFile = mstrcpy( string );

    printf( "enter root output filename: " );
    gets( string );
    nameRootOutFile = mstrcpy( string );

    printf( "enter check SL flag (0=do not check for spurious SL matches; 1=check for spurious SL matches): " );
    gets( string );
    bundle.checkSLFlag = atoi( string );

    printf( "enter single/paired end flag (1=single end only; 2=paired end only; 3=both single and paired end): " );
    gets( string );
    sepeFlag = atoi( string );
  }

  if( bundle.checkSLFlag < 0 || bundle.checkSLFlag > 1 )
  {
    sprintf( _msg_, "bad check SL flag value (%d): must be 0 or 1", bundle.checkSLFlag );
    EMSG( _msg_ );
    return( -1 );
  }

  if( sepeFlag < 1 || sepeFlag > 3 )
  {
    sprintf( _msg_, "bad single/paired end flag value (%d): must be between 1 and 3", sepeFlag );
    EMSG( _msg_ );
    return( -1 );
  }

  bundle.nameBamFile      = mstrcpy( nameBamFile );
  sprintf( string, "%s.log", nameRootOutFile );
  lfp = fopen( string, "w+" );
  if( lfp == NULL )
  {
    sprintf( _msg_, "unable to open file %s", string );
    EMSG( _msg_ );
    return( -1 );
  }

  timeofday( starttime );
  bundle.starttime = starttime;

  for( i = 0; i < argc; ++i )
  {
    fprintf( stderr, " %s", argv[i] );
  }
  fprintf( stderr, "\n" );
  fprintf( stderr, "Program version:           %s\n", VERSION );
  fprintf( stderr, "Runtime:                   %s\n", starttime );
  fprintf( stderr, "BAM file:                  %s\n", nameBamFile );
  fprintf( stderr, "Output root filename:      %s\n", nameRootOutFile );
  fprintf( stderr, "Check SL flag:             %d\n", bundle.checkSLFlag );
  fprintf( stderr, "SE/PE flag:                %d\n", sepeFlag );
  fprintf( stderr, "SE bit test:               %d%d%d%d\n", SE_BIT_TEST & BF2_B08 ? 1 : 0, SE_BIT_TEST & BF2_B09 ? 1 : 0, SE_BIT_TEST & BF2_B10 ? 1 : 0, SE_BIT_TEST & BF2_B11 ? 1 : 0 );
  fprintf( stderr, "PE bit test:               %d%d%d%d\n", PE_BIT_TEST & BF2_B08 ? 1 : 0, PE_BIT_TEST & BF2_B09 ? 1 : 0, PE_BIT_TEST & BF2_B10 ? 1 : 0, PE_BIT_TEST & BF2_B11 ? 1 : 0 );
#ifdef FASTER_PROC_REG
  fprintf( stderr, "xprocRegAlign modified for efficiency\n" );
#else
  fprintf( stderr, "xprocRegAlign not modified for efficiency\n" );
#endif

  fprintf( lfp, "Command line:  " );
  for( i = 0; i < argc; ++i )
  {
    fprintf( lfp, " %s", argv[i] );
  }
  fprintf( lfp, "\n" );
  fprintf( lfp, "Program version:           %s\n", VERSION );
  fprintf( lfp, "Runtime:                   %s\n", starttime );
  fprintf( lfp, "BAM file:                  %s\n", nameBamFile );
  fprintf( lfp, "Output root filename:      %s\n", nameRootOutFile );
  fprintf( lfp, "Check SL flag:             %d\n", bundle.checkSLFlag );
  fprintf( lfp, "SE/PE flag:                %d\n", sepeFlag );

  for( i = 0; i < MXCOUNTER; ++i )
  {
    counts[0][i] = (int64_t)0;
    counts[1][i] = (int64_t)0;
    counts[2][i] = (int64_t)0;
  }

  for( i = 0; i < MXHIST; ++i )
  {
    hist[0][i] = (int64_t)0;
    hist[1][i] = (int64_t)0;
  }

  /*
  ** Open BAM file.
  */
  fprintf( stderr, "Open BAM file...\n" );
  fp = bgzf_open( nameBamFile, "r" );
  if( fp == NULL )
  {
    sprintf( _msg_, "unable to open BAM file %s", nameBamFile );
    EMSG( _msg_ );
    return( -1 );
  }

  /*
  ** Read BAM file header.
  */
  fprintf( stderr, "Read BAM file header...\n" );
  bamReadHeader( fp, &bamHeader, &status );
  IF_STATUS_MAIN( "bad status: bamReadHeader" );

  /*
  ** Report BAM file information.
  */
  fprintf( stderr, "BAM header program information\n" );
  fprintf( stderr, "  name:        %s\n", bamHeader.program.name != NULL ? bamHeader.program.name : "N/A" );
  fprintf( stderr, "  command:     %s\n", bamHeader.program.command != NULL ? bamHeader.program.command : "N/A" );
  fprintf( stderr, "  description: %s\n", bamHeader.program.description != NULL ? bamHeader.program.description : "N/A" );
  fprintf( stderr, "  version:     %s\n", bamHeader.program.version != NULL ? bamHeader.program.version : "N/A" );
  fprintf( stderr, "\n" );

  fprintf( lfp, "BAM header program information\n" );
  fprintf( lfp, "  name:        %s\n", bamHeader.program.name != NULL ? bamHeader.program.name : "N/A" );
  fprintf( lfp, "  command:     %s\n", bamHeader.program.command != NULL ? bamHeader.program.command : "N/A" );
  fprintf( lfp, "  description: %s\n", bamHeader.program.description != NULL ? bamHeader.program.description : "N/A" );
  fprintf( lfp, "  version:     %s\n", bamHeader.program.version != NULL ? bamHeader.program.version : "N/A" );
  fprintf( lfp, "\n" );

  /*
  ** Allocate chromosome structures.
  */
  fprintf( stderr, "Allocate chromosome structures...\n" );
  xallocGenome( &bamHeader, &genome, &status );
  IF_STATUS_MAIN( "bad status: xallocGenome" );

  /*
  ** Process BAM file: read alignments, filter, and count coverage.
  */
  fprintf( stderr, "Read alignments and count coverage...\n" );
  xprocBam( fp, &bundle, &bamHeader, &alData, &numAlData, counts[0], &status );
  IF_STATUS_MAIN( "bad status: xprocBam" );

  /*
  ** Close BAM file.
  */
  fprintf( stderr, "Close BAM file...\n" );
  bgzf_close( fp );


  if( sepeFlag & 0x2 )
  {
    /*
    ** Initialize genome counters/flags for paired end duplicate removal.
    */
    fprintf( stderr, "Initialize genome counters and flags...\n" );
    xinitGenome( &genome, &status );
    IF_STATUS_MAIN( "bad status: xinitGenome" );

    fprintf( stderr, "Find paired end read starts and duplicates...\n" );
    xfindStartReadPE( alData, numAlData, &genome, &numAlignPE, counts[2], &status );
    IF_STATUS_MAIN( "bad status: xfindStartReadPE" );
    fprintf( stderr, "%d paired end alignments processed\n", numAlignPE );
  
    if( numAlignPE > 0 )
    {
      /*
      ** Mark coverage regions.
      */
      fprintf( stderr, "Mark region coverage set...\n" );
      xmarkCovRegSet( &genome, RS_PAIRED, counts[2], &status );
      IF_STATUS_MAIN( "bad status: xmarkCovRegSet" );

/*
      xdumpStartRead( "", &genome, &status );
      IF_STATUS_MAIN( "bad status: xdumpStartRead" );
*/

      /*
      ** Estimate PCR duplication rate.
      */
      fprintf( stderr, "Estimate PCR duplication rate...\n" );
      xestimatePCRDupRate( &bundle, &bamHeader, &genome, RS_PAIRED, &pcrDupRatePE, hist[0], counts[2], &status );
      IF_STATUS_MAIN( "bad status: xestimatePCRDupRate" );
    
      fprintf( stderr, "Identify duplicate reads for removal...\n" );
      xremoveDuplicateReads( &bundle, &bamHeader, alData, numAlignPE, &genome, RS_PAIRED, pcrDupRatePE, counts[2], &status );
      IF_STATUS_MAIN( "bad status: xremoveDuplicateReads" );

/*
      xdumpCovReg( &genome, &status );
      IF_STATUS_MAIN( "bad status: xdumpCovReg" );
    
      xdumpCovNonReg( &genome, &status );
      IF_STATUS_MAIN( "bad status: xdumpCovNonReg" );
*/

      /*
      ** Final consistency tests.
      */
      fprintf( stderr, "Consistency check...\n" );
      xcheckAlign( &bundle, &bamHeader, alData, numAlignPE, &genome, counts[2], &status );
      IF_STATUS_MAIN( "bad status: xcheckAlign" );
    } /* if numAlignPE */
  } /* if sepeFlag */


  if( sepeFlag & 0x1 )
  {
    /*
    ** Initialize genome counters/flags for single end duplicate removal.
    */
    fprintf( stderr, "Initialize genome counters and flags...\n" );
    xinitGenome( &genome, &status );
    IF_STATUS_MAIN( "bad status: xinitGenome" );
  
    fprintf( stderr, "Find single end read starts and duplicates...\n" );
    xfindStartReadSE( alData, numAlData, &genome, &numAlignSE, counts[1], &status );
    IF_STATUS_MAIN( "bad status: xfindStartReadSE" );
    fprintf( stderr, "%d single end alignments processed\n", numAlignSE );
  
    if( numAlignSE > 0 )
    {
      /*
      ** Mark coverage regions.
      */
      fprintf( stderr, "Mark region coverage set...\n" );
      xmarkCovRegSet( &genome, RS_SINGLE, counts[1], &status );
      IF_STATUS_MAIN( "bad status: xmarkCovRegSet" );

/*
** File names do not distinguish between SE and PE read sets.
       xdumpStartRead( "", &genome, &status );
       IF_STATUS_MAIN( "bad status: xdumpStartRead" );
*/

      /*
      ** Estimate PCR duplication rate.
      */
      fprintf( stderr, "Estimate PCR duplication rate...\n" );
      xestimatePCRDupRate( &bundle, &bamHeader, &genome, RS_SINGLE, &pcrDupRateSE, hist[1], counts[1], &status );
      IF_STATUS_MAIN( "bad status: xestimatePCRDupRate" );
    
      fprintf( stderr, "Identify duplicate reads for removal...\n" );
      xremoveDuplicateReads( &bundle, &bamHeader, alData, numAlignSE, &genome, RS_SINGLE, pcrDupRateSE, counts[1], &status );
      IF_STATUS_MAIN( "bad status: xremoveDuplicateReads" );

/*
      xdumpCovReg( &genome, &status );
      IF_STATUS_MAIN( "bad status: xdumpCovReg" );
    
      xdumpCovNonReg( &genome, &status );
      IF_STATUS_MAIN( "bad status: xdumpCovNonReg" );
*/

      /*
      ** Final consistency tests.
      */
      fprintf( stderr, "Consistency check...\n" );
      xcheckAlign( &bundle, &bamHeader, alData, numAlignSE, &genome, counts[1], &status );
      IF_STATUS_MAIN( "bad status: xcheckAlign" );
    }
  }

  /*
  ** Write reads in high coverage regions.
  */
  fprintf( stderr, "Write duplicates selected for removal...\n" );
  sprintf( string, "%s.rmdup.lst", nameRootOutFile );
  rfp = fopen( string, "w+" );
  if( rfp == NULL )
  {
    sprintf( _msg_, "unable to open file %s", string );
    EMSG( _msg_ );
    return( -1 );
  }
  xwriteReads( rfp, &bundle, &bamHeader, alData, numAlData, &genome, &status );
  IF_STATUS_MAIN( "bad status: xwriteReads" );
  fclose( rfp );

  /*
  ** Report counts.
  */
  fprintf( stderr, "Report counts...\n" );
  xreportCounts( stderr, &bundle, pcrDupRateSE, pcrDupRatePE, hist, counts, &status );
  IF_STATUS_MAIN( "bad status: xreportCounts" );
  xreportCounts( lfp, &bundle, pcrDupRateSE, pcrDupRatePE, hist, counts, &status );
  IF_STATUS_MAIN( "bad status: xreportCounts" );

  fprintf( stderr, "finished\n\n" );

  /*
  ** Close log file.
  */
  fclose( lfp );

  return( 0 );
}


static int xallocGenome( BamHeader *bamHeader, Genome *genome, int *fstatus )
{
  int i;
  int ichr;
  int nchr;
  int32_t len;

  nchr = bamHeader->numRefSeq;

  genome->numChr = nchr;
  genome->chr    = (Chromosome *)malloc( nchr * sizeof( Chromosome ) );
  if( genome->chr == NULL )
  {
    EMSG( "unable to allocate memory" );
    *fstatus = -1;
    return( 0 );
  }

  for( ichr = 0; ichr < nchr; ++ichr )
  {
    len = bamHeader->refSeq[ichr].length;

    genome->chr[ichr].name = mstrcpy( bamHeader->refSeq[ichr].name );
    genome->chr[ichr].len  = (int64_t)len;
    for( i = 0; i < 4; ++i )
    {
      genome->chr[ichr].cover[i] = (int64_t *)calloc( (size_t)len, sizeof( int64_t ) );
      if( genome->chr[ichr].cover[i] == NULL )
      {
        EMSG( "unable to allocate memory" );
        *fstatus = -1;
        return( 0 );
      }
    }

    genome->chr[ichr].flag = (uint16_t *)calloc( (size_t)len, sizeof( uint16_t ) );
    if( genome->chr[ichr].flag == NULL )
    {
      EMSG( "unable allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
  }

  *fstatus = 0;

  return( 0 );
}


static int xinitGenome( Genome *genome, int *fstatus )
{
  int i;
  int ichr, nchr;
  int len;

  Chromosome *chr;

  nchr = genome->numChr;
  for( ichr = 0; ichr < nchr; ++ichr )
  {
    chr = &(genome->chr[ichr]);
    len = chr->len;
    for( i = 0; i < 4; ++i )
    {
      memset( chr->cover[i], 0, len * sizeof( int64_t ) );
    }
    memset( chr->flag, 0, len * sizeof( uint16_t ) );
  }

  *fstatus = 0;

  return( 0 );
}


static int xcmpSortAlDataSE1( const void *fa, const void *fb )
{
  int ist;
  int ia, ib;

  AlData *a;
  AlData *b;

  a = (AlData *)fa;
  b = (AlData *)fb;

  ia = ( a->flag & SF_PAIREDEND ) ? 1 : 0;
  ib = ( b->flag & SF_PAIREDEND ) ? 1 : 0;

  ist = ia - ib;
  if( ist )
  {
    return( ist );
  }

  ist = a->refid - b->refid;
  if( ist )
  {
    return( ist );
  }

  return( a->pos5p - b->pos5p );
}


static int xcmpSortAlDataPE1( const void *fa, const void *fb )
{
  int ist;
  int ia, ib;

  AlData *a;
  AlData *b;

  a = (AlData *)fa;
  b = (AlData *)fb;

  ia = ( a->flag & SF_PAIREDEND ) ? 1 : 0;
  ib = ( b->flag & SF_PAIREDEND ) ? 1 : 0;

  ist = ib - ia;
  if( ist )
  {
    return( ist );
  }

  ist = a->refid - b->refid;
  if( ist )
  {
    return( ist );
  }

  return( a->pos5p - b->pos5p );
}


static int xcmpSortAlData2( const void *fa, const void *fb )
{
  int ist;
  int ia, ib;

  AlData *a;
  AlData *b;

  a = (AlData *)fa;
  b = (AlData *)fb;

  ist = strcmp( a->read_name, b->read_name );
  if( ist )
  {
    return( ist );
  }

  ia = ( a->flag & SF_PROPER_PAIR ) ? 1 : 0;
  ib = ( b->flag & SF_PROPER_PAIR ) ? 1 : 0;

  ist = ib - ia;
  if( ist )
  {
    return( ist );
  }
  
  return( a->pos - b->pos );
}


static int xcmpSortAlData3( const void *fa, const void *fb )
{
  AlData *a;
  AlData *b;

  a = (AlData *)fa;
  b = (AlData *)fb;

  return( a->index - b->index );
}


static int xcmpSortRedPair( const void *fa, const void *fb )
{
  RedPair *a;
  RedPair *b;
  int32_t amapq;
  int32_t bmapq;
  int32_t abasq;
  int32_t bbasq;
  int32_t dif;

  a = (RedPair *)fa;
  b = (RedPair *)fb;

  dif = a->refid[0] - b->refid[0];
  if( dif )
  {
    return( dif );
  }

  dif = a->pos5p[0] - b->pos5p[0];
  if( dif )
  {
    return( dif );
  }

  dif = a->refid[1] - b->refid[1];
  if( dif )
  {
    return( dif );
  }

  dif = a->pos5p[1] - b->pos5p[1];
  if( dif )
  {
    return( dif );
  }

  amapq = a->mapq[0] + a->mapq[1];
  bmapq = b->mapq[0] + b->mapq[1];
  dif = bmapq - amapq;
  if( dif )
  {
    return( dif );
  }

  abasq = a->sbasq[0] + a->sbasq[1];
  bbasq = b->sbasq[0] + b->sbasq[1];
  return( bbasq - abasq );
}


static int xcmpSortRedOne( const void *fa, const void *fb )
{
  RedOne *a;
  RedOne *b;
  int32_t dif;

  a = (RedOne *)fa;
  b = (RedOne *)fb;

  dif = a->refid - b->refid;
  if( dif )
  {
    return( dif );
  }

  dif = a->pos5p - b->pos5p;
  if( dif )
  {
    return( dif );
  }

  dif = b->mapq - a->mapq;
  if( dif )
  {
    return( dif );
  }

  return( b->sbasq - a->sbasq );
}


static int xprocBam( BGZF *fp, Bundle *bundle, BamHeader *bamHeader, AlData **falData, int *fnumAlData, int64_t counts[MXCOUNTER], int *fstatus )
{
  int ibas;
  int status;
  int cigarFlag;
  int mdFlag;
  int ndsc;
  int ialign, nalign;
  int uniqMapFlag;
  int allocBamDsc;

  int clip5[2];
  int clip3[2];

  int32_t  sbasq;
  uint32_t flag;
  int64_t  fp0;

  int64_t  bpos, epos;

  BamAlign    align;
  BamDsc     *dsc;
  AlData     *alData;

  fprintf( stderr, "Test BAM file...\n" );
  bamTestCigar( fp, &cigarFlag, &mdFlag, &status );
  IF_STATUS_ZERO( "bad status: bamTestCigar" );

  fprintf( stderr, "Init BAM align structure...\n" );
  bamInitAlign( &align, &status );
  IF_STATUS_ZERO( "bad status: bamInitAlign" );

  /*
  ** Get pointer to file start.
  */
  fp0 = bgzf_tell( fp );

  /*
  ** Count alignments in BAM file.
  */
  nalign = 0;
  fprintf( stderr, "Count alignments...\n" );
  while( bamReadAlign( fp, &align, &status ) && status == 0 )
  {
    flag = align.flag;

    /*
    ** Skip unaligned entries.
    */
    if( flag & SF_UNMAPPED )
    {
      bamFreeAlign( &align );
      continue;
    }

    /*
    ** Skip chimeric alignments.
    */
    if( flag & SF_SUP_ALN )
    {
      bamFreeAlign( &align );
      continue;
    }

    /*
    ** Skip alternate alignments.
    */
    if( flag & SF_SEC_ALN )
    {
      bamFreeAlign( &align );
      continue;
    }

#ifdef UNIQ_MAP
    xtestUniqMap( &align, &uniqMapFlag, &status );
    IF_STATUS_ZERO( "bad status: xtestUniqMap" );
    if( !uniqMapFlag )
    {
fprintf( stderr, "multimap: %s %u\n", align.read_name, align.flag );
      bamFreeAlign( &align );
      continue;
    }
fprintf( stderr, "uniqmap: %s %u\n", align.read_name, align.flag );
#endif

    ++nalign;

    if( nalign == INT_MAX )
    {
      EMSG( "too many alignments in BAM file" );
      *fstatus = -1;
      return( 0 );
    }

    bamFreeAlign( &align );
  }

  fprintf( stderr, "%d alignments counted\n", nalign );

  /*
  ** Allocate AlData.
  */
  alData = (AlData *)calloc( nalign, sizeof( AlData ) );
  if( alData == NULL )
  {
    EMSG( "unable to allocate memory" );
    *fstatus = -1;
    return( 0 );
  }
  
  dsc         = NULL;
  allocBamDsc = 0;

  bgzf_seek( fp, fp0, SEEK_SET );

  /*
  ** Process alignments in BAM file.
  ** Notes:
  **   o  skip alignments
  **        o  not aligned
  **        o  not chimeric
  **        o  not alternate
  **        o  not uniquely mapped (optional)
  */
  bpos   = 0;
  epos   = 0;
  ialign = 0;
  fprintf( stderr, "Read alignments...\n" );
  while( bamReadAlign( fp, &align, &status ) && status == 0 )
  {
    ++counts[C_NUM_ENTRY_BAM];

    flag = align.flag;

    /*
    ** Skip unaligned entries.
    */
    if( flag & SF_UNMAPPED )
    {
      ++counts[C_UNALIGNED_BAM];
      bamFreeAlign( &align );
      continue;
    }

    ++counts[C_ALIGNED_BAM];

    /*
    ** Skip chimeric alignments.
    */
    if( flag & SF_SUP_ALN )
    {
      ++counts[C_CHIMERIC_BAM];
      bamFreeAlign( &align );
      continue;
    }

    /*
    ** Count and skip alternate alignments.
    */
    if( flag & SF_SEC_ALN )
    {
      ++counts[C_SECONDARY_BAM];
      bamFreeAlign( &align );
      continue;
    }

#ifdef UNIQ_MAP
    xtestUniqMap( &align, &uniqMapFlag, &status );
    IF_STATUS_ZERO( "bad status: xtestUniqMap" );
    if( !uniqMapFlag )
    {
      bamFreeAlign( &align );
      continue;
    }
#endif

    ++counts[C_ACCEPTED_BAM];

    if( align.flag & SF_PROPER_PAIR )
    {
      ++counts[C_PROPERPAIR_BAM];
    }
    else
    {
      ++counts[C_NOT_PROPERPAIR_BAM];
    }

    /*
    ** Get discrepancies.
    */
    bamGetDsc( &align, cigarFlag, mdFlag, &dsc, &ndsc, &allocBamDsc, clip5, clip3, &status );
    IF_STATUS_ZERO( "bad status: bamGetDsc" );

    /*
    ** Get posBeg and posEnd.
    */
    xbamCalcEndsRead( bamHeader, &align, dsc, ndsc, clip5, clip3, bundle->checkSLFlag, &bpos, &epos, counts, &status );
    IF_STATUS_ZERO( "bad status: xbamCalcEndsRead" );

    sbasq = 0;
    for( ibas = 0; ibas < align.l_seq; ++ibas )
    {
      sbasq += align.qual[ibas];
    }

    /*
    **
    */
    alData[ialign].refid     = align.refid;
    alData[ialign].read_name = mstrcpy( align.read_name );
    alData[ialign].pos       = align.pos;
    alData[ialign].mapq      = align.mapq;
    alData[ialign].flag      = align.flag;
    alData[ialign].pos5p     = !(flag & SF_REV_CMP) ? bpos : epos;
    alData[ialign].tlen      = align.tlen;
    alData[ialign].sbasq     = sbasq;
    alData[ialign].xflag     = (uint8_t)0;
    alData[ialign].mateindex = (int32_t)-1;
    alData[ialign].hcrindex  = (int32_t)-1;

    ++ialign;

    bamFreeAlign( &align );
  }

  if( ialign != nalign )
  {
    sprintf( _msg_, "inconsistent alignment counts %d %d", ialign, nalign );
    EMSG( _msg_ );
    *fstatus = -1;
    return( 0 );
  }

  *falData    = alData;
  *fnumAlData = nalign;

  *fstatus = 0;

  return( 0 );
}


static int xfindStartReadPE( AlData *alData, int numAlData, Genome *genome, int *fnumAlignPE, int64_t counts[MXCOUNTER], int *fstatus )
{
  int ialign, nalign;
  int ialign0, ialign1;  
  int ipair, jpair, npair;
  int ione, jone, none;
  int iend;
  int ichr, ibas;

  RedPair *redPair;
  RedOne  *redOne;
  Chromosome *chr;

  nalign = numAlData;

  /*
  ** Sort AlData by paired/single end, genome position, and set index.
  ** Notes:
  **   o  conditions
  **        o  paired/single end
  **        o  refid
  **        o  pos5p
  */
  qsort( alData, nalign, sizeof( AlData ), xcmpSortAlDataPE1 );

  for( ialign = 0; ialign < nalign; ++ialign )
  {
    if( !( alData[ialign].flag & SF_PAIREDEND ) )
    {
      nalign = ialign;
      break;
    }
    alData[ialign].index = ialign;
  }

  /*
  ** Sort AlData by read name.
  ** Notes:
  **   o  conditions
  **        o  read name
  **        o  proper pair flag (reverse order)
  **        o  position in genome
  */
  qsort( alData, nalign, sizeof( AlData ), xcmpSortAlData2 );

  /*
  ** Count distinct proper pair read names.
  **   So, find first read of the first proper pair.
  */
  for( ialign = 0; ialign < nalign; ++ialign )
  {
    if( alData[ialign].flag & SF_PROPER_PAIR )
    {
      break;
    }
  }

  npair = 1;
  for( ialign += 1; ialign < nalign; ++ialign )
  {
    if( !(alData[ialign].flag & SF_PROPER_PAIR ) )
    {
      continue;
    }

    if( strcmp( alData[ialign].read_name, alData[ialign-1].read_name ) != 0 )
    {
      ++npair;
    }
  }

  fprintf( stderr, "%d proper read pairs counted\n", npair );

  /*
  ** Allocate redPair.
  */
  redPair = (RedPair *)calloc( npair, sizeof( RedPair ) );
  if( redPair == NULL )
  {
    EMSG( "unable to allocate memory" );
    *fstatus = -1;
    return( 0 );
  }

  for( ipair = 0; ipair < npair; ++ipair )
  {
    redPair[ipair].index[0] = -1;
    redPair[ipair].index[1] = -1;
  }

  /*
  ** Assign alignments to redPair.
  */
  for( ialign = 0; ialign < nalign; ++ialign )
  {
    if( alData[ialign].flag & SF_PROPER_PAIR )
    {
      break;
    }
  }

  ipair  = 0;
  iend   = !(alData[ialign].flag & SF_REV_CMP ) ? 0 : 1;
  redPair[ipair].index[iend] = alData[ialign].index;
  redPair[ipair].read_name   = alData[ialign].read_name;
  redPair[ipair].refid[iend] = alData[ialign].refid;
  redPair[ipair].pos5p[iend] = alData[ialign].pos5p;
  redPair[ipair].pos[iend]   = alData[ialign].pos;
  redPair[ipair].mapq[iend]  = alData[ialign].mapq;
  redPair[ipair].sbasq[iend] = alData[ialign].sbasq;
  redPair[ipair].xflag       = (uint8_t)0;
  ++redPair[ipair].nalign;
  if( alData[ialign].tlen > 0 )
  {
    redPair[ipair].tlen  = alData[ialign].tlen;
  }

  for( ialign += 1; ialign < nalign; ++ialign )
  {
    if( !( alData[ialign].flag & SF_PROPER_PAIR ) )
    {
      continue;
    }

    if( strcmp( alData[ialign].read_name, alData[ialign-1].read_name ) != 0 )
    {
      ++ipair;
    }

    iend = !(alData[ialign].flag & SF_REV_CMP ) ? 0 : 1;
    redPair[ipair].index[iend] = alData[ialign].index;
    redPair[ipair].refid[iend] = alData[ialign].refid;
    redPair[ipair].pos5p[iend] = alData[ialign].pos5p;
    redPair[ipair].pos[iend]   = alData[ialign].pos;
    redPair[ipair].mapq[iend]  = alData[ialign].mapq;
    redPair[ipair].sbasq[iend] = alData[ialign].sbasq;
    redPair[ipair].xflag       = (uint8_t)0;
    ++redPair[ipair].nalign;
    if( redPair[ipair].read_name == NULL )
    {
      redPair[ipair].read_name = alData[ialign].read_name;
    }
    if( alData[ialign].tlen > 0 )
    {
      redPair[ipair].tlen = alData[ialign].tlen;
    }
  }

  if( ( ipair + 1 ) != npair )
  {
    sprintf( _msg_, "inconsistent read pair counts %d %d", ipair + 1, npair );
    EMSG( _msg_ );
    *fstatus = -1;
    return( 0 );
  }

  /*
  ** Sort alData back to original order.
  ** Notes:
  **   o  sort conditions
  **        o  alData index
  */
  qsort( alData, nalign, sizeof( AlData ), xcmpSortAlData3 );

  /*
  ** Check for missing mates.
  */
  for( ipair = 0; ipair < npair; ++ipair )
  {
    if( redPair[ipair].nalign < 2 )
    {
      sprintf( _msg_, "mateless proper paired end read %s", redPair[ipair].read_name );
      IMSG( _msg_ );
    }
  }

  /*
  ** Sort redPair to pos5p[0] and pos5p[1].
  ** Notes:
  **   o  sort conditions
  **        o  read 1 reference index
  **        o  read 1 start
  **        o  read 2 reference index
  **        o  read 2 start
  **        o  mapping quality of two reads (reverse order)
  **        o  base quality sum (reverse order)
  */
  qsort( redPair, npair, sizeof( RedPair ), xcmpSortRedPair );

  /*
  ** Mark duplicate pairs.
  ** Note: don't mark first copy of duplicates.
  */
  jpair = 0;
  for( ipair = 1; ipair < npair; ++ipair )
  {
    if( redPair[ipair].index[0] >= 0 && redPair[jpair].index[0] >= 0 &&
        redPair[ipair].index[1] >= 0 && redPair[jpair].index[1] >= 0 &&
        redPair[ipair].refid[0] == redPair[jpair].refid[0] &&
        redPair[ipair].pos5p[0] == redPair[jpair].pos5p[0] &&
        redPair[ipair].refid[1] == redPair[jpair].refid[1] &&
        redPair[ipair].pos5p[1] == redPair[jpair].pos5p[1] )
    {
      redPair[ipair].xflag |= BF1_B01;
      ++counts[C_DUP_PAIR_READ];
    }
    else
    {
      jpair = ipair;
    }
  }

  /*
  ** Mark alData entries:
  **   o  proper pair xflag bit and mate index, as required
  **   o  duplicate read pair xflag bit, as required
  */
  for( ipair = 0; ipair < npair; ++ipair )
  {
    ialign0 = redPair[ipair].index[0];
    ialign1 = redPair[ipair].index[1];
    if( ialign0 >= 0 && ialign1 >= 0 )
    {
      alData[ialign0].xflag |= BF1_B00;
      alData[ialign0].mateindex = redPair[ipair].index[1];
      if( redPair[ipair].xflag & BF1_B01 )
      {
        alData[ialign0].xflag |= BF1_B01;
      }

      alData[ialign1].xflag |= BF1_B00;
      alData[ialign1].mateindex = redPair[ipair].index[0];
      if( redPair[ipair].xflag & BF1_B01 )
      {
        alData[ialign1].xflag |= BF1_B01;
      }
    }
  }

  /*
  ** Count paired read starts and duplicate read starts at chromosome bases.
  */
  for( ipair = 0; ipair < npair; ++ipair )
  {
    if( redPair[ipair].index[0] >= 0 && redPair[ipair].index[1] >= 0 )
    {
      ichr = redPair[ipair].refid[0];
      ibas = redPair[ipair].pos5p[0];
      chr  = &(genome->chr[ichr]);
      ++(chr->cover[0][ibas-1]);
      ++counts[C_START_PAIRED_READ];
      if( redPair[ipair].xflag & BF1_B01 )
      {
        ++(chr->cover[1][ibas-1]);
        ++counts[C_DUP_START_PAIRED_READ];
      }

      ichr = redPair[ipair].refid[1];
      ibas = redPair[ipair].pos5p[1];
      chr  = &(genome->chr[ichr]);
      ++(chr->cover[0][ibas-1]);
      ++counts[C_START_PAIRED_READ];
      if( redPair[ipair].xflag & BF1_B01 )
      {
        ++(chr->cover[1][ibas-1]);
        ++counts[C_DUP_START_PAIRED_READ];
      }
    }
  }

  free( redPair );

  /*
  ** Count unpaired read starts.
  */
  none = 0;
  for( ialign = 0; ialign < nalign; ++ialign )
  {
    if( alData[ialign].xflag & BF1_B00 )
    {
      continue;
    }

    ++none;
  }

  /*
  ** Allocate redOne.
  */
  redOne = (RedOne *)calloc( none, sizeof( RedOne ) );
  if( redOne == NULL )
  {
    EMSG( "unable to allocate memory" );
    *fstatus = -1;
    return( 0 );
  }

  /*
  ** Assign unpaired read alignments.
  */
  ione = 0;
  for( ialign = 0; ialign < nalign; ++ialign )
  {
    if( alData[ialign].xflag & BF1_B00 )
    {
      continue;
    }

    redOne[ione].index     = alData[ialign].index;
    redOne[ione].read_name = alData[ialign].read_name;
    redOne[ione].refid     = alData[ialign].refid;
    redOne[ione].pos5p     = alData[ialign].pos5p;
    redOne[ione].pos       = alData[ialign].pos;
    redOne[ione].mapq      = alData[ialign].mapq;
    redOne[ione].sbasq     = alData[ialign].sbasq;
    redOne[ione].xflag     = (uint8_t)0;

    ++ione;
  }

  if( ione != none )
  {
    EMSG( "inconsistent counts" );
    *fstatus = -1;
    return( 0 );
  }

  /*
  ** Sort redOne by 
  **   o  refid
  **   o  pos5p
  **   o  mapq (reverse order)
  **   o  base quality sum (reverse order)
  */
  qsort( redOne, none, sizeof( RedOne ), xcmpSortRedOne );

  /*
  ** Mark unpaired duplicates.
  */
  jone = 0;
  for( ione = 1; ione < none; ++ione )
  {
    if( redOne[ione].refid == redOne[jone].refid &&
        redOne[ione].pos5p == redOne[jone].pos5p )
    {
      redOne[ione].xflag |= BF1_B01;
    }
    else
    {
      jone = ione;
    }
  }

  /*
  ** Mark alData entry:
  **   o  duplicate read xflag bit, as required
  */
  for( ione = 0; ione < none; ++ione )
  {
    ialign = redOne[ione].index;
    if( redOne[ione].xflag & BF1_B01 )
    {
      alData[ialign].xflag |= BF1_B01;
    }
  }

  /*
  ** Count unpaired read starts and duplicate read starts at chromosome bases.
  */
  for( ione = 0; ione < none; ++ione )
  {
    ichr = redOne[ione].refid;
    ibas = redOne[ione].pos5p;
    chr  = &(genome->chr[ichr]);
    ++(chr->cover[2][ibas-1]);
    ++counts[C_START_UNPAIRED_READ];
    if( redOne[ione].xflag & BF1_B01 )
    {
      ++(chr->cover[3][ibas-1]);
      ++counts[C_DUP_START_UNPAIRED_READ];
    }
  }

  free( redOne );

  *fnumAlignPE = nalign;

  *fstatus = 0;

  return( 0 );
}


static int xfindStartReadSE( AlData *alData, int numAlData, Genome *genome, int *fnumAlignSE, int64_t counts[MXCOUNTER], int *fstatus )
{
  int ialign, nalign;
  int ichr, ibas;
  int ione, jone, none;

  RedOne *redOne;
  Chromosome *chr;

  nalign = numAlData;

  /*
  ** Sort AlData by paired/single end, genome position, and set index.
  ** Notes:
  **   o  conditions
  **        o  paired/single end
  **        o  refid
  **        o  pos5p
  */
  qsort( alData, nalign, sizeof( AlData ), xcmpSortAlDataSE1 );

  for( ialign = 0; ialign < nalign; ++ialign )
  {
    if( alData[ialign].flag & SF_PAIREDEND )
    {
      nalign = ialign;
      break;
    }
    alData[ialign].index = ialign;
  }

  /*
  ** Sort AlData by read name.
  ** Notes:
  **   o  conditions
  **        o  read name
  **        o  proper pair flag (reverse order)
  **        o  position in genome
  */
  qsort( alData, nalign, sizeof( AlData ), xcmpSortAlData2 );

  none = nalign;

  /*
  ** 
  */
  redOne = (RedOne *)calloc( none, sizeof( RedOne ) );
  if( redOne == NULL )
  {
    EMSG( "unable to allocate memory" );
    *fstatus = -1;
    return( 0 );
  }

  /*
  ** Assign unpaired read alignments.
  */
  ione = 0;
  for( ialign = 0; ialign < nalign; ++ialign )
  {
    redOne[ione].index     = alData[ialign].index;
    redOne[ione].read_name = alData[ialign].read_name;
    redOne[ione].refid     = alData[ialign].refid;
    redOne[ione].pos5p     = alData[ialign].pos5p;
    redOne[ione].pos       = alData[ialign].pos;
    redOne[ione].mapq      = alData[ialign].mapq;
    redOne[ione].sbasq     = alData[ialign].sbasq;
    redOne[ione].xflag     = (uint8_t)0;

    ++ione;
  }

  if( ione != none )
  {
    EMSG( "inconsistent counts" );
    *fstatus = -1;
    return( 0 );
  }

  /*
  ** Sort alData back to original order.
  ** Notes:
  **   o  sort conditions
  **        o  alData index
  */
  qsort( alData, nalign, sizeof( AlData ), xcmpSortAlData3 );

  /*
  ** Sort redOne by 
  **   o  refid
  **   o  pos5p
  **   o  mapq
  */
  qsort( redOne, none, sizeof( RedOne ), xcmpSortRedOne );

  /*
  ** Mark unpaired duplicates.
  */
  jone = 0;
  for( ione = 1; ione < none; ++ione )
  {
    if( redOne[ione].refid == redOne[jone].refid &&
        redOne[ione].pos5p == redOne[jone].pos5p )
    {
      redOne[ione].xflag |= BF1_B01;
    }
    else
    {
      jone = ione;
    }
  }

  /*
  ** Mark alData entry:
  **   o  duplicate read xflag bit, as required
  */
  for( ione = 0; ione < none; ++ione )
  {
    ialign = redOne[ione].index;
    if( redOne[ione].xflag & BF1_B01 )
    {
      alData[ialign].xflag |= BF1_B01;
    }
  }

  /*
  ** Count unpaired read starts and duplicate read starts at chromosome bases.
  */
  for( ione = 0; ione < none; ++ione )
  {
    ichr = redOne[ione].refid;
    ibas = redOne[ione].pos5p;
    chr  = &(genome->chr[ichr]);
    ++(chr->cover[2][ibas-1]);
    ++counts[C_START_UNPAIRED_READ];
    if( redOne[ione].xflag & BF1_B01 )
    {
      ++(chr->cover[3][ibas-1]);
      ++counts[C_DUP_START_UNPAIRED_READ];
    }
  }

  free( redOne );

  *fnumAlignSE = nalign;

  *fstatus = 0;

  return( 0 );
}


static int xdumpStartRead( char *nameRoot, Genome *genome, int *fstatus )
{
  int ichr;
  int ibas;

  char ofn[8192];

  FILE *ofp;

  for( ichr = 0; ichr < genome->numChr; ++ichr )
  {
    sprintf( ofn, "%s%s.txt", nameRoot, genome->chr[ichr].name );
    ofp = fopen( ofn, "w+" );
    if( ofp == NULL )
    {
      sprintf( _msg_, "unable to open file %s", ofn );
      EMSG( _msg_ );
      *fstatus = -1;
      return( 0 );
    }

    for( ibas = 0; ibas < genome->chr[ichr].len; ++ibas )
    {
      fprintf( ofp, "%d %ld %ld %u %d  %d %d %d %d\n",
               ibas + 1,
               genome->chr[ichr].cover[0][ibas] + genome->chr[ichr].cover[2][ibas],
               genome->chr[ichr].cover[1][ibas] + genome->chr[ichr].cover[3][ibas],
               (unsigned int)genome->chr[ichr].flag[ibas],
               ( (unsigned int)( genome->chr[ichr].flag[ibas] & BF2_B04 ) ) ? 1 : 0,
               ( (unsigned int)( genome->chr[ichr].flag[ibas] & BF2_B08 ) ) ? 1 : 0,
               ( (unsigned int)( genome->chr[ichr].flag[ibas] & BF2_B09 ) ) ? 1 : 0,
               ( (unsigned int)( genome->chr[ichr].flag[ibas] & BF2_B10 ) ) ? 1 : 0,
               ( (unsigned int)( genome->chr[ichr].flag[ibas] & BF2_B11 ) ) ? 1 : 0 );
    }

    fclose( ofp );
  }

  *fstatus = 0;

  return( 0 );
}


static int xtestUniqMap( BamAlign *align, int *funiqMapFlag, int *fstatus )
{
  int ifld, nfld;
  int uniqMapFlag;

  BamOptField *field;


  /*
  ** Check that the NH optional field is '1'.
  */
  uniqMapFlag = -1;  
  nfld  = align->numOptField;
  field = align->field;
  for( ifld = 0; ifld < nfld; ++ifld )
  {
    if( field[ifld].tag[0] == 'N' && field[ifld].tag[1] == 'H' )
    {
      if( *((uint8_t *)field[ifld].value) == 1 )
      {
        uniqMapFlag = 1;
      }
      else
      {
        uniqMapFlag = 0;
      }
      break;
    }
  }  

  if( uniqMapFlag == -1 )
  {
    EMSG( "missing NH 'optional' field" );
    *fstatus = -1;
    return( 0 );
  }

  *funiqMapFlag = uniqMapFlag;

  *fstatus = 0;

  return( 0 );
}


/*
** Set read implied start and end in genomic coordinates.
** Notes:
**   o  returns ends of the implied read in the reference;
**      that is, it extends the genomic coordinates to the
**      first and last read base when the alignment is clipped
*/
static int xbamCalcEndsRead( BamHeader *bamHeader, BamAlign *bamAlign, BamDsc *bamDsc, int numBamDsc, int clip5[2], int clip3[2], int checkSLFlag, int64_t *fbpos, int64_t *fepos, int64_t counts[MXCOUNTER], int *fstatus )
{
  int n;
  int idsc;
  int slFlag;
  int shortFlag;
  int status;

  int64_t adj;
  int64_t iadj, jadj;
  int64_t bpos, epos;

  if( bamAlign->flag & SF_UNMAPPED )
  {
    *fbpos = 0;
    *fepos = 0;

    *fstatus = 0;
    return( 0 );
  }

  /*
  ** Check for spurious SL match.
  */
  slFlag = 0;
  iadj   = 0;
  if( checkSLFlag )
  {
    xcheckSLMatch( bamAlign, bamDsc, numBamDsc, clip5, clip3, &slFlag, &iadj, counts, &status );
    IF_STATUS_ZERO( "bad status: xcheckSLMatch" );
  }

  if( !slFlag )
  {
    ++counts[C_NOT_SLMATCH];
  }
  else
  {
    ++counts[C_SLMATCH];
  }

/*
** This appears to have essentially no effect.
**
  xcheckShortMatch( bamAlign, bamDsc, numBamDsc, clip5, clip3, &shortFlag, &jadj, counts, &status );
  IF_STATUS_ZERO( "bad status: xcheckShortMatch" )

  if( !shortFlag )
  {
    ++counts[C_NOT_SHORTMATCH];
  }
  else
  {
    ++counts[C_SHORTMATCH];
  }

  iadj += jadj;
*/

  adj = 0L;
  for( idsc = 0; idsc < numBamDsc; ++idsc )
  {
    if( bamDsc[idsc].type == 2 || bamDsc[idsc].type == 4 )
    {
      /*
      ** Deletion.
      */
      adj += bamDsc[idsc].len;
    }
    else
    if( bamDsc[idsc].type == 3 )
    {
      /*
      ** Insertion.
      */
      adj -= bamDsc[idsc].len;
    }
  }

  if( bamAlign->flag & SF_REV_CMP )
  {
    bpos = bamAlign->pos - (int64_t)( clip5[0] + clip5[1] );
  }
  else
  {
    bpos = bamAlign->pos - (int64_t)( clip5[0] + clip5[1] ) + iadj;
  }
  epos = bpos + ( bamAlign->l_seq + clip5[1] + clip3[1] ) - 1 + ( adj - iadj );

  n = bamHeader->refSeq[bamAlign->refid].length;
  if( bpos < 1 ) bpos = 1;
  if( bpos > n ) bpos = n;
  if( epos < 1 ) epos = 1;
  if( epos > n ) epos = n;

  *fbpos = bpos;
  *fepos = epos;

  return( 0 );
}


/*
** Look for spurious SL match at read start.
*/
static int xcheckSLMatch( BamAlign *bamAlign, BamDsc *bamDsc, int numBamDsc, int clip5[2], int clip3[2], int *fslFlag, int64_t *fiadj, int64_t counts[MXCOUNTER], int *fstatus )
{
  int iop, nop;
  int isl;
  int ioff;
  int ibas, ibas1;
  int status;
  int nmis, nmat;

  char    cop[8];
  int32_t lop[8];

  static int       sizSeq = 0;
  static char     *seq    = NULL;
  static int       luop   = 0;
  static uint32_t *uop    = NULL;

#define BSL_A	( 1 << 0 )
#define BSL_C	( 1 << 1 )
#define BSL_G	( 1 << 2 )
#define BSL_T	( 1 << 3 )

  static int initFlag = 0;

  static uint8_t bsl[32] = {
                              BSL_G,    // 1
                              BSL_A,
                              ( BSL_A | BSL_G ),
                              ( BSL_C | BSL_T ),
                              ( BSL_A | BSL_C | BSL_T ),
                              ( BSL_A | BSL_C | BSL_T ), // 6
                              ( BSL_A | BSL_G | BSL_T ),
                              ( BSL_A | BSL_T ),
                              ( BSL_A | BSL_G | BSL_T ),
                              ( BSL_A | BSL_C | BSL_G | BSL_T ),
                              ( BSL_A | BSL_C | BSL_G ), // 11
                              ( BSL_A | BSL_C ),
                              ( BSL_A | BSL_C ),
                              ( BSL_A | BSL_C | BSL_T ),
                              ( BSL_A | BSL_C | BSL_T ),
                              ( BSL_A | BSL_T ),  // 16
                              ( BSL_A | BSL_T ),
                              BSL_T,
                              BSL_T,
                              ( BSL_G | BSL_T )
                            };

  static int bbas[128];
  static int koff[16];
  static int mxkoff = 15;
  static int mxmis  = 2;
  static int mnmat1  = 2;
  static int mnmat2  = 10;
  
  static char *slSeq[14] = {
                             "GGTTTAATTACCCAAGTTTGAG",
                             "GGTTTTAACCCAGTTACTCAAG",
                             "GGTTTTAACCCAGTTAACCAAG",
                             "GGTTTTAACCCAGTTTAACCAAG",
                             "GGTTTTAACCCAGTTACCAAG",
                             "GGTTTAAAACCCAGTTACCAAG",
                             "GGTTTTAACCCAGTTAATTGAG",
                             "GGTTTTTACCCAGTTAACCAAG",
                             "GGTTTATACCCAGTTAACCAAG",
                             "GGTTTTAACCCAAGTTAACCAAG",
                             "GGTTTTAACCAGTTAACTAAG",
                             "GGTTTTAACCCATATAACCAAG",
                             "GGTTTTAACCCAGTTAACTAAG",
                             "GGTTTAAACCCAGTTAACAAG"
                           };

  static char *rslSeq[14];

                            
  if( !initFlag )
  {
    int i, n;
    char string[1024];

    memset( bbas, 0, 128 * sizeof( int ) );

    bbas[(int)'a'] = BSL_A; bbas[(int)'A'] = BSL_A;
    bbas[(int)'c'] = BSL_C; bbas[(int)'C'] = BSL_C;
    bbas[(int)'g'] = BSL_G; bbas[(int)'G'] = BSL_G;
    bbas[(int)'t'] = BSL_T; bbas[(int)'T'] = BSL_T;

    koff[0] = 0;
    koff[1] = -1;
    koff[2] = 1;
    koff[3] = -2;
    koff[4] = 2;
    koff[5] = -3;
    koff[6] = 3;
    koff[7] = -4;
    koff[8] = 4;
    koff[9] = -5;
    koff[10] = 5;
    koff[11] = -6;
    koff[12] = 6;
    koff[13] = -7;
    koff[14] = 7;

    for( isl = 0; isl < 14; ++isl )
    {
      n = strlen( slSeq[isl] );
      for( i = 0; i < n; ++i )
      {
        string[i] = slSeq[isl][n-i-1];
      }
      string[n] = '\0';
      rslSeq[isl] = mstrcpy( string );
    }

    initFlag = 1;
  }

  nop = bamAlign->n_cigar_op;

  /*
  ** Need at least
  **   o  (short) match
  **   o  (too long) intron
  **   o  match
  */
  if( nop < 3 )
  {
    *fslFlag = 0;
    *fiadj   = (int64_t)0;
    *fstatus = 0;
    return( 0 );
  }

  /*
  ** Get CIGAR ops into uop in the correct order.
  */
  if( nop > luop )
  {
    luop += nop + 128;
    uop = (uint32_t *)realloc( uop, luop * sizeof( uint32_t ) );
    if( uop == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fslFlag = 0;
      *fiadj   = (int64_t)0;
      *fstatus = -1;
      return( 0 );
    }
  }

  if( bamAlign->flag & SF_REV_CMP )
  {
    for( iop = 0; iop < nop; ++iop )
    {
      uop[iop] = bamAlign->cigar[nop-iop-1];
    }
  }
  else
  {
    memcpy( uop, bamAlign->cigar, nop * sizeof( uint32_t ) );
  }

  /*
  ** Check for possible soft clipping preceding possible match/intron/match.
  */
  bamCigarOp( uop[0], &(cop[0]), &(lop[0]), &status );
  IF_STATUS_ZERO( "bad status: bamCigarOp" );

  if( cop[0] == 'S' && nop < 4 )
  {
    *fslFlag = 0;
    *fiadj   = (int64_t)0;
    *fstatus = 0;
    return( 0 );
  }

  if( cop[0] != 'S' )
  {
    bamCigarOp( uop[1], &(cop[1]), &(lop[1]), &status );
    IF_STATUS_ZERO( "bad status: bamCigarOp" );
    ibas1 = lop[0] - 1;
  }
  else
  {
    lop[3] = lop[0];
    bamCigarOp( uop[1], &(cop[0]), &(lop[0]), &status );
    IF_STATUS_ZERO( "bad status: bamCigarOp" );

    bamCigarOp( uop[2], &(cop[1]), &(lop[1]), &status );
    IF_STATUS_ZERO( "bad status: bamCigarOp" );
    ibas1 = lop[0] + lop[3] - 1;
  }

  /*
  ** Skip if not short match followed by an intron.
  */
  if( cop[0] != 'M' || cop[1] != 'N' || lop[0] < mnmat1 || lop[0] > 20 )
  {
    *fslFlag = 0;
    *fiadj   = (int64_t)0;
    *fstatus = 0;
    return( 0 );
  }

  /*
  ** Unpack read sequence and put in original orientation.
  */
  bamUnpackSeq( bamAlign, &seq, &sizSeq, &status );
  IF_STATUS_ZERO( "bad status: bamUnpackSeq" );

  if( bamAlign->flag & SF_REV_CMP )
  {
    rcseq2( seq, &status );
    IF_STATUS_ZERO( "bad status rcseq2" );
  }

  /*
  ** Long intron?
  ** If yes, check all SL sequences and allow two mismatches.
  */
  if( lop[0] >= ( mnmat2 - mxmis ) && lop[1] >= 3000 )
  {
    /*
    ** Allow 2 mismatches 
    */
    for( ioff = 0; ioff < mxkoff; ++ioff )
    {
      nmat = 0;
      nmis = 0;
      for( ibas = 0; ibas < 20 && ibas1 - ibas + koff[ioff] >= 0; ++ibas )
      {
        if( !( bbas[(int)seq[ibas1-ibas+koff[ioff]]] & bsl[ibas] ) )
        {
          ++nmis;
          if( nmis > mxmis )
          {
            break;
          }
        }
        else
        {
          ++nmat;
        }
      }

      if( nmat >= mnmat2 )
      {
        if( nmis <= mxmis )
        {
          for( isl = 0; isl < 14; ++isl )
          {
            nmat = 0;
            nmis = 0;
            for( ibas = 0; ibas < strlen( rslSeq[isl] ) && ibas1 - ibas + koff[ioff] >= 0; ++ibas )
            {
              if( !( seq[ibas1-ibas+koff[ioff]] == rslSeq[isl][ibas] ) )
              {
                ++nmis;
                if( nmis > mxmis )
                {
                  break;
                }
              }
              else
              {
                ++nmat;
              }
            }

            if( nmat >= mnmat2 && nmis <= mxmis )
            {
#ifdef SL_REP
fprintf( stdout, "SL: A 1  lops: %d %d  isl: %d  koff: %d  nmat: %d  nmis: %d  bpos: %d (%d)  seq: %s  isl: %d cl: %d %d  rn: %s\n",
         lop[0], lop[1], isl, koff[ioff], nmat, nmis, bamAlign->pos - clip5[0] - clip5[1] + lop[1], bamAlign->pos, seq, isl, clip5[0], clip5[1], bamAlign->read_name );
#endif
              ++counts[C_LONGINTRON_SLMATCH];
              *fslFlag = 1;
              *fiadj   = (int64_t)lop[1];
              *fstatus = 0;
              return( 0 );
            }
          }
        }
      }
    }
  }
  else
  if( lop[1] >= 20 )
  {
    /*
    ** Allow no mismatches 
    */
    for( ioff = 0; ioff < mxkoff; ++ioff )
    {
      for( isl = 0; isl < 2; ++isl )
      {
        nmat = 0;
        nmis = 0;
        for( ibas = 0; ibas < strlen( rslSeq[isl] ) && ibas1 - ibas + koff[ioff] >= 0; ++ibas )
        {
          if( !( seq[ibas1-ibas+koff[ioff]] == rslSeq[isl][ibas] ) )
          {
            ++nmis;
            break;
          }
          else
          {
            ++nmat;
          }
        }

        if( nmat >= mnmat1 && nmis == 0 )
        {
#ifdef SL_REP
fprintf( stdout, "SL: B 1  lops: %d %d  isl: %d  koff: %d  nmat: %d  nmis: %d  bpos: %d (%d)  seq: %s  isl: %d cl: %d %d  rn: %s\n",
         lop[0], lop[1], isl, koff[ioff], nmat, nmis, bamAlign->pos - clip5[0] - clip5[1] + lop[1], bamAlign->pos, seq, isl, clip5[0], clip5[1], bamAlign->read_name );
#endif
          ++counts[C_SHORTINTRON_SLMATCH];
          *fslFlag = 1;
          *fiadj   = (int64_t)lop[1];
          *fstatus = 0;
          return( 0 );
        }
      }
    }
  }

  *fslFlag = 0;
  *fiadj   = (int64_t)0;

  *fstatus = 0;

  return( 0 );
}


/*
** Check for short matches at read start and set to soft-clipping.
*/
static int xcheckShortMatch( BamAlign *bamAlign, BamDsc *bamDsc, int numBamDsc, int clip5[2], int clip3[2], int *fshortFlag, int64_t *fiadj, int64_t counts[MXCOUNTER], int *fstatus )
{
  int iop, nop;
  int mnmat1;
  int status;

  char    cop[8];
  int32_t lop[8];

  static int       luop   = 0;
  static uint32_t *uop    = NULL;

  mnmat1 = 5;

  nop = bamAlign->n_cigar_op;

  /*
  ** Need at least
  **   o  (short) match
  **   o  (too long) intron
  **   o  match
  */
  if( nop < 3 )
  {
    *fshortFlag = 0;
    *fiadj   = (int64_t)0;
    *fstatus = 0;
    return( 0 );
  }

  /*
  ** Get CIGAR ops into uop in the correct order.
  */
  if( nop > luop )
  {
    luop += nop + 128;
    uop = (uint32_t *)realloc( uop, luop * sizeof( uint32_t ) );
    if( uop == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fshortFlag = 0;
      *fiadj   = (int64_t)0;
      *fstatus = -1;
      return( 0 );
    }
  }

  if( bamAlign->flag & SF_REV_CMP )
  {
    for( iop = 0; iop < nop; ++iop )
    {
      uop[iop] = bamAlign->cigar[nop-iop-1];
    }
  }
  else
  {
    memcpy( uop, bamAlign->cigar, nop * sizeof( uint32_t ) );
  }

  /*
  ** Check for possible soft clipping preceding possible match/intron/match.
  */
  bamCigarOp( uop[0], &(cop[0]), &(lop[0]), &status );
  IF_STATUS_ZERO( "bad status: bamCigarOp" );

  if( cop[0] != 'S' )
  {
    bamCigarOp( uop[1], &(cop[1]), &(lop[1]), &status );
    IF_STATUS_ZERO( "bad status: bamCigarOp" );
  }
  else
  {
    lop[3] = lop[0];
    bamCigarOp( uop[1], &(cop[0]), &(lop[0]), &status );
    IF_STATUS_ZERO( "bad status: bamCigarOp" );

    bamCigarOp( uop[2], &(cop[1]), &(lop[1]), &status );
    IF_STATUS_ZERO( "bad status: bamCigarOp" );
  }

  /*
  ** Skip if not short match followed by an intron.
  */
  if( cop[0] != 'M' || cop[1] != 'N' || lop[0] >= mnmat1 )
  {
    *fshortFlag = 0;
    *fiadj      = (int64_t)0;
    *fstatus    = 0;
    return( 0 );
  }

  *fshortFlag = 1;
  *fiadj      = (int64_t)lop[1];

  *fstatus = 0;

  return( 0 );
}


/*
** Identify and mark high coverage regions using D-segments.
*/
static int xmarkCovReg( RegPar *regPar, Genome *genome, int rsFlag, int64_t counts[MXCOUNTER], int *fstatus )
{
  int ichr, nchr;
  int ibas, jbas, nbas;
  int sfalse;
  int strue;
  int ival;
  int lmin;
  int dscore;
  int index;

  int threshold;
  int score, sx;
  int i0, i1;
  int nstart;

  double rsc;
  double rmin;

  uint16_t bflag;

  Chromosome *chr;

  sfalse    = regPar->sfalse;
  strue     = regPar->strue;
  dscore    = regPar->dscore;
  lmin      = regPar->lmin;
  rmin      = regPar->rmin;
  bflag     = regPar->bflag;
  threshold = 1;

  if( rsFlag == RS_SINGLE )
  {
    index = 2;
  }
  else
  if( rsFlag == RS_PAIRED )
  {
    index = 0;
  }
  else
  {
    EMSG( "unrecognized read type value" );
    *fstatus = -1;
    return( 0 );
  }

  nchr = genome->numChr;

  for( ichr = 0; ichr < nchr; ++ichr )
  {
    chr    = &(genome->chr[ichr]);
    nbas   = chr->len;

    i0     = 0;
    i1     = 0;
    sx     = 0;
    score  = 0;
    for( ibas = 0; ibas < nbas; ++ibas )
    {
      /*
      ** Score data using read start counts.
      */
      ival    = chr->cover[index][ibas];
      score  += ( ival < threshold ) ? sfalse : strue;

      if( score >= sx )
      {
        sx = score;
        i1 = ibas;
      }

      if( score <= 0 || score <= sx + dscore || ibas == nbas - 1 )
      {
        nstart = 0;
        for( jbas = i0; jbas <= i1; ++jbas )
        {
          nstart += chr->cover[index][jbas];
        }
        rsc = (double)nstart / (double)( i1 - i0 + 1 );
        if( sx > strue && rsc >= rmin && i1 - i0 + 1 >= lmin )
        {
          for( jbas = i0; jbas <= i1; ++jbas )
          {
            chr->flag[jbas] |= bflag;
          }
        }

        score  = 0;
        sx     = 0;
        i0     = ibas + 1;
        i1     = ibas + 1;
      }
    }
  }

  *fstatus = 0;

  return( 0 );
}


static int xmarkCovRegSet( Genome *genome, int rsFlag, int64_t counts[MXCOUNTER], int *fstatus )
{
  int status;

  RegPar regPar;

  regPar.sfalse = -1;
  regPar.lmin   = 10;

  /*
  ** Mark regions with read starts per base >= .01
  */
  regPar.strue  = 500;
  regPar.dscore = -500;
  regPar.rmin   = 0.01;
  regPar.bflag  = BF2_B08;
  xmarkCovReg( &regPar, genome, rsFlag, counts, &status );
  IF_STATUS_ZERO( "bad status: xmarkCovReg" );

  /*
  ** Mark regions with read starts per base >= .1
  */
  regPar.strue  = 50;
  regPar.dscore = -50;
  regPar.rmin   = 0.1;
  regPar.bflag  = BF2_B09;
  xmarkCovReg( &regPar, genome, rsFlag, counts, &status );
  IF_STATUS_ZERO( "bad status: xmarkCovReg" );

  /*
  ** Mark regions with read starts per base >= 1
  */
  regPar.strue  = 5;
  regPar.dscore = -5;
  regPar.rmin   = 1.0;
  regPar.bflag  = BF2_B10;
  xmarkCovReg( &regPar, genome, rsFlag, counts, &status );
  IF_STATUS_ZERO( "bad status: xmarkCovReg" );

  /*
  ** Mark regions with read starts per base >= 10
  */
  regPar.strue  = 2;
  regPar.dscore = -2;
  regPar.rmin   = 10.0;
  regPar.bflag  = BF2_B11;
  xmarkCovReg( &regPar, genome, rsFlag, counts, &status );
  IF_STATUS_ZERO( "bad status: xmarkCovReg" );

  *fstatus = 0;

  return( 0 );
}


static int xcountCovRegion( Genome *genome, uint16_t bit_true, uint16_t bit_false, int minSpc, int minLen, int *fnstart, int *fndup, int64_t counts[MXCOUNTER], int *fstatus )
{
  int ichr, nchr;
  int ibas, jbas, nbas;
  int ibas0, ibas1;
  int nh0, nh1;
  int ndup[2];
  int nstart[2];

  Chromosome *chr;

  nchr = genome->numChr;

  ndup[0]   = 0;
  ndup[1]   = 0;
  nstart[0] = 0;
  nstart[1] = 0;

  for( ichr = 0; ichr < nchr; ++ichr )
  {
    chr  = &(genome->chr[ichr]);
    nbas = chr->len;

    for( ibas = 0; ibas < nbas; ++ibas )
    {
      /*
      ** Find region start.
      */
      if( ( chr->flag[ibas] & bit_true ) && !( chr->flag[ibas] & bit_false ) )
      {
        /*
        ** Find region end.
        */
        ibas0 = ibas;
        ibas1 = ibas;
        for( ; ibas < nbas; ++ibas )
        {
          if( !( chr->flag[ibas] & bit_true ) || ( chr->flag[ibas] & bit_false ) )
          {
            ibas1 = ibas - 1;
            break;
          }
        }

        if( ibas1 - ibas0 + 1 < minLen )
        {
          continue;
        }

        /*
        ** Count region read starts and duplicates, allowing
        ** for 'edge effects'.
        */
        nh0 = 0;
        for( jbas = ibas0 - 1; jbas >= 0; --jbas )
        {
          if( chr->flag[jbas] & bit_false || nh0 == minSpc )
          {
            break;
          }
          ++nh0;
        }
        ibas0 += minSpc - nh0;

        nh1 = 0;
        for( jbas = ibas1 + 1; jbas < nbas; ++jbas )
        {
          if( chr->flag[jbas] & bit_false || nh1 == minSpc )
          {
            break;
          }
          ++nh1;
        }
        ibas1 -= minSpc - nh1;

        if( ibas1 - ibas0 + 1 >= minLen )
        {
          for( jbas = ibas0; jbas <= ibas1; ++jbas )
          {
            nstart[0] += chr->cover[0][jbas];   /* PE read start */
            ndup[0]   += chr->cover[1][jbas];   /* PE duplicate */
            nstart[1] += chr->cover[2][jbas];   /* SE read start */
            ndup[1]   += chr->cover[3][jbas];   /* SE duplicate */
          }
        } 
      } /* if chr->flag[ibas] ... */
    } /* for ibas */
  } /* for ichr */

  fndup[0]   = ndup[0];
  fndup[1]   = ndup[1];
  fnstart[0] = nstart[0];
  fnstart[1] = nstart[1];

  *fstatus = 0;

  return( 0 );
}


static int xcountLoCovRegion( Genome *genome, uint16_t bit_false, int minSpc, int minLen, int *fnstart, int *fndup, int32_t rsFlag, int64_t hist[MXHIST], int64_t counts[MXCOUNTER], int *fstatus )
{
  int nd;
  int ichr, nchr;
  int ibas, jbas, nbas;
  int ibas0, ibas1;
  int nh0, nh1;
  int ndup[2];
  int nstart[2];
  int irs;

  Chromosome *chr;

  nchr = genome->numChr;

  ndup[0]   = 0;
  ndup[1]   = 0;
  nstart[0] = 0;
  nstart[1] = 0;

  if( rsFlag == RS_SINGLE )
  {
    irs = 3;
  }
  else
  if( rsFlag == RS_PAIRED )
  {
    irs = 1;
  }
  else
  {
    EMSG( "unrecognized rsFlag value" );
    *fstatus = -1;
    return( 0 );
  }

  for( ichr = 0; ichr < nchr; ++ichr )
  {
    chr  = &(genome->chr[ichr]);
    nbas = chr->len;

    for( ibas = 0; ibas < nbas; ++ibas )
    {
      /*
      ** Find region start.
      */
      if( !( chr->flag[ibas] & bit_false ) )
      {
        /*
        ** Find region end.
        */
        ibas0 = ibas;
        ibas1 = ibas;
        for( ; ibas < nbas; ++ibas )
        {
          if( chr->flag[ibas] & bit_false )
          {
            ibas1 = ibas - 1;
            break;
          }
        }

        if( ibas1 - ibas0 + 1 < minLen )
        {
          continue;
        }

        /*
        ** Count region read starts and duplicates, allowing
        ** for 'edge effects'.
        */
        nh0 = 0;
        for( jbas = ibas0 - 1; jbas >= 0; --jbas )
        {
          if( chr->flag[jbas] & bit_false || nh0 == minSpc )
          {
            break;
          }
          ++nh0;
        }
        ibas0 += minSpc - nh0;

        nh1 = 0;
        for( jbas = ibas1 + 1; jbas < nbas; ++jbas )
        {
          if( chr->flag[jbas] & bit_false || nh1 == minSpc )
          {
            break;
          }
          ++nh1;
        }
        ibas1 -= minSpc - nh1;

        if( ibas1 - ibas0 + 1 >= minLen )
        {
          for( jbas = ibas0; jbas <= ibas1; ++jbas )
          {
            nstart[0] += chr->cover[0][jbas];
            ndup[0]   += chr->cover[1][jbas];
            nstart[1] += chr->cover[2][jbas];
            ndup[1]   += chr->cover[3][jbas];
            nd         = chr->cover[irs][jbas];
            if( nd >= MXHIST ) nd = MXHIST - 1;
            ++hist[nd];
          }
        } 
      } /* if chr->flag[ibas] ... */
    } /* for ibas */
  } /* for ichr */

  fndup[0]   = ndup[0];
  fndup[1]   = ndup[1];
  fnstart[0] = nstart[0];
  fnstart[1] = nstart[1];

  *fstatus = 0;

  return( 0 );
}


static int xestimatePCRDupRate( Bundle *bundle, BamHeader *bamHeader, Genome *genome, int rsFlag, double *fpcrDupRate, int64_t hist[MXHIST], int64_t counts[MXCOUNTER], int *fstatus )
{
  int minSpc;
  int minLen;
  int ndup[8][2];
  int nstart[8][2];
  int status;

  uint16_t bit_true;
  uint16_t bit_false;

  char *fptr;

  minSpc = 50;
  minLen = 200;

  /*
  ** Calculate detected duplication rate for regions with >= .01 and < .1 read start rate.
  */
  bit_true  = BF2_B08;
  bit_false = BF2_B09 | BF2_B10 | BF2_B11;
  xcountCovRegion( genome, bit_true, bit_false, minSpc, minLen, nstart[0], ndup[0], counts, &status );
  IF_STATUS_ZERO( "bad status: xcountCovRegion" );

  /*
  ** Calculate detected duplication rate for regions with >= .1 and < 1 read start rate.
  */
  bit_true  = BF2_B09;
  bit_false = BF2_B10 | BF2_B11;
  xcountCovRegion( genome, bit_true, bit_false, minSpc, minLen, nstart[1], ndup[1], counts, &status );
  IF_STATUS_ZERO( "bad status: xcountCovRegion" );

  /*
  ** Calculate detected duplication rate for regions with >= 1.0 and < 10 read start rate.
  */
  bit_true  = BF2_B10;
  bit_false = BF2_B11;
  xcountCovRegion( genome, bit_true, bit_false, minSpc, minLen, nstart[2], ndup[2], counts, &status );
  IF_STATUS_ZERO( "bad status: xcountCovRegion" );

  /*
  ** Calculate detected duplication rate for regions with >= 10.0 read start rate.
  */
  bit_true  = BF2_B11;
  bit_false = (uint16_t)0;
  xcountCovRegion( genome, bit_true, bit_false, minSpc, minLen, nstart[3], ndup[3], counts, &status );
  IF_STATUS_ZERO( "bad status: xcountCovRegion" );

  if( rsFlag == RS_SINGLE )
  {
    bit_false = SE_BIT_TEST;
  }
  else
  if( rsFlag == RS_PAIRED )
  {
    bit_false = PE_BIT_TEST;
  }
  else
  {
    EMSG( "unrecognized read set value" );
    *fstatus = -1;
    return( 0 );
  }
  xcountLoCovRegion( genome, bit_false, minSpc, minLen, nstart[4], ndup[4], rsFlag, hist, counts, &status );
  IF_STATUS_ZERO( "bad status: xcountLoCovRegion" );

  fptr = strrchr( bundle->nameBamFile, '/' );
  fptr = fptr != NULL ? fptr + 1 : bundle->nameBamFile;

  /*
  **
  */
  fprintf( stderr, "DURAT: %.4f %.4f %.4f %.4f  %.4f %.4f %.4f %.4f  %d %d %d %d  %d %d %d %d  %d %d %d %d  %d %d %d %d  %.4f %.4f  %d %d  %d %d %s\n",
           nstart[0][0] > 0 ? (double)ndup[0][0] / (double)nstart[0][0] : -1,
           nstart[1][0] > 0 ? (double)ndup[1][0] / (double)nstart[1][0] : -1,
           nstart[2][0] > 0 ? (double)ndup[2][0] / (double)nstart[2][0] : -1,
           nstart[3][0] > 0 ? (double)ndup[3][0] / (double)nstart[3][0] : -1,

           nstart[0][1] > 0 ? (double)ndup[0][1] / (double)nstart[0][1] : -1,
           nstart[1][1] > 0 ? (double)ndup[1][1] / (double)nstart[1][1] : -1,
           nstart[2][1] > 0 ? (double)ndup[2][1] / (double)nstart[2][1] : -1,
           nstart[3][1] > 0 ? (double)ndup[3][1] / (double)nstart[3][1] : -1,

           nstart[0][0], nstart[1][0], nstart[2][0], nstart[3][0],
           ndup[0][0],   ndup[1][0],   ndup[2][0],   ndup[3][0],

           nstart[0][1], nstart[1][1], nstart[2][1], nstart[3][1],
           ndup[0][1],   ndup[1][1],   ndup[2][1],   ndup[3][1],

           nstart[4][0] > 0 ? (double)ndup[4][0] / (double)nstart[4][0] : -1,
           nstart[4][1] > 0 ? (double)ndup[4][1] / (double)nstart[4][1] : -1,

           nstart[4][0], nstart[4][1],
           ndup[4][0],   ndup[4][1],

           fptr );

  if( rsFlag == RS_SINGLE )
  {
    *fpcrDupRate = nstart[4][1] > 0 ? (double)ndup[4][1] / (double)nstart[4][1] : -1;
  }
  else
  if( rsFlag == RS_PAIRED )
  {
    *fpcrDupRate = nstart[4][0] > 0 ? (double)ndup[4][0] / (double)nstart[4][0] : -1;
  }
  else
  {
    EMSG( "unrecognized read set value" );
    *fstatus = -1;
    return( 0 );
  }

  *fstatus = 0;

  return( 0 );
}


static int xinitXorShift1024( void )
{
  int i;
  uint64_t seed[16];

  /*
  ** 8 random bytes from http://www.random.org/
  */
  seed[0]  = 0x3426822a9273d29c;
  seed[1]  = 0x7ae1e06e43ea4419;
  seed[2]  = 0x11f20b1635dc29be;
  seed[3]  = 0x6f11870aed507e92;
  seed[4]  = 0xadb155fcc94d4d71;
  seed[5]  = 0x1d29120c00897afa;
  seed[6]  = 0x2f1810ea691c6628;
  seed[7]  = 0x26a6d78fa6d97f1f;
  seed[8]  = 0x2f888f10354265d7;
  seed[9]  = 0xa649e3e7c91a1b99;
  seed[10] = 0x9cb155cdfc6ab1c3;
  seed[11] = 0x663056cc26c262dc;
  seed[12] = 0xe41182067eed1d39;
  seed[13] = 0xc791834e4dfdf7e6;
  seed[14] = 0xfa3e14468a029959;
  seed[15] = 0x2c47105282429488;

  for( i = 0; i < 16; ++i )
  {
    xrs[i] = seed[i];
  }

  return( 0 );
}


#ifdef FASTER_PROC_REG

static int xprocRegAlign( Genome *genome, HiCovReg *hiCovReg, AlData *alData, int ialign0, int ialign1, int nstart, int ndup, double pcrDupRate, int *fnrem, int *fstatus )
{
  int ireg;
  int nrem;
  int ialign;
  int idup, jdup, kdup, mdup, tdup;
  int status;

  DupSet tdupSet;

  static int mxDupSet = 0;
  static DupSet *dupSet = NULL;

  static int mxa2d = 0;
  static int *a2d  = NULL;

  ireg = alData[ialign0].hcrindex;

  /*
  ** Consistency tests.
  */
  xcheckRegAlign( genome, alData, ialign0, ialign1, hiCovReg, nstart, ndup, &status );
  IF_STATUS_ZERO( "bad status: xcheckRegAlign" );

  /*
  ** Re-allocate duplicate array, if necessary.
  */
  if( ndup > mxDupSet )
  {
    mxDupSet += ndup + 1024;
    dupSet = (DupSet *)realloc( dupSet, mxDupSet * sizeof( DupSet ) );
    if( dupSet == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
  }

  /*
  ** Calculate number of duplicates to remove.
  */
  mdup = (int)( pcrDupRate * (double)nstart + 0.5 );

  /*
  ** Do we remove all duplicates?
  */
  if( mdup >= ndup )
  {
    /*
    ** Mark all duplicates (and their mates) for removal and
    ** return.
    */
    nrem = 0;
    for( ialign = ialign0; ialign <= ialign1; ++ialign )
    {
      if( alData[ialign].xflag & BF1_B01 && !( alData[ialign].xflag & BF1_B03 ) )
      {
        alData[ialign].xflag |= BF1_B03;
        ++nrem;
      }

      if( alData[ialign].xflag & BF1_B00 && alData[alData[ialign].mateindex].xflag & BF1_B01 && !( alData[alData[ialign].mateindex].xflag & BF1_B03 ) )
      {
        alData[alData[ialign].mateindex].xflag |= BF1_B03;
        ++nrem;
      }
    }

    *fnrem   += nrem;
    *fstatus  = 0;
    return( 0 );
  }

  if( ialign1 - ialign0 + 1 > mxa2d )
  {
    mxa2d += ialign1 - ialign0 + 1;
    a2d = (int *)realloc( a2d, mxa2d * sizeof( int ) );
  }

  /*
  ** Copy duplicate indices.
  ** Note: watch and account for alignments
  **       removed previously as mates.
  */
  idup = 0;
  for( ialign = ialign0; ialign <= ialign1; ++ialign )
  {
    if( alData[ialign].xflag & BF1_B03 )
    {
      --ndup;
      --mdup;
      continue;
    }

    if( alData[ialign].xflag & BF1_B01 )
    {
      dupSet[idup].index    = ialign;
      dupSet[idup].flag     = 0;
      dupSet[idup].mateidup = -1;
      a2d[ialign-ialign0]   = idup;
      ++idup;
      continue;
    }

    a2d[ialign-ialign0] = -1;
  }

  if( idup == 0 )
  {
    *fstatus = 0;
    return( 0 );
  }

  if( idup != ndup )
  {
    sprintf( _msg_, "inconsistent counts for idup/ndup: %d != %d", idup, ndup );
    EMSG( _msg_ );
    *fstatus = -1;
    return( 0 );
  }

  /*
  ** Set mate location in dupSet.
  */
  for( idup = 0; idup < ndup; ++idup )
  {
    if( alData[dupSet[idup].index].xflag & BF1_B00 && alData[alData[dupSet[idup].index].mateindex].hcrindex == ireg )
    {
      if( alData[dupSet[idup].index].mateindex < ialign0 || alData[dupSet[idup].index].mateindex > ialign1 )
      {
         EMSG( "unexpected condition" );
         *fstatus = -1;
         return( 0 );
      }
      dupSet[idup].mateidup = a2d[alData[dupSet[idup].index].mateindex-ialign0];
    }
  }

#ifdef VERIFY_LONG
  /*
  ** Consistency test: check that mates in the region are in dupSet.
  */
  xcheckDupSet( alData, dupSet, ndup, hiCovReg, ireg, &status );
  IF_STATUS_ZERO( "bad status: xcheckDupSet" );
#endif

  /*
  ** Shuffle randomly dupSet elements.
  */
  kdup = ndup;
  while( kdup > 0 )
  {
    jdup = xorShift1024() % kdup;

    /*
    ** Adjust mateidups... we exchange the kdup-1 and jdup dupSet[] elements
    */
    tdup = dupSet[jdup].mateidup;
    if( dupSet[kdup-1].mateidup >= 0 ) dupSet[dupSet[kdup-1].mateidup].mateidup = jdup;
    if( tdup >= 0 ) dupSet[tdup].mateidup = kdup - 1;

    /* 
    ** Now exchange the dupSet[] elements.
    */
    memcpy( &tdupSet, &(dupSet[kdup-1]), sizeof( DupSet ) );
    memcpy( &(dupSet[kdup-1]), &(dupSet[jdup]), sizeof( DupSet ) );
    memcpy( &(dupSet[jdup]), &tdupSet, sizeof( DupSet ) );
    --kdup;
  }

  if( mdup > 1000000 )
  {
    fprintf( stderr, "(ndup: %d  mdup: %d)", ndup, mdup );
  }

  /*
  ** The dupSet elements are shuffled randomly so select
  ** from idup = 0 to idup = mdup-1, in order.
  */
  nrem = 0;
  for( idup = 0; idup < mdup; ++idup )
  {
    /*
    ** Check for previously removed duplicates/mates.
    */
    if( alData[dupSet[idup].index].xflag & BF1_B03 ||
        ( alData[dupSet[idup].index].xflag & BF1_B00 && alData[alData[dupSet[idup].index].mateindex].xflag & BF1_B03 ) )
    {
      EMSG( "unexpected condition" );
fprintf( stderr, "ireg: %d  rn: %s  index: %d  xflag: %u  mateindex: %d\n",
ireg, alData[dupSet[idup].index].read_name, dupSet[idup].index, alData[dupSet[idup].index].xflag, alData[dupSet[idup].index].mateindex );
fprintf( stderr, "idup: %d  ndup: %d  mdup: %d  dupSet[idup].index: %d  dupSet[idup].mateidup: %d  dupSet[dupSet[idup].mateidup].index: %d\n",
idup, ndup, mdup, dupSet[idup].index, dupSet[idup].mateidup, dupSet[dupSet[idup].mateidup].index );
      *fstatus = -1;
      return( 0 );
    }

    /*
    ** Mark alignment and mate for removal.
    */
    alData[dupSet[idup].index].xflag |= BF1_B03;
    ++nrem;

    if( alData[dupSet[idup].index].xflag & BF1_B00 && alData[alData[dupSet[idup].index].mateindex].xflag & BF1_B02 )
    {
      if( alData[alData[dupSet[idup].index].mateindex].xflag & BF1_B03 )
      {
        EMSG( "unexpected condition" );
        *fstatus = -1;
        return( 0 );
      }

      alData[alData[dupSet[idup].index].mateindex].xflag |= BF1_B03;
      ++nrem;

      /*
      ** If mate is in ireg, increment idup and replace it
      ** in dupSet, if necessary.
      */
      if( alData[alData[dupSet[idup].index].mateindex].hcrindex == ireg )
      {
/*
        for( kdup = idup + 1; kdup < mdup; ++kdup )
        {
          if( dupSet[kdup].index == alData[dupSet[idup].index].mateindex )
          {
            memcpy( &(dupSet[kdup]), &(dupSet[idup+1]), sizeof( DupSet ) );
            break;
          }
        }
*/

        if( alData[dupSet[idup].index].mateindex != dupSet[dupSet[idup].mateidup].index )
        {
          EMSG( "inconsistent dupSet[].index and dupSet[].mateidup values" );
fprintf( stderr, "ireg: %d  rn: %s  index: %d  xflag: %u  mateindex: %d\n",
ireg, alData[dupSet[idup].index].read_name, dupSet[idup].index, alData[dupSet[idup].index].xflag, alData[dupSet[idup].index].mateindex );
fprintf( stderr, "idup: %d  ndup: %d  mdup: %d  dupSet[idup].index: %d  dupSet[idup].mateidup: %d  dupSet[dupSet[idup].mateidup].index: %d\n",
idup, ndup, mdup, dupSet[idup].index, dupSet[idup].mateidup, dupSet[dupSet[idup].mateidup].index );
          *fstatus = -1;
          return( 0 );
        }

        kdup = dupSet[idup].mateidup;
        if( dupSet[idup+1].mateidup >= 0 ) dupSet[dupSet[idup+1].mateidup].mateidup = kdup;
        memcpy( &(dupSet[kdup]), &(dupSet[idup+1]), sizeof( DupSet ) );

        ++idup;
      }
    }
  }

  *fnrem   += nrem;
  *fstatus  = 0;

  return( 0 );
}

#else

static int xprocRegAlign( Genome *genome, HiCovReg *hiCovReg, AlData *alData, int ialign0, int ialign1, int nstart, int ndup, double pcrDupRate, int *fnrem, int *fstatus )
{
  int ireg;
  int nrem;
  int ialign;
  int idup, jdup, kdup, mdup;
  int status;

  DupSet tdupSet;

  static int mxDupSet = 0;
  static DupSet *dupSet = NULL;

  ireg = alData[ialign0].hcrindex;

  /*
  ** Consistency tests.
  */
  xcheckRegAlign( genome, alData, ialign0, ialign1, hiCovReg, nstart, ndup, &status );
  IF_STATUS_ZERO( "bad status: xcheckRegAlign" );

  /*
  ** Re-allocate duplicate array, if necessary.
  */
  if( ndup > mxDupSet )
  {
    mxDupSet += ndup + 1024;
    dupSet = (DupSet *)realloc( dupSet, mxDupSet * sizeof( DupSet ) );
    if( dupSet == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
  }

  /*
  ** Calculate number of duplicates to remove.
  */
  mdup = (int)( pcrDupRate * (double)nstart + 0.5 );

  /*
  ** Do we remove all duplicates?
  */
  if( mdup >= ndup )
  {
    /*
    ** Mark all duplicates (and their mates) for removal and
    ** return.
    */
    nrem = 0;
    for( ialign = ialign0; ialign <= ialign1; ++ialign )
    {
      if( alData[ialign].xflag & BF1_B01 && !( alData[ialign].xflag & BF1_B03 ) )
      {
        alData[ialign].xflag |= BF1_B03;
        ++nrem;
      }

      if( alData[ialign].xflag & BF1_B00 && alData[alData[ialign].mateindex].xflag & BF1_B01 && !( alData[alData[ialign].mateindex].xflag & BF1_B03 ) )
      {
        alData[alData[ialign].mateindex].xflag |= BF1_B03;
        ++nrem;
      }
    }

    *fnrem   += nrem;
    *fstatus  = 0;
    return( 0 );
  }

  /*
  ** Copy duplicate indices.
  ** Note: watch and account for alignments
  **       removed previously as mates.
  */
  idup = 0;
  for( ialign = ialign0; ialign <= ialign1; ++ialign )
  {
    if( alData[ialign].xflag & BF1_B03 )
    {
      --ndup;
      --mdup;
      continue;
    }

    if( alData[ialign].xflag & BF1_B01 )
    {
      dupSet[idup].index = ialign;
      dupSet[idup].flag  = 0;
      ++idup;
    }
  }

  if( idup == 0 )
  {
    *fstatus = 0;
    return( 0 );
  }

  if( idup != ndup )
  {
    sprintf( _msg_, "inconsistent counts for idup/ndup: %d != %d", idup, ndup );
    EMSG( _msg_ );
    *fstatus = -1;
    return( 0 );
  }

#ifdef VERIFY_LONG
  /*
  ** Consistency test: check that mates in the region are in dupSet.
  */
  xcheckDupSet( alData, dupSet, ndup, hiCovReg, ireg, &status );
  IF_STATUS_ZERO( "bad status: xcheckDupSet" );
#endif

  /*
  ** Shuffle randomly dupSet elements.
  */
  kdup = ndup;
  while( kdup > 0 )
  {
    jdup = xorShift1024() % kdup;

    memcpy( &tdupSet, &(dupSet[kdup-1]), sizeof( DupSet ) );
    memcpy( &(dupSet[kdup-1]), &(dupSet[jdup]), sizeof( DupSet ) );
    memcpy( &(dupSet[jdup]), &tdupSet, sizeof( DupSet ) );
    --kdup;
  }

  if( mdup > 1000000 )
  {
    fprintf( stderr, "(ndup: %d  mdup: %d)", ndup, mdup );
  }

  /*
  ** The dupSet elements are shuffled randomly so select
  ** from idup = 0 to idup = mdup-1, in order.
  */
  nrem = 0;
  for( idup = 0; idup < mdup; ++idup )
  {
    /*
    ** Check for previously removed duplicates/mates.
    */
    if( alData[dupSet[idup].index].xflag & BF1_B03 ||
        ( alData[dupSet[idup].index].xflag & BF1_B00 && alData[alData[dupSet[idup].index].mateindex].xflag & BF1_B03 ) )
    {
      EMSG( "unexpected condition" );
      *fstatus = -1;
      return( 0 );
    }

    /*
    ** Mark alignment and mate for removal.
    */
    alData[dupSet[idup].index].xflag |= BF1_B03;
    ++nrem;

    if( alData[dupSet[idup].index].xflag & BF1_B00 && alData[alData[dupSet[idup].index].mateindex].xflag & BF1_B02 )
    {
      if( alData[alData[dupSet[idup].index].mateindex].xflag & BF1_B03 )
      {
        EMSG( "unexpected condition" );
        *fstatus = -1;
        return( 0 );
      }

      alData[alData[dupSet[idup].index].mateindex].xflag |= BF1_B03;
      ++nrem;

      /*
      ** If mate is in ireg, increment idup and replace it
      ** in dupSet, if necessary.
      */
      if( alData[alData[dupSet[idup].index].mateindex].hcrindex == ireg )
      {
        for( kdup = idup + 1; kdup < mdup; ++kdup )
        {
          if( dupSet[kdup].index == alData[dupSet[idup].index].mateindex )
          {
            memcpy( &(dupSet[kdup]), &(dupSet[idup+1]), sizeof( DupSet ) );
            break;
          }
        }
        ++idup;
      }
    }
  }

  *fnrem   += nrem;
  *fstatus  = 0;

  return( 0 );
}

#endif


static int xcheckRegAlign( Genome *genome, AlData *alData, int ialign0, int ialign1, HiCovReg *hiCovReg, int nstart, int ndup, int *fstatus )
{
  int ireg;
  int ialign;

  ireg = alData[ialign0].hcrindex;

  if( hiCovReg[ireg].nstart != nstart )
  {
    EMSG( "inconsistent read start count for region" );
    sprintf( _msg_, "ireg: %d  refid: %d  nstart: %d  hiCovReg.nstart: %d", ireg, hiCovReg[ireg].refid, nstart, hiCovReg[ireg].nstart );
    EMSG( _msg_ );
    *fstatus = -1;
    return( 0 );
  }

  if( hiCovReg[ireg].ndup != ndup )
  {
    int ichr;
    int ibas;
    sprintf( _msg_, "inconsistent detected duplicate read count for region: %d != %d", ndup, hiCovReg[ireg].ndup );
    EMSG( _msg_ );
    fprintf( stderr, "ireg: %d [%d->%d] starts: %d  ialign: %d->%d\n",
             ireg,
             hiCovReg[ireg].beg,
             hiCovReg[ireg].end,
             nstart,
             ialign0,
             ialign1 );
    for( ialign = ialign0-2; ialign <= ialign1+2; ++ialign )
    {
      fprintf( stderr, "  ialign: %d  refid: %d  pos5p: %d  xflag: %d %d %d %d  flag: %u  rn: %s\n",
               ialign,
               alData[ialign].refid,
               alData[ialign].pos5p,
               alData[ialign].xflag & BF1_B00 ? 1 : 0,
               alData[ialign].xflag & BF1_B01 ? 1 : 0,
               alData[ialign].xflag & BF1_B02 ? 1 : 0,
               alData[ialign].xflag & BF1_B03 ? 1 : 0,
               alData[ialign].flag,
               alData[ialign].read_name );
    }
    ichr = hiCovReg[ireg].refid;
    fprintf( stderr, "\n" );
    fprintf( stderr, "cover\n" );
    for( ibas = hiCovReg[ireg].beg - 2; ibas <= hiCovReg[ireg].end + 2; ++ibas )
    {
      fprintf( stderr, "BAS: %d  %ld  %ld  %ld %ld %ld %ld  %d %d %d %d\n",
               ibas,
               genome->chr[ichr].cover[0][ibas-1] + genome->chr[ichr].cover[2][ibas-1],
               genome->chr[ichr].cover[1][ibas-1] + genome->chr[ichr].cover[3][ibas-1],
               genome->chr[ichr].cover[0][ibas-1],
               genome->chr[ichr].cover[1][ibas-1],
               genome->chr[ichr].cover[2][ibas-1],
               genome->chr[ichr].cover[3][ibas-1],
               genome->chr[ichr].flag[ibas-1] & BF2_B08 ? 1 : 0,
               genome->chr[ichr].flag[ibas-1] & BF2_B09 ? 1 : 0,
               genome->chr[ichr].flag[ibas-1] & BF2_B10 ? 1 : 0,
               genome->chr[ichr].flag[ibas-1] & BF2_B11 ? 1 : 0 );
    }

    *fstatus = -1;
    return( 0 );
  }

  for( ialign = ialign0; ialign <= ialign1; ++ialign )
  {
    /*
    ** Must be on same chromosome.
    */
    if( ialign > ialign0 && alData[ialign].refid != alData[ialign0].refid )
    {
      sprintf( _msg_, "reads align to inconsistent chromosomes: %d != %d", alData[ialign].refid, alData[ialign0].refid );
      EMSG( _msg_ );
      *fstatus = -1;
      return( 0 );
    }

    /*
    ** pos5p must be in high coverage region.
    */
    if( alData[ialign].pos5p < hiCovReg[ireg].beg || alData[ialign].pos5p > hiCovReg[ireg].end )
    {
      sprintf( _msg_, "reads alignment start not in high coverage region: %d not between %d and %d  ialign: [%d,%d] %d",
               alData[ialign].pos5p, hiCovReg[ireg].beg, hiCovReg[ireg].end, ialign0, ialign1, ialign );
      EMSG( _msg_ );
      *fstatus = -1;
      return( 0 );
    }

    /*
    ** Must have high coverage region bit set.
    */
    if( !( alData[ialign].xflag & BF1_B02 ) )
    {
      sprintf( _msg_, "reads alignment high coverage region bit not set" );
      EMSG( _msg_ );
      *fstatus = -1;
      return( 0 );
    }

    /*
    ** Consistent high coverage region bit set for proper pair reads.
    */
    if( alData[ialign].xflag & BF1_B00 && ( ( alData[ialign].xflag & BF1_B01 ) != ( alData[alData[ialign].mateindex].xflag & BF1_B01 ) ) )
    {
      sprintf( _msg_, "inconsistent xflag bit set for properly paired reads: %u %u", alData[ialign].xflag & BF1_B01, alData[alData[ialign].mateindex].xflag & BF1_B01 );
      EMSG( _msg_ );
      *fstatus = -1;
      return( 0 );
    }

    /*
    ** Must be in same high coverage region.
    */
    if( ialign > ialign0 && alData[ialign].hcrindex != alData[ialign0].hcrindex )
    {
      sprintf( _msg_, "reads align to inconsistent high coverage regions: %d != %d", alData[ialign].hcrindex, alData[ialign0].hcrindex );
      EMSG( _msg_ );
      *fstatus = -1;
      return( 0 );
    }

    /*
    ** Check for consistent mate region information.
    */
    if( alData[ialign].xflag & BF1_B00 &&
        alData[alData[ialign].mateindex].xflag & BF1_B02 &&
        alData[alData[ialign].mateindex].hcrindex == ireg &&
        alData[alData[ialign].mateindex].refid == hiCovReg[ireg].refid &&
        ( alData[alData[ialign].mateindex].pos5p < hiCovReg[ireg].beg || alData[alData[ialign].mateindex].pos5p > hiCovReg[ireg].end ) )
    {
      sprintf( _msg_, "inconsistent mate information" );
      EMSG( _msg_ );
      sprintf( _msg_, "ialign: %d  ialign0: %d  ialign1: %d\n", ialign, ialign0, ialign1 );
      EMSG( _msg_ );
      sprintf( _msg_, "ireg: %d  beg: %d  end: %d\n", ireg, hiCovReg[ireg].beg, hiCovReg[ireg].end );
      EMSG( _msg_ );
      sprintf( _msg_, "ialign: %d  hcrindex: %d  pos5p: %d  mate: index: %d  hcrindex: %d  pos5p: %d\n",
               ialign,
               alData[ialign].hcrindex,
               alData[ialign].pos5p,
               alData[alData[ialign].mateindex].index,
               alData[alData[ialign].mateindex].hcrindex,
               alData[alData[ialign].mateindex].pos5p );
      EMSG( _msg_ );
      *fstatus = -1;
      return( 0 );
    }
  }

  *fstatus = 0;

  return( 0 );
}


static int xcheckDupSet( AlData *alData, DupSet *dupSet, int ndup, HiCovReg *hiCovReg, int ireg, int *fstatus )
{
  int idup, jdup;
  int flag;

  if( ndup > 100000 )
  {
    fprintf( stderr, "(ndup: %d)", ndup );
  }
  for( idup = 0; idup < ndup; ++idup )
  {
    if( ndup > 100000 && idup % 10000 == 0 )
    {
      fprintf( stderr, ":" );
    }

    if( alData[dupSet[idup].index].xflag & BF1_B00 &&
        alData[alData[dupSet[idup].index].mateindex].xflag & BF1_B02 &&
        alData[alData[dupSet[idup].index].mateindex].hcrindex == ireg &&
        alData[dupSet[idup].index].pos5p < alData[alData[dupSet[idup].index].mateindex].pos5p )
    {
      flag = 0;
      for( jdup = idup + 1; jdup < ndup; ++jdup )
      {
        if( dupSet[jdup].index == alData[alData[dupSet[idup].index].mateindex].index )
        { 
          flag = 1;
          break;
        } 
      }

      if( !flag )
      {
        EMSG( "missing mate in dupSet" );
        fprintf( stderr, "ndup: %d", ndup );
        fprintf( stderr, "idup: %d  ireg: %d  iref: %d  pos5p: %d  hcrindex: %d  dup.iref: %d  dup.pos5p: %d  dup.hcrindex: %d\n",
                 idup,
                 ireg,
                 alData[dupSet[idup].index].refid,
                 alData[dupSet[idup].index].pos5p,
                 alData[dupSet[idup].index].hcrindex,
                 alData[alData[dupSet[idup].index].mateindex].refid,
                 alData[alData[dupSet[idup].index].mateindex].pos5p,
                 alData[alData[dupSet[idup].index].mateindex].hcrindex );
        fprintf( stderr, "hiCovReg: %d to %d\n", hiCovReg[alData[dupSet[idup].index].hcrindex].beg, hiCovReg[alData[dupSet[idup].index].hcrindex].end );
        *fstatus = -1;
        return( 0 );
      }
    }
  }

  *fstatus = 0;

  return( 0 );
}


static int xreportHiCovRegDist( HiCovReg *hiCovReg, int nreg, int *fstatus )
{
  int i;
  int binwid;

  int64_t len, ibin;
  int64_t mxlen;
  int64_t hist[1024];
  int64_t sum, tot;

  binwid = 50;

  for( i = 0; i < 1024; ++i )
  {
    hist[i] = 0;
  }

  tot = 0;
  mxlen = 0L;
  for( i = 0; i < nreg; ++i )
  {
    len = hiCovReg[i].end - hiCovReg[i].beg + 1;
    if( len > mxlen )
    {
      mxlen = len;
    }
    ibin = len / binwid;
    ibin = ibin < 1000 ? ibin : 999;
    ++hist[ibin];
    ++tot;
  }

  fprintf( stderr, "High coverage region length distribution\n" );
  fprintf( stderr, "Region   num. of   rev. cumul.     fraction of rev. cumul.\n" );
  fprintf( stderr, "range    regions   num. regions    num. regions\n" );
  for( i = 0; i < 1000; ++i )
  {
    if( hist[i] > 0 )
    {
      sum += hist[i];
      fprintf( stderr, "%d-%d  %ld  %ld  %.3f\n",
               i * binwid,
               ( i + 1 ) * binwid - 1,
               hist[i],
               tot - sum,
               (double)( tot - sum ) / (double)tot );
    }
  }
  fprintf( stderr, "maximum length region: %ld\n", mxlen );
  fprintf( stderr, "\n" );

  *fstatus = 0;

  return( 0 );
}


static int xremoveDuplicateReads( Bundle *bundle, BamHeader *bamHeader, AlData *alData, int numAlData, Genome *genome, int rsFlag, double pcrDupRate, int64_t counts[MXCOUNTER], int *fstatus )
{
  int flag;
  int ichr, nchr;
  int ibas, nbas;
  int ibas0, ibas1;
  int ireg, mreg, nreg;
  int ialign, nalign;
  int ialign0, ialign1;
  int ndup, nstart;
  int tbas[2];
  int tdup[2], tstart[2];
  int mrem, nrem;
  int status;

  uint16_t btest;

  Chromosome *chr;
  HiCovReg   *hiCovReg;

  if( rsFlag == RS_SINGLE )
  {
    btest = SE_BIT_TEST;
  }
  else
  if( rsFlag == RS_PAIRED )
  {
    btest = PE_BIT_TEST;
  }
  else
  {
    EMSG( "unrecognized read type value" );
    *fstatus = -1;
    return( 0 );
  }

  nchr   = genome->numChr;
  nalign = numAlData;

  /*
  ** Low coverage region counts for statistics.
  */
  tbas[1]   = 0;
  tstart[1] = 0;
  tdup[1]   = 0;
  mreg      = 0;
  nstart    = 0;
  ndup      = 0;
  for( ichr = 0; ichr < nchr; ++ichr )
  {
    chr  = &(genome->chr[ichr]);
    nbas = chr->len;

    ibas0 = 0;
    ibas1 = 0;
    flag  = 0;
    for( ibas = 0; ibas < nbas; ++ibas )
    {
      if( !flag && !( chr->flag[ibas] & btest ) )
      {
        ibas0  = ibas;
        nstart = chr->cover[0][ibas] + chr->cover[2][ibas];
        ndup   = chr->cover[1][ibas] + chr->cover[3][ibas];
        flag   = 1;
      }
      else
      if( flag )
      {
        if( chr->flag[ibas] & btest )
        {
          ibas1     = ibas - 1;
          tbas[1]   += ibas1 - ibas0 + 1;
          tstart[1] += nstart;
          tdup[1]   += ndup;
          flag      = 0;
          counts[C_NOT_HCOV_BASE] += ibas1 - ibas0 + 1;
          ++mreg;
          continue;
        }

        nstart += chr->cover[0][ibas] + chr->cover[2][ibas];
        ndup   += chr->cover[1][ibas] + chr->cover[3][ibas];
      }
    } /* for ibas */
    if( flag )
    {
      ibas1      = ibas - 1;
      tbas[1]   += ibas1 - ibas0 + 1;
      tstart[1] += nstart;
      tdup[1]   += ndup;
      counts[C_NOT_HCOV_BASE] += ibas1 - ibas0 + 1;
      ++mreg;
    }
  }

  /*
  ** Count high coverage regions and allocate region storage.
  */
  nreg = 0;
  for( ichr = 0; ichr < nchr; ++ichr )
  {
    chr  = &(genome->chr[ichr]);
    nbas = chr->len;

    flag = 0;
    for( ibas = 0; ibas < nbas; ++ibas )
    {
      if( !flag && chr->flag[ibas] & btest )
      {
        flag = 1;
      }
      else
      if( flag && !( chr->flag[ibas] & btest ) )
      {
        flag = 0;
        ++nreg;
      }
    } /* for ibas */

    if( flag )
    {
      ++nreg;
    }
  } /* for ichr */

  hiCovReg = (HiCovReg *)calloc( nreg, sizeof( HiCovReg ) );
  if( hiCovReg == NULL )
  {
    EMSG( "unable to allocate memory" );
    *fstatus = -1;
    return( 0 );
  }
  bundle->hiCovReg    = hiCovReg;
  bundle->numHiCovReg = nreg;

  /*
  ** Record high coverage regions.
  */
  tbas[0]   = 0;
  tstart[0] = 0;
  tdup[0]   = 0;
  ireg = 0;
  for( ichr = 0; ichr < nchr; ++ichr )
  {
    chr  = &(genome->chr[ichr]);
    nbas = chr->len;

    ibas0 = 0;
    ibas1 = 0;
    flag  = 0;
    for( ibas = 0; ibas < nbas; ++ibas )
    {
      if( !flag && chr->flag[ibas] & btest )
      {
        ibas0  = ibas;
        nstart = chr->cover[0][ibas] + chr->cover[2][ibas];
        ndup   = chr->cover[1][ibas] + chr->cover[3][ibas];
        flag   = 1;
      }
      else
      if( flag )
      {
        if( !( chr->flag[ibas] & btest ) )
        {
          ibas1 = ibas - 1;
  
          hiCovReg[ireg].refid  = ichr;
          hiCovReg[ireg].beg    = ibas0 + 1;
          hiCovReg[ireg].end    = ibas1 + 1;
          hiCovReg[ireg].nstart = nstart; 
          hiCovReg[ireg].ndup   = ndup; 
          tbas[0]   += ibas1 - ibas0 + 1;
          tstart[0] += nstart;
          tdup[0]   += ndup;
          flag       = 0;
          counts[C_HCOV_BASE] += ibas1 - ibas0 + 1;

          ++ireg;
          continue;
        }

        nstart += chr->cover[0][ibas] + chr->cover[2][ibas];
        ndup   += chr->cover[1][ibas] + chr->cover[3][ibas];
      }
    } /* for ibas */

    if( flag )
    {
      ibas1 = ibas - 1;

      hiCovReg[ireg].refid  = ichr;
      hiCovReg[ireg].beg    = ibas0 + 1;
      hiCovReg[ireg].end    = ibas1 + 1;
      hiCovReg[ireg].nstart = nstart; 
      hiCovReg[ireg].ndup   = ndup; 
      tbas[0]   += ibas1 - ibas0 + 1;
      tstart[0] += nstart;
      tdup[0]   += ndup;
      counts[C_HCOV_BASE] += ibas1 - ibas0 + 1;
      ++ireg;
    }
  } /* for ichr */

  if( ireg != nreg )
  {
    sprintf( _msg_, "inconsistent region counts: %d != %d", ireg, nreg );
    EMSG( _msg_ );
    *fstatus = -1;
    return( 0 );
  }

  /*
  ** Report high coverage region length distribution.
  */
  xreportHiCovRegDist( hiCovReg, nreg, &status );
  IF_STATUS_ZERO( "bad status: xreportHiCovRegDist" )

  /*
  ** Check that the alignments are sorted correctly.
  */
  for( ialign = 1; ialign < nalign; ++ialign )
  {
    if( alData[ialign].refid < alData[ialign-1].refid )
    {
      sprintf( _msg_, "mis-sorted alignment data: refids: %d < %d for ialign %d", alData[ialign].refid, alData[ialign-1].refid, ialign );
      EMSG( _msg_ );
      *fstatus = -1;
      return( 0 );
    }
    else
    if( alData[ialign].refid == alData[ialign-1].refid &&
        alData[ialign].pos5p <  alData[ialign-1].pos5p )
    {
      sprintf( _msg_, "mis-sorted alignment data: pos5p: %d < %d for ialign %d", alData[ialign].pos5p, alData[ialign-1].pos5p, ialign );
      EMSG( _msg_ );
      *fstatus = -1;
      return( 0 );
    }
  }

  /*
  ** Set high coverage region information in AlData.
  ** Note: set for both properly paired and not
  **       properly paired alignments.
  */
  ialign = 0;
  for( ireg = 0; ireg < nreg; ++ireg )
  {
    ichr  = hiCovReg[ireg].refid;
    ibas0 = hiCovReg[ireg].beg;
    ibas1 = hiCovReg[ireg].end;

    while( ialign < nalign && alData[ialign].refid < ichr )
    {
      ++ialign;
    }

    while( ialign < nalign && alData[ialign].refid == ichr && alData[ialign].pos5p < ibas0 )
    {
      ++ialign;
    }

    if( alData[ialign].refid > ichr )
    {
      continue;
    }

    for( ; ialign < nalign; ++ialign )
    {
      if( alData[ialign].pos5p > ibas1 || alData[ialign].refid != ichr )
      {
        break;
      }

      alData[ialign].xflag    |= BF1_B02;
      alData[ialign].hcrindex  = ireg;
    }
  }

  /*
  ** Initialize random number generator.
  */
  xinitXorShift1024();

  /*
  ** Select detected duplicate reads for removal in high coverage regions.
  */
  ireg     = 0;
  flag     = 0;
  nstart   = 0;
  ndup     = 0;
  ialign0  = 0;
  ialign1  = 0;
  nrem     = 0;
  for( ialign = 0; ialign < nalign; ++ialign )
  {
    if( ialign % 100000 == 0 )
    {
      fprintf( stderr, "." );
    }

    if( flag == 0 && alData[ialign].xflag & BF1_B02 )
    {
      /*
      ** High coverage region start.
      */
      nstart  = 1;
      ndup    = ( alData[ialign].xflag & BF1_B01 ) ? 1 : 0;
      ialign0 = ialign;
      flag    = 1;
    }
    else
    if( flag == 1 )
    {
      /*
      ** High coverage region end (previous alignment).
      */
      if( !( alData[ialign].xflag & BF1_B02 ) || ( ialign > 0 && alData[ialign].hcrindex != alData[ialign-1].hcrindex ) )
      {
        ialign1 = ialign - 1;

        xprocRegAlign( genome, hiCovReg, alData, ialign0, ialign1, nstart, ndup, pcrDupRate, &nrem, &status );
        IF_STATUS_ZERO( "bad status: xprocRegAlign" );

        flag = 0;
        ++ireg;
        --ialign;
        continue;
      }

      ++nstart;
      if( alData[ialign].xflag & BF1_B01 )
      {
        ++ndup;
      }
    }
  }

  if( flag )
  {
    ialign1 = ialign - 1;

    xprocRegAlign( genome, hiCovReg, alData, ialign0, ialign1, nstart, ndup, pcrDupRate, &nrem, &status );
    IF_STATUS_ZERO( "bad status: xprocRegAlign" );
    ++ireg;
  }
  fprintf( stderr, "\n" );

  if( ireg != nreg )
  {
    sprintf( _msg_, "inconsistent regions counts: %d != %d", ireg, nreg );
    EMSG( _msg_ );
    *fstatus = -1;
    return( 0 );
  }

  /*
  ** Check for consistent read removal counts.
  */
  mrem = 0;
  for( ialign = 0; ialign < nalign; ++ialign )
  {
    if( alData[ialign].xflag & BF1_B03 )
    {
      ++mrem;
    }
  }

  if( mrem != nrem )
  {
    sprintf( _msg_, "inconsistent read removal counts: %d != %d", mrem, nrem );
    EMSG( _msg_ );
    *fstatus = -1;
    return( 0 );
  }

/*
  for( ireg = 0; ireg < nreg; ++ireg )
  {
    fprintf( stderr, "HCR: %d %d %d %d %d %d %.4f %.4f\n",
             hiCovReg[ireg].refid,
             hiCovReg[ireg].beg,
             hiCovReg[ireg].end,
             hiCovReg[ireg].nstart,
             hiCovReg[ireg].ndup,
             hiCovReg[ireg].end - hiCovReg[ireg].beg + 1,
             (double)hiCovReg[ireg].nstart / (double)( hiCovReg[ireg].end - hiCovReg[ireg].beg + 1 ),
             hiCovReg[ireg].nstart > 0 ? (double)hiCovReg[ireg].ndup / (double)hiCovReg[ireg].nstart : 0 );
  }
*/

  fprintf( stderr, "SUMHCR: %d %d %d %d %.4f %.4f\n",
           nreg,
           tbas[0],
           tstart[0],
           tdup[0],
           (double)tstart[0] / (double)tbas[0],
           tstart[0] > 0 ? (double)tdup[0] / (double)tstart[0] : 0.0 );

  fprintf( stderr, "SUMLCR: %d %d %d %d %.4f %.4f\n",
           mreg,
           tbas[1],
           tstart[1],
           tdup[1],
           (double)tstart[1] / (double)tbas[1],
           tstart[1] > 0 ? (double)tdup[1] / (double)tstart[1] : 0.0 );

  /*
  ** Mark duplicate reads for removal in low coverage regions.
  */
  for( ialign = 0; ialign < nalign; ++ialign )
  {
    if( !( alData[ialign].xflag & BF1_B02 ) )
    {
      ++counts[C_NOT_HCOV_READ];
      if( alData[ialign].xflag & BF1_B01 )
      {
        alData[ialign].xflag |= BF1_B03;
        ++counts[C_NOT_HCOV_DUP_DETECT];
        ++counts[C_NOT_HCOV_DUP_REMOVE];
      }
    }
    else
    {
      ++counts[C_HCOV_READ];
      if( alData[ialign].xflag & BF1_B01 )
      {
        ++counts[C_HCOV_DUP_DETECT];
        if( alData[ialign].xflag & BF1_B03 )
        {
          ++counts[C_HCOV_DUP_REMOVE];
        }
      }
    }
  }

  *fstatus = 0;

  return( 0 );
}


static int xcheckAlign( Bundle *bundle, BamHeader *bamHeader, AlData *alData, int numAlData, Genome *genome, int64_t counts[MXCOUNTER], int *fstatus )
{
  int ireg, nreg;
  int ichr, nchr;
  int ibas, nbas;
  int ialign, nalign;

  HiCovReg   *hiCovReg;
  Chromosome *chr;

  hiCovReg = bundle->hiCovReg;
  nreg     = bundle->numHiCovReg;
  nalign   = numAlData;

  static int  mreg = 0;
  static int *reg  = NULL;

  for( ialign = 0; ialign < nalign; ++ialign )
  {
    if( !( alData[ialign].xflag & BF1_B01 ) &&
           alData[ialign].xflag & BF1_B03 )
    {
      *fstatus = -1;
      return( 0 );
      EMSG( "unexpected condition" );
    }

    if( alData[ialign].xflag & BF1_B00 &&
        alData[ialign].xflag & BF1_B01 &&
        !( alData[alData[ialign].mateindex].xflag & BF1_B01 ) )
    {
      EMSG( "unexpected condition" );
      *fstatus = -1;
      return( 0 );
    }

    if( !( alData[ialign].xflag & BF1_B01 ) && alData[ialign].xflag & BF1_B03 )
    {
      EMSG( "non-detected duplicate marked for removal" );
      *fstatus = -1;
      return( 0 );
    }

    /*
    ** Is detected duplicate read?
    */
    if( alData[ialign].xflag & BF1_B01 )
    {
      /*
      ** Is in high coverage region?
      */
      if( alData[ialign].xflag & BF1_B02 )
      {
        ireg = alData[ialign].hcrindex;
  
        if( alData[ialign].refid != hiCovReg[ireg].refid ||
            alData[ialign].pos5p < hiCovReg[ireg].beg ||
            alData[ialign].pos5p > hiCovReg[ireg].end )
        {
          EMSG( "unexpected condition" );
          *fstatus = -1;
          return( 0 );
        }
      }
    }
  }

  /*
  ** Check for alignments in high coverage regions that are
  ** not marked as in high coverage region.
  */
#ifdef VERIFY_LONG
  nchr = genome->numChr;
  for( ichr = 0; ichr < nchr; ++ichr )
  {
    chr = &(genome->chr[ichr]);
    nbas = chr->len;
    if( nbas > mreg )
    {
      mreg = nbas * 2;
      reg = (int *)realloc( reg, mreg * sizeof( int ) );
      if( reg == NULL )
      {
        EMSG( "unable to allocate memory" );
        *fstatus = -1;
        return( 0 );
      }
    }
    memset( reg, 0, nbas * sizeof( int ) );

    nreg = bundle->numHiCovReg;
    for( ireg = 0; ireg < nreg; ++ireg )
    {
      if( hiCovReg[ireg].refid < ichr )
      {
        continue;
      }
      else
      if( hiCovReg[ireg].refid > ichr )
      {
        break;
      }

      for( ibas = hiCovReg[ireg].beg; ibas <= hiCovReg[ireg].end; ++ibas )
      {
        reg[ibas-1] = ireg + 1;
      }
    }

    for( ialign = 0; ialign < nalign; ++ialign )
    {
      if( alData[ialign].refid < ichr )
      {
        continue;
      }

      if( alData[ialign].refid > ichr )
      {
        break;
      }

      if( reg[alData[ialign].pos5p-1] && !( alData[ialign].xflag & BF1_B02 ) )
      {
        EMSG( "alignment not marked in high coverage region" );
        *fstatus = -1;
        return( 0 );
      }
    }
  }

  free( reg );
#endif

  *fstatus = 0;

  return( 0 );
}


static int xwriteReads( FILE *ofp, Bundle *bundle, BamHeader *bamHeader, AlData *alData, int numAlData,  Genome *genome, int *fstatus )
{
  int nrem;
  int ialign, nalign;

  nalign   = numAlData;

  nrem = 0;
  for( ialign = 0; ialign < nalign; ++ialign )
  {
    if( alData[ialign].xflag & BF1_B03 )
    {
      ++nrem;
    }
  }
  fprintf( ofp, "HEADER: %d | %s | %s\n",
           nrem,
           bundle->starttime,
           bundle->nameBamFile );

  for( ialign = 0; ialign < nalign; ++ialign )
  {
    if( !( alData[ialign].xflag & BF1_B03 ) )
    {
      continue;
    }

    fprintf( ofp, "ALIGN: %d %d %d %d %d %s %d %u %d %s\n",
             alData[ialign].xflag & BF1_B02 ? alData[ialign].hcrindex : -1,
             alData[ialign].xflag & BF1_B00 ? 1 : 0,   /* properly paired read */
             alData[ialign].xflag & BF1_B01 ? 1 : 0,   /* detected duplicate read */
             alData[ialign].xflag & BF1_B02 ? 1 : 0,   /* high coverage region read */
             alData[ialign].xflag & BF1_B03 ? 1 : 0,   /* read selected for removal */
             bamHeader->refSeq[alData[ialign].refid].name,
             alData[ialign].pos5p,
             alData[ialign].flag,
             alData[ialign].pos,
             alData[ialign].read_name );
  }

  *fstatus = 0;

  return( 0 );
}


static int xreportCounts( FILE *ofp, Bundle *bundle, double pcrDupRateSE, double pcrDupRatePE, int64_t hist[2][MXHIST], int64_t counts[3][MXCOUNTER], int *fstatus )
{
  int i;
  int widFld;

  int64_t sum, tot;

  char *cptr;

  widFld = 9;

  fprintf( ofp, "BAM file entries\n" );
  fprintf( ofp, "  %*ld entries read\n", widFld, counts[0][C_NUM_ENTRY_BAM] );
  fprintf( ofp, "  %*ld unaligned skipped\n", widFld, counts[0][C_UNALIGNED_BAM] );
  fprintf( ofp, "  %*ld chimeric skipped\n", widFld, counts[0][C_CHIMERIC_BAM] );
  fprintf( ofp, "  %*ld secondary skipped\n", widFld, counts[0][C_SECONDARY_BAM] );
  fprintf( ofp, "  %*ld not paired end skipped\n", widFld, counts[0][C_NOT_PAIREDEND_BAM] );
  fprintf( ofp, "  %*ld accepted\n", widFld, counts[0][C_ACCEPTED_BAM] );

  fprintf( ofp, "SL matches\n" );
  fprintf( ofp, "  %*ld alignments with SL match\n", widFld, counts[0][C_SLMATCH] );
  fprintf( ofp, "  %*ld alignments with SL match and short intron\n", widFld, counts[0][C_SHORTINTRON_SLMATCH] );
  fprintf( ofp, "  %*ld alignments with SL match and long intron\n", widFld, counts[0][C_LONGINTRON_SLMATCH] );
  fprintf( ofp, "  %*ld alignments without SL match\n", widFld, counts[0][C_NOT_SLMATCH] );

  fprintf( ofp, "\n" );

/*
  fprintf( ofp, "Short exon match at alignment start\n" );
  fprintf( ofp, "  %*ld alignments with short match\n", widFld, counts[0][C_SHORTMATCH] );
  fprintf( ofp, "  %*ld alignments without short match\n", widFld, counts[0][C_NOT_SHORTMATCH] );

  fprintf( ofp, "\n" );
*/

  fprintf( ofp, "Paired end reads\n" );

  fprintf( ofp, "  Proper paired reads\n" );
  fprintf( ofp, "    %*ld proper paired reads\n", widFld, counts[2][C_START_PAIRED_READ] );
  fprintf( ofp, "    %*ld duplicate proper paired reads\n", widFld, counts[2][C_DUP_START_PAIRED_READ] );
  fprintf( ofp, "    %*.4f overall duplicate rate\n", widFld, counts[2][C_START_PAIRED_READ] > 0 ? (double)counts[2][C_DUP_START_PAIRED_READ] / (double)counts[2][C_START_PAIRED_READ] : 0.0 );

  fprintf( ofp, "  Not proper paired reads\n" );
  fprintf( ofp, "    %*ld not proper paired reads\n", widFld, counts[2][C_START_UNPAIRED_READ] );
  fprintf( ofp, "    %*ld duplicate not proper paired reads\n", widFld, counts[2][C_DUP_START_UNPAIRED_READ] );
  fprintf( ofp, "    %*.4f overall duplicate rate\n", widFld, counts[2][C_START_UNPAIRED_READ] > 0 ? (double)counts[2][C_DUP_START_UNPAIRED_READ] / (double)counts[2][C_START_UNPAIRED_READ] : 0.0 );

  fprintf( ofp, "  Duplicate read removal\n" );
  fprintf( ofp, "    %*.4f estimated PCR duplicate rate\n", widFld, pcrDupRatePE );
  fprintf( ofp, "    %*ld duplicates detected\n", widFld, counts[2][C_NOT_HCOV_DUP_DETECT] + counts[2][C_HCOV_DUP_DETECT] );
  fprintf( ofp, "    %*ld duplicates marked for removal\n", widFld, counts[2][C_NOT_HCOV_DUP_REMOVE] + counts[2][C_HCOV_DUP_REMOVE] );
  fprintf( ofp, "    %*ld low coverage region bases\n", widFld, counts[2][C_NOT_HCOV_BASE] );
  fprintf( ofp, "    %*ld low coverage region reads\n", widFld, counts[2][C_NOT_HCOV_READ] );
  fprintf( ofp, "    %*ld low coverage region duplicates detected\n", widFld, counts[2][C_NOT_HCOV_DUP_DETECT] );
  fprintf( ofp, "    %*ld low coverage region duplicates marked for removal\n", widFld, counts[2][C_NOT_HCOV_DUP_REMOVE] );
  fprintf( ofp, "    %*ld high coverage region bases\n", widFld, counts[2][C_HCOV_BASE] );
  fprintf( ofp, "    %*ld high coverage region reads\n", widFld, counts[2][C_HCOV_READ] );
  fprintf( ofp, "    %*ld high coverage region duplicates detected\n", widFld, counts[2][C_HCOV_DUP_DETECT] );
  fprintf( ofp, "    %*ld high coverage region duplicates marked for removal\n", widFld, counts[2][C_HCOV_DUP_REMOVE] );

  fprintf( ofp, "\n" );

  fprintf( ofp, "  Number of duplicate read starts at low coverage region bases\n" );
  fprintf( ofp, "  (Note: for example, for two copies of a read, we count one\n" );
  fprintf( ofp, "  as a duplicate)\n" );
  fprintf( ofp, "\n" );
  fprintf( ofp, "  number of     number of low coverage region bases\n" );
  fprintf( ofp, "  duplicates    (and sum and fractions)\n" );
  tot = (int64_t)0;
  for( i = 1; i < MXHIST; ++i )
  {
    tot += hist[0][i];
  }
  sum = (int64_t)0;
  for( i = 1; i < MXHIST; ++i )
  {
    if( hist[0][i] > (int64_t)0 )
    {
      sum += hist[0][i];
      fprintf( ofp, "  %3d            %*ld    %.3f    %*ld    %.3f\n",
               i,
               widFld, (long int)hist[0][i],
               (double)hist[0][i] / (double)tot,
               widFld, sum,
               (double)sum / (double)tot );
    }
  }

  fprintf( ofp, "\n" );


  fprintf( ofp, "Single end reads\n" );

  fprintf( ofp, "  %*ld single end reads\n", widFld, counts[1][C_START_UNPAIRED_READ] );
  fprintf( ofp, "  %*ld duplicate single end reads\n", widFld, counts[1][C_DUP_START_UNPAIRED_READ] );
  fprintf( ofp, "  %*.4f overall duplicate rate\n", widFld, counts[1][C_START_UNPAIRED_READ] > 0 ? (double)counts[1][C_DUP_START_UNPAIRED_READ] / (double)counts[1][C_START_UNPAIRED_READ] : 0.0 );

  fprintf( ofp, "  Duplicate read removal\n" );
  fprintf( ofp, "    %*.4f estimated PCR duplicate rate\n", widFld, pcrDupRateSE );
  fprintf( ofp, "    %*ld duplicates detected\n", widFld, counts[1][C_NOT_HCOV_DUP_DETECT] + counts[1][C_HCOV_DUP_DETECT] );
  fprintf( ofp, "    %*ld duplicates marked for removal\n", widFld, counts[1][C_NOT_HCOV_DUP_REMOVE] + counts[1][C_HCOV_DUP_REMOVE] );
  fprintf( ofp, "    %*ld low coverage region bases\n", widFld, counts[1][C_NOT_HCOV_BASE] );
  fprintf( ofp, "    %*ld low coverage region reads\n", widFld, counts[1][C_NOT_HCOV_READ] );
  fprintf( ofp, "    %*ld low coverage region duplicates detected\n", widFld, counts[1][C_NOT_HCOV_DUP_DETECT] );
  fprintf( ofp, "    %*ld low coverage region duplicates marked for removal\n", widFld, counts[1][C_NOT_HCOV_DUP_REMOVE] );
  fprintf( ofp, "    %*ld high coverage region bases\n", widFld, counts[1][C_HCOV_BASE] );
  fprintf( ofp, "    %*ld high coverage region reads\n", widFld, counts[1][C_HCOV_READ] );
  fprintf( ofp, "    %*ld high coverage region duplicates detected\n", widFld, counts[1][C_HCOV_DUP_DETECT] );
  fprintf( ofp, "    %*ld high coverage region duplicates marked for removal\n", widFld, counts[1][C_HCOV_DUP_REMOVE] );

  fprintf( ofp, "\n" );

  fprintf( ofp, "  Number of duplicate read starts at low coverage region bases\n" );
  fprintf( ofp, "  (Note: for example, for two copies of a read, we count one\n" );
  fprintf( ofp, "  as a duplicate)\n" );
  fprintf( ofp, "\n" );
  fprintf( ofp, "  number of     number of low coverage\n" );
  fprintf( ofp, "  duplicates    region bases\n" );

  for( i = 1; i < MXHIST; ++i )
  {
    tot += hist[1][i];
  }
  sum = (int64_t)0;
  for( i = 1; i < MXHIST; ++i )
  {
    if( hist[1][i] > (int64_t)0 )
    {
      sum += hist[1][i];
      fprintf( ofp, "  %3d            %*ld    %.3f    %*ld    %.3f\n",
               i,
               widFld, (long int)hist[1][i],
               (double)hist[1][i] / (double)tot,
               widFld, sum,
               (double)sum / (double)tot );
    }
  }

  fprintf( ofp, "\n" );

  fprintf( ofp, "HEADER_SUMMARY: \
bam_file_name \
#_bam_entries \
#_unaligned_skipped \
#_chimeric_skipped \
#_secondary_skipped \
#_not_paired_end_skipped \
#_accepted_alignments \
\
#_sl_matches \
#_sl_match_with_short_intron \
#_sl_match_with_long_intron \
#_without_sl_match \
\
pe_#_proper_paired_reads \
pe_#_duplicate_proper_paired_reads \
pe_overall_proper_paired_read_duplicate_rate \
\
pe_#_not_proper_paired_reads \
pe_#_duplicate_not_proper_paired_reads \
pe_overall_not_proper_paired_read_duplicate_rate \
\
pe_estimated_pcr_duplicate_rate \
pe_#_duplicates_detected \
pe_#_duplicates_marked_for_removal \
pe_#_low_coverage_region_bases \
pe_#_low_coverage_region_reads \
pe_#_low_coverage_region_duplicates_detected \
pe_#_low_coverage_region_duplicates_marked_for_removal \
pe_#_high_coverage_region_bases \
pe_#_high_coverage_region_reads \
pe_#_high_coverage_region_duplicates_detected \
pe_#_high_coverage_region_duplicates_marked_for_removal \
\
se_#_single_end_reads \
se_#_duplicate_single_end_reads \
se_overall_single_end_read_duplicate_rate \
\
se_estimated_pcr_duplicate_rate \
se_#_duplicates_detected \
se_#_duplicates_marked_for_removal \
se_#_low_coverage_region_bases \
se_#_low_coverage_region_reads \
se_#_low_coverage_region_duplicates_detected \
se_#_low_coverage_region_duplicates_marked_for_removal \
se_#_high_coverage_region_bases \
se_#_high_coverage_region_reads \
se_#_high_coverage_region_duplicates_detected \
se_#_high_coverage_region_duplicates_marked_for_removal \
\
bam_file_path\n" );

  if( strrchr( bundle->nameBamFile, '/' ) != NULL )
  {
    cptr = strrchr( bundle->nameBamFile, '/' ) + 1;
  }
  else
  {
    cptr = bundle->nameBamFile;
  }

  fprintf( ofp, "SUMMARY: %s %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %.6f %ld %ld %.6f %.4f %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %.6f %.4f %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %s\n",
           cptr,
           counts[0][C_NUM_ENTRY_BAM],
           counts[0][C_UNALIGNED_BAM],
           counts[0][C_CHIMERIC_BAM],
           counts[0][C_SECONDARY_BAM],
           counts[0][C_NOT_PAIREDEND_BAM],
           counts[0][C_ACCEPTED_BAM],

           counts[0][C_SLMATCH],
           counts[0][C_SHORTINTRON_SLMATCH],
           counts[0][C_LONGINTRON_SLMATCH],
           counts[0][C_NOT_SLMATCH],

           counts[2][C_START_PAIRED_READ],
           counts[2][C_DUP_START_PAIRED_READ],
           counts[2][C_START_PAIRED_READ] > 0 ? (double)counts[2][C_DUP_START_PAIRED_READ] / (double)counts[2][C_START_PAIRED_READ] : 0.0,

           counts[2][C_START_UNPAIRED_READ],
           counts[2][C_DUP_START_UNPAIRED_READ],
           counts[2][C_START_UNPAIRED_READ] > 0 ? (double)counts[2][C_DUP_START_UNPAIRED_READ] / (double)counts[2][C_START_UNPAIRED_READ] : 0.0, 

           pcrDupRatePE,
           counts[2][C_NOT_HCOV_DUP_DETECT] + counts[2][C_HCOV_DUP_DETECT],
           counts[2][C_NOT_HCOV_DUP_REMOVE] + counts[2][C_HCOV_DUP_REMOVE],
           counts[2][C_NOT_HCOV_BASE],
           counts[2][C_NOT_HCOV_READ],
           counts[2][C_NOT_HCOV_DUP_DETECT],
           counts[2][C_NOT_HCOV_DUP_REMOVE],
           counts[2][C_HCOV_BASE],
           counts[2][C_HCOV_READ],
           counts[2][C_HCOV_DUP_DETECT],
           counts[2][C_HCOV_DUP_REMOVE],

           counts[1][C_START_UNPAIRED_READ],
           counts[1][C_DUP_START_UNPAIRED_READ],
           counts[1][C_START_UNPAIRED_READ] > 0 ? (double)counts[1][C_DUP_START_UNPAIRED_READ] / (double)counts[1][C_START_UNPAIRED_READ] : 0.0,

           pcrDupRateSE,
           counts[1][C_NOT_HCOV_DUP_DETECT] + counts[1][C_HCOV_DUP_DETECT],
           counts[1][C_NOT_HCOV_DUP_REMOVE] + counts[1][C_HCOV_DUP_REMOVE],
           counts[1][C_NOT_HCOV_BASE],
           counts[1][C_NOT_HCOV_READ],
           counts[1][C_NOT_HCOV_DUP_DETECT],
           counts[1][C_NOT_HCOV_DUP_REMOVE],
           counts[1][C_HCOV_BASE],
           counts[1][C_HCOV_READ],
           counts[1][C_HCOV_DUP_DETECT],
           counts[1][C_HCOV_DUP_REMOVE],

           bundle->nameBamFile );
               
  *fstatus = 0;

  return( 0 );
}


static int xdumpCovReg( Genome *genome, int *fstatus )
{
  int ichr, nchr;
  int ibas, nbas;
  int ireg0, ireg1;
  int flag;
  int nstart, ndup;

  uint16_t t1, t2;

  Chromosome *chr;

  nchr = genome->numChr;

  for( ichr = 0; ichr < nchr; ++ichr )
  {
    chr  = &(genome->chr[ichr]);
    nbas = chr->len;
    flag = 0;
    for( ibas = 0; ibas < nbas; ++ibas )
    {
      t1 = chr->flag[ibas] & BF2_B08;
      t2 = chr->flag[ibas] & ( BF2_B09 | BF2_B10 | BF2_B11 );

      if( !flag && t1 && !t2 )
      {
        ireg0 = ibas + 1;
        nstart = chr->cover[0][ibas] + chr->cover[2][ibas];
        ndup   = chr->cover[1][ibas] + chr->cover[3][ibas];
        flag  = 1;
      }
      else
      if( flag && ( !t1 || t2 ) )
      {
        ireg1 = ibas;
        fprintf( stderr, "LOCOREG: %s %d %d  %d %d  %.4f\n", chr->name, ireg0, ireg1, nstart, ndup, nstart > 0 ? (double)ndup / (double)nstart : 0 );
        flag = 0;
      }
      else
      {
        nstart += chr->cover[0][ibas] + chr->cover[2][ibas];
        ndup   += chr->cover[1][ibas] + chr->cover[3][ibas];
      }
    }

    if( flag )
    {
      ireg1 = ibas;
      fprintf( stderr, "LOCOREG: %s %d %d  %d %d  %.4f\n", chr->name, ireg0, ireg1, nstart, ndup, nstart > 0 ? (double)ndup / (double)nstart : 0 );
    }
  }

  *fstatus = 0;

  return( 0 );
}


static int xdumpCovNonReg( Genome *genome, int *fstatus )
{
  int ichr, nchr;
  int ibas, nbas;
  int ireg0, ireg1;
  int flag;
  int nstart, ndup;

  uint16_t t1;

  Chromosome *chr;

  nchr = genome->numChr;

  for( ichr = 0; ichr < nchr; ++ichr )
  {
    chr  = &(genome->chr[ichr]);
    nbas = chr->len;
    flag = 0;
    for( ibas = 0; ibas < nbas; ++ibas )
    {
      t1 = chr->flag[ibas] & ( BF2_B08 | BF2_B09 | BF2_B10 | BF2_B11 );

      if( !flag && !t1 )
      {
        ireg0 = ibas + 1;
        nstart = chr->cover[0][ibas] + chr->cover[2][ibas];
        ndup   = chr->cover[1][ibas] + chr->cover[3][ibas];
        flag  = 1;
      }
      else
      if( flag && t1 )
      {
        ireg1 = ibas;
        fprintf( stderr, "NOCOREG: %s %d %d  %d %d  %.4f\n", chr->name, ireg0, ireg1, nstart, ndup, nstart > 0 ? (double)ndup / (double)nstart : 0.0 );
        flag = 0;
      }
      else
      {
        nstart += chr->cover[0][ibas] + chr->cover[2][ibas];
        ndup   += chr->cover[1][ibas] + chr->cover[3][ibas];
      }
    }

    if( flag )
    {
      ireg1 = ibas;
      fprintf( stderr, "NOCOREG: %s %d %d  %d %d  %.4f\n", chr->name, ireg0, ireg1, nstart, ndup, nstart > 0 ? (double)ndup / (double)nstart : 0.0 );
    }
  }

  *fstatus = 0;

  return( 0 );
}


