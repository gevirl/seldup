/** seldup.h **/

/*
*|   File: seldup.h                                                           |*
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
** Version: 20161117
*/

#define SF_PAIREDEND	0x1
#define SF_PROPER_PAIR	0x2
#define SF_UNMAPPED	0x4
#define SF_REV_CMP	0x10
#define SF_SEC_ALN	0x100
#define SF_SUP_ALN	0x800


#define RS_SINGLE	1
#define RS_PAIRED	2


#define C_NUM_ENTRY_BAM				0
#define C_UNALIGNED_BAM				1
#define C_ALIGNED_BAM				2
#define C_CHIMERIC_BAM				3
#define C_SECONDARY_BAM				4
#define C_ACCEPTED_BAM				5
#define C_PAIREDEND_BAM				6
#define C_NOT_PAIREDEND_BAM			7
#define C_DUP_PAIR_READ				8
#define C_START_PAIRED_READ			9
#define C_DUP_START_PAIRED_READ			10
#define C_START_UNPAIRED_READ			11
#define C_DUP_START_UNPAIRED_READ		12

#define C_PROPERPAIR_BAM			13
#define C_NOT_PROPERPAIR_BAM			14
#define C_SLMATCH				15
#define C_NOT_SLMATCH				16
#define C_SHORTMATCH				17
#define C_NOT_SHORTMATCH			18
#define C_LONGINTRON_SLMATCH			19
#define C_SHORTINTRON_SLMATCH			20

#define C_HCOV_BASE				21
#define C_HCOV_READ				22
#define C_HCOV_DUP_DETECT			23
#define C_HCOV_DUP_REMOVE			24
#define C_NOT_HCOV_BASE				25
#define C_NOT_HCOV_READ				26
#define C_NOT_HCOV_DUP_DETECT			27
#define C_NOT_HCOV_DUP_REMOVE			28


typedef struct
{
  char    *name;
  int64_t  len;
  uint16_t *flag; /* see above */
  int64_t *cover[4];
} Chromosome;


typedef struct
{
  int         numChr;
  Chromosome *chr;
} Genome;


/*
** flag bits
**   1: 0=not rDNA, 1=is rDNA
** score[i][j] indices
**   i: 0=top strand; 1=bottom strand
**   j: 0=low representation value base; 1=high representation value base
** numHiRvBase[i] indices
**   i: 0=top strand; 1=bottom strand; 2=both strands
*/
typedef struct
{
  int64_t  indexExon;   /* GTF file index */
  int64_t  indexChr;    /* BAM chromosome index */
  uint8_t  flag;
  char    *nameExon;
  char    *nameTranscript;
  char    *nameGene;
  char    *nameChr;
  uint8_t  strand;      /* 0=top strand; 1=bottom strand */
  int64_t  beg;
  int64_t  end;
} Exon;


typedef struct
{
  int64_t  indexExon;
  char    *nameChr;
  uint8_t  strand;
  int64_t  beg;
  int64_t  end;

  char    *nameTranscript;
  char    *nameGene;
  int      nexon;

  int64_t  sumTopCov;
  int64_t  sumBotCov;
  int64_t  sumCombCov;

  int64_t  sumTopLoRv;
  int64_t  sumBotLoRv;
  int64_t  sumCombLoRv;

  double   dcpmTop;
  double   dcpmBot;
  double   dcpmComb;
} Transcript
;

typedef struct
{
  int64_t  indexExon;
  char    *nameChr;
  uint8_t  strand;
  int64_t  beg;
  int64_t  end;

  char    *nameGene;
  int      nexon;

  int64_t  sumTopCov;
  int64_t  sumBotCov;
  int64_t  sumCombCov;

  int64_t  sumTopLoRv;
  int64_t  sumBotLoRv;
  int64_t  sumCombLoRv;

  int64_t  sumOverlap;   /* number of overlapping exon bases in gene */
  int64_t  sumNoOverlap; /* number of non-overlapping exon bases in gene */
  int64_t  sumHiRepVal;  /* number of high representation value bases in gene */

  double   dcpmTop;
  double   dcpmBot;
  double   dcpmComb;
} Gene;


