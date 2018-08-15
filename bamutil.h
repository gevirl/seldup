/** bamutil.h **/

/*
*|   File: bamutil.h                                                          |*
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
** Version: 20160802
*/

/*
** The master file location is
**
**   whim:/users/bge/src/bamutil/bam2fast.h
**
*/

#define BAM_ADDITIONAL
#define BAM_DSC_CHR_POS
#define INTRON_AS_DISCREPANCY
// #define BAM_DSC_NO_SORT

#ifndef BAMUTIL_H
#define BAMUTIL_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include "bgzf.h"
#include "util.h"


typedef struct
{
  char *name;
  int32_t length;
} BamRefSeq;


typedef struct
{
  char *name;
  char *command;
  char *description;
  char *version;
} BamProgram;


typedef struct
{
  char      *description;  // bamReadHeader stores a NULL-terminated string
  int        numRefSeq;
  BamRefSeq *refSeq;
  BamProgram program;
} BamHeader;


typedef struct
{
  char  tag[2];
  char  type[2];  /* type[1] is B-type array type */
  int   alen;     /* B-type array length */
  void *value;
} BamOptField;


typedef struct
{
  char  type;  /* 0=undefined; 1=int; 2=sub; 3=del */
  int   len;   /* length of tok, if defined */
  int   mxlen; /* allocated length of tok */
  char *tok;
} BamMdTok;


/*
** Defined for BamDsc.base (below).
*/
#define BAM_DSC_REF_A	( 1 << 0 )
#define BAM_DSC_REF_C	( 1 << 1 )
#define BAM_DSC_REF_G	( 1 << 2 )
#define BAM_DSC_REF_T	( 1 << 3 )
#define BAM_DSC_REF_N	( 1 << 4 )
#define BAM_DSC_RED_N	( 1 << 5 )

typedef struct
{
  char     type;    /* 1=sub; 2=del; 3=ins; 4=intron (optional) */
  int32_t  pos;     /* discrepancy start location in read coordinates (deletion: immediately to 5' end of deletion start); 1-based */
  int32_t  len;     /* indel length */
  uint8_t  base;    /* reference base(s) unknown: 0; A: bit 0; C: bit 1; G: bit2; T: bit 3; ref N: bit 4; read N: bit 5 Note: bamGetDsc does not check for Ns not marked as dscs */
#ifdef BAM_DSC_CHR_POS
  int32_t  posChr;  /* pos in the SAM format BAM file is int32_t; 1-based */
#endif
} BamDsc;


typedef struct
{
  int      mx_read_name;
  char    *read_name;  // original read name (with original suffix)
  int      dscore0;    // alignment score
  int      dscore;     // alignment score
  int      pscore;     // alignment score
  int      black_list; // 0=not aligned to black-listed region; 1=aligned to black-listed region
  int      clip5[2];   // soft: clip5[0]; hard: clip5[1]
  int      clip3[2];   // soft: clip3[0]; hard: clip3[1]
  int8_t   hclip;
  int64_t  begRead;    // read implied start in genomic coordinates
  int64_t  endRead;    // read implied end in genomic coordinates
  int64_t  begAlign;   // alignment start in genomic coordinates (same as bamAlign.pos)
  int64_t  endAlign;   // alignment end in genomic coordinates
  int      allocDsc;
  int      numDsc;
  uint8_t  bitFlag;
  int32_t  mx_l_seq;
  char    *seq;
  int8_t   template_end; // template end of read for paired end read set 0=fwd; 1=rev; -1=unknown
  BamDsc  *dsc;
  int32_t  nclpDsc;
  int32_t  nfltDsc;
} BamAdditional;


/*
** Notes:
**   o  the mx_* records give the number of units
**      of allocated memory for the corresponding
**      variable. For example, mx_cigar_op has the
**      number of allocated uint32_t values for
**      *cigar. The point is to allow reusing a BamAlign
**      structure to minimize allocation/freeing.
**
**   o  bamReadAlign expands (reallocs) read_name, cigar, seq, and qual
**      to lengths mx_read_name, mx_cigar_op, and mx_l_seq.
**   o  bamReadAlign expands (reallocs) 'field' but not values stored in
**      each field structure (it mallocs these)
*/
typedef struct
{
  int32_t   refid;
  int32_t   pos;         /* this is POS; that is, starting at base 1 */
  uint32_t  bin;
  int32_t   mapq;        /* mapping quality */
  uint32_t  flag;
  uint32_t  n_cigar_op;  /* number of CIGAR operations */
  int32_t   l_seq;       /* sequence length */
  int32_t   next_refid;
  int32_t   next_pos;
  int32_t   tlen;        /* template length */
  int32_t   mx_read_name;
  char     *read_name;
  int32_t   mx_cigar_op;
  uint32_t *cigar;
  int32_t   mx_l_seq;
  char     *seq;         /* bases - in BAM orientation and packed two bases/byte*/
  char     *qual;        /* quality values as char (not ascii-encoded) */
  int32_t   mx_field;
  int32_t   numOptField;
  BamOptField *field;
#ifdef BAM_ADDITIONAL
  BamAdditional addl;
#endif
} BamAlign;


int bamReadHeader( BGZF *fp, BamHeader *header, int *fstatus );
int bamInitAlign( BamAlign *bamAlign, int *fstatus );
int bamTestCigar( BGZF *fp, int *fcigarFlag, int *fmdFlag, int *fstatus );
int bamReadAlign( BGZF *fp, BamAlign *align, int *fstatus );
int bamDumpSamHeader( BamHeader *header );
int bamDumpSamAlign( BamHeader *header, BamAlign *align );
int bamWriteSamHeader( FILE *fp, BamHeader *header, int *fstatus );
int bamWriteSamAlign( FILE *fp, BamHeader *header, BamAlign *align, int *fstatus );
int bamCigarOp( uint32_t uop, char *op, int32_t *len, int *fstatus );
int bamCigarUop( char cop, int32_t lop, uint32_t *fuop, int *fstatus );
int bamUnpackSeq( BamAlign *align, char **sbuf, int *lbuf, int *fstatus );
int bamTestCigar( BGZF *fp, int *fcigarFlag, int *fmdFlag, int *fstatus );
int bamCountDsc( BamAlign *bamAlign, int cigarFlag, int mdFlag, int *fnumSub, int *fnumDel, int *fnumIns, int *fstatus );
int bamGetDsc( BamAlign *bamAlign, int cigarFlag, int mdFlag, BamDsc **fbamDsc, int *fnumBamDsc, int *fallocBamDsc, int clip5[2], int clip3[2], int *fstatus );
int bamDelNDsc( BamDsc **fbamDsc, int *fnumBamDsc, int *numDelN, int *fstatus );
int bamCalcEndsRead( BamAlign *bamAlign, BamDsc *bamDsc, int numBamDsc, int clip5[2], int clip3[2], int64_t *fbpos, int64_t *fepos, int *fstatus );
int bamCalcEndsAlign( BamAlign *bamAlign, BamDsc *bamDsc, int numBamDsc, int clip5[2], int clip3[2], int64_t *fbpos, int64_t *fepos, int *fstatus );
int bamFreeAlign( BamAlign *align );
int bamMoveAlign( BamAlign *dst, BamAlign *src, int *fstatus );
int getAlignSet( char *nameRead, BGZF *fp, BamAlign **fbamAlign, int *fmxBamAlign, int *fnumBamAlign, int *faflag, int *fstatus );
int strnum_cmp(const void *_a, const void *_b);
int bamReportAlign( char *label, BamAlign *bamAlign, int numBamAlign, char *refSeq, int64_t lenRefSeq, char *nameRefSeq, FILE *afp, int *fstatus );
int bamSetBufAlign( BamAlign *bamAlign, char *refSeq, int64_t lenRefSeq, char **frbuf, char **fqbuf, char **fdbuf, char **fsbuf, char **ftbuf, int *flbuf, int64_t *fbegExt, int *fstatus );
int bamSetRefBase( BamDsc *bamDsc, int numBamDsc, char *refSeq, int64_t lenRefSeq, int *fstatus );
int bamGetMdTok( char *string, BamMdTok **fmdTok, int *fmxTok, int *fntok, int *fstatus );
#endif
