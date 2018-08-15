/** rmdupsam.c **/

/*
*|   File: rmdupsam.c                                                         |*
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


#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"


#define MXTOK	1024


typedef struct
{
  char     *nameRead;
  char     *refSeq;
  uint32_t  flag;
  int32_t   pos;
  char     *key;
} RmDup;


static int xreadRmDupFile( char *nameFile, RmDup **frmDup, int *fnumRmDup, int *fstatus );
static int xcmpSortRmDup( const void *fa, const void *fb );
static int xgetEnd( uint32_t flag, int *fiend, int *fstatus );
static int xprocBam( char *nameFile, RmDup *rmDup, int numRmDup, int *fstatus );


int main( int argc, char **argv )
{
  int status;
  int numRmDup;

  char nameRmDupFile[8192];
  char nameBamFile[8192];

  RmDup *rmDup;

  if( argc == 3 )
  {
    strcpy( nameBamFile, argv[1] );
    strcpy( nameRmDupFile, argv[2] );

  }
  else
  {
    printf( "rmdupsam removes duplicate alignments from a SAM format\n" );
    printf( "input file or from SAM format alignments piped into stdin\n" );
    printf( "Suggested usage:\n" );
    printf( "  samtools view -h <input BAM filename> | rmdupsam - <seldup *.rmdup.lst filename> | samtools view -b -o <output BAM filename>\n" );
    printf( "\n" );

    printf( "enter BAM filename: " );
    gets( nameBamFile );
    printf( "enter rmdup filename: " );
    gets( nameRmDupFile );
  }

  /*
  ** Read rmdup file.
  */
  xreadRmDupFile( nameRmDupFile, &rmDup, &numRmDup, &status );
  IF_STATUS_MAIN( "bad status: xreadRmDupFile" );
  fprintf( stderr, "%d alignments marked for removal\n", numRmDup );

  /*
  ** Sort rmDup by key.
  */
  qsort( rmDup, numRmDup, sizeof( RmDup ), xcmpSortRmDup );

  /*
  ** Process BAM file.
  */
  xprocBam( nameBamFile, rmDup, numRmDup, &status );
  IF_STATUS_MAIN( "bad status: xprocBam" );

  return( 0 );
}


static int xgetEnd( uint32_t flag, int *fiend, int *fstatus )
{
  int iend;

  if( flag & 0x4 && flag & 0x8 )
  {
    iend = 1;
  }
  else
  if( !( flag & 0x4 ) && !( flag & 0x8 ) )
  {
    iend = 0;
  }
  else
  if( flag & 0x4 )
  {
    iend = 1;
  }
  else
  {
    iend = 2;
  }

  *fiend = iend;

  *fstatus = 0;

  return( 0 );
}


static int xcmpSortRmDup( const void *fa, const void *fb )
{
  RmDup *a;
  RmDup *b;

  a = (RmDup *)fa;
  b = (RmDup *)fb;

  return( strcmp( a->key, b->key ) );
}


static int xcmpSearchRmDup( const void *fa, const void *fb )
{
  char  *a;
  RmDup *b;

  a = (char *)fa;
  b = (RmDup *)fb;

  return( strcmp( a, b->key ) );
}


static int xreadRmDupFile( char *nameFile, RmDup **frmDup, int *fnumRmDup, int *fstatus )
{
  int itok;
  int irm, nrm;
  int iend;
  int status;

  size_t mline;

  uint32_t flag;

  char *pline;
  char *cptr;
  char *eptr;
  char *stok[MXTOK];
  char  string[8192];

  RmDup *rmDup;

  FILE *fp;

  fp = fopen( nameFile, "r" );
  if( fp == NULL )
  {
    sprintf( _msg_, "unable to open file %s", nameFile );
    EMSG( _msg_ );
    *fstatus = -1;
    return( 0 );
  }

  mline = (size_t)0;
  pline = NULL;

  rmDup = NULL;
  irm   = 0;
  nrm   = -1;
  while( getline( &pline, &mline, fp ) >= 0 )
  {
    cptr = pline;
    itok = 0;
    while( itok < MXTOK && ( stok[itok] = strtok( cptr, " \n" ) ) != NULL )
    {
      ++itok;
      cptr = NULL;
    }

    if( strcmp( stok[0], "ALIGN:" ) == 0 )
    {
      if( nrm == -1 )
      {
        EMSG( "unexpected condition" );
        *fstatus = -1;
        return( 0 );
      }

      flag = strtoul( stok[8], &eptr, 10 );
      rmDup[irm].nameRead = mstrcpy( stok[10] );
      rmDup[irm].flag     = flag;
      rmDup[irm].refSeq   = mstrcpy( stok[6] );
      rmDup[irm].pos      = atoi( stok[9] );

      xgetEnd( flag, &iend, &status );
      IF_STATUS_ZERO( "bad status: xgetEnd" );
     
      sprintf( string, "%s_%d", rmDup[irm].nameRead, iend );
      rmDup[irm].key = mstrcpy( string );
      ++irm;
    }
    else
    if( strcmp( stok[0], "HEADER:" ) == 0 )
    {
      nrm = atoi( stok[1] );
      rmDup = (RmDup *)calloc( nrm, sizeof( RmDup ) );
      if( rmDup == NULL )
      {
        EMSG( "unable to allocate memory" );
        *fstatus = -1;
        return( 0 );
      }
    }
  }

  fclose( fp );

  if( irm != nrm )
  {
    EMSG( "inconsistent rmdup file entry count" );
    *fstatus = -1;
    return( 0 );
  }

  *frmDup    = rmDup;
  *fnumRmDup = nrm;

  *fstatus = 0;

  return( 0 );
}


static int xprocBam( char *nameFile, RmDup *rmDup, int numRmDup, int *fstatus )
{
  int i;
  int itok;
  int iend;
  int nrm;
  int nline;
  int status;

  char *pline;
  char *cptr;
  char *eptr;
  char *stok[MXTOK];
  char  string[8192];

  uint32_t flag;

  size_t mline;

  RmDup *pdup;

  FILE *fp;

  if( nameFile[0] == '-' && nameFile[1] == '\0' )
  {
    fp = stdin;
  }
  else
  {
     fp = fopen( nameFile, "r" );
     if( fp == NULL )
     {
       sprintf( _msg_, "unable to open file %s", nameFile );
       *fstatus = -1;
       return( 0 );
     }
  }


  nrm   = 0;
  nline = 0;
  mline = (size_t)0;
  pline = NULL;
  while( getline( &pline, &mline, fp ) >= 0 )
  {
    ++nline;
    if( pline[0] == '@' )
    {
      fprintf( stdout, "%s", pline );
      continue;
    }

    cptr = pline;
    itok = 0;
    while( itok < MXTOK && ( stok[itok] = strtok( cptr, "\t\n" ) ) != NULL )
    {
      ++itok;
      cptr = NULL;
    }

    if( itok < 11 )
    {
      sprintf( _msg_, "too few columns %d < 11 in line %d", itok, nline );
      EMSG( _msg_ );
      for( i = 0; i < itok; ++i )
      {
        if( i > 0 ) fprintf( stderr, "\t" );
        fprintf( stderr, "%s", stok[i] );
      }
      fprintf( stderr, "\n" );
      fprintf( stderr, "skip\n" );
      continue;
    }

    flag = strtoul( stok[1], &eptr, 10 );
    xgetEnd( flag, &iend, &status );
    IF_STATUS_ZERO( "bad status: xgetEnd" );

    sprintf( string, "%s_%d", stok[0], iend );
    pdup = bsearch( string, rmDup, numRmDup, sizeof( RmDup ), xcmpSearchRmDup );
    if( pdup )
    {
      ++nrm;
      continue;
    }

    fprintf( stdout, "%s", stok[0] );
    for( i = 1; i < itok; ++i )
    {
      fprintf( stdout, "\t%s", stok[i] );
    }
    fprintf( stdout, "\n" );
  }

  if( fp != stdout )
  {
    fclose( fp );
  }

  fprintf( stderr, "%d lines read from SAM file\n", nline );
  fprintf( stderr, "%d alignments removed from BAM file\n", nrm );

  *fstatus = 0;

  return( 0 );
}

