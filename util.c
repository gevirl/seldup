/** util.c **/

/*
*|   File: util.c                                                             |*
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
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <fenv.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/wait.h>
#include <time.h>
#include <regex.h>
#include <ctype.h>
#include "util.h"



/*******************************************************************************
**
** Program reports
**
*******************************************************************************/

int timeofday( char *string )
{
  time_t now;

  time( &now );
  time2string( now, 0, 1, string );

  return( 0 );
}


int time2string( time_t intime, int flagFormat, int flagTime, char *string )
{
  struct tm *ptm;

  ptm = localtime( &intime );

  if( flagFormat == 0 )
  {
    if( flagTime )
    {
      sprintf( string, "%02d%02d%02d:%02d%02d%02d",
             ptm->tm_year + 1900,
             ptm->tm_mon  + 1,
             ptm->tm_mday,
             ptm->tm_hour,
             ptm->tm_min,
             ptm->tm_sec );
    }
    else
    {
      sprintf( string, "%02d%02d%02d",
             ptm->tm_year + 1900,
             ptm->tm_mon  + 1,
             ptm->tm_mday );
    }
  }
  else
  if( flagFormat == 1 )
  {
    if( flagTime )
    {
      sprintf( string, "%04d-%02d-%02d %02d:%02d:%02d",
             ptm->tm_year + 1900,
             ptm->tm_mon  + 1,
             ptm->tm_mday,
             ptm->tm_hour, 
             ptm->tm_min,
             ptm->tm_sec );
    }
    else
    {
      sprintf( string, "%04d-%02d-%02d",
             ptm->tm_year + 1900,
             ptm->tm_mon  + 1,
             ptm->tm_mday );
    }
  }
  else
  {
    strcpy( string, ctime( &intime ) );
  }

  return( 0 );
}


/*******************************************************************************
**
** Strings
**
*******************************************************************************/

char *mstrcpy( char *string )
{
  int i;
  int numByte;
  char *cptr;

  if( string == NULL )
  {
    return( NULL );
  }

  i = 0;
  while( string[i] != '\0' )
  {
    ++i;
  }

  numByte = ( i + 1 ) * sizeof( char );
  cptr = ( char *)malloc( numByte );
  if( cptr == NULL )
  {
    fprintf( stderr,
             "mstrcpy: unable to allocate memory\n" );
    return( NULL );
  }

  memcpy( cptr, string, numByte );

  return( cptr );
}


/*******************************************************************************
**
** Bioinformatics
**
*******************************************************************************/


char setcmptab( char *cmp )
{
  int i;

  for( i = 0; i < 256; ++i )
  {
    cmp[i] = (char)i;
  }

  /*
  ** Encode the relative peak height with the character 'case'.
  */
  cmp['a'] = 't';
  cmp['c'] = 'g';
  cmp['g'] = 'c';
  cmp['t'] = 'a';
  cmp['u'] = 'a'; /* U = T */

  cmp['m'] = 'K'; /* aC    */
  cmp['k'] = 'M'; /* gT    */
  cmp['r'] = 'Y'; /* aG    */
  cmp['y'] = 'R'; /* cT    */
  cmp['s'] = 'S'; /* cG    */
  cmp['w'] = 'W'; /* aT    */

  cmp['x'] = 'x'; /* ACGT  */
  cmp['n'] = 'n'; /* ACGT  */

  cmp['A'] = 'T';
  cmp['C'] = 'G';
  cmp['G'] = 'C';
  cmp['T'] = 'A';
  cmp['U'] = 'A'; /* U = T */

  cmp['M'] = 'k'; /* Ac    */
  cmp['K'] = 'm'; /* Gt    */
  cmp['R'] = 'y'; /* Ag    */
  cmp['Y'] = 'r'; /* Ct    */
  cmp['S'] = 's'; /* Cg    */
  cmp['W'] = 'w'; /* At    */

  cmp['X'] = 'X'; /* ACGT  */
  cmp['N'] = 'N'; /* ACGT  */

  cmp['-'] = '-';

  return( 0 );
}

/*
** Reverse complement DNA sequence in place.
*/
int rcseq2( char *seq, int *fstatus )
{
  int   i;
  int   len;

  static int initFlag = 0;
  static char cmp[256];

  static int   lbuf = 0;
  static char *cbuf = NULL;

  if( initFlag == 0 )
  {
    setcmptab( cmp );
    initFlag = 1;
  }

  len = strlen( seq );
  if( len >= lbuf )
  {
    lbuf += len + 8192;
    cbuf  = (char *)realloc( cbuf, lbuf * sizeof( char ) );
    if( cbuf == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
  }

  for( i = 0; i < len; ++i )
  {
    cbuf[len-i-1] = cmp[(int)seq[i]];
  }
  memcpy( seq, cbuf, len * sizeof( char ) );
  seq[len] = '\0';

  *fstatus = 0;

  return( 0 );
}


/*******************************************************************************
**
** Memory allocation
**
*******************************************************************************/

int testFree( void *ptr )
{
  if( ptr != NULL )
  {
    free( ptr );
  }

  return( 0 );
}


