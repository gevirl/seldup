/** util.h **/

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

#ifndef UTIL_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <fenv.h>
#include <stdint.h>
#include <sys/times.h>
#include <sys/types.h>
#include <time.h>
#include <signal.h>
#include <dirent.h>
#include <unistd.h>
#include <asm/param.h>
#include <math.h>
#include <errno.h>
#include <execinfo.h>


#define ___MXMSG___	8192

static char _msg_[___MXMSG___];


typedef struct
{
  int  pid;        /* process id */
  char name[1024]; /* executable filename in parentheses */
  char state;      /* execution state: "RSDZTW" */
  int  ppid;       /* parent process id */
  int  gpid;       /* process group id */
  int64_t minorFault; /* number of minor faults */
  int64_t majorFault; /* number of major faults */
  int64_t utime;      /* user time (jiffies) */
  int64_t stime;      /* system time (jiffies) */
  int64_t nthread;    /* number of threads */
  int64_t vsiz;       /* process virtual memory size (bytes) */
  int64_t rsiz;       /* process resident memory size (bytes) */
  int64_t iodelay;    /* aggregated block I/O delays (centiseconds) */
} ProcStat;


/*
** MFREE
*/
#define MFREE( a )      { if( (a) != NULL ) { free( a ); } }

#ifndef VEMSG
#define MSG( a, b )     ___XMSG___( a, b, __FILE__, __func__, __LINE__ )
#define EMSG( b )       ___XMSG___( "ERROR", b, __FILE__, __func__, __LINE__ )
#define WMSG( b )       ___XMSG___( "WARN", b, __FILE__, __func__, __LINE__ )
#define IMSG( b )       ___XMSG___( "INFO", b, __FILE__, __func__, __LINE__ )
#define DMSG( b )       ___XMSG___( "DIAG", b, __FILE__, __func__, __LINE__ )
#define RMSG( b )       ___XMSG___( "REPT", b, __FILE__, __func__, __LINE__ )
#else
#define MSG( a, b )     ___XVMSG___( a, b, __FILE__, __func__, __LINE__ )
#define EMSG( b )       ___XVMSG___( "ERROR", b, __FILE__, __func__, __LINE__ )
#define WMSG( b )       ___XVMSG___( "WARN", b, __FILE__, __func__, __LINE__ )
#define IMSG( b )       ___XVMSG___( "INFO", b, __FILE__, __func__, __LINE__ )
#define DMSG( b )       ___XVMSG___( "DIAG", b, __FILE__, __func__, __LINE__ )
#define RMSG( b )       ___XVMSG___( "REPT", b, __FILE__, __func__, __LINE__ )
#endif


static int ___XMSG___( char *typeError, char *messageError, const char nameFile[], const char nameFunc[], int numLine )
{
  fprintf( stderr,
           "%s: %s (%d): %s\n",
           typeError,
           ( strcmp( nameFunc, "main" ) != 0 ) ? nameFunc : nameFile,
           numLine,
           messageError );

  return( 0 );
}

static int ___XVMSG___( char *typeError, char *messageError, const char nameFile[], const char nameFunc[], int numLine )
{
  fprintf( stderr,
           "%s: %s in %s (%d): %s\n",
           typeError,
           nameFunc,
           nameFile,
           numLine,
           messageError );

  return( 0 );
}

#define IF_STATUS_MAIN(...) \
  if( status != 0 ) { \
    int _ist_; \
    char ___estr___[___MXMSG___]; \
    _ist_ = snprintf( ___estr___, ___MXMSG___, __VA_ARGS__ ); \
    if( _ist_ == ___MXMSG___ ) \
    { \
       EMSG( "string exceeds buffer size" ); \
       exit( -1 ); \
    } \
    fprintf( stderr, \
             "ERROR: %s (%d): %s\n", \
             __FILE__, \
             __LINE__, \
             ___estr___ ); \
    return( -1 ); }

#define IF_STATUS_ZERO(...) \
  if( status != 0 ) { \
    int _ist_; \
    char ___estr___[___MXMSG___]; \
    _ist_ = snprintf( ___estr___, ___MXMSG___, __VA_ARGS__ ); \
    if( _ist_ == ___MXMSG___ ) \
    { \
       EMSG( "string exceeds buffer size" ); \
       exit( -1 ); \
    } \
    fprintf( stderr, \
             "ERROR: %s (%d): %s\n", \
             __func__, \
             __LINE__, \
             ___estr___ ); \
    *fstatus = -1; \
    return( 0 ); }

#define IF_STATUS_NULL(...) \
  if( status != 0 ) { \
    int _ist_; \
    char ___estr___[___MXMSG___]; \
    _ist_ = snprintf( ___estr___, ___MXMSG___, __VA_ARGS__ ); \
    if( _ist_ == ___MXMSG___ ) \
    { \
       EMSG( "string exceeds buffer size" ); \
       exit( -1 ); \
    } \
    fprintf( stderr, \
             "ERROR: %s (%d): %s\n", \
             ( strcmp( __func__, "main" ) == 0 ) ? __FILE__ : __func__, \
             __LINE__, \
             ___estr___ ); \
    *fstatus = -1; \
    return( NULL ); }

/*
** Utility functions declarations.
*/

int timeofday( char *string );
int time2string( time_t intime, int flagFormat, int flagTime, char *string );
int runtime( char *format );
int mprocstat( int fpid, ProcStat *fprocStat, int *fstatus );
char *mstrcpy( char *string );
char setcmptab( char *cmp );
int rcseq2( char *seq, int *fstatus );
int testFree( void *ptr );


#endif

#define UTIL_H
