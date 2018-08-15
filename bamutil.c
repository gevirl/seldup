/** bamutil.c **/

/*
*|   File: bamutil.c                                                          |*
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
** Notes:
**   o  strings that depend on the read orientation
**      in the BAM file are stored as is; that is,
**      if the read aligns revcmp; then the sequence
**      remains reverse complemented, and the quality
**      values (in both the QUAL field and the OQ
**      optional field), the CIGAR operations, and
**      the MD field entries remain reversed.
**
*/

#include "bamutil.h"


#define MAX_OPT_FIELD	4096


static int xbamParseHeader( BamHeader *bamHeader, int *fstatus );
static float unpackFloat32( const uint8_t *buffer );
#ifndef BAM_DSC_NO_SORT
static int xcmpSortBamDsc( const void *fa, const void *fb );
#endif


typedef union
{
  int16_t  i16;
  int32_t  i32;
  uint16_t u16;
  uint32_t u32;
  float    f32;
} ValBuf;


static float unpackFloat32( const uint8_t *buffer )
{
  ValBuf valBuf;
  valBuf.u32 = buffer[0] | (buffer[1] << 8) | (buffer[2] << 16) | (buffer[3] << 24);
  return ( valBuf.f32 );
}


int bamReadHeader( BGZF *fp, BamHeader *header, int *fstatus )
{
  int status;

  int32_t iseq;
  int32_t bsiz, rsiz;
  size_t sizBuf;
  uint8_t *buf;

  sizBuf = 16384;
  buf    = (uint8_t *)calloc( sizBuf, sizeof( uint8_t ) );
  if( buf == NULL )
  {
    EMSG( "unable to allocate memory" );
    *fstatus = -1;
    return( 0 );
  }

  /*
  ** Read header information.
  */
  bsiz = 4 * sizeof( char ) + sizeof( int32_t );
  rsiz = bgzf_read( fp, buf, bsiz );
  if( rsiz != bsiz )
  {
    sprintf( _msg_, "bad read: bgzf_read: bsiz: %d != rsiz: %d", bsiz, rsiz );
    EMSG( _msg_ );
    *fstatus = -1;
    return( 0 );
  }

  /*
  ** Check magic string.
  */
  if( buf[0] != 'B' || buf[1] != 'A' || buf[2] != 'M' || buf[3] != 1 )
  {
    EMSG( "not a BAM file" );
    *fstatus = -1;
    return( 0 );
  }

  bsiz = (int32_t)unpackInt32( &(buf[4]) );

  if( bsiz >= sizBuf )
  {
    sizBuf += bsiz;
    buf = (uint8_t *)realloc( buf, sizBuf * sizeof( uint8_t ) );
    if( buf == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
  }

  rsiz = bgzf_read( fp, buf, bsiz );
  if( rsiz != bsiz )
  {
    EMSG( "bad read: bgzf_read" );
    *fstatus = -1;
    return( 0 );
  }

  header->description = calloc( bsiz + 1, sizeof( char ) );
  if( header->description == NULL )
  {
    EMSG( "unable to allocate memory" );
    *fstatus = -1;
    return( 0 );
  }

  memcpy( header->description, buf, bsiz * sizeof( char ) );
  header->description[bsiz] = '\0';

  /*
  ** Read reference information.
  */
  bsiz = sizeof( int32_t );
  rsiz = bgzf_read( fp, buf, bsiz );
  if( rsiz != bsiz )
  {
    EMSG( "bad read: bgzf_read" );
    *fstatus = -1;
    return( 0 );
  }

  header->numRefSeq = (int32_t)unpackInt32( buf );
  header->refSeq    = (BamRefSeq *)calloc( header->numRefSeq, sizeof( BamRefSeq ) );
  if( header->refSeq == NULL )
  {
    EMSG( "unable to allocate memory" );
    *fstatus = -1;
    return( 0 );
  }

  for( iseq = 0; iseq < header->numRefSeq; ++iseq )
  {
    bsiz = sizeof( int32_t );
    rsiz = bgzf_read( fp, buf, bsiz );
    if( rsiz != bsiz )
    {
      EMSG( "bad read: bgzf_read" );
      *fstatus = -1;
      return( 0 );
    }

    bsiz = (int32_t)unpackInt32( buf );    
    if( bsiz >= sizBuf )
    {
      sizBuf += bsiz;
      buf = (uint8_t *)realloc( buf, sizBuf * sizeof( uint8_t ) );
      if( buf == NULL )
      {
        EMSG( "unable to allocate memory" );
        *fstatus = -1;
        return( 0 );
      }
    }
    rsiz = bgzf_read( fp, buf, bsiz );
    if( rsiz != bsiz )
    {
      EMSG( "bad read: bgzf_read" );
      *fstatus = -1;
      return( 0 );
    }
    header->refSeq[iseq].name = mstrcpy( (char *)buf );
    bsiz = sizeof( int32_t );
    rsiz = bgzf_read( fp, buf, bsiz );
    if( rsiz != bsiz )
    {
      EMSG( "bad read: bgzf_read" );
      *fstatus = -1;
      return( 0 );
    }

    header->refSeq[iseq].length = (int32_t)unpackInt32( buf );
  }

  free( buf );

  xbamParseHeader( header, &status );
  IF_STATUS_ZERO( "bad status: xbamParseHeader" );

  *fstatus = 0;

  return( 0 );
}


static int xbamParseHeader( BamHeader *bamHeader, int *fstatus )
{
  int i, j;
  int lbuf;
  int irec;

  char *pbuf;
  char *vbuf;

  bamHeader->program.name        = NULL;
  bamHeader->program.command     = NULL;
  bamHeader->program.description = NULL;
  bamHeader->program.version     = NULL;

  lbuf = strlen( bamHeader->description ) + 1;
  pbuf = (char *)malloc( lbuf * sizeof( char ) );
  if( pbuf == NULL )
  {
    EMSG( "unable to allocate memory" );
    *fstatus = -1;
    return( 0 );
  }
  memcpy( pbuf, bamHeader->description, lbuf );

  vbuf = (char *)malloc( lbuf * sizeof( char ) );
  if( vbuf == NULL )
  {
    EMSG( "unable to allocate memory" );
    *fstatus = -1;
    return( 0 );
  }

  irec = 0;
  for( i = 0; pbuf[i] != '\0'; ++i )
  {
    if( pbuf[i] == '@' )
    {
      if( strncmp( &(pbuf[i]), "@PG", 3 ) == 0 )
      {
        irec = 1;
        i += 2;
      }
      else
      {
        for( ; pbuf[i] != '\0' && pbuf[i] != '\n'; ++i );
        if( pbuf[i] == '\0' ) goto finished;
        irec = 0;
      }
    }
    else
    if( pbuf[i] == '\t' )
    {
      j = 0;
      for( ++i; pbuf[i] != '\0' && pbuf[i] != '\t' && pbuf[i] != '\n'; ++i  )
      {
        vbuf[j] = pbuf[i];
        ++j;
      }
      vbuf[j] = '\0';

      if( pbuf[i] == '\t' ) --i;

      if( irec == 1 )
      {
        if( strncmp( vbuf, "PN:", 3 ) == 0 )
        {
          bamHeader->program.name = mstrcpy( &(vbuf[3]) );
        }
        else
        if( strncmp( vbuf, "CL:", 3 ) == 0 )
        {
          bamHeader->program.command = mstrcpy( &(vbuf[3]) );
        }
        else
        if( strncmp( vbuf, "DS:", 3 ) == 0 )
        {
          bamHeader->program.description = mstrcpy( &(vbuf[3]) );
        }
        if( strncmp( vbuf, "VN:", 3 ) == 0 )
        {
          bamHeader->program.version = mstrcpy( &(vbuf[3]) );
        }
        if( pbuf[i] == '\0' ) goto finished;
      }
    }
  }

  finished:

  free( pbuf );

  *fstatus = 0;

  return( 0 );
}


/*
** bamTestCigar tests whether or not the CIGAR string encodes mismatches
** and whether or not the alignments appear to contain consistently the
** MD field for describing mismatches. Sets these flags by reading the
** first 10000 alignments.
**
**   cigarFlag     1        CIGAR codes have [=X]
**   cigarFlag     0        CIGAR codes have [M]
**   cigarFlag    -1        CIGAR codes have [=XM] (contradictory)
**   mdFlag        1        MD flags appear in all tested entries
**   mdFlag        0        MD flags appear in none of the tested entries
**   mdFlag       -1        MD flags appear in some of the tested entries but not all
*/
int bamTestCigar( BGZF *fp, int *fcigarFlag, int *fmdFlag, int *fstatus )
{
  int status;
  int istat;
  int itst, ntst;
  int cigarFlag;
  int mdFlag;
  int iop, nop;
  int ifld, nfld;
  int nm, ne, nx, nmd;
  int op;

  uint32_t uval;

  int64_t fpos;
  int64_t lstat;

  BamAlign align;
  BamOptField *field;

  ntst = 10000;

  fpos = bgzf_tell( fp );
  if( fpos < 0 )
  {
    EMSG( "bad status: bgzf_tell" );
    *fstatus = -1;
    return( 0 );
  }

  bamInitAlign( &align, &status );
  IF_STATUS_ZERO( "bad status: bamInitAlign" );

  nm   = 0;
  ne   = 0;
  nx   = 0;
  nmd  = 0;
  for( itst = 0; itst < ntst; )
  {
    istat = bamReadAlign( fp, &align, &status );
    IF_STATUS_ZERO( "bad status: bamReadAlign" );

    if( istat != 1 )
    {
      break;
    }

    if( align.flag & 0x4 )
    {
      continue;
    }

    nop = align.n_cigar_op;
    for( iop = 0; iop < nop; ++iop )
    {
      uval = align.cigar[iop];
      op  = (int)( uval & 0xf );
      if( op == 0 ) ++nm;  /* M */
      else
      if( op == 7 ) ++ne;  /* = */
      else
      if( op == 8 ) ++nx;  /* X */
    }

    nfld  = align.numOptField;
    field = align.field;
    for( ifld = 0; ifld < nfld; ++ifld )
    {
      if( field[ifld].tag[0] == 'M' && field[ifld].tag[1] == 'D' )
      {
        ++nmd;
      }
    }
    ++itst;
  }

  if( nm && !ne && !nx )
  {
    cigarFlag = 0;
  }
  else
  if( !nm && ( ne || nx ) )
  {
    cigarFlag = 1;
  }
  else
  {
    cigarFlag = -1;
  }

  if( nmd == 0 )
  {
    mdFlag = 0;
  }
  else
  if( nmd && nmd == itst )
  {
    mdFlag = 1;
  }
  else
  {
    mdFlag = -1;
  }

  bamFreeAlign( &align );

  testFree( align.read_name );
  testFree( align.cigar );
  testFree( align.seq );
  testFree( align.qual );

  lstat = bgzf_seek( fp, fpos, SEEK_SET );
  if( lstat < 0 )
  {
    EMSG( "bad status: bgzf_tell" );
    *fstatus = -1;
    return( 0 );
  }

  *fcigarFlag = cigarFlag;
  *fmdFlag    = mdFlag;

  *fstatus = 0;

  return( 0 );
}


int bamReadAlign( BGZF *fp, BamAlign *align, int *fstatus )
{
  int i;
  int32_t numOptField;
  int32_t bsiz, rsiz, msiz;
  int32_t len, lenNam;

  uint32_t ubuf;
  uint16_t usbuf;

  uint8_t *pbuf;
  char     tag[2];
  char     val_type;

  static int sizBuf   = 0;
  static uint8_t *buf = NULL;
  static BamOptField optField[MAX_OPT_FIELD];


  bsiz = 1024;
  if( bsiz >= sizBuf )
  {
    sizBuf += bsiz;
    buf = (uint8_t *)realloc( buf, sizBuf * sizeof( uint8_t ) );
    if( buf == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
  }

  /*
  ** Read alignment record size.
  ** Note: bgzf_read returns 0 on EOF.
  */
  bsiz = sizeof( int32_t );
  rsiz = bgzf_read( fp, buf, bsiz );

  if( rsiz == 0 )
  {
    *fstatus = 0;
    return( 0 );
  }

  if( rsiz != bsiz )
  {
    EMSG( "bad read: bgzf_read" );
    *fstatus = -1;
    return( 0 );
  }

  bsiz = (int32_t)unpackInt32( buf );

  msiz = bsiz;
  if( bsiz >= sizBuf )
  {
    sizBuf += bsiz;
    buf = (uint8_t *)realloc( buf, sizBuf * sizeof( uint8_t ) );
    if( buf == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
  }

  /*
  ** Read alignment record.
  */
  rsiz = bgzf_read( fp, buf, bsiz );
  if( rsiz != bsiz ) 
  {
    EMSG( "bad read: bgzf_read" );
    *fstatus = -1;
    return( 0 );
  } 

  pbuf = buf;

  align->refid = (int32_t)unpackInt32( pbuf );

  /*
  ** Note: in the BAM file, pos is 0-based, but I store
  **       it as 1-based.
  */
  pbuf += 4;
  align->pos = (int32_t)unpackInt32( pbuf ) + 1;

  pbuf += 4;
  ubuf = (uint32_t)unpackUInt32( pbuf );
  usbuf = ( ubuf >> 16 ) & 0xffff;
  align->bin  = (uint32_t)usbuf;
  align->mapq = (int32_t)( ( ubuf >>  8 ) & 0xff );
  lenNam      = (int32_t)( ubuf & 0xff );

  pbuf += 4;
  ubuf = (uint32_t)unpackUInt32( pbuf );
  usbuf = ( ubuf >> 16 ) & 0xffff;
  align->flag = (uint32_t)usbuf;
  usbuf = ubuf & 0xffff;
  align->n_cigar_op = (uint32_t)usbuf;

  pbuf += 4;
  align->l_seq = (int32_t)unpackInt32( pbuf );

  pbuf += 4;
  align->next_refid = (int32_t)unpackInt32( pbuf );

  /*
  ** Note: in the BAM file, next_pos is 0-based and I store
  **       it as 1-based.
  */
  pbuf += 4;
  align->next_pos = (int32_t)unpackInt32( pbuf ) + 1;

  pbuf += 4;
  align->tlen = (int32_t)unpackInt32( pbuf );

  if( lenNam > align->mx_read_name )
  {
    align->mx_read_name += lenNam;
    align->read_name = (char *)realloc( align->read_name, align->mx_read_name * sizeof( char ) );
    if( align->read_name == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
  }
  pbuf += 4;

  memcpy( align->read_name, pbuf, lenNam * sizeof( char ) );

  pbuf += lenNam * sizeof( char );

  if( align->n_cigar_op > align->mx_cigar_op )
  {
    align->mx_cigar_op += align->n_cigar_op;
    align->cigar = (uint32_t *)realloc( align->cigar, align->mx_cigar_op * sizeof( uint32_t ) );
    if( align->cigar == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
  }

  if( align->n_cigar_op > 0 )
  {
    int iop;
    for( iop = 0; iop < align->n_cigar_op; ++iop )
    {
      align->cigar[iop] = unpackUInt32( &(pbuf[iop*4]) );
    }
    pbuf += align->n_cigar_op * sizeof( uint32_t );
  }

  if( align->l_seq >= align->mx_l_seq )
  {
    align->mx_l_seq += align->l_seq + 1;
    align->seq = (char *)realloc( align->seq, align->mx_l_seq * sizeof( char ) );
    if( align->seq == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }

    align->qual = (char *)realloc( align->qual, align->mx_l_seq * sizeof( char ) );
    if( align->qual == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
  }

  /*
  ** Note: the sequence is packed two bases to a byte.
  */
  memcpy( align->seq, pbuf, ( align->l_seq + 1 ) / 2 * sizeof( uint8_t ) );
  pbuf += ( align->l_seq + 1 ) / 2 * sizeof( uint8_t );

  /*
  ** Note: the quality values are chars rather than encoded as alphabetics
  **       so one must add 33 to them in order to store them as ASCII chars.
  */
  memcpy( align->qual, pbuf, align->l_seq * sizeof( char ) );
  pbuf += align->l_seq * sizeof( char );

  /*
  ** Extract optional fields here.
  ** Notes:
  **   o  MD Z String for mismatching positions. Regex : [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)* (SAM)
  */
  numOptField = 0;
  while( pbuf < buf + msiz )
  {
    if( ( numOptField + 1 ) == MAX_OPT_FIELD )
    {
      EMSG( "optional field buffer too small" );
      *fstatus = -1;
      return( 0 );
    }
    memcpy( tag, pbuf, 2 * sizeof( char ) );
    pbuf += 2;
    val_type = (char)*pbuf;
    ++pbuf;
 // fprintf( stdout, "fop: tag: %.2s  type: %c\n", tag, val_type );
    switch( val_type )
    {
      case 'A':
      /* printable character: [!-~] */
        memcpy( &(optField[numOptField].tag), &tag, 2 * sizeof( char ) );
        optField[numOptField].type[0] = 'A';
        optField[numOptField].value   = malloc( sizeof( char ) );
        *((char *)optField[numOptField].value) = *((char *)pbuf);
        ++numOptField;
        ++pbuf;
      break;

      case 'c':
      /* int8_t */
        memcpy( &(optField[numOptField].tag), &tag, 2 * sizeof( char ) );
        optField[numOptField].type[0] = 'c';
        optField[numOptField].value   = malloc( sizeof( int8_t ) );
        *((int8_t *)optField[numOptField].value) = *((int8_t *)pbuf);
        ++numOptField;
        ++pbuf;
      break;

      case 'C':
      /* uint8_t */
        memcpy( &(optField[numOptField].tag), &tag, 2 * sizeof( char ) );
        optField[numOptField].type[0] = 'C';
        optField[numOptField].value   = malloc( sizeof( uint8_t ) );
        *((uint8_t *)optField[numOptField].value) = *((uint8_t *)pbuf);
        ++numOptField;
        ++pbuf;
      break;

      case 's':
      /* int16_t */
        memcpy( &(optField[numOptField].tag), &tag, 2 * sizeof( char ) );
        optField[numOptField].type[0] = 's';
        optField[numOptField].value   = malloc( sizeof( int16_t ) );
        *((int16_t *)optField[numOptField].value) = unpackInt16( pbuf );
        ++numOptField;
        pbuf += 2;
      break;

      case 'S':
      /* uint16_t */
        memcpy( &(optField[numOptField].tag), &tag, 2 * sizeof( char ) );
        optField[numOptField].type[0] = 'S';
        optField[numOptField].value   = malloc( sizeof( uint16_t ) );
        *((uint16_t *)optField[numOptField].value) = unpackUInt16( pbuf );
        ++numOptField;
        pbuf += 2;
      break;

      case 'i':
      /* int32_t */
        memcpy( &(optField[numOptField].tag), &tag, 2 * sizeof( char ) );
        optField[numOptField].type[0] = 'i';
        optField[numOptField].value   = malloc( sizeof( int32_t ) );
        *((int32_t *)optField[numOptField].value) = unpackInt32( pbuf );
        ++numOptField;
        pbuf += 4;
      break;

      case 'I':
      /* uint32_t */
        memcpy( &(optField[numOptField].tag), &tag, 2 * sizeof( char ) );
        optField[numOptField].type[0] = 'I';
        optField[numOptField].value   = malloc( sizeof( uint32_t ) );
        *((uint32_t *)optField[numOptField].value) = unpackUInt32( pbuf );
        ++numOptField;
        pbuf += 4;
      break;

      case 'f':
      /* single precision float: [-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)? */
        memcpy( &(optField[numOptField].tag), &tag, 2 * sizeof( char ) );
        optField[numOptField].type[0] = 'f';
        optField[numOptField].value   = malloc( sizeof( float ) );
        *((float *)optField[numOptField].value) = unpackFloat32( pbuf );
        ++numOptField;
        pbuf += 4;
      break;

      case 'H':
      /* byte array in Hex format: [0-9A-F]+ */
        memcpy( &(optField[numOptField].tag), &tag, 2 * sizeof( char ) );
        optField[numOptField].type[0] = 'H';
        optField[numOptField].value   = (void *)mstrcpy( (char *)pbuf );
        for( ; *pbuf != '\0'; ++pbuf );
// fprintf( stdout, "  string: %s", (char *)optField[numOptField].value );
        ++pbuf;
        ++numOptField;
      break;

      case 'Z':
      /* printable string, including space: [ !-~]+ (NULL terminated - no length specifier */
        memcpy( &(optField[numOptField].tag), &tag, 2 * sizeof( char ) );
        optField[numOptField].type[0] = 'Z';
        optField[numOptField].value   = (void *)mstrcpy( (char *)pbuf );
        for( ; *pbuf != '\0'; ++pbuf );
// fprintf( stdout, "  string: %s", (char *)optField[numOptField].value );
        ++pbuf;
        ++numOptField;
      break;

      case 'B':
      /* integer or numeric array: [cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+ */
        memcpy( &(optField[numOptField].tag), &tag, 2 * sizeof( char ) );
        optField[numOptField].type[0] = 'B';
        val_type = *pbuf;
        optField[numOptField].type[1] = val_type;
        ++pbuf;
        len = (int32_t)unpackInt32( pbuf );
        optField[numOptField].alen = len;
        pbuf += 4;

        /*
        ** A 'B'-typed (array) tag-value pair is stored as follows. The
        ** first two bytes keep the two-character tag. The 3rd byte is
        ** always 'B'. The 4th byte, matching /^[cCsSiIf]$/, indicates
        ** the type of an element in the array. Bytes from 5 to 8
        ** encode a little-endian 32-bit integer which gives the number
        ** of elements in the array. Bytes starting from the 9th store
        ** the array in the little-endian byte order; the number of
        ** these bytes is determined by the type and the length of the
        ** array.
        ** An integer may be stored as one of `cCsSiI' in BAM,
        ** representing int8_t, uint8_t, int16_t, uint16_t, int32_t,
        ** and uint32_t, respectively. In SAM, all single (i.e.,
        ** non-array) integer types are stored as `i', regardless of
        ** magnitude.
        */

        switch( val_type )
        {
          case 'c':
            optField[numOptField].value = malloc( len * sizeof( int8_t ) );
            for( i = 0; i < len; ++i )
              ((int8_t *)optField[numOptField].value)[i] = ((int8_t *)pbuf)[i];
            pbuf += len;
          break;
  
          case 'C':
            optField[numOptField].value = malloc( len * sizeof( uint8_t ) );
            for( i = 0; i < len; ++i )
              ((uint8_t *)optField[numOptField].value)[i] = ((uint8_t *)pbuf)[i];
            pbuf += len;
          break;
  
          case 's':
            optField[numOptField].value = malloc( len * sizeof( int16_t ) );
            for( i = 0; i < len; ++i )
              ((int16_t *)optField[numOptField].value)[i] = unpackInt16( &(pbuf[i*2]) );
            pbuf += len * 2;
          break;

          case 'S':
            optField[numOptField].value = malloc( len * sizeof( uint16_t ) );
            for( i = 0; i < len; ++i )
              ((uint16_t *)optField[numOptField].value)[i] = unpackUInt16( &(pbuf[i*2]) );
            pbuf += len * 2;
          break;
  
          case 'i':
            optField[numOptField].value = malloc( len * sizeof( int32_t ) );
            for( i = 0; i < len; ++i )
              ((int32_t *)optField[numOptField].value)[i] = unpackInt32( &(pbuf[i*4]) );
            pbuf += len * 4;
          break;
  
          case 'I':
            optField[numOptField].value = malloc( len * sizeof( uint32_t ) );
            for( i = 0; i < len; ++i )
              ((uint32_t *)optField[numOptField].value)[i] = unpackUInt32( &(pbuf[i*4]) );
            pbuf += len * 4;
          break;
  
          case 'f':
            optField[numOptField].value = malloc( len * sizeof( float ) );
            for( i = 0; i < len; ++i )
              ((float *)optField[numOptField].value)[i] = unpackFloat32( &(pbuf[i*4]) );
            pbuf += len * 4;
          break;
  
          default:
            sprintf( _msg_, "unrecognized optional field type '%c'", val_type );
            EMSG( _msg_ );
            *fstatus = -1;
            return( 0 );
        }

        ++numOptField;
      break;

      default:
        sprintf( _msg_, "unrecognized optional field type '%c'", val_type );
        EMSG( _msg_ );
        *fstatus = -1;
        return( 0 );
    }
// fputc( (int)'\n', stdout );
  }

  if( numOptField > align->mx_field )
  {
    align->field = (BamOptField *)realloc( align->field, numOptField * sizeof( BamOptField ) );
    if( align->field == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
    align->mx_field = numOptField;
  }
  memcpy( align->field, optField, numOptField * sizeof( BamOptField ) );
  align->numOptField = numOptField;

 *fstatus = 0;

  return( 1 );
}


int bamInitAlign( BamAlign *align, int *fstatus )
{

  align->mx_read_name = 0;
  align->read_name    = NULL;
  align->mx_cigar_op  = 0;
  align->cigar        = NULL;
  align->mx_l_seq     = 0;
  align->seq          = NULL;
  align->qual         = NULL;
  align->mx_field     = 0;
  align->field        = NULL;

  align->refid        = 0;
  align->pos          = 0;
  align->bin          = 0;
  align->mapq         = 0;
  align->flag         = 0;
  align->n_cigar_op   = 0;
  align->l_seq        = 0;
  align->next_refid   = 0;
  align->next_pos     = 0;
  align->tlen         = 0;
  align->l_seq        = 0;
  align->l_seq        = 0;
  align->l_seq        = 0;
  align->numOptField  = 0;

#ifdef BAM_ADDITIONAL
  align->addl.mx_read_name  = 0;
  align->addl.read_name     = NULL;
  align->addl.pscore        = 0;
  align->addl.dscore        = 0;
  align->addl.black_list    = 0;
  align->addl.clip5[0]      = 0;
  align->addl.clip5[1]      = 0;
  align->addl.clip3[0]      = 0;
  align->addl.clip3[1]      = 0;
  align->addl.hclip         = 0;
  align->addl.begRead       = 0;
  align->addl.endRead       = 0;
  align->addl.allocDsc      = 0;
  align->addl.numDsc        = 0;
  align->addl.dsc           = NULL;
  align->addl.bitFlag       = 0;
  align->addl.template_end  = -1;
  align->addl.mx_l_seq      = 0;
  align->addl.seq           = NULL;
#endif

  *fstatus = 0;

  return( 0 );
}


int bamFreeAlign( BamAlign *align )
{
  int32_t i;

/*
** These are realloced from call-to-call of bamReadAlign.
**
  testFree( align->read_name );
  testFree( align->cigar );
  testFree( align->seq );
  testFree( align->qual );
*/

  for( i = 0; i < align->numOptField; ++i )
  {
    testFree( align->field[i].value );
    align->field[i].value = NULL;
  }

  align->numOptField = 0;

  return( 0 );
}


int bamCigarOp( uint32_t uop, char *op, int32_t *len, int *fstatus )
{
  static char cop[] = "MIDNSHP=X";

  *op  = cop[(int)(uop & 0xf)];
  *len = (int32_t)( ( uop >> 4 ) & 0xfffffff );

  *fstatus = 0;

  return( 0 );
}


int bamCigarUop( char cop, int32_t lop, uint32_t *fuop, int *fstatus )
{
  static int initFlag = 0;
  static uint32_t umap[256];

  if( !initFlag )
  {
    memset( umap, 0, 256 * sizeof( uint32_t ) );
    umap[(int)'M'] = (uint32_t)0;
    umap[(int)'I'] = (uint32_t)1;
    umap[(int)'D'] = (uint32_t)2;
    umap[(int)'N'] = (uint32_t)3;
    umap[(int)'S'] = (uint32_t)4;
    umap[(int)'H'] = (uint32_t)5;
    umap[(int)'P'] = (uint32_t)6;
    umap[(int)'='] = (uint32_t)7;
    umap[(int)'X'] = (uint32_t)8;
    initFlag = 1;
  }

  /*
  ** cigar CIGAR: op len<<4|op. `MIDNSHP=X'!`012345678' uint32 t[n cigar op]
  */
  *fuop = ( (uint32_t)lop << 4 ) | umap[(int)cop];

  *fstatus = 0;

  return( 0 );
}


int bamUnpackSeq( BamAlign *align, char **sbuf, int *lbuf, int *fstatus )
{
  int i, j;

  /*
  ** Is the buffer sufficiently long?
  */
  if( *sbuf == NULL || ( align->l_seq + 1 ) > *lbuf )
  {
    *lbuf = align->l_seq + 1024;
    *sbuf = (char *)realloc( *sbuf, (size_t)*lbuf * sizeof( char ) );
    if( *sbuf == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
  }

  /*
  ** Unpack sequence into buffer.
  */
  for( i = 0, j = 0; i < ( align->l_seq + 1 ) / 2; ++i )
  {
    (*sbuf)[j] = "=ACMGRSVTWYHKDBN"[(int)( ( (align->seq[i]) >> 4 ) & 0xf )];
    ++j;
    if( j >= align->l_seq )
    {
      break;
    }
    (*sbuf)[j] = "=ACMGRSVTWYHKDBN"[(int)(align->seq[i] & 0xf)];
    ++j;
  }

  (*sbuf)[align->l_seq] = '\0';

  *fstatus = 0;

  return( 0 );
}


int bamCountDsc( BamAlign *bamAlign, int cigarFlag, int mdFlag, int *fnumSub, int *fnumDel, int *fnumIns, int *fstatus )
{
  int i;
  int iop, nop;
  int len, dlen;
  int flag;
  int idel;
  int nsub, ndel, nins;
  int status;

  char  cop;
  char *sptr;

  static int  mdel = 0;
  static int *dlidx;

  if( bamAlign->flag & 0x4 )
  {
    *fnumSub = 0;
    *fnumDel = 0;
    *fnumIns = 0;

    *fstatus = 0;
    return( 0 );
  }

  if( mdel == 0 )
  {
    mdel = bamAlign->l_seq;
    dlidx = (int *)malloc( mdel * sizeof( int ) );
    if( dlidx == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
  }

  /*
  ** Count deletions and insertions and record deletion indices.
  */
  nsub = 0;
  ndel = 0;
  nins = 0;
  nop  = bamAlign->n_cigar_op;
  for( iop = 0; iop < nop; ++iop )
  {
    bamCigarOp( bamAlign->cigar[iop], &cop, &len, &status );
    if( cop == 'D' )
    {
      if( ndel + 1 >= mdel )
      {
        mdel += bamAlign->l_seq;
        dlidx = (int *)realloc( dlidx, mdel * sizeof( int ) );
        if( dlidx == NULL )
        {
          EMSG( "unable to allocate memory" );
          *fstatus = -1;
          return( 0 );
        }
      }
      dlidx[ndel] = iop;
      ++ndel;
    }
    else
    if( cop == 'I' )
    {
      ++nins;
    }
    else
    if( cop == 'X' )
    {
      nsub += len;
    }
  }

  /*
  ** Count substitutions if there is an MD field.
  */
  if( mdFlag )
  {
    flag = 0;
    nop  = bamAlign->numOptField;
    for( iop = 0; iop < nop; ++iop )
    {
      if( bamAlign->field[iop].tag[0] == 'M' && bamAlign->field[iop].tag[1] == 'D' )
      {
        sptr = (char *)( bamAlign->field[iop].value );
        flag = 1;
      }
    }

    if( flag == 0 )
    {
      sprintf( _msg_, "missing MD field for BAM file entry %s with flag %u continuing...", bamAlign->read_name, bamAlign->flag );
      EMSG( _msg_ );
      *fnumSub = nsub;
      *fnumDel = ndel;
      *fnumIns = nins;
      *fstatus = 0;
      return( 0 );
    }

    idel = 0;
    len  = strlen( sptr );
    for( i = 0; i < len; ++i )
    {
      if( sptr[i] >= (int)'A' && sptr[i] <= (int)'Z' )
      {
        ++nsub;
      }
      else
      if( sptr[i] == '^' )
      {
        bamCigarOp( bamAlign->cigar[dlidx[idel]], &cop, &dlen, &status );
        i += dlen;
        ++idel;
      }
    }
  }

  *fnumSub = nsub;
  *fnumDel = ndel;
  *fnumIns = nins;

  *fstatus = 0;

  return( 0 );
}


/*
**   position in read:
**     o  deletion: upstream base
**     o  insertion: base
**     o  substitution: base
**
** Notes:
**   by default, the discrepancies are sorted (qsort) by
**     primary: indel, then substitutions
**     secondary: position in read
**   (this is the original behavior before allowing for CIGAR strings with '='
**   and without MD field)
**   turn off sorting by defining preprocessor definition 'BAM_DSC_NO_SORT'
**
**   from the sam format document
**
**   MD Z String for mismatching positions. Regex : [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*6
**   The MD field aims to achieve SNP/indel calling without looking at the reference.
**   For example, a string `10A5^AC6' means from the leftmost reference base in the
**   alignment, there are 10 matches followed by an A on the reference which differs
**   from the aligned read base; the next 5 reference bases are matches followed by
**   a 2bp deletion from the reference; the deleted sequence is AC; the last 6 bases
**   are matches. The MD field ought to match the CIGAR string.
**
**   commentators observe
**     o  the MD string cannot represent unambiguously a deletion followed by a
**        substitution
**     o  some notice zeros included between adjacent substitutions and between
**        deletions followed by substitutions. They refer to the bwa source code.
**        URL: http://seqanswers.com/forums/showthread.php?t=8978
**          Should have looked before I posted my question but nevertheless this
**          post could still be useful for someone else.
**
**            30M8D6M 27T2^ATGCATTT0G3T1
**
**          There is a 0 separating the 8 deletions and the single mismatch after
**          the deletions.
**
**          I haven't found why they are necessary, but sometimes it helps to have
**          them visually. They generally occur between SNPs, or between a deletion
**          then a SNP.
**
**          For example, "5^AC0C5" with a cigar "5M2D6M", or "5A0C5" with a cigar
**          "12M". In the former it is easy to see where the deletion ends (the 0)
**          and the next base (a C SNP) starts.
**
**          The SAMtools code puts them in, so others follow the same lead. You
**          could ask the samtools help list.
**     o  some say that clipping is not recorded in the MD field
**        URL: http://www.biostars.org/p/7331/
**          another caveat is that with most short read aligners we are still
**          talking about a local alignment, so what are ostensibly mismatches
**          (and matches) at either end of the query sequence are instead soft
**          clipped positions (S). This may or may not be how you want to score
**          it, so rebuilding the alignment is often a necessary step.
**     o  some say that one needs to use both the CIGAR string and the MD field
**        in order to correctly identify the alignment discrepancies
**   my observations
**     o  do not count clipped bases as errors
**
** Warnings:
**   o  discrepancy positions are in read coordinates in alignment orientation!
**   o  clipped bases are not counted as discrepancies
**   o  indel sequence is not stored in the BamDsc->base (only substitutions are stored)
**
** Tested:
** d: OK
**   DSC: [53=1D48=] (3)  50 2:53:x:1  INDEL 95803963
** i: OK
**   DSC: [14M1I35M] (3)  49 3:15:x:1  INDEL 105355457/1
** d->s: OK
**   DSC: [13=1D56=1X31=] (5)  50 2:13:x:1 1:70:x:1  INDEL 109634975
** i->s: OK
**   DSC: [1=1I58=1X7=1X3=1X5=1D23=] (11)  50 3:2:x:1 1:61:x:1 1:69:x:1 1:73:x:1 2:78:x:1  INDEL 107308237
** d->d: OK
**   DSC: [78=1D1=1D3=1X18=] (7)  50 2:78:x:1 2:79:x:1 1:83:x:1  INDEL 111323453
** i->i: OK
**   DSC: [25=1I3=1I71=] (5)  50 3:26:x:1 3:30:x:1  INDEL 129522037
** d->i: OK
**   DSC: [16=1D47=2I1=2I33=] (7)  50 2:16:x:1 3:64:x:2 3:67:x:2  INDEL 5952221
** i->d: OK
**   DSC: [8M1I19M1D22M] (5)  27^G22 3:9:x:1 2:28:x:1  INDEL 69538007/1
** 
*/
int bamGetDsc( BamAlign *bamAlign, int cigarFlag, int mdFlag, BamDsc **fbamDsc, int *fnumBamDsc, int *fallocBamDsc, int clip5[2], int clip3[2], int *fstatus )
{
  int i;
  int ired, mchr, ochr; /* Note: ired gives 0-based read position; ochr gives (0-based) offset from alignment start in ref. without introns; mchr gives (0-based) offset from alignment start in ref. with introns */
  int idsc, jdsc;
  int clen;
  int icig, ncig;
  int ifld, nfld;
  int indelFlag;
  int mxlen;
  int allocBamDsc;
  int status;
  uint32_t *pcig;

  char  cop;
  char *cptr, *nptr;

  BamDsc       *bamDsc;
  BamOptField  *field;

  static int      initFlag     = 0;
  static int32_t  sizMap2red   = 0;
  static int32_t  *map2red     = NULL; /* Note: map2red gives 0-based read position. */
  static int32_t  *map2chr     = NULL; /* Note: map2chr gives 0-based chromosome position. */
  static int       sizSeq      = 0;
  static char     *seq         = NULL;

  if( initFlag == 0 )
  {
    sizMap2red = 10 * bamAlign->l_seq;
    map2red = (int32_t *)malloc( sizMap2red * sizeof( int32_t ) );
    if( map2red == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }

    map2chr = (int32_t *)malloc( sizMap2red * sizeof( int32_t ) );
    if( map2chr == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }

    initFlag = 1;
  }

  /*
  ** Calculate a safe size of the map2red array.
  */
  ncig  = bamAlign->n_cigar_op;
  pcig  = bamAlign->cigar;
  mxlen = bamAlign->l_seq;
  for( icig = 0; icig < ncig; ++icig )
  {
    bamCigarOp( pcig[icig], &cop, &clen, &status );
    if( cop == 'D' )
    {
      mxlen += clen;
    }
  }
  mxlen += 8192;  /* Reduce reallocation activity. */

  if( mxlen > sizMap2red )
  {
    sizMap2red += mxlen;
    map2red = (int32_t *)realloc( map2red, sizMap2red * sizeof( int32_t ) );
    if( map2red == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }

    map2chr = (int32_t *)malloc( sizMap2red * sizeof( int32_t ) );
    if( map2chr == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
  }

  bamDsc      = *fbamDsc;
  allocBamDsc = *fallocBamDsc;

  if( bamDsc == NULL || allocBamDsc == 0 )
  {  
    allocBamDsc = 1024;
    bamDsc = (BamDsc *)malloc( allocBamDsc * sizeof( BamDsc ) );
    if( bamDsc == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
  }

  if( bamAlign->l_seq * 2 >= allocBamDsc )
  {
    allocBamDsc += bamAlign->l_seq * 2;
    bamDsc = (BamDsc *)malloc( allocBamDsc * sizeof( BamDsc ) );
    if( bamDsc == NULL )
    {
      EMSG( "unable to allocate memory" );
      *fstatus = -1;
      return( 0 );
    }
  }

  /*
  ** Skip unmapped entries.
  */
  if( bamAlign->flag & 0x4 )
  {
    *fbamDsc      = bamDsc;
    *fallocBamDsc = allocBamDsc;
    *fnumBamDsc   = 0;
    clip5[0]      = 0;
    clip3[0]      = 0;
    clip5[1]      = 0;
    clip3[1]      = 0;

    *fstatus = 1;
    return( 0 );
  }


  /*
  ** Expand read sequence.
  */
  bamUnpackSeq( bamAlign, &seq, &sizSeq, &status );
  IF_STATUS_ZERO( "bad status: bamUnpackSeq" );


  /*
  **  CIGAR string
  **    6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
  **
  **  Op  BAM  Description
  **  M   0    alignment match (can be a sequence match or mismatch)
  **  I   1    insertion to the reference
  **  D   2    deletion from the reference
  **  N   3    skipped region from the reference
  **  S   4    soft clipping (clipped sequences present in SEQ)
  **  H   5    hard clipping (clipped sequences NOT present in SEQ)
  **  P   6    padding (silent deletion from padded reference)
  **  =   7    sequence match
  **  X   8    sequence mismatch
  **  
  **  o  H can only be present as the first and/or last operation.
  **  o  S may only have H operations between them and the ends of the CIGAR string.
  **  o  For mRNA-to-genome alignment, an N operation represents an intron. For other types of
  **  o  alignments, the interpretation of N is not defined.
  **  o  Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
  **
  **
  **  Problems:
  **    o  CIGAR strings represent partial read edits, excluding (typically) substitutions
  **    o  MD strings represent partial reference edits, excluding bases inserted in read
  **
  **  Notes:
  **    o  set offset[] for match/mismatch bases, which will appear in the MD string.
  **    o  increment offset map on insertions in read
  **    o  exclude clipped bases
  **
  */
  /*
  ** Set CIGAR string-based values.
  */
  ired       = 0;   /* location in read */
  ochr       = 0;   /* offset from alignment start in reference without introns */
  mchr       = 0;   /* offset from alignment start in reference with introns */
  pcig       = bamAlign->cigar;
  ncig       = bamAlign->n_cigar_op;
  clip5[0]   = 0;
  clip3[0]   = 0;
  clip5[1]   = 0;
  clip3[1]   = 0;
  idsc       = 0;
  indelFlag  = 0;

  for( icig = 0; icig < ncig; ++icig )
  {
    bamCigarOp( pcig[icig], &cop, &clen, &status );
    IF_STATUS_ZERO( "bad status: bamCigarOp" );

/*
{
  char tstr[8192];
  sprintf( tstr, "%d%c", clen, cop );
  strcat( cigar, tstr );
}
*/

    if( cop == 'M' || cop == '=' )
    {
      if( ochr + clen + 1024 >= sizMap2red )
      {
        sizMap2red += ochr + clen + 1024;
        map2red = (int32_t *)realloc( map2red, sizMap2red * sizeof( int32_t ) );
        if( map2red == NULL )
        {
          EMSG( "unable to allocate memory" );
          *fstatus = -1;
          return( 0 );
        }

        map2chr = (int32_t *)malloc( sizMap2red * sizeof( int32_t ) );
        if( map2chr == NULL )
        {
          EMSG( "unable to allocate memory" );
          *fstatus = -1;
          return( 0 );
        }
      }

      /*
      ** There may be substitutions in this 'matching' region
      ** so set map2red values.
      */
      for( i = 0; i < clen; ++i )
      {
        map2red[ochr+i] = ired + i;
        map2chr[ochr+i] = mchr + i;
      }
      ired += clen;
      ochr += (int64_t)clen;
      mchr += (int64_t)clen;
    }
    else
    if( cop == 'X' )
    {
      if( cigarFlag == 0 )
      {
        EMSG( "unexpected condition" );
        *fstatus = -1;
        return( 0 );
      }

// fprintf( stderr, "cigar: X clen: %d\n", clen );

      if( ochr + clen + 1024 >= sizMap2red )
      {
        sizMap2red += ochr + clen + 1024;
        map2red = (int32_t *)realloc( map2red, sizMap2red * sizeof( int32_t ) );
        if( map2red == NULL )
        {
          EMSG( "unable to allocate memory" );
          *fstatus = -1;
          return( 0 );
        }

        map2chr = (int32_t *)malloc( sizMap2red * sizeof( int32_t ) );
        if( map2chr == NULL )
        {
          EMSG( "unable to allocate memory" );
          *fstatus = -1;
          return( 0 );
        }
      }

      /*
      ** Record substitution.
      ** Note: we cannot assign the correct ref base value to
      **       bamDsc[idsc].base so it becomes unknown.
      **
      ** Allow for <n>X where <n> > 1 by 'expanding' into single
      ** substitutions.
      */
      for( i = 0; i < clen; ++i )
      {
        bamDsc[idsc].type = 1;
        bamDsc[idsc].pos  = ired + i + 1;
        bamDsc[idsc].len  = 1;
        bamDsc[idsc].base = (uint8_t)0;
        if( seq[ired+i] == 'N' )
        {
          bamDsc[idsc].base = BAM_DSC_RED_N;
        }
#ifdef BAM_DSC_CHR_POS
        if( cigarFlag == 1 ) /* cigar string has 'X's? */
        {
          bamDsc[idsc].posChr = bamAlign->pos + mchr + i;
        }
#endif
        ++idsc;

        /*
        ** We may want (someday/for some reason) to get this map
        ** value again so set map2red values.
        */
        map2red[ochr+i] = ired + i;
        map2chr[ochr+i] = mchr + i;
      }

// fprintf( stderr, "cigar: add %d substitutions: idsc %d\n", clen, idsc );

      ired += clen;
      ochr += (int64_t)clen;
      mchr += (int64_t)clen;
    }
    else
    if( cop == 'D' )
    {
      /*
      ** Record deletion.
      */
      bamDsc[idsc].type = 2;
      bamDsc[idsc].pos  = ired; /* point to (last) read base upstream of deletion */
      bamDsc[idsc].len  = clen;
      bamDsc[idsc].base = 0;
#ifdef BAM_DSC_CHR_POS
      bamDsc[idsc].posChr = bamAlign->pos + mchr; /* point to first deleted ref base */
#endif
      ++idsc;

      /*
      ** There cannot be substitutions in this reference
      ** region so don't set map2red values but move ochr
      ** to deletion end.
      */
      ochr += (int64_t)clen;
      mchr += (int64_t)clen;

      indelFlag = 1;
    }
    else
    if( cop == 'I' )
    {
      /*
      ** Record insertion discrepancy.
      */
      bamDsc[idsc].type = 3;
      bamDsc[idsc].pos  = ired + 1;
      bamDsc[idsc].len  = clen;
      bamDsc[idsc].base = 0;
#ifdef BAM_DSC_CHR_POS
      bamDsc[idsc].posChr = bamAlign->pos + mchr - 1; /* point to (last) ref base upstream of insertion */
#endif
      ++idsc;

      /*
      ** Increment read base counter.
      */
      ired += clen;

      indelFlag = 1;
    }
    else
    if( cop == 'H' )
    {
      /*
      ** Do not increment read or reference base
      ** counters because the positions are relative
      ** to the first base in the BAM file read
      ** sequence, and hard clipped bases are not in
      ** the read sequence.
      */
      if( icig == 0 )
      {
        clip5[1] += clen;
      }
      else
      {
        clip3[1] += clen;
      }
    }
    else
    if( cop == 'S' )
    {
      /*
      ** Increment read base counter but not
      ** reference position counter because
      ** the reference position is the first
      ** aligned base.
      */
      if( icig == 0 )
      {
        clip5[0] += clen;
      }
      else
      {
        clip3[0] += clen;
      }

      ired += clen;
    }
    else
    if( cop == 'N' )
    {
#ifdef INTRON_AS_DISCREPANCY
      /*
      ** Record deletion.
      */
      bamDsc[idsc].type = 4;
      bamDsc[idsc].pos  = ired; /* point to base upstream of intron (deletion) */
      bamDsc[idsc].len  = clen;
      bamDsc[idsc].base = 0;
#ifdef BAM_DSC_CHR_POS
      bamDsc[idsc].posChr = bamAlign->pos + mchr;
#endif
      ++idsc;
#endif

      /*
      ** Assume that the MD string has an '^' at the
      ** intron location. If not, our substitution
      ** bases will be wrong.
      ** Watch for missing '^' and warn user. 
      if( intronFlag == 0 )
      {
        EMSG( "check intron handling with valgrind" );
        *fstatus = -1;
        return( 0 );
        intronFlag = 1;
      }
      */

      /*
      ** There cannot be substitutions in this reference
      ** region so don't set map2red values but move ochr
      ** to intron end.
      */
      mchr += (int64_t)clen;

      indelFlag = 1;
    }
    else
    if( cop == 'P' )
    {
      sprintf( _msg_, "unsupported CIGAR operation %c", cop );
      *fstatus = -1;
      return( 0 );
    }
    else
    {
      sprintf( _msg_, "unknown CIGAR operation '%c'", cop );
      EMSG( _msg_ );
      *fstatus = -1;
      return( 0 );
    }
  }

  /*
  ** Use MD field information, if it's expected to exist.
  ** We look for substitution information
  **   o  if CIGAR has Xs, we get substitution reference bases
  **   o  if CIGAR does not have Xs, we get/add substitutions
  ** Note: use map2red to get read base positions (MD field
  **       positions are in the reference sequence.
  */
  if( mdFlag != 0 ) /* MD field appeared to be present in at least some entries? */
  {
    /*
    ** Find MD string, if it's present in this entry.
    */
    nfld  = bamAlign->numOptField;
    field = bamAlign->field;
    cptr  = NULL;

    for( ifld = 0; ifld < nfld; ++ifld )
    {
      if( field[ifld].tag[0] == 'M' && field[ifld].tag[1] == 'D' )
      {
        cptr = (char *)field[ifld].value;
        break;
      }
    }
  
    if( mdFlag == 1 && cptr == NULL )
    {
      /*
      ** MD field appeared to be present in all entries but it's
      ** not in this one.
      */
      sprintf( _msg_, "missing 'MD' field in BAM entry: read: %s flag: %d", bamAlign->read_name, bamAlign->flag );
      EMSG( _msg_ );
  
      *fbamDsc      = bamDsc;
      *fallocBamDsc = allocBamDsc;
      *fnumBamDsc   = idsc;
  
      *fstatus = 0;
      return( 0 );
    }

    ochr = 0;
    while( cptr != NULL && *cptr != '\0' )
    {
      clen = strtol( cptr, &nptr, 10 );
      if( clen == 0 && nptr == cptr )
      {
        EMSG( "unexpected condition" );
        *fstatus = -1;
        return( 0 );
      }

      if( *nptr == 0 )
      {
        break;
      }
  
      ochr += clen;
 
      jdsc = 0; 
      cptr = nptr;
      for( ; cptr != NULL && *cptr != '\0' && !( isdigit( *cptr ) ); ++cptr )
      {
        if( *cptr != '^' )
        {
          if( cigarFlag != 1 )
          {
            /*
            ** There are no 'X' (substitutions) marked in CIGAR codes.
            ** Note: initialize bamDsc[].base bit flags
            */
            bamDsc[idsc].type = 1;
            bamDsc[idsc].pos  = map2red[ochr] + 1;
            bamDsc[idsc].len  = 1;

            if( *cptr == 'A' )
              bamDsc[idsc].base = BAM_DSC_REF_A;
            else
            if( *cptr == 'C' )
              bamDsc[idsc].base = BAM_DSC_REF_C;
            else
            if( *cptr == 'G' )
              bamDsc[idsc].base = BAM_DSC_REF_G;
            else
            if( *cptr == 'T' )
              bamDsc[idsc].base = BAM_DSC_REF_T;
            else
            if( *cptr == 'N' )  /* N in chromosome */
            {
              bamDsc[idsc].base = BAM_DSC_REF_N;
            }
            else
            {
              bamDsc[idsc].base = 0;
            }
  
            if( toupper( seq[map2red[ochr]] ) == 'N' ) /* N in read */
            {
              bamDsc[idsc].base |= BAM_DSC_RED_N;
            }

#ifdef BAM_DSC_CHR_POS
            bamDsc[idsc].posChr = bamAlign->pos + map2chr[ochr];
#endif

            ++idsc;
          }
          else
          {
            /*
            ** There are 'X's (substitutions) marked in CIGAR codes.
            **
            ** Here we must find the existing substitution
            ** discrepancy by matching the CIGAR-based
            ** reference position to ochr here.
            **
            ** The substitution positions in both the CIGAR and
            ** MD sources are ordered from alignment start to end.
            **
            ** Note: we set the bamDsc[].base bit flag for a read 'N'
            **       so do not initialize bamDsc[].base here.
            **
            */
            for( ; jdsc < idsc; ++jdsc )
            {
              if( bamDsc[jdsc].type == 1 && 
                  bamDsc[jdsc].pos == ( map2red[ochr] + 1 ) )
              {
                if( *cptr == 'A' )
                  bamDsc[jdsc].base |= BAM_DSC_REF_A;
                else
                if( *cptr == 'C' )
                  bamDsc[jdsc].base |= BAM_DSC_REF_C;
                else
                if( *cptr == 'G' )
                  bamDsc[jdsc].base |= BAM_DSC_REF_G;
                else
                if( *cptr == 'T' )
                  bamDsc[jdsc].base |= BAM_DSC_REF_T;
                else
                if( *cptr == 'N' )  /* N in chromosome */
                {
                  bamDsc[jdsc].base |= BAM_DSC_REF_N;
                }
              }

              if( bamDsc[jdsc].pos >= ( map2red[ochr] + 1 ) )
              {
                break;
              }
            }
          }

          ++ochr;
        }
        else
        {
          /*
          ** Increment 'ochr' over deleted base string.
          ** SAM format: MD Z String for mismatching positions. Regex : [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
          */
          ++cptr;
          for( ; *cptr >= 'A' && *cptr <= 'Z'; ++cptr )
          {
            ++ochr;
          }
          --cptr;
        }
      }
    }
  }

  /*
  ** Unless BAM_DSC_NO_SORT is defined,
  ** sort dsc entries by (a) indel, then substitution,
  ** and (b) position in read. This was the behavior
  ** before allowing no MD field.
  */
#ifndef BAM_DSC_NO_SORT
  if( idsc > 1 && indelFlag )
  {
    qsort( bamDsc, idsc, sizeof( BamDsc ), xcmpSortBamDsc );
  }
#endif

  *fbamDsc      = bamDsc;
  *fallocBamDsc = allocBamDsc;
  *fnumBamDsc   = idsc;

  *fstatus = 0;

  return( 0 );
}


#ifndef BAM_DSC_NO_SORT
static int xcmpSortBamDsc( const void *fa, const void *fb )
{
  int istat;

  BamDsc *a;
  BamDsc *b;

  a = (BamDsc *)fa;
  b = (BamDsc *)fb;

  istat = ( a->type < 2 ) - ( b->type < 2 );
  if( istat )
  {
    return( istat );
  }

  if( a->pos < b->pos )
  {
    return( -1 );
  }
  else
  if( a->pos > b->pos )
  {
    return( 1 );
  }

  return( 0 );
}
#endif


/*
** The function 'strnum_cmp' is from Heng Li's samtools.
*/
int strnum_cmp(const void *_a, const void *_b)
{
  const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
  const unsigned char *pa = a, *pb = b;
  while (*pa && *pb) {
    if (isdigit(*pa) && isdigit(*pb)) {
      while (*pa == '0') ++pa;
      while (*pb == '0') ++pb;
      while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
      if (isdigit(*pa) && isdigit(*pb)) {
        int i = 0;
        while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
        return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
      } else if (isdigit(*pa)) return 1;
      else if (isdigit(*pb)) return -1;
      else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
    } else {
      if (*pa != *pb) return (int)*pa - (int)*pb;
      ++pa; ++pb;
    }
  }
  return *pa? 1 : *pb? -1 : 0;
}


/*
** Set read implied start and end in genomic coordinates.
** Notes:
**   o  returns ends of the implied read in the reference;
**      that is, it extends the genomic coordinates to the
**      first and last read base when the alignment is clipped
*/
int bamCalcEndsRead( BamAlign *bamAlign, BamDsc *bamDsc, int numBamDsc, int clip5[2], int clip3[2], int64_t *fbpos, int64_t *fepos, int *fstatus )
{
  int idsc;

  int64_t adj;
  int64_t bpos, epos;

  if( bamAlign->flag & 0x4 )
  {
    *fbpos = 0;
    *fepos = 0;

    *fstatus = 0;
    return( 0 );
  }

  adj = 0L;
  for( idsc = 0; idsc < numBamDsc; ++idsc )
  {
    if( bamDsc[idsc].type == 2 )
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

  bpos = bamAlign->pos - (int64_t)( clip5[0] + clip5[1] );
  epos = bpos + ( bamAlign->l_seq + clip5[1] + clip3[1] ) - 1 + adj;

  *fbpos = bpos;
  *fepos = epos;

  return( 0 );
}

