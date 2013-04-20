
#include <algorithm>
#include <iostream>
#include <vector>

#include "merge.h"


using std::sort;
using std::cerr;
using std::endl;
using std::vector;


#define MERGE_SIZE 128

#define MAX( a, b ) (((a) > (b)) ? (a) : (b))
#define MIN( a, b ) (((a) < (b)) ? (a) : (b))

enum cmp_t { LT, GT, EQ };
enum res_t { SUCCESS, FAILURE, ERROR };


char nuc2bits( const char nuc )
{
    switch ( nuc ) {
    case 'A': return 1;
    case 'C': return 2;
    case 'G': return 4;
    case 'T': return 8;
    case 'M': return 1 | 2;
    case 'R': return 1 | 4;
    case 'W': return 1 | 8;
    case 'S': return 2 | 4;
    case 'Y': return 2 | 8;
    case 'K': return 4 | 8;
    case 'V': return 1 | 2 | 4;
    case 'H': return 1 | 2 | 8;
    case 'D': return 1 | 4 | 8;
    case 'B': return 2 | 4 | 8;
    default:  return 1 | 2 | 4 | 8;
    }
}

char bits2nuc( const char bits )
{
    switch ( bits ) {
    case 1:           return 'A';
    case 2:           return 'C';
    case 4:           return 'G';
    case 8:           return 'T';
    case (1 | 2):     return 'M';
    case (1 | 4):     return 'R';
    case (1 | 8):     return 'W';
    case (2 | 4):     return 'S';
    case (2 | 8):     return 'Y';
    case (4 | 8):     return 'K';
    case (1 | 2 | 4): return 'V';
    case (1 | 2 | 8): return 'H';
    case (1 | 4 | 8): return 'D';
    case (2 | 4 | 8): return 'B';
    default:          return 'N';
    }
}

inline
cmp_t pos_cmp( const triple_t & x, const triple_t & y )
{
    // test the column, then the insertion, then equal
    if ( x.col > y.col )
        return GT;
    else if ( y.col > x.col )
        return LT;
    else if ( x.ins > y.ins )
        return GT;
    else if ( y.ins < x.ins )
        return LT;
    else
        return EQ;
}


bool ncontrib_cmp( const aligned_t & x, const aligned_t & y )
{
    return x.ncontrib > y.ncontrib;
}


inline
void cerr_triple( const triple_t & x, bool end=true )
{
    cerr << x.col << " " << x.ins << " " << bits2nuc( x.nuc );
    if ( end )
        cerr << endl;
    else
        cerr << ", ";
}


res_t merge_two(
    const aligned_t & xs,
    const aligned_t & ys,
    const opts_t & opts,
    aligned_t & merged
    )
{
    int overlap = 0;
    int xidx = 0;
    int yidx = 0;
    int midx = 0;
    int mlen = 0;
    int cmp = 0;

    if ( !xs.len || !ys.len )
        return FAILURE;

    // if there is absolutely no hope to reach min_overlap, skip
    if ( xs.rpos < ys.lpos + opts.min_overlap &&
            ys.rpos < xs.lpos + opts.min_overlap )
        return FAILURE;

    cmp = pos_cmp( xs.data[ xidx ], ys.data[ yidx ] );

    // cerr_triple( xs.data[ xidx ], false );
    // cerr_triple( xs.data[ yidx ] );

    // disregard overhangs
    if ( cmp == LT ) {
        for ( xidx = 1; cmp == LT && xidx < xs.len; ++xidx )
            cmp = pos_cmp( xs.data[ xidx ], ys.data[ yidx ] );
        // if its not a match, it's a gap, then backup one in xs
        if ( cmp != EQ )
            --xidx;
        // mlen is now the # of xs already visited
        mlen = xidx + 1;
    }
    else if ( cmp == GT ) {
        for ( yidx = 1; cmp == GT && yidx < ys.len; ++yidx )
            cmp = pos_cmp( xs.data[ xidx ], ys.data[ yidx ] );
        // if its not a match, it's a gap, then backup one in ys
        if ( cmp != EQ )
            --yidx;
        // mlen is now the # of ys already visited
        mlen = yidx + 1;
    }

    // compute the amount of overlap
    for ( ; xidx < xs.len && yidx < ys.len; ) {
        cmp = pos_cmp( xs.data[ xidx ], ys.data[ yidx ] );
        if ( cmp == LT ) {
            if ( !opts.tol_gaps )
                return FAILURE;
            ++xidx;
            ++mlen;
        }
        else if ( cmp == GT ) {
            if ( !opts.tol_gaps )
                return FAILURE;
            ++yidx;
            ++mlen;
        }
        // if the nucleotides match, move ahead
        else if ( xs.data[ xidx ].nuc == ys.data[ yidx ].nuc ||
                    ( opts.tol_ambigs && xs.data[ xidx ].nuc & ys.data[ yidx ].nuc ) ) {
                ++overlap;
                ++xidx;
                ++yidx;
                ++mlen;
        }
        // nucleotides do not match, abort early
        else
            return FAILURE;
    }

    if ( overlap < opts.min_overlap )
        return FAILURE;

    // get the remainder of either xs or ys, whichever remains
    if ( xidx < xs.len )
        mlen += xs.len - xidx;
    else if ( yidx < ys.len )
        mlen += ys.len - yidx;

    merged.len = mlen;
    merged.data = new triple_t[ merged.len ];
    merged.lpos = MIN( xs.lpos, ys.lpos );
    merged.rpos = MAX( xs.rpos, ys.rpos );
    merged.ncontrib = xs.ncontrib + ys.ncontrib;

    if ( !merged.data )
        return FAILURE;

    xidx = 0;
    yidx = 0;
    
    cerr << "leading and matching part" << endl;

    for ( ; xidx < xs.len && yidx < ys.len; ) {
        cmp = pos_cmp( xs.data[ xidx ], ys.data[ yidx ] );
        if ( cmp == LT )
            merged.data[ midx++ ] = xs.data[ xidx++ ];
        else if ( cmp == GT )
            merged.data[ midx++ ] = ys.data[ yidx++ ];
        else {
            merged.data[ midx ] = xs.data[ xidx ];
            merged.data[ midx ].nuc = MIN( xs.data[ xidx ].nuc, ys.data[ yidx ].nuc );
            ++xidx;
            ++yidx;
            ++midx;
        }
    }

    cerr << "trailing part" << endl;

    if ( xidx < xs.len )
        for ( ; xidx < xs.len; ++xidx )
            merged.data[ midx++ ] = xs.data[ xidx ];
    else if ( yidx < ys.len )
        for ( ; yidx < ys.len; ++yidx )
            merged.data[ midx++ ] = ys.data[ yidx ];

    if ( midx != mlen )
        cerr << "error: failed to fill 'merged' data" << endl;

    return SUCCESS;
}

int merge_clusters(
    vector< aligned_t > & clusters,
    const opts_t & opts
    )
{
    int i, j, nclusters;
    aligned_t merged;
    res_t res = FAILURE; 

    sort( clusters.begin(), clusters.end(), ncontrib_cmp );

begin:
    nclusters = 0;

    for ( i = 0; i < clusters.size(); ++i ) {
        for ( j = i + 1; j < clusters.size(); ++j ) {
            res = merge_two( clusters[ i ], clusters[ j ], opts, merged );
            if ( res == SUCCESS ) {
                aligned_destroy( clusters[ i ] );
                aligned_destroy( clusters[ j ] );
                // replace i and remove j
                clusters[ i ] = merged;
                clusters.erase( clusters.begin() + j );
                goto begin;
            }
            else if ( res == ERROR )
                return -1;
        if ( clusters[ i ].ncontrib >= opts.min_reads )
            ++nclusters;
        }
    }

    return nclusters;
}

aligned_t * merge__(
    const int nreads,
    const aligned_t * const reads,
    const opts_t * const opts,
    int * nclusters
    )
{
    int i, j, merge_size = MERGE_SIZE;
    vector< aligned_t > clusters;
    aligned_t * clusters_;
    aligned_t merged;
    res_t res;

    clusters.push_back( reads[ 0 ] );

    for ( i = 1; i < nreads; ) {
        for ( j = 0; j < clusters.size(); ++j ) {
            res = merge_two( reads[ i ], clusters[ j ], *opts, merged );
            if ( res == SUCCESS ) {
                // destroy our clusters
                if ( clusters[ j ].ncontrib > 1 )
                    aligned_destroy( clusters[ j ] );
                clusters[ j ] = merged;
                if ( clusters.size() > merge_size ) {
                    if ( merge_clusters( clusters, *opts ) < 0 )
                        goto error;
                    merge_size *= 2;
                }
                goto next;
            }
            else if ( res == ERROR )
                goto error;
        }
        clusters.push_back( reads[ i ] );
next:
    ++i;
    }

    *nclusters = merge_clusters( clusters, *opts );
    
    if ( *nclusters < 0 )
        goto error;
    
    clusters_ = new aligned_t[ *nclusters ];

    if ( !clusters_ )
        goto error;

    for ( i = 0, j = 0; i < clusters.size(); ++i ) {
        if ( clusters[ i ].ncontrib >= opts->min_reads )
            clusters_[ j++ ] = clusters[ i ];
        // only free clusters we allocated
        else if ( clusters[ i ].ncontrib > 1 )
            aligned_destroy( clusters[ i ] );
    }
    
    goto end;

error:
    *nclusters = 0;
    return NULL;

end:
    return clusters_;
}

void aligned_destroy( aligned_t read )
{
    delete[] read.data;
    read.data = NULL;
    read.len = 0;
}
