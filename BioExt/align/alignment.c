
#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alignment.h"

//____________________________________________________________________________________

// match or skip whole codons
#define HY_111_111 0
#define HY_111_000 1
#define HY_000_111 2

// match 3 in the ref to 1 in the query
#define HY_111_100 3
#define HY_111_010 4
#define HY_111_001 5

#define HY_3X1_START 3
#define HY_3X1_COUNT 3

// match 3 in the ref to 2 in the query
#define HY_111_110 6
#define HY_111_101 7
#define HY_111_011 8

#define HY_3X2_START 6
#define HY_3X2_COUNT 3

// match 3 in the ref to 4 in the query
#define HY_1110_1111 9
#define HY_1101_1111 10
#define HY_1011_1111 11
#define HY_0111_1111 12

#define HY_3X4_START 9
#define HY_3X4_COUNT 4

// match 5 in the ref to 5 in the query
#define HY_11100_11111 13
#define HY_11010_11111 14
#define HY_11001_11111 15
#define HY_10110_11111 16
#define HY_10101_11111 17
#define HY_10011_11111 18
#define HY_01110_11111 19
#define HY_01101_11111 20
#define HY_01011_11111 21
#define HY_00111_11111 22

#define HY_3X5_START   13
#define HY_3X5_COUNT   10

#define HY_ALIGNMENT_TYPES_COUNT 23

//____________________________________________________________________________________

#define ALLOCA( type, len ) ( (type *) calloc( ( len ), sizeof( type ) ) )
#define ISNULL( var ) ( ( var ) == NULL )

#define MAX(a, b) ( (( a ) > ( b )) ? ( a ) : ( b ) )

#define FALSE 0
#define TRUE  1

#define A_LARGE_NUMBER 1.e100

//____________________________________________________________________________________

void print_score_matrix( FILE * const file
                       , const double * const score_matrix
                       , const long score_rows
                       , const long score_cols
                       )
{
    int i, j;

    for ( i = 0; i < score_rows; ++i ) {
        if ( i == 0 )
            fprintf( file, "[" );
        else
            fprintf( file, " " );
        
        for ( j = 0; j < score_cols; ++j ) {
            const double v = score_matrix[ i * score_cols + j ];
            if ( j == 0 )
                fprintf( file, "[ % 4.1f", v );
            else
                fprintf( file, ", % 4.1f", v );
        }
        
        if ( i == ( score_rows - 1 ) )
            fprintf( file, " ]]\n" );
        else
            fprintf( file, " ],\n" );
    }
}

//____________________________________________________________________________________

static long CodonAlignStringsStep( double * const score_matrix
                                 , double * const insertion_matrix
                                 , double * const deletion_matrix
                                 , long * const reference
                                 , long * const query
                                 , const long r
                                 , const long q
                                 , const long score_cols
                                 , const long char_count
                                 , const double miscall_cost
                                 , const double open_insertion
                                 , const double open_deletion
                                 , const double extend_insertion
                                 , const double extend_deletion
                                 , const double * const cost_matrix
                                 , const long cost_stride
                                 , const double * const codon3x5
                                 , const double * const codon3x4
                                 , const double * const codon3x2
                                 , const double * const codon3x1
                                 )
{
    /**
     * r is CODON position in the reference,
     * q is NUC position in the query,
     * curr is a pointer to the current position in the scoring matrix,
     * prev is a pointer to the previous CODON in the scoring matrix
     * NOTE: moving by score_cols in the scoring matrix changes the CODON
     *       position in the scoring matrix, as we're only interested in CODON
     *       alignments to the reference
     * rpos is the position of r in the reference
     */
    const long curr = ( r - 0 ) * score_cols + q, // where we currently are
               prev = ( r - 1 ) * score_cols + q, // up a codon in the reference
               offset3x5 = HY_3X5_COUNT * char_count * char_count * char_count, // both 3x5 and 3x4 are
               offset3x4 = HY_3X4_COUNT * char_count * char_count * char_count, // full codons
               offset3x2 = HY_3X2_COUNT * char_count * char_count,
               offset3x1 = HY_3X1_COUNT * char_count,
               rpos = r * 3; // the real position in R
    // mutable vars
    long r_codon = -1,
         q_codon = -1,
         best_choice = 0,
         i, choice, partial_codons[ 10 ]; // we need to multiply by 3 to get the NUC position
    // 3x5 codon specifications (negative indices)
    const long codon_spec_3x5[ 10 ][ 3 ] = {
        { 5, 4, 3 }, // 11100
        { 5, 4, 2 }, // 11010
        { 5, 4, 1 }, // 11001
        { 5, 3, 2 }, // 10110
        { 5, 3, 1 }, // 10101
        { 5, 2, 1 }, // 10011
        { 4, 3, 2 }, // 01110
        { 4, 3, 1 }, // 01101
        { 4, 2, 1 }, // 01011
        { 3, 2, 1 }  // 00111
    };
    // 3x4 codon specifications (negative indices)
    const long codon_spec_3x4[ 4 ][ 3 ] = {
        { 4, 3, 2 }, // 1110
        { 4, 3, 1 }, // 1101
        { 4, 2, 1 }, // 1011
        { 3, 2, 1 }  // 0111
    };
    double choices[ HY_ALIGNMENT_TYPES_COUNT ],
           max_score = -A_LARGE_NUMBER,
           penalty;

    // store the scores of our choices in choices,
    // pre-initialize to -infinity
    for ( i = 0; i < HY_ALIGNMENT_TYPES_COUNT; i++ ) {
        choices[ i ] = -A_LARGE_NUMBER;
    }

    // if we're at least a CODON away from the edge...
    // (psst, r is CODONs remember?)
    if ( r >= 1 ) {
        // if we're doing affine gaps (deletions)
        if ( deletion_matrix ) {
            choices[ HY_111_000 ] = MAX( score_matrix[ prev ] - open_deletion,
                                         deletion_matrix[ prev ] - ( r > 1 ? extend_deletion : open_deletion ) );
            deletion_matrix[ curr ] = choices[ HY_111_000 ];
        }
        else {
            choices[ HY_111_000 ] = score_matrix[ prev ] - open_deletion;
        }

        r_codon = ( reference[ rpos - 3 ]   * char_count
                  + reference[ rpos - 2 ] ) * char_count
                  + reference[ rpos - 1 ] ;

        if ( r_codon < 0 ) {
            r_codon = cost_stride - 1;
        }
    }

    // if we're at least 1 codon away from the edge
    if ( q >= 3 ) {
        // if we're doing affine gaps (insertions)
        if ( insertion_matrix ) {
            choices[ HY_000_111 ] = MAX( score_matrix[ curr - 3 ] - open_insertion,
                                         insertion_matrix[ curr - 3 ] - ( q > 3 ? extend_insertion : open_insertion ) );
            insertion_matrix[ curr ] = choices[ HY_000_111 ];
        }
        else {
            choices[ HY_000_111 ] = score_matrix[ curr - 3 ] - open_insertion;
        }

        q_codon = ( query[ q - 3 ]   * char_count
                  + query[ q - 2 ] ) * char_count
                  + query[ q - 1 ] ;

        if ( q_codon < 0 ) {
            q_codon = cost_stride - 1;
        }
    }

    // if q_codon and r_codon both exist, set the score equal to match
    if ( q_codon >= 0 ) {
        if ( r_codon >= 0 ) {
            choices[ HY_111_111 ] = score_matrix[ prev - 3 ]
                                  + cost_matrix[ r_codon * cost_stride + q_codon ];
        }
    }

    // we disallow partial moves in the reference, so those use to be here but are now gone

    // HERE BE DRAGONS!!!!

    // miscall matches, starting with 3x5, then 3x4, then 3x2, finally 3x1
    if ( r_codon >= 0 ) {
        // 3x5 partial codons
        if ( q >= 5 ) {
            // fill in the partial codons table
            // use a 10x1 array to allow for load hoisting,
            // we don't want to be bouncing cachelines in this critical inner loop
            for ( i = 0; i < HY_3X5_COUNT; ++i ) {
                partial_codons[ i ] = ( query[ q - codon_spec_3x5[ i ][ 0 ] ]   * char_count
                                      + query[ q - codon_spec_3x5[ i ][ 1 ] ] ) * char_count
                                      + query[ q - codon_spec_3x5[ i ][ 2 ] ] ;
            }

            // go over each choice, fill it in
            for ( i = 0; i < HY_3X5_COUNT; ++i ) {
                if ( partial_codons[ i ] >= 0 ) {
                    choice = HY_3X5_START + i;

                    // if we have a double ragged edge, don't penalize
                    if ( ( q == 5              && choice == HY_00111_11111 )
                      || ( q == score_cols - 1 && choice == HY_11100_11111 ) )
                        penalty = 0.;
                    // if we have a single ragged edge, penalize by a single miscall
                    // we don't have to worry about specifying each case here,
                    // as the 00111_11111 case takes preference above,
                    // so we don't have to explicitly avoid it
                    else if ( q == 5 && choice >= HY_01110_11111 )
                        penalty = miscall_cost;
                    // if we have a single ragged edge, penalize by a single miscall
                    // unfortunately these cases are spread out,
                    // so we have to enumerate them explicitly here
                    else if ( q == score_cols - 1
                           && ( choice == HY_11010_11111
                             || choice == HY_10110_11111
                             || choice == HY_01110_11111 ) )
                        penalty = miscall_cost;
                    // if we don't fit into any of these special cases,
                    // the miscall penalty is double (as we're matching 3 to 5)
                    else
                        penalty = 2. * miscall_cost;

                    choices[ choice ] = score_matrix[ prev - 5 ] - penalty
                                      + codon3x5[ r_codon * offset3x5 + HY_3X5_COUNT * partial_codons[ i ] + i ];
                }
            }
        }

        // 3x4 partial codons
        if ( q >= 4 ) {
            // fill in partial codons table
            for ( i = 0; i < HY_3X4_COUNT; ++i ) {
                partial_codons[ i ] = ( query[ q - codon_spec_3x4[ i ][ 0 ] ]   * char_count
                                      + query[ q - codon_spec_3x4[ i ][ 1 ] ] ) * char_count
                                      + query[ q - codon_spec_3x4[ i ][ 2 ] ] ;
            }

            // fill in choices
            for ( i = 0; i < HY_3X4_COUNT; ++i ) {
                if ( partial_codons[ i ] >= 0 ) {
                    choice = HY_3X4_START + i;

                    // if we have a ragged edge,
                    // penalize it not at all
                    if ( ( q == 4              && choice == HY_0111_1111 )
                      || ( q == score_cols - 1 && choice == HY_1110_1111 ) )
                        penalty = 0.;
                    // otherwise it's just a single miscall penalty
                    else
                        penalty = miscall_cost;

                    choices[ choice ] = score_matrix[ prev - 4 ] - penalty
                                      + codon3x4[ r_codon * offset3x4 + HY_3X4_COUNT * partial_codons[ i ] + i ];
                }
            }
        }

        // 3x2
        if ( q >= 2 ) {
            // only a single partial codon
            partial_codons[ 0 ] = query[ q - 2 ] * char_count
                                + query[ q - 1 ] ;

            // fill in choices
            if ( partial_codons[ 0 ] >= 0 ) {
                for ( i = 0; i < HY_3X2_COUNT; ++i ) {
                    choice = HY_3X2_START + i;

                    // if we have a ragged edge at the beginning or end,
                    // respectively, don't penalize it
                    if ( ( q == 2              && choice == HY_111_011 )
                      || ( q == score_cols - 1 && choice == HY_111_110 ) )
                        penalty = 0.;
                    // otherwise it's just a single miscall penalty
                    else
                        penalty = miscall_cost;

                    choices[ choice ] = score_matrix[ prev - 2 ] - penalty
                                      + codon3x2[ r_codon * offset3x2 + HY_3X2_COUNT * partial_codons[ 0 ] + i ];
                }
            }
        }

        // 3x1
        if ( q >= 1 ) {
            // only a single partial codon
            partial_codons[ 0 ] = query[ q - 1 ];

            // fill in choices
            if ( partial_codons[ 0 ] >= 0 ) {
                for ( i = 0; i < HY_3X1_COUNT; ++i ) {
                    choice = HY_3X1_START + i;

                    // if we have a double ragged edge,
                    // don't enforce a miscall penalty
                    if ( ( q == 1              && choice == HY_111_001 )
                      || ( q == score_cols - 1 && choice == HY_111_100 ) )
                        penalty = 0.;
                    // if we have a single ragged edge,
                    // enforce only a single miscall penalty
                    else if ( ( q == 1              && choice == HY_111_010 )
                           || ( q == score_cols - 1 && choice == HY_111_010 ) )
                        penalty = miscall_cost;
                    // otherwise we need a double miscall penalty,
                    // for the two positions we're inserting
                    else
                        penalty = 2. * miscall_cost;

                    choices[ choice ] = score_matrix[ prev - 1 ] - penalty
                                      + codon3x1[ r_codon * offset3x1 + HY_3X1_COUNT * partial_codons[ 0 ] + i ];
                }
            }
        }
    }

    // find the best possible choice
    for ( i = 0; i < HY_ALIGNMENT_TYPES_COUNT; ++i ) {
        /* if ( i > 0 )
         fprintf( stderr, ", " );
         fprintf( stderr, "( %ld, %.3g )", i, choices[ i ] ); */
        if ( choices[ i ] > max_score ) {
            best_choice = i;
            max_score = choices[ i ];
        }
    }

    /* fprintf( stderr, "\nscore: %.3g best: %ld\n", max_score, best_choice ); */
    // assign the score to the current position
    score_matrix[ curr ] = max_score;
    return best_choice;
}

//____________________________________________________________________________________

static inline void BacktrackAlign( signed char * const edit_ops
                                 , long * edit_ptr
                                 , long * r
                                 , long * q
                                 , double deletion
                                 , double insertion
                                 , double match
                                 )
{
    if ( match >= deletion && match >= insertion ) {
        --(*r);
        --(*q);
        edit_ops[ (*edit_ptr)++ ] = 0;
    }
    else if ( deletion >= insertion ) {
        --(*r);
        edit_ops[ (*edit_ptr)++ ] = -1;
    }
    else {
        --(*q);
        edit_ops[ (*edit_ptr)++ ] = 1;
    }
}

//____________________________________________________________________________________

static inline void BacktrackAlignCodon( signed char * const edit_ops
                                      , long * edit_ptr
                                      , long * r
                                      , long * q
                                      , const long code
                                      )
{
    long r_str[ 5 ] = { 0, 0, 0, 0, 0 },
         q_str[ 5 ] = { 0, 0, 0, 0, 0 },
         idx         = 2;
    long frameshift = TRUE;

    switch ( code ) {
        // match
    case HY_111_111:
        r_str[0] = r_str[1] = r_str[2] = 1;
        q_str[0] = q_str[1] = q_str[2] = 1;
        frameshift = FALSE;
        break;

        // deletion
    case HY_111_000:
        r_str[0] = r_str[1] = r_str[2] = 1;
        frameshift = FALSE;
        break;

        // insertion
    case HY_000_111:
        q_str[0] = q_str[1] = q_str[2] = 1;
        frameshift = FALSE;
        break;

        // 3x2
    case HY_111_110:
        r_str[0] = r_str[1] = r_str[2] = 1;
        q_str[0] = q_str[1] = 1;
        break;

    case HY_111_101:
        r_str[0] = r_str[1] = r_str[2] = 1;
        q_str[0] = q_str[2] = 1;
        break;

    case HY_111_011:
        r_str[0] = r_str[1] = r_str[2] = 1;
        q_str[1] = q_str[2] = 1;
        break;

        // 3x1
    case HY_111_100:
        r_str[0] = r_str[1] = r_str[2] = 1;
        q_str[0] = 1;
        break;

    case HY_111_010:
        r_str[0] = r_str[1] = r_str[2] = 1;
        q_str[1] = 1;
        break;

    case HY_111_001:
        r_str[0] = r_str[1] = r_str[2] = 1;
        q_str[2] = 1;
        break;

        // 3x4
    case HY_1110_1111:
        q_str[0] = q_str[1] = q_str[2] = q_str[3] = 1;
        r_str[0] = r_str[1] = r_str[2] = 1;
        idx = 3;
        break;

    case HY_1101_1111:
        q_str[0] = q_str[1] = q_str[2] = q_str[3] = 1;
        r_str[0] = r_str[1] = r_str[3] = 1;
        idx = 3;
        break;

    case HY_1011_1111:
        q_str[0] = q_str[1] = q_str[2] = q_str[3] = 1;
        r_str[0] = r_str[2] = r_str[3] = 1;
        idx = 3;
        break;

    case HY_0111_1111:
        q_str[0] = q_str[1] = q_str[2] = q_str[3] = 1;
        r_str[1] = r_str[2] = r_str[3] = 1;
        idx = 3;
        break;

        // 3x5
    case HY_11100_11111:
        q_str[0] = q_str[1] = q_str[2] = q_str[3] = q_str[4] = 1;
        r_str[0] = r_str[1] = r_str[2] = 1;
        idx = 4;
        break;

    case HY_11010_11111:
        q_str[0] = q_str[1] = q_str[2] = q_str[3] = q_str[4] = 1;
        r_str[0] = r_str[1] = r_str[3] = 1;
        idx = 4;
        break;

    case HY_11001_11111:
        q_str[0] = q_str[1] = q_str[2] = q_str[3] = q_str[4] = 1;
        r_str[0] = r_str[1] = r_str[4] = 1;
        idx = 4;
        break;

    case HY_10110_11111:
        q_str[0] = q_str[1] = q_str[2] = q_str[3] = q_str[4] = 1;
        r_str[0] = r_str[2] = r_str[3] = 1;
        idx = 4;
        break;

    case HY_10101_11111:
        q_str[0] = q_str[1] = q_str[2] = q_str[3] = q_str[4] = 1;
        r_str[0] = r_str[2] = r_str[4] = 1;
        idx = 4;
        break;

    case HY_10011_11111:
        q_str[0] = q_str[1] = q_str[2] = q_str[3] = q_str[4] = 1;
        r_str[0] = r_str[3] = r_str[4] = 1;
        idx = 4;
        break;

    case HY_01110_11111:
        q_str[0] = q_str[1] = q_str[2] = q_str[3] = q_str[4] = 1;
        r_str[1] = r_str[2] = r_str[3] = 1;
        idx = 4;
        break;

    case HY_01101_11111:
        q_str[0] = q_str[1] = q_str[2] = q_str[3] = q_str[4] = 1;
        r_str[1] = r_str[2] = r_str[4] = 1;
        idx = 4;
        break;

    case HY_01011_11111:
        q_str[0] = q_str[1] = q_str[2] = q_str[3] = q_str[4] = 1;
        r_str[1] = r_str[3] = r_str[4] = 1;
        idx = 4;
        break;

    case HY_00111_11111:
        q_str[0] = q_str[1] = q_str[2] = q_str[3] = q_str[4] = 1;
        r_str[2] = r_str[3] = r_str[4] = 1;
        idx = 4;
        break;

    default:
        assert( 0 );
    }

    for ( ; idx >= 0 ; --idx ) {
        if ( r_str[ idx ] ) {
            if ( q_str[ idx ] ) {
                --(*r);
                --(*q);
                edit_ops[ (*edit_ptr)++ ] = 0;
            }
            else {
                --(*r);
                edit_ops[ (*edit_ptr)++ ] = -( frameshift ? 2 : 1 );
            }
        }
        else {
            --(*q);
            edit_ops[ (*edit_ptr)++ ] = ( frameshift ? 2 : 1 );
        }
    }
}


//____________________________________________________________________________________

static inline void MatchScore( char * r_str
                             , char * q_str
                             , const long r
                             , const long q
                             , const long * const char_map
                             , const double * const cost_matrix
                             , const long cost_stride
                             , double * const score
                             )
{
    const long r_char = char_map[ (int) r_str[ r - 1 ] ];

    if ( r_char >= 0 ) {
        const long q_char = char_map[ (int) q_str[ q - 1 ] ];
        if ( q_char >= 0 )
            (*score) += cost_matrix[ r_char * cost_stride + q_char ];
    }
}

//____________________________________________________________________________________

double AlignStrings( char * const r_str
                   , char * const q_str
                   , char ** r_res
                   , char ** q_res
                   , const long char_count
                   , const long * const char_map
                   , const double * const cost_matrix
                   , const long cost_stride
                   , const char gap
                   , const double open_insertion
                   , const double extend_insertion
                   , const double open_deletion
                   , const double extend_deletion
                   , const double miscall_cost
                   , const long do_local
                   , const long do_affine
                   , const long do_codon
                   , const double * const codon3x5
                   , const double * const codon3x4
                   , const double * const codon3x2
                   , const double * const codon3x1
                   , double * const score_matrix
                   , double * const deletion_matrix
                   , double * const insertion_matrix
                   )
{
    const unsigned long r_len = strlen( r_str ),
                        q_len = strlen( q_str ),
                        ref_stride = ( do_codon ? 3 : 1 ),
                        score_rows = r_len / ref_stride + 1,
                        score_cols = q_len + 1;
    long i, j, k;
    double score = 0.;

    if ( do_codon && ( r_len % 3 != 0 ) ) {
        *r_res = NULL;
        *q_res = NULL;
        return -A_LARGE_NUMBER;
    }

    // handle some corner cases,
    // return early if possible
    if ( score_rows <= 1 ) {
        if ( score_cols > 1 ) {
            *r_res = ALLOCA( char, q_len + 1 );
            *q_res = ALLOCA( char, q_len + 1 );

            if ( ISNULL( *r_res ) || ISNULL( *q_res ) ) {
                free( *r_res );
                free( *q_res );
                *r_res = NULL;
                *q_res = NULL;
                return 0.;
            }

            // no ref, just query, which remains untouched
            memcpy( *q_res, q_str, q_len + 1 );
            // ref full of gaps
            memset( *r_res, gap, sizeof( char ) * q_len );
            // null terminate
            r_res[ q_len ] = '\0';
            q_res[ q_len ] = '\0';

            // compute score
            if ( ! do_local ) {
                if ( do_affine )
                    score = -open_insertion - ( q_len - 1 ) * extend_insertion;
                else
                    score = -open_insertion * q_len;
            }
        }
    }
    else {
        if ( score_rows <= 1 ) {
            *r_res = ALLOCA( char, r_len + 1 );
            *q_res = ALLOCA( char, r_len + 1 );

            if ( ISNULL( *r_res ) || ISNULL( *q_res ) ) {
                free( *r_res );
                free( *q_res );
                *r_res = NULL;
                *q_res = NULL;
                return 0.;
            }

            // no query, just ref, which remains untouched
            memcpy( *r_res, r_str, r_len + 1 );
            // ref full of gaps
            memset( *q_res, gap, sizeof( char ) * r_len );
            // null terminate
            r_res[ r_len ] = '\0';
            q_res[ r_len ] = '\0';

            // if do local, score is 0
            if ( ! do_local ) {
                if ( do_affine )
                    score = -open_deletion - ( r_len - 1 ) * extend_deletion;
                else
                    score = -open_deletion * r_len;
            }
        }
        else {
            long edit_ptr = 0;
            // don't forget the optional termination character
            signed char * const edit_ops = ALLOCA( signed char, r_len + q_len );
#if 0
            double * const score_matrix = ALLOCA( double, score_rows * score_cols ),
                   * const deletion_matrix  = do_affine ? ALLOCA( double, score_rows * score_cols ) : NULL,
                   * const insertion_matrix = do_affine ? ALLOCA( double, score_rows * score_cols ) : NULL;
#endif
            // encode each string using the character map (char_map)
            long * const r_enc = ALLOCA( long, r_len ),
                 * const q_enc = ALLOCA( long, q_len );

            if ( ISNULL( edit_ops )
#if 0
              || ISNULL( score_matrix )
              || ( do_affine && ( ISNULL( insertion_matrix )
                               || ISNULL( deletion_matrix  ) ) )
#endif
              || ISNULL( r_enc )
              || ISNULL( q_enc ) ) {
                *r_res = NULL;
                *q_res = NULL;
                goto end;
            }

#if 1
            // if this is set, then 8 0-bytes is equivalent to a 0. double,
            // which is done because we calloc'd
#ifdef __STDC_IEC_559__
            memset( score_matrix, 0, sizeof( double ) * score_rows * score_cols );
            if ( do_affine ) {
                memset( deletion_matrix, 0, sizeof( double ) * score_rows * score_cols );
                memset( insertion_matrix, 0, sizeof( double ) * score_rows * score_cols );
            }
#else
            for ( i = 0; i < score_rows * score_cols; ++i ) {
                score_matrix[ i ] = 0.;
            }
            if ( do_affine )
                for ( i = 0; i < score_rows * score_cols; ++i ) {
                    deletion_matrix[ i ] = 0.;
                    insertion_matrix[ i ] = 0.;
                }
#endif
#endif
            if ( do_codon ) {
                for ( i = 0; i < r_len; ++i )
                    r_enc[ i ] = char_map[ (int) r_str[ i ] ];

                for ( i = 0; i < q_len; ++i )
                    q_enc[ i ] = char_map[ (int) q_str[ i ] ];
            }

            // pre-initialize the values in the various matrices
            if ( ! do_local ) {
                double cost;

                // initialize gap costs in first column and first row
                // they are 0 for local alignments, so ignore
                if ( do_affine ) {
                    // first handle insertions
                    cost = -open_insertion;
                    insertion_matrix[ 0 ] = cost;

                    for ( i = 1; i < score_cols; ++i, cost -= extend_insertion ) {
                        score_matrix[ i ] = cost;
                        insertion_matrix[ i ] = cost;
                        deletion_matrix[ i ] = cost;
                    }

                    // then deletions
                    cost = -open_deletion;
                    deletion_matrix[ 0 ] = cost;

                    for ( i = score_cols; i < score_rows * score_cols; i += score_cols, cost -= extend_deletion ) {
                        score_matrix[ i ] = cost;
                        insertion_matrix[ i ] = cost;
                        deletion_matrix[ i ] = cost;
                    }
                }
                else {
                    // handle the do_local, regular (non codon-alignment) case
                    if ( ! do_codon ) {
                        cost = -open_insertion;

                        for ( i = 1; i < score_cols; ++i, cost -= open_insertion )
                            score_matrix[ i ] = cost;

                        cost = -open_deletion;

                        for ( i = score_cols; i < score_rows * score_cols; i += score_cols, cost -= open_deletion )
                            score_matrix[ i ] = cost;

                        // handle the do_local, do_codon case
                    }
                    else {
                        cost = -open_insertion;

                        for ( i = 1; i < score_cols; ++i, cost -= open_insertion )
                            score_matrix[ i ] = cost - ( i % 3 != 1 ? miscall_cost : 0 );

                        cost = -open_deletion;

                        for ( i = score_cols, j = 0; i < score_rows * score_cols; i += score_cols, cost -= open_insertion, ++j )
                            score_matrix[ i ] = cost - ( j % 3 != 0 ? miscall_cost : 0 );
                    }
                }

                // if we're doing a local alignment,
                // the costs of opening a deletion or an insertion
                // remain the same no matter how far down the ref or query
                // we've traveled, respectively
            }
            else {
                if ( do_affine ) {
                    if ( do_codon ) {
                        // XXX: should we be including the frameshift penalty here? I think not
                        // fill in the first row of the affine deletion matrix
                        // with the deletion cost plus the miscall penalty
                        for ( i = 1; i < score_cols; ++i )
                            deletion_matrix[ i ] = -open_deletion - ( i % 3 != 1 ? miscall_cost : 0 );

                        // fill in the first column of the affine insertion matrix
                        // with the insertion cost plus the miscall penalty
                        for ( i = score_cols, j = 0; i < score_rows * score_cols; i += score_cols, ++j )
                            insertion_matrix[ i ] = -open_insertion - ( j % 3 != 0 ? miscall_cost : 0 );
                    }
                    else {
                        // fill in the first row of the affine deletion matrix
                        // with the deletion cost
                        for ( i = 1; i < score_cols; ++i )
                            deletion_matrix[ i ] = -open_deletion;

                        // fill in the first column of the affine insertion matrix
                        // with the insertion cost
                        for ( i = score_cols; i < score_rows * score_cols; i += score_cols )
                            insertion_matrix[ i ] = -open_insertion;
                    }
                }
            }

            if ( do_codon ) {
                for ( i = 1; i < score_rows; ++i )
                    for ( j = 1; j < score_cols; ++j )
                        CodonAlignStringsStep( score_matrix
                                             , insertion_matrix
                                             , deletion_matrix
                                             , r_enc
                                             , q_enc
                                             , i
                                             , j
                                             , score_cols
                                             , char_count
                                             , miscall_cost
                                             , open_insertion
                                             , open_deletion
                                             , extend_insertion
                                             , extend_deletion
                                             , cost_matrix
                                             , cost_stride
                                             , codon3x5
                                             , codon3x4
                                             , codon3x2
                                             , codon3x1
                                             );

                // not doing codon alignment
            }
            else {
                for ( i = 1; i < score_rows; ++i ) {
                    const long r_char = char_map[ (int) r_str[ i - 1 ] ];

                    for ( j = 1; j < score_cols; ++j ) {
                        const long curr = ( i - 0 ) * score_cols + j,
                                   prev = ( i - 1 ) * score_cols + j;
                        // ref but not query is deletion
                        // query but not ref is insertion
                        double deletion  = score_matrix[ prev ] - open_deletion,
                               insertion = score_matrix[ curr - 1 ] - open_insertion,
                               match     = score_matrix[ prev - 1 ];

                        // if there is a match bonus or penalty, add it in
                        if ( r_char >= 0 ) {
                            const long q_char = char_map[ (int) q_str[ j - 1 ] ];

                            if ( q_char >= 0 ) {
                                match += cost_matrix[ r_char * cost_stride + q_char ];
                            }
                        }

                        // if we're doing affine gaps,
                        // look up potential moves in the affine gap matrices
                        if ( do_affine ) {
                            deletion  = MAX( deletion,
                                             deletion_matrix[ prev ] - ( i > 1 ? extend_deletion : open_deletion ) );
                            insertion = MAX( insertion,
                                             insertion_matrix[ curr - 1 ] - ( j > 1 ? extend_insertion : open_insertion ) );
                            // store the values back in the gap matrices
                            deletion_matrix[ curr ] = deletion;
                            insertion_matrix[ curr ] = insertion;
                        }

                        score_matrix[ curr ] = MAX( match, MAX( deletion, insertion ) );
                    }
                }
            }

            // set these indices to point at the ends
            // of the ref and query, respectively
            i = r_len;
            j = q_len;
            // grab maximum score from the last entry in the table
            score = score_matrix[ score_rows * score_cols - 1 ];

            // if we're doing a local alignment,
            // find the best score in the last row and column of the scoring matrix
            // and start backtracking from there ( if it's better than the score
            // we've already found, that is )
            if ( do_local ) {
                // grab the best score from the last column of the score matrix,
                // skipping the very last entry ( we already checked it )
                for ( k = score_cols - 1; k < score_rows * score_cols - 1; k += score_cols )
                    if ( score_matrix[ k ] > score ) {
                        score = score_matrix[ k ];
                        // if do_codon, k / score_cols indexes into the codon space
                        // of the reference, which is resolved by multiplication
                        // by ref_stride ( which is 3 ), otherwise this
                        // directly indexes into the reference
                        i = ref_stride * ( k / score_cols );
                    }

                // grab the best score from the last row of the score matrix,
                // skipping the very last entry ( we already checked it )
                for ( k = ( score_rows - 1 ) * score_cols; k < score_rows * score_cols - 1; ++k )
                    if ( score_matrix[ k ] > score ) {
                        score = score_matrix[ k ];
                        // if we've found a better score here,
                        // don't forget to reset the ref index
                        i = r_len;
                        // remove the initial value!
                        j = k - ( score_rows - 1 ) * score_cols;
                    }

                // fill in the edit_ops with the difference
                // between r_len and i
                for ( k = i; k < r_len; ++k )
                    edit_ops[ edit_ptr++ ] = -1;

                // fill in the edit_ops with the difference
                // between q_len and j
                for ( k = j; k < q_len; ++k )
                    edit_ops[ edit_ptr++ ] = 1;
            }

            // backtrack now

            /*
            // prints the score matrix
            for ( long m = 0; m < score_rows; ++m ) {
               for ( long n = 0; n < score_cols; ++n ) {
                   if ( n > 0 )
                       fprintf( stderr, "," );
                   fprintf( stderr, "% 3.3g", score_matrix[ m * score_cols + n ] );
               }
               fprintf( stderr, "\n" );
            }
            fprintf( stderr, "\n" );
            */

            if ( do_codon ) {
                // if either index hits 0, we're done
                // or if both indices fall below 3, we're done
                while ( i && j && ( i >= 3 || j >= 3 ) ) {
                    // perform a step
                    const long code = CodonAlignStringsStep(
                                        score_matrix
                                      , insertion_matrix
                                      , deletion_matrix
                                      , r_enc
                                      , q_enc
                                      // divide by 3 to index into codon space
                                      , ( i / 3 )
                                      , j
                                      , score_cols
                                      , char_count
                                      , miscall_cost
                                      , open_insertion
                                      , open_deletion
                                      , extend_insertion
                                      , extend_deletion
                                      , cost_matrix
                                      , cost_stride
                                      , codon3x5
                                      , codon3x4
                                      , codon3x2
                                      , codon3x1
                                      );
                    // alter edit_ops and decrement i and j
                    // according to the step k we took
                    BacktrackAlignCodon( edit_ops, &edit_ptr, &i, &j, code );

                    // if anything drops below 0, something bad happened
                    if ( i < 0 || j < 0 ) {
                        *r_res = NULL;
                        *q_res = NULL;
                        score = -A_LARGE_NUMBER;
                        goto end;
                    }

                    // handle the affine cases
                    if ( do_affine ) {
                        // divide by 3 to index into codon space
                        k = ( i / 3 ) * score_cols + j;

                        // reference matched but not query, a deletion
                        if ( code == HY_111_000 ) {
                            // while deletion is preferential to match
                            while ( i >= 3
                                    && score_matrix[ k ] - open_deletion
                                    <= deletion_matrix[ k ] - extend_deletion ) {
                                // take a codon out of the reference
                                i -= 3;
                                edit_ops[ edit_ptr++ ] = -1;
                                edit_ops[ edit_ptr++ ] = -1;
                                edit_ops[ edit_ptr++ ] = -1;
                                // move up a row in the score_matrix
                                // which is a codon in the reference
                                k -= score_cols;
                            }

                            // query matched but not reference, insertion
                        }
                        else if ( code == HY_000_111 ) {
                            // while insertion is preferential to match
                            while ( j >= 3
                                    && score_matrix[ k ] - open_insertion
                                    <= insertion_matrix[ k ] - extend_insertion ) {
                                // take a codon out of the query
                                j -= 3;
                                edit_ops[ edit_ptr++ ] = 1;
                                edit_ops[ edit_ptr++ ] = 1;
                                edit_ops[ edit_ptr++ ] = 1;
                                // move up 3 in the score_matrix
                                // which is a codon in the query
                                k -= 3;
                            }
                        }
                    }
                }
            }
            else {
                if ( do_affine ) {
                    while ( i && j ) {
                        long curr = ( i - 0 ) * score_cols + j,
                             prev = ( i - 1 ) * score_cols + j,
                             best_choice = 0;
                        // check the current affine scores and the match score
                        double scores[ 3 ] = {
                            deletion_matrix[ curr ],
                            insertion_matrix[ curr ],
                            score_matrix[ prev - 1 ]
                        }, max_score = scores[ best_choice ];
                        MatchScore( r_str, q_str, i, j, char_map, cost_matrix, cost_stride, &scores[2] );

                        // look at choice other than 0
                        for ( k = 1; k < 3; ++k )
                            if ( scores[ k ] > max_score ) {
                                max_score = scores[ k ];
                                best_choice = k;
                            }

                        switch ( best_choice ) {
                        case 0:
                            // we have at least 1 deletion
                            --i;
                            edit_ops[ edit_ptr++ ] = -1;

                            // deletion is travel in the reference but not query,
                            // look at scores back in the reference,
                            // and while they are better for the deletion case,
                            // move backwards in the reference
                            while ( i
                                    && score_matrix[ curr - score_cols ] - open_deletion
                                    <= deletion_matrix[ curr - score_cols ] - extend_deletion
                                  ) {
                                --i;
                                edit_ops[ edit_ptr++ ] = -1;
                                curr -= score_cols;
                            }

                            break;

                        case 1:
                            // we have at least 1 insertion
                            --j;
                            edit_ops[ edit_ptr++ ] = 1;

                            // insertion is travel in the query but not the reference,
                            // look at scores back in the query,
                            // and while they are better than for the insertion case,
                            // move backwards in the query
                            while ( j
                                    && score_matrix[ curr - 1 ] - open_insertion
                                    <= insertion_matrix[ curr - 1 ] - extend_insertion
                                  ) {
                                --j;
                                edit_ops[ edit_ptr++ ] = 1;
                                --curr;
                            }

                            break;

                        case 2:
                            // it's a match! move back in both
                            --i;
                            --j;
                            edit_ops[ edit_ptr++ ] = 0;
                            break;
                        }
                    }

                    // no affine gaps, no codons
                }
                else {
                    while ( i && j ) {
                        const long curr = ( i - 0 ) * score_cols + j,
                                   prev = ( i - 1 ) * score_cols + j;
                        double deletion  = score_matrix[ prev ] - open_deletion,
                               insertion = score_matrix[ curr - 1 ] - open_insertion,
                               match     = score_matrix[ prev - 1 ];
                        MatchScore( r_str, q_str, i, j, char_map, cost_matrix, cost_stride, &match );
                        BacktrackAlign( edit_ops, &edit_ptr, &i, &j, deletion, insertion, match );
                    }
                }
            }

            // for anything that remains,
            // don't forget it!!!

            // reference
            while ( --i >= 0 )
                edit_ops[ edit_ptr++ ] = -1;

            // then query
            while ( --j >= 0 )
                edit_ops[ edit_ptr++ ] = 1;

            if ( edit_ptr > 0 ) {
                // reset indices to 0
                i = j = 0;
                // rebuild the strings from the edit_ops
                // with room for the null terminator
                *r_res = ALLOCA( char, edit_ptr + 1 );
                *q_res = ALLOCA( char, edit_ptr + 1 );

                if ( ISNULL( *r_res ) || ISNULL( *q_res ) ) {
                    free( *r_res );
                    free( *q_res );
                    *r_res = NULL;
                    *q_res = NULL;
                    score = 0.;
                    goto end;
                }

                for ( --edit_ptr, k = 0; edit_ptr >= 0; --edit_ptr, ++k ) {
                    switch ( edit_ops[ edit_ptr ] ) {
                        // match! include characters from both strings
                    case 0:
                        (*r_res)[ k ] = r_str[ i++ ];
                        (*q_res)[ k ] = q_str[ j++ ];
                        break;

                        // insertion!
                    case 1:
                        (*r_res)[ k ] = gap;
                        (*q_res)[ k ] = tolower( q_str[ j++ ] );
                        break;

                    case 2:
                        (*r_res)[ k ] = gap;
                        (*q_res)[ k ] = tolower( q_str[ j++ ] );
                        break;

                    case -1:
                        (*r_res)[ k ] = tolower( r_str[ i++ ] );
                        (*q_res)[ k ] = gap;
                        break;

                    case -2:
                        (*r_res)[ k ] = tolower( r_str[ i++ ] );
                        (*q_res)[ k ] = gap;
                        break;
                    }
                }

                // make sure to null-terminate
                (*r_res)[ k ] = '\0';
                (*q_res)[ k ] = '\0';
            }

end:
#if 0
            {
                double max_score = -A_LARGE_NUMBER;

                if ( do_codon ) {
                    long l, m, n;

                    for ( i = 0; i < char_count; ++i ) {
                        for ( j = 0; j < char_count; ++j ) {
                            for ( k = 0; k < char_count; ++k ) {
                                const long r_codon = ( i * char_count + j ) * char_count + k;
                                for ( l = 0; l < char_count; ++l ) {
                                    for ( m = 0; m < char_count; ++m ) {
                                        for ( n = 0; n < char_count; ++n ) {
                                            const long q_codon = ( l * char_count + m ) * char_count + n;
                                            const double v = cost_matrix[ r_codon * cost_stride + q_codon ];
                                            if ( v > max_score )
                                                max_score = v;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                else {
                    for ( i = 0; i < char_count; ++i ) {
                        const long r_char = char_map[ i ];
                        for ( j = 0; j < char_count; ++j ) {
                            const long q_char = char_map[ j ];
                            const double v = cost_matrix[ r_char * cost_stride + q_char ];
                            if ( v > max_score )
                                max_score = v;
                        }
                    }
                }
   
                fprintf( stderr, "max_score: %g\n", max_score );

                const double scorep = score / ( 1.0 * q_len / ref_stride );
    
                if ( scorep > max_score ) {
                    fprintf( stderr, "\nscore: %g, seq: %s\n", score / ( 1.0 * q_len / ref_stride ), *q_res );
                    print_score_matrix( stderr, score_matrix, score_rows, score_cols );
                }
            }
#endif

            free( edit_ops );
#if 0
            free( score_matrix );
            free( deletion_matrix );
            free( insertion_matrix );
#endif
            free( r_enc );
            free( q_enc );
        }
    }

    return score;
}
