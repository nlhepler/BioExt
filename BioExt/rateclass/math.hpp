
#include <cmath>
#include <vector>
#include <utility>

#ifndef MATH_H
#define MATH_H

namespace math
{
    inline
    double lg_choose( const int n, const int k )
    {
        double rv = 0.0;

        for ( int i = 1; i <= k; ++i )
            rv += std::log( n - k + i ) - std::log( i );

        return rv;
    }


    inline
    double prob_background( const double lg_bg, const double lg_invbg, const int cov, const int k )
    {
        double rv = std::exp( cov * lg_invbg );

        for ( int i = 1; i < k; ++i )
            rv += std::exp( lg_choose( cov, i ) + i * lg_bg + ( cov - i ) * lg_invbg );

        rv = 1.0 - rv;

        return ( rv < 0.0 ) ? 0.0 : rv;
    }


    inline
    double weighted_harmonic_mean( std::vector< std::pair< double, double > > & xs )
    {
        double num = 0.0, den = 0.0;
        std::vector< std::pair< double, double > >::const_iterator it;

        for ( it = xs.begin(); it != xs.end(); ++it ) {
            num += it->first;
            den += it->first / it->second;
        }

        return num / den;
    }
}

#endif // MATH_H
