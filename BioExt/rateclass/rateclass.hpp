
#include <utility>
#include <vector>

#ifndef RATECLASS_H
#define RATECLASS_H

namespace rateclass
{
    void params_json_dump(
            FILE * const file,
            const double lg_L,
            const double aicc,
            const std::vector< std::pair< double, double > > & params,
            const double bg = 0.0
            );

    class rateclass_t
    {
    private:
        const std::vector< std::pair< int, int > > & data;
        const int factor;

    public:
        rateclass_t( const std::vector< std::pair< int, int > > & data, const int factor = 1 );
        void learn(
            double & lg_L,
            double & aicc,
            std::vector< std::pair< double, double > > & params,
            const int nrestart = 50
            ) const;
    };
}

#endif // RATECLASS_H
