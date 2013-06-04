
#ifndef BAMITER_H
#define BAMITER_H

namespace util
{
    template < class T, class U, class V >
    class triple
    {
    public:
        T first;
        U second;
        V third;
    
        triple()
        {
        }

        triple( const triple< T, U, V > & tpl ) :
            first( tpl.first ),
            second( tpl.second ),
            third( tpl.third )
        {
        }

        triple( const T & first, const U & second, const V & third ) :
            first( first ),
            second( second ),
            third( third )
        {
        }
    };

    template < class T, class U, class V >
    triple< T, U, V > make_triple(
        const T & first,
        const U & second,
        const V & third
        )
    {
        return triple< T, U, V >( first, second, third );
    }
}

#endif // BAMITER_H
