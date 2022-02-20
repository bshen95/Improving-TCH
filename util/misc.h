/*
 * katch/util/misc.h
 *
 *
 * This file is part of
 *
 * KaTCH -- Karlsruhe Time-Dependent Contraction Hierarchies
 *
 * Copyright (C) 2015
 *
 * Institut fuer Theroretische Informatik,
 * Karlsruher Institut fuer Technology (KIT),
 * 76131 Karlsruhe, Germany
 *
 * Author: Gernot Veit Batz
 *
 *
 *
 * KaTCH is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaTCH is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with Contraction Hierarchies; see the file COPYING;
 * if not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef KATCH_MISC_H_
#define KATCH_MISC_H_

#include <algorithm>
#include <assert.h>
#include <cstdlib>
#include <chrono>
#include <functional>
#include <limits>
#include <vector>
#include <cmath>
#include "datastr/graph/basic.h"

#define KATCH_IMPLIES(a, b) ((!(a)) || (b))
#define KATCH_EQUIV(a, b) (KATCH_IMPLIES((a),(b)), KATCH_IMPLIES((b),(a)))

#define KATCH_BE_VERBOSE true

#define KATCH_STATUS(s) if (KATCH_BE_VERBOSE) std::cout << "STATUS [" << __FILE__ << ":" << __LINE__  << "] " << s << std::flush
#define KATCH_WARNING(s) if (KATCH_BE_VERBOSE) std::cerr << "WARNING [" << __FILE__ << ":" << __LINE__  << "] " << s << std::flush
#define KATCH_ERROR(s) if (KATCH_BE_VERBOSE) std::cerr << "ERROR [" << __FILE__ << ":" << __LINE__  << "] " << s << std::flush

#define KATCH_CONTINUE_STATUS(s) if (KATCH_BE_VERBOSE) std::cout << s << std::flush
#define KATCH_CONTINUE_ERROR(s) if (KATCH_BE_VERBOSE) std::cerr << s << std::flush
#define KATCH_CONTINUE_WARNING(s) if (KATCH_BE_VERBOSE) std::cerr << s << std::flush

namespace katch
{
namespace util
{

// Check whether a range (of "less" and "equal_to" comparable elements) contains duplicates
template
<
    typename Iterator,
    typename LtFn = std::less<typename std::iterator_traits<Iterator>::value_type>,
    typename EqFn = std::equal_to<typename std::iterator_traits<Iterator>::value_type>
>
bool has_duplicates
(
        const Iterator& begin,
        const Iterator& end,
        const LtFn& lt = std::less<typename std::iterator_traits<Iterator>::value_type>(),
        const EqFn& eq = std::equal_to<typename std::iterator_traits<Iterator>::value_type>()
)
{
    std::vector<typename std::iterator_traits<Iterator>::value_type> vec(begin, end);

    std::sort(vec.begin(), vec.end(), lt);
    return std::adjacent_find(vec.begin(), vec.end(), eq) != vec.end();
}

// Check whether a vector (of "less" and "equal_to" comparable elements) contains duplicates
template <typename T>
bool has_duplicates(const std::vector<T>& vec)
{
    return has_duplicates(vec.begin(), vec.end());
}


// Remove duplicates from a vector of "less" and "equal_to" comparable elements
template <typename T, typename LtFn = std::less<T>, typename EqFn = std::equal_to<T>>
void remove_duplicates(std::vector<T>& vec, const LtFn& lt = std::less<T>(), const EqFn& eq = std::equal_to<T>())
{
    std::sort(vec.begin(), vec.end(), lt);
    auto new_end = std::unique(vec.begin(), vec.end(), eq);

    vec.resize( std::distance(vec.begin(), new_end) );
    assert( ! has_duplicates(vec) );
}

template <typename Graph, typename Iterator>
std::vector<NodeIterator> neighborhood(const Graph& graph, Iterator begin, Iterator end)
{
    std::vector<NodeIterator> result;

    for ( auto it = begin ; it != end ; ++it )
    {
        const auto& x = *it;
        std::vector<NodeIterator> neighbors_of_x;

        for ( auto e = graph.out_edges_begin(x) ; e != graph.out_edges_end(x) ; ++e )
            neighbors_of_x.push_back(graph.get_target(e));

        for ( auto e = graph.in_edges_begin(x) ; e != graph.in_edges_end(x) ; ++e )
            neighbors_of_x.push_back(graph.get_source(e));

        remove_duplicates(neighbors_of_x);
        result.insert(result.end(), neighbors_of_x.begin(), neighbors_of_x.end());
    }

    return result;
}

// Compute the position of the most significant bit (where the right most position is 0)
inline uint32_t msb(const uint32_t x)
{
    int digit = std::numeric_limits<uint32_t>::digits - 1;
    while( (x & (1 << digit)) == 0 ) --digit;

    assert( digit >= 0 );
    return digit;
}

// Compute x mod m
inline double modulo(const double x, const double m)
{
    if (x >= 0)
    {
        if ( x < m ) return x;
        if ( x < m+m ) return x - m;
    }

    double result = fmod(x, m);
    if (result < 0) result += m;

    assert( 0 <= result );
    assert( result < m );
    assert( fabs(result - (x - m * floor(x/m))) < 0.000001 );

    return result;
}

// compute x mod m
inline int modulo(const int x, const int m)
{
    int result = x % m;
    if (result < 0) result += m;

    assert(result >= 0);
    assert( result == (x - m * floor((double)x/(double)m)) );

    return result;
}


double random(const double a, const double b)
{
    assert( a < b );
    double result =  (double(rand()) / double(RAND_MAX)) * (b - a) + a;

    assert( a <= result + 0.000001 );
    assert( result <= b + 0.000001 );
    return result;
}

using TimePoint = std::chrono::high_resolution_clock::time_point;

inline TimePoint time_stamp()
{
    return std::chrono::high_resolution_clock::now();
}

inline double get_duration_in_seconds(const TimePoint& from, const TimePoint& to)
{
    assert( from <= to );
    return std::chrono::duration_cast<std::chrono::duration<double>>(to - from).count();
}


inline double geo_dist(double lat_a, double lon_a, double lat_b, double lon_b){

// The formula used in this function was derived as following:
//
// lat_a = input
// lon_a = input
// lat_b = input
// lon_b = input
//
// ----- Step 1: Convert from degree to radian
//
// lat_a /= 180
// lat_a *= pi
// lon_a /= 180
// lon_a *= pi
// lat_b /= 180
// lat_b *= pi
// lon_b /= 180
// lon_b *= pi
//
// ----- Step 2: Convert to 3D unit vectors
//
// x_a = cos(lon_a) * cos(lat_a)
// y_a = sin(lon_a) * cos(lat_a)
// z_a = sin(lat_a)
//
// x_b = cos(lon_b) * cos(lat_b)
// y_b = sin(lon_b) * cos(lat_b)
// z_b = sin(lat_b)
//
// ----- Step 3: Compute scalar product
//
// c = x_a * x_b + y_a * y_b + z_a * z_b
//
// ----- Step 4: acos(c) is length in radian, scaling up by earth radius gives distance
//
// length = earth_radius*acos(c)
// Output length
//
// ----- Steps 2 to 4 merged in one formula is
//
// length = R*acos(
//	cos(lon_a) * cos(lat_a) * cos(lon_b) * cos(lat_b) +
//	sin(lon_a) * cos(lat_a) * sin(lon_b) * cos(lat_b) +
//	sin(lat_a) * sin(lat_b)
// )
//
// ----- It can be simplified as follows using standard trigonometric identities and basic math operations
// See https://en.wikipedia.org/wiki/List_of_trigonometric_identities
//
// = R*acos(
//	cos(lat_a) * cos(lat_b) * (
//		cos(lon_a) * cos(lon_b) +
//		sin(lon_a) * sin(lon_b)
//	)+
//	sin(lat_a) * sin(lat_b)
// ) = R*acos(
//	cos(lat_a) * cos(lat_b) * (
//		1/2 * (cos(lon_a - lon_b) + cos(lon_a + lon_b)) +
//		1/2 * (cos(lon_a - lon_b) - cos(lon_a + lon_b))
//	)+
//	sin(lat_a) * sin(lat_b)
// ) = R*acos(
//	1/2 * cos(lat_a) * cos(lat_b) * (
//		cos(lon_a - lon_b) + cos(lon_a + lon_b) +
//		cos(lon_a - lon_b) - cos(lon_a + lon_b)
//	)+
//	sin(lat_a) * sin(lat_b)
// ) = R*acos(
//	1/2 * cos(lat_a) * cos(lat_b) * (
//		cos(lon_a - lon_b) +
//		cos(lon_a - lon_b)
//	)+
//	sin(lat_a) * sin(lat_b)
// ) = R*acos(
//	cos(lat_a) * cos(lat_b) * cos(lon_a - lon_b) +
//	sin(lat_a) * sin(lat_b)
// ) = R*acos(
//	1/2 * (cos(lat_a - lat_b) + cos(lat_a + lat_b)) * cos(lon_a - lon_b) +
//	1/2 * (cos(lat_a - lat_b) - cos(lat_a + lat_b))
// ) = R*acos(
//	1/2 * (
//		(cos(lat_a - lat_b) + cos(lat_a + lat_b)) * cos(lon_a - lon_b) +
//		(cos(lat_a - lat_b) - cos(lat_a + lat_b))
//	)
// )

        const double pi_div_180 = 3.14159265359/180.0;
        const double earth_radius = 6371000.785; // in meter

        // To help the auto-vectorizer figure out that using SIMD is a good idea here, we use a local array with a constant size instead of variables.
        // GCC 5's auto vectorizer requires that the array has a power of two size.
        // We thus add a fourth dummy value.
        //
        // If the code was vectorized then the compiler would have to insert this dummy value implicitly. We are helping it by doing it explicitly.
        // If SIMD is disabled, the array should be converted by the optimizer to 4 variables. The dummy variable should subsequently be eliminated using dead code elimination.
        double vec[4] = {
                lat_a - lat_b,
                lat_a + lat_b,
                lon_a - lon_b,
                0
        };

        // It is important that this loop has four iterations, as otherwise the auto-vectorizer does not do its job.
        for(unsigned i=0; i<4; ++i)
            vec[i] = cos(vec[i]*pi_div_180);

        double len = (vec[0] + vec[1]) * vec[2] + vec[0] - vec[1];
        len *= 0.5;
        len  = acos(len);
        len *= earth_radius;

        // len is in meter

        return len;
    }




} /* namespace util */
} /* namespace katch */

#endif /* KATCH_MISC_H_ */
