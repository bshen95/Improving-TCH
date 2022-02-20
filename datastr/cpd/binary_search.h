//
// Created by Bojie Shen on 8/8/20.
//

#ifndef KATCH_SEARCH_H
#define KATCH_SEARCH_H

#include <cassert>

template<class Iter, class Pred>
Iter binary_find_first_true(Iter begin, Iter end, Pred p){
    if(begin == end)
        return end;
    if(!p(*(end-1)))
        return end;

    while(end - begin > 1){
        Iter mid = begin + (end-begin-1)/2;

        if(p(*mid))
            end = mid+1;
        else
            begin = mid+1;
    }
    return begin;
}


template<class Iter, class Pred>
Iter binary_find_last_true(Iter begin, Iter end, Pred p){
    assert(begin != end);
    assert(p(*begin));

    while(end - begin > 1){
        Iter mid = begin + (end-begin)/2;

        if(p(*mid))
            begin = mid;
        else
            end = mid;
    }
    return begin;
}



#endif //KATCH_BINARY_SEARCH_H
