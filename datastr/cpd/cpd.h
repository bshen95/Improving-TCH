//
// Created by Bojie Shen on 8/8/20.
//



#ifndef  KATCH_CPD_H
#define  KATCH_CPD_H

#pragma once
#include <vector>
#include <set>
#include <algorithm>
#include <string>
#include <iostream>
#include "binary_search.h"
#include "range.h"
#include "io/vec_io.h"
#include <unordered_set>
#include <map>

using namespace std;
//! Compressed Path database. Allows to quickly query the first out arc id of
//! any shortest source-target-path. There may be at most 15 outgoing arcs for
//! any node.
namespace katch {
    class CPD {
    public:
        CPD() : begin{0} {}

        //! Adds a new node s to the CPD. first_move should be an array that
        //! maps every target node onto a 15-bit bitfield that has a bit set
        //! for every valid first move. get_first_move is free to return any of
        //! them.

        void append_row(unsigned source_node, const vector<unsigned short> &allowed_first_move);

        void append_compressed_row(const vector<unsigned>& compressed_row);

        std::vector<unsigned> dbg_compress_row( const vector<unsigned short> &allowed_first_move);

        void append_row_with_landmark(unsigned source_node, const vector<unsigned short> &allowed_first_move);

        void append_row_using_mapper(unsigned source_node, const vector<unsigned short> &allowed_first_move);

        void
        append_row_with_landmark_using_mapper(unsigned source_node, const vector<unsigned short> &allowed_first_move);

        void append_rows(const CPD &other);

        void append_rows(const CPD &other, unsigned row_id);

        //! Get the first move.
        //! An ID of 0xF means that there is no path.
        //! If source_node == target_node then return value is undefined.
        unsigned char get_first_move(unsigned source_node, unsigned target_node) const {
            target_node <<= 8;
            target_node |= 0xFF;
            return *binary_find_last_true(
                    entry.begin() + begin[source_node],
                    entry.begin() + begin[source_node + 1],
                    [=](unsigned x) { return x <= target_node; }
            ) & 0xFF;
        }

        unsigned node_count() const {
            // get the number of rows
            // in inverse centroid cpd, this returns the number of centroids.
            return begin.size() - 1;
        }

        unsigned entry_count() const {
            return entry.size();
        }

        friend bool operator==(const CPD &l, const CPD &r) {
            return l.begin == r.begin && l.entry == r.entry;
        }

        friend bool operator!=(const CPD &l, const CPD &r) {
            return !(l == r);
        }

        void save(std::FILE *f) const {
            save_vector(f, begin);
            save_vector(f, entry);
        }

        void load(std::FILE *f) {
            begin = load_vector<unsigned>(f);
            entry = load_vector<unsigned>(f);
//            get_row_distribution(0);
//            get_row_distribution(100);
//            get_row_distribution(1000);
//            get_row_distribution(10000);

//            print_symbol_distribution();
        }

        unsigned long long get_entry_size() const {
            return entry.size();
        }

        unsigned long long get_row_size(unsigned row_id) const {
            return begin[row_id + 1] - begin[row_id];
        }


        const vector<unsigned> &get_entry(unsigned row_number) const {
            entry.begin() + begin[row_number];
            entry.begin() + begin[row_number + 1];

            return entry;
        }

        const vector<unsigned> &get_entry() const {
            return entry;
        }

        const vector<unsigned> &get_begin() const {
            return begin;
        }

        void print_symbol_distribution(){
            std::map<unsigned short, unsigned> symbol_map;
            for(unsigned i = 0; i < begin.size()-1; i ++){
                vector<unsigned short> d_row = decode_row(i);
                for(const auto& symbol : d_row) {
                    if (symbol_map.find(symbol) != symbol_map.end()) {
                        symbol_map[symbol]++;
                    } else {
                        symbol_map.insert({symbol, 0});
                    }
                }
            }
            vector<pair<unsigned short, unsigned> > symbol_vec;
            for (auto& it :symbol_map) {
                symbol_vec.push_back(it);
            }

            std::sort(
                    symbol_vec.begin(), symbol_vec.end(),
                    [](const std::pair<unsigned short, unsigned> &lhs,
                       const std::pair<unsigned short, unsigned> &rhs) { return lhs.second > rhs.second; }
            );
            for(auto it = symbol_vec.begin(); it != symbol_vec.end(); ++it) {
                std::cout <<"Symbols: "<< it->first << " Number: " << it->second <<"\n";
            }
        }

        vector<pair<unsigned short, unsigned>>  get_row_distribution(unsigned row_number){
            vector<unsigned short> d_row = decode_row(row_number);
            std::map<unsigned short, unsigned> symbol_map;

            for(const auto& symbol : d_row){
                if(symbol_map.find(symbol) != symbol_map.end()){
                    symbol_map[symbol] ++;
                }else{
                    symbol_map.insert({ symbol, 1 });
                }
            }
            vector<pair<unsigned short, unsigned> > symbol_vec;
            for (auto& it :symbol_map) {
                symbol_vec.push_back(it);
            }

            std::sort(
                    symbol_vec.begin(), symbol_vec.end(),
                    [](const std::pair<unsigned short, unsigned> &lhs,
                            const std::pair<unsigned short, unsigned> &rhs) { return lhs.second > rhs.second; }
            );
            for(auto it = symbol_vec.begin(); it != symbol_vec.end(); ++it) {
                std::cout <<"Symbols: "<< it->first << " Number: " << it->second <<"\n";
            }

            return symbol_vec;



        }


        vector<unsigned short> decode_row(unsigned row_number){
            vector<unsigned short>decoded_row;
            for(unsigned b = begin[row_number]; b< begin[row_number+1]; ++b ) {
                unsigned cur_start = entry[b] >> 8;
                unsigned short cur_symbol = entry[b] &  0xFF;
                unsigned next_start = b < begin[row_number+1]-1 ? entry[b+1] >> 8 : node_count();
                for (unsigned s_index = cur_start; s_index < next_start; s_index++) {
                    decoded_row.push_back(cur_symbol);
                }
            }
            assert(decoded_row.size() == node_count());
//            vector<unsigned> compressed_row = dbg_compress_row(decoded_row);
//            for(int i = 0; i < compressed_row.size(); i++){
//                if(compressed_row[i] != entry[begin[row_number]+i]){
//                    std::cout<<"error"<<std::endl;
//                }
//            }
//            assert(compressed_row.back() == entry[begin[row_number+1]-1]);
            return decoded_row;
        }

        vector<unsigned> get_compressed_runs(unsigned row_number){
            vector<unsigned >compressed_runs;
            for(unsigned b = begin[row_number]; b< begin[row_number+1]; ++b ) {
                compressed_runs.push_back(entry[b]);
            }
            return compressed_runs;
        }


        vector<unsigned short> get_compressed_symbol(unsigned row_number){
            vector<unsigned short>compressed_symbol;
            for(unsigned b = begin[row_number]; b< begin[row_number+1]; ++b ) {
                unsigned cur_start = entry[b] >> 8;
                unsigned short cur_symbol = entry[b] &  0xFF;
                compressed_symbol.push_back(cur_symbol);
            }
            return compressed_symbol;
        }

        vector<unsigned short> getIntersection( const std::vector<unsigned short> &s1,  const std::vector<unsigned short>& s2){
            // assume s1 and s2 only contains unique elements.
            vector<unsigned short> intersection;
            for(const auto& e_s1 : s1){
                for(const auto& e_s2 : s2){
                    if(e_s1 == e_s2) intersection.push_back(e_s1);
                }
            }

            return intersection;
        }



        pair<int, set<unsigned short>> getIntersection(set<unsigned short> &s1, const set<unsigned short> &s2) {
            set<unsigned short> intersect = set<unsigned short>();
            if (s2.empty() && s1.empty()) {
                return pair<int, set<unsigned short>>(1, intersect);
            } else if (s2.empty()) {
                return pair<int, set<unsigned short>>(1, s1);
            } else if (s1.empty()) {
                return pair<int, set<unsigned short>>(1, s2);
            }
            set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                             std::inserter(intersect, intersect.begin()));
            return pair<int, set<unsigned short>>(0, intersect);
        }


        void append_row_multiple_symbols(unsigned source_node, const vector<set<unsigned short>> &set);

        void append_row_multiple_symbols(unsigned source_node, const vector<unsigned short> & cpd_result);

        void push_back_row(const vector<unsigned> &compressed_row) {
            entry.insert(entry.end(), compressed_row.begin(), compressed_row.end());
            begin.push_back(entry.size());
        }

    protected:
        std::vector<unsigned> begin;
        std::vector<unsigned> entry;

        const unsigned &get_allowed(unsigned x, unsigned s, const vector<unsigned> &fmoves) const;

        const set<unsigned> &
        get_allowed_multiple_row(unsigned x, unsigned s, const vector<set<unsigned>> &fmoves) const;


    };
}

#endif //KATCH_CPD_H
