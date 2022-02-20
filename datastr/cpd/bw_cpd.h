//
// Created by Bojie Shen on 17/9/21.
//

#ifndef KATCH_BW_CPD_H
#define KATCH_BW_CPD_H
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
    class BW_CPD {
    public:
        BW_CPD() : begin{0} {}

        //! Adds a new node s to the CPD. first_move should be an array that
        //! maps every target node onto a 15-bit bitfield that has a bit set
        //! for every valid first move. get_first_move is free to return any of
        //! them.

        void append_row(unsigned source_node, const vector<bool> &allowed_first_move);

        void append_row(unsigned source_node, const vector<unsigned short> &allowed_first_move);

        void append_compressed_row(const vector<unsigned>& compressed_row);

        void append_rows(const BW_CPD&other);

        void append_rows(const BW_CPD &other, unsigned row_id);

        bool get_reachability_linear_search(unsigned source_node, unsigned target_node) const {
            target_node <<= 1;
            target_node |= 0001;
            for(auto it = entry.begin() + begin[source_node];  it != entry.begin() + begin[source_node + 1];++it){
                if(*it > target_node){
                    return *(it-1) & 0001;
                }
            }
            return  *(entry.begin() + begin[source_node + 1] - 1 )  & 0001;
        }

        //! Get the first move.
        //! An ID of 0xF means that there is no path.
        //! If source_node == target_node then return value is undefined.
        bool get_reachability(unsigned source_node, unsigned target_node) const {
            target_node <<= 1;
            target_node |= 0001;
            return *binary_find_last_true(
                    entry.begin() + begin[source_node],
                    entry.begin() + begin[source_node + 1],
                    [=](unsigned x) { return x <= target_node; }
            ) & 0001;
        }



        unsigned node_count() const {
            // get the number of rows
            // in inverse centroid cpd, this returns the number of centroids.
            return begin.size() - 1;
        }

        unsigned entry_count() const {
            return entry.size();
        }

        friend bool operator==(const BW_CPD &l, const BW_CPD &r) {
            return l.begin == r.begin && l.entry == r.entry;
        }

        friend bool operator!=(const BW_CPD &l, const BW_CPD &r) {
            return !(l == r);
        }

        void save(std::FILE *f) const {
            save_vector(f, begin);
            save_vector(f, entry);
        }

        void load(std::FILE *f) {
            begin = load_vector<unsigned>(f);
            entry = load_vector<unsigned>(f);
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

        void push_back_row(const vector<unsigned> &compressed_row) {
            entry.insert(entry.end(), compressed_row.begin(), compressed_row.end());
            begin.push_back(entry.size());
        }

    protected:
        std::vector<unsigned> begin;
        std::vector<unsigned> entry;


    };
}

#endif //KATCH_BW_CPD_H
