//
// Created by Bojie Shen on 24/10/21.
//


#include "rev_bw_cpd.h"

// compile with -O3 -DNDEBUG
namespace  katch {

    void  REV_BW_CPD::append_row(
            unsigned source_node, const vector<bool> &allowed_first_move) {

        unsigned node_begin = 0;
        bool allowed_up_to_now = allowed_first_move[0];
        if(source_node == 0 ){
            allowed_up_to_now = allowed_first_move[1];
        }
        for (unsigned i = 1; i < allowed_first_move.size(); ++i) {
            bool local = allowed_first_move[i];
            if (allowed_up_to_now != allowed_first_move[i]) {
                if (i == source_node) {
                    continue;
                }else{
                    entry.push_back((node_begin << 1) | allowed_up_to_now);
                    node_begin = i;
                    allowed_up_to_now = local;
                }
            }
        }
        entry.push_back((node_begin << 1) | allowed_up_to_now);

        begin.push_back(entry.size());
    }



    void REV_BW_CPD::append_row(
            unsigned source_node, const vector<unsigned short> &allowed_first_move) {
        unsigned node_begin = 0;
        unsigned short allowed_up_to_now = allowed_first_move[0];

        unsigned cur_it = 0;
        while(allowed_up_to_now == 2){
            allowed_up_to_now = allowed_first_move[++cur_it];
        }
        for (; cur_it < allowed_first_move.size(); ++cur_it) {
            unsigned short local = allowed_first_move[cur_it];
            if (allowed_up_to_now != allowed_first_move[cur_it]) {
                if (local == 2) {
                    continue;
                }else{
                    entry.push_back((node_begin << 1) | (bool)allowed_up_to_now);
                    node_begin = cur_it;
                    allowed_up_to_now = local;
                }
            }
        }
        entry.push_back((node_begin << 1) | (bool)allowed_up_to_now);
        begin.push_back(entry.size());
    }


    void  REV_BW_CPD::append_rows(const  REV_BW_CPD &other) {
        unsigned offset = begin.back();
        for (auto x:make_range(other.begin.begin() + 1, other.begin.end()))
            begin.push_back(x + offset);
        std::copy(other.entry.begin(), other.entry.end(), back_inserter(entry));
    }

    void  REV_BW_CPD::append_rows(const  REV_BW_CPD &other, unsigned row_id) {
        unsigned start = other.begin[row_id];
        unsigned end = other.begin[row_id + 1];
        std::copy(other.entry.begin() + start, other.entry.begin() + end, back_inserter(entry));
        begin.push_back(entry.size());
    }



    void  REV_BW_CPD::append_compressed_row(const vector<unsigned>& compressed_row){
        std::copy(compressed_row.begin(), compressed_row.end(), back_inserter(entry));
        begin.push_back(entry.size());
    }

}
