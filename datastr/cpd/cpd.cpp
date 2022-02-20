//
// Created by Bojie Shen on 8/8/20.
//


#include "cpd.h"

// compile with -O3 -DNDEBUG
namespace  katch {
    void CPD::append_row_multiple_symbols(
            unsigned source_node, const vector<set<unsigned short>> &allowed_first_move) {

        auto get_allowed_local = [&](unsigned x) {
            return allowed_first_move[x];
        };

        unsigned int node_begin = 0;
        set<unsigned short> allowed_up_to_now = get_allowed_local(0);

        for (unsigned i = 1; i < allowed_first_move.size(); ++i) {
            pair<int, set<unsigned short>> allowed_next = getIntersection(allowed_up_to_now, get_allowed_local(i));
            if (allowed_next.first == 0 && allowed_next.second.empty()) {

                entry.push_back((node_begin << 8) | *allowed_up_to_now.begin());
                node_begin = i;
                allowed_up_to_now = get_allowed_local(i);
            } else
                allowed_up_to_now = allowed_next.second;
        }
        entry.push_back((node_begin << 8) | *allowed_up_to_now.begin());

        begin.push_back(entry.size());
    }

    void CPD::append_row_multiple_symbols(
            unsigned source_node, const vector<unsigned short> &allowed_first_move) {

        auto get_allowed_local = [&](unsigned x) {
            return allowed_first_move[x*2+1] == 0xFC ?
            std::vector<unsigned short>{allowed_first_move[x*2]} :
            std::vector<unsigned short>{allowed_first_move[x*2],allowed_first_move[x*2+1]};
        };

        unsigned int node_begin = 0;
        std::vector<unsigned short> allowed_up_to_now =  get_allowed_local(0);
        if (allowed_up_to_now[0] == 0xFE) {
            allowed_up_to_now = get_allowed_local(1);
        }

        for (unsigned i = 1; i < allowed_first_move.size()/2; ++i) {
            std::vector<unsigned short>  allowed_next = get_allowed_local(i);
            if(allowed_next[0] == 0xFE ){
                continue;
            }
            std::vector<unsigned short>  intersection = getIntersection(allowed_up_to_now, allowed_next);
            if(intersection.empty()){
                //no intersection, always take first element, try not take 0xFC symbol.
                entry.push_back((node_begin << 8) | allowed_up_to_now[0]);
                node_begin = i;
                allowed_up_to_now = allowed_next;
            }else{
                allowed_up_to_now = intersection ;
            }
        }
        entry.push_back((node_begin << 8) | allowed_up_to_now[0]);

        begin.push_back(entry.size());
    }


//    void CPD::append_row(
//            unsigned source_node, const vector<unsigned short> &allowed_first_move) {
//
//        auto get_allowed_local = [&](unsigned x) {
//            return allowed_first_move[x];
//        };
//
//        unsigned node_begin = 0;
//        unsigned short allowed_up_to_now = get_allowed_local(0);
//
//        if (allowed_up_to_now == 0xFE) {
//            allowed_up_to_now = get_allowed_local(1);
//        }
//        for (unsigned i = 1; i < allowed_first_move.size(); ++i) {
//            unsigned short local = get_allowed_local(i);
//            if (allowed_up_to_now != get_allowed_local(i)) {
//                if (local == 0xFE) {
//                    continue;
//                } else {
//                    entry.push_back((node_begin << 8) | allowed_up_to_now);
//                    node_begin = i;
//                    allowed_up_to_now = local;
//                }
//
//            }
//        }
//        entry.push_back((node_begin << 8) | allowed_up_to_now);
//
//        begin.push_back(entry.size());
//    }



    void CPD::append_row(
            unsigned source_node, const vector<unsigned short> &allowed_first_move) {
        unsigned node_begin = 0;
        unsigned short allowed_up_to_now = allowed_first_move[0];

        unsigned cur_it = 0;
        while(allowed_up_to_now == 0xFE){
            allowed_up_to_now = allowed_first_move[++cur_it];
        }
        for (; cur_it < allowed_first_move.size(); ++cur_it) {
            unsigned short local = allowed_first_move[cur_it];
            if (allowed_up_to_now != allowed_first_move[cur_it]) {
                if (local ==0xFE) {
                    continue;
                }else{
                    entry.push_back((node_begin << 8) | allowed_up_to_now);
                    node_begin = cur_it;
                    allowed_up_to_now = local;
                }
            }
        }
        entry.push_back((node_begin << 8) | allowed_up_to_now);
        begin.push_back(entry.size());
    }

    std::vector<unsigned> CPD::dbg_compress_row( const vector<unsigned short> &allowed_first_move) {
        vector<unsigned> compressed_runs;
        unsigned node_begin = 0;
        unsigned short allowed_up_to_now = allowed_first_move[0];

        if (allowed_up_to_now == 0xFE) {
            allowed_up_to_now = allowed_first_move[1];
        }
        for (unsigned i = 1; i < allowed_first_move.size(); ++i) {
            unsigned short local = allowed_first_move[i];
            if (allowed_up_to_now != allowed_first_move[i]) {
                if (local == 0xFE) {
                    continue;
                } else {
                    compressed_runs.push_back((node_begin << 8) | allowed_up_to_now);
                    node_begin = i;
                    allowed_up_to_now = local;
                }

            }
        }
        compressed_runs.push_back((node_begin << 8) | allowed_up_to_now);
        return compressed_runs;
    }

    void CPD::append_rows(const CPD &other) {
        unsigned offset = begin.back();
        for (auto x:make_range(other.begin.begin() + 1, other.begin.end()))
            begin.push_back(x + offset);
        std::copy(other.entry.begin(), other.entry.end(), back_inserter(entry));
    }

    void CPD::append_rows(const CPD &other, unsigned row_id) {
//        int offset = begin.back();
//
//        for (auto x:make_range(other.begin.begin() + 1, other.begin.end()))
//            begin.push_back(x + offset);
//
        unsigned start = other.begin[row_id];
        unsigned end = other.begin[row_id + 1];
        std::copy(other.entry.begin() + start, other.entry.begin() + end, back_inserter(entry));
        begin.push_back(entry.size());
    }



    void CPD::append_compressed_row(const vector<unsigned>& compressed_row){
        std::copy(compressed_row.begin(), compressed_row.end(), back_inserter(entry));
        begin.push_back(entry.size());
    }


    const set<unsigned> &
    CPD::get_allowed_multiple_row(unsigned x, unsigned s, const vector<set<unsigned>> &fmoves) const {
        return fmoves[x];
    }


    const unsigned &CPD::get_allowed(unsigned x, unsigned s, const vector<unsigned> &fmoves) const {
        return fmoves[x];
    }

}
