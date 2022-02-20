//
// Created by Bojie Shen on 22/10/21.
//

#ifndef KATCH_REV_CPD_H
#define KATCH_REV_CPD_H
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
#include <sys/mman.h>
#include <fcntl.h>

using namespace std;
//! Compressed Path database. Allows to quickly query the first out arc id of
//! any shortest source-target-path. There may be at most 15 outgoing arcs for
//! any node.
namespace katch {
    class REV_CPD {
    public:
        REV_CPD () : number_of_items_in_row {0}{}

        //! Adds a new node s to the CPD. first_move should be an array that
        //! maps every target node onto a 15-bit bitfield that has a bit set
        //! for every valid first move. get_first_move is free to return any of
        //! them.



        void append_row(const vector<unsigned char>& _first_move){
            std::copy(_first_move.begin(), _first_move.end(), back_inserter(entry));
            number_of_items_in_row = _first_move.size();
        }


        void append_rows(const REV_CPD &other) {
            std::copy(other.entry.begin(), other.entry.end(), back_inserter(entry));
            number_of_items_in_row  = other.number_of_items_in_row;
        }

        unsigned char get_first_move(unsigned source_node) const {
//            return entry[row_index + source_node];
//            return vm_entry[row_index + source_node];
//            return *(current_row + source_node);
            return current_row_entry[source_node];
        }

        void fetch_target_row(unsigned target_row_index){
            row_index = number_of_items_in_row * target_row_index;
            current_row = vm_entry+row_index;
        }

        void load_to_cache(unsigned target_row_index){
            row_index = number_of_items_in_row * target_row_index;
            current_row = vm_entry+row_index;
            for(unsigned i = 0; i < number_of_items_in_row; i++){
//                const auto& a  = *(current_row + i);
                current_row_entry[i] = *(current_row + i);
            }
        }


        unsigned node_count() const {
            return 10000;
        }

        unsigned entry_count() const {
            return entry.size();
        }

        friend bool operator==(const REV_CPD &l, const REV_CPD  &r) {
            return l.entry == r.entry;
        }

        friend bool operator!=(const REV_CPD  &l, const REV_CPD  &r) {
            return !(l == r);
        }

        void save(std::FILE *f) const {
            save_vector_and_row_size(f,number_of_items_in_row, entry);
        }

        void load(std::FILE *f) {
            entry =load_vector_and_row_size<unsigned char>(f,number_of_items_in_row);
        }
        // load to vm, preventing oom
        void load_in_vm(const string& file_name){
            int fd  = open(file_name.c_str(), O_RDONLY, S_IRUSR|S_IWUSR);
            read(fd, &number_of_items_in_row, sizeof(number_of_items_in_row));
            size_t s;
            read(fd, &s, sizeof(s));

            vm_entry = (unsigned char*) mmap(nullptr, s + 16,  PROT_READ, MAP_PRIVATE, fd, 0);
            if ( vm_entry == MAP_FAILED) {
                perror("mmap");
                return ;
            }
            //row_size: 8 , vector_size 8;
            vm_entry = vm_entry  + 16;

            current_row_entry.resize(number_of_items_in_row);
        }

        unsigned long long get_entry_size() const {
            return entry.size();
        }



        const vector<unsigned char> &get_entry() const {
            return entry;
        }


    protected:
        std::vector<unsigned char> entry;
        std::vector<unsigned char> current_row_entry;
        unsigned long number_of_items_in_row;
        unsigned long row_index;
        unsigned char* vm_entry;
        unsigned char * current_row;
    };

}

#endif //KATCH_REV_CPD_H