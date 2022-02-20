//
// Created by Bojie Shen on 8/8/20.
//

#ifndef KATCH_VEC_IO_H
#define KATCH_VEC_IO_H
#pragma once
#include <cstdio>
#include <stdexcept>
#include <string>
#include<fstream>
template<class T>
inline void save_vector(std::FILE*file, const std::vector<T>&v){
    int s = v.size();
    if(std::fwrite(&s, sizeof(s), 1, file) != 1)
        throw std::runtime_error("std::fwrite failed");
    if(std::fwrite(&v[0], sizeof(T), v.size(), file) != v.size())
        throw std::runtime_error("std::fwrite failed");
}

template<class T>
inline void save_vector_and_row_size(std::FILE*file, unsigned long row_size, const std::vector<T>&v){
    if(std::fwrite(&row_size, sizeof(row_size), 1, file) != 1)
        throw std::runtime_error("std::fwrite failed");

    size_t s = v.size();
    if(std::fwrite(&s, sizeof(s), 1, file) != 1)
        throw std::runtime_error("std::fwrite failed");
    if(std::fwrite(&v[0], sizeof(T), v.size(), file) != v.size())
        throw std::runtime_error("std::fwrite failed");
}


template<class T>
inline void save_vector(const std::string& file_name, const std::vector<T>&v){
    FILE*file = fopen(file_name.c_str(), "wb");
    int s = v.size();
    if(std::fwrite(&s, sizeof(s), 1, file) != 1)
        throw std::runtime_error("std::fwrite failed");
    if(std::fwrite(&v[0], sizeof(T), v.size(), file) != v.size())
        throw std::runtime_error("std::fwrite failed");
    fclose(file);
}



template<class T>
inline std::vector<T>load_vector(std::FILE*file){
    int s;
    if(std::fread(&s, sizeof(s), 1, file) != 1)
        throw std::runtime_error("std::fread failed");
    std::vector<T>v(s);

    if((int)std::fread(&v[0], sizeof(T), s, file) != s)
        throw std::runtime_error("std::fread failed");

    return v; // NVRO
}

template<class T>
inline std::vector<T>load_vector_and_row_size(std::FILE*file, unsigned long & row_size){
    if(std::fread(& row_size, sizeof( row_size), 1, file) != 1)
        throw std::runtime_error("std::fread failed");
    size_t s;
    if(std::fread(&s, sizeof(s), 1, file) != 1)
        throw std::runtime_error("std::fread failed");
    std::vector<T>v(s);

    if((int)std::fread(&v[0], sizeof(T), s, file) != s)
        throw std::runtime_error("std::fread failed");

    return v; // NVRO
}


template<class T>
inline std::vector<T>load_vector(const std::string& file_name){
    FILE*file = fopen(file_name.c_str(), "r");
    int s;
    if(std::fread(&s, sizeof(s), 1, file) != 1)
        throw std::runtime_error("std::fread failed");
    std::vector<T>v(s);

    if((int)std::fread(&v[0], sizeof(T), s, file) != s)
        throw std::runtime_error("std::fread failed");
    fclose(file);
    return v; // NVRO
}


inline void save_string(std::FILE* f, const std::string& s) {
    size_t len = s.size();
    if (std::fwrite(&len, sizeof(size_t), 1, f) != 1)
        throw std::runtime_error("std::fwrite failed");
    if (std::fwrite(&s[0], sizeof(char), s.size(), f) != len)
        throw std::runtime_error("std::fwrite failed");
}

inline std::string load_string(std::FILE* f) {
    size_t len;
    if (std::fread(&len, sizeof(size_t),1, f) != 1)
        throw std::runtime_error("std::fread failed");

    char* temp = new char[len+1];
    if (std::fread(temp, sizeof(char), len, f) != len)
        throw std::runtime_error("std::fread failed");

    temp[len] = '\0';
    std::string res = std::string(temp);
    return res;
}


template<class T>
void RT_save_vector(const std::string&file_name, const std::vector<T>&vec){
    static_assert(std::is_pod<T>::value, "Cannot find non-trivial serialization code for this type, maybe a header is missing or serialization is simply not available");
    std::ofstream out(file_name, std::ios::binary);
    if(!out)
        throw std::runtime_error("Can not open \""+file_name+"\" for writing.");
    out.write(reinterpret_cast<const char*>(&vec[0]), vec.size()*sizeof(T));
}

template<class T>
std::vector<T>RT_load_vector(const std::string&file_name){
    static_assert(std::is_pod<T>::value, "Cannot find non-trivial serialization code for this type, maybe a header is missing or serialization is simply not available");
    std::ifstream in(file_name, std::ios::binary);
    if(!in)
        throw std::runtime_error("Can not open \""+file_name+"\" for reading.");
    in.seekg(0, std::ios::end);
    unsigned long long file_size = in.tellg();
    if(file_size % sizeof(T) != 0)
        throw std::runtime_error("File \""+file_name+"\" can not be a vector of the requested type because it's size is no multiple of the element type's size.");
    in.seekg(0, std::ios::beg);
    std::vector<T>vec(file_size / sizeof(T));
    in.read(reinterpret_cast<char*>(&vec[0]), file_size);
    return vec; // NVRO
}

inline bool is_permutation(const std::vector<unsigned int>&p){
    std::vector<bool>found(p.size(), false);
    for(unsigned x:p){
        if(x >= p.size())
            return false;
        if(found[x])
            return false;
        found[x] = true;
    }
    return true;
}



inline std::vector<unsigned int> invert_permutation(const std::vector<unsigned int>&p){
    assert(is_permutation(p) && "p must be a permutation");

    std::vector<unsigned int> inv_p(p.size());
    for(int i=0; i<p.size(); ++i)
        inv_p[p[i]] = i;

    return inv_p; // NVRO
}

#endif //KATCH_VEC_IO_H
