//
// Created by Bojie Shen on 18/5/21.
//

#ifndef KATCH_TPGR_FORMAT_NOISE_H
#define KATCH_TPGR_FORMAT_NOISE_H

#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <random>
#include <datastr/base/pwl_ttf.h>

#include "io/edge_info.h"
#include "datastr/base/point.h"
#include "datastr/base/ttf_wrapper.h"
#include "util/misc.h"

namespace katch {
    static constexpr double PERIOD = 864000.0;
    static constexpr int BUCKET_SIZE = 8;

    using TTFImpl = PwlTTF<static_cast<int>(PERIOD), BUCKET_SIZE>;
    using TTFRef = TTFReference<TTFImpl>;
    using TTF = TTFWrapper<TTFImpl>;

    struct Edge{
        unsigned int src;
        unsigned int tgt;
        std::vector<Point> sample_points;
    };

    std::vector<Edge> read_edges(const std::string &input_file_name) {
        std::vector<Edge> result;
        size_t check_n_edges = 0;
        size_t check_n_points = 0;
        unsigned edge_larger_than_10_times = 0;
        unsigned total_edges = 0;
        if (input_file_name == "")
        {
            KATCH_ERROR("Empty input file name given.\n");
            return result;
        }

        KATCH_STATUS("Reading TPGR file '" << input_file_name << "'...");

        std::ifstream input_tpgr_file(input_file_name);

        if ( ! input_tpgr_file.is_open() )
        {
            KATCH_CONTINUE_STATUS(" ABORT\n");
            KATCH_ERROR("Unable to open file '" << input_file_name << "'\n");
            return result;
        }

        unsigned int n_nodes;
        input_tpgr_file >> n_nodes;

        unsigned int n_edges;
        input_tpgr_file >> n_edges;

        unsigned int n_points;
        input_tpgr_file >> n_points;

        unsigned int period;
        input_tpgr_file >> period;

        if ( period != 864000 )
        {
            KATCH_CONTINUE_STATUS(" ABORT\n");
            KATCH_ERROR("Period of time functions is " << period << " (expected 864000).\n");
            return result;
        }

        unsigned int edge_counter = 0;
        while (!input_tpgr_file.eof() && edge_counter < n_edges) {
            ++check_n_edges;
            ++edge_counter;

            unsigned int src;
            input_tpgr_file >> src;

            unsigned int tgt;
            input_tpgr_file >> tgt;

            unsigned int sample_size;
            input_tpgr_file >> sample_size;

            if (sample_size == 0) {
                KATCH_CONTINUE_STATUS(" ABORT\n");
                KATCH_ERROR("TPGR file corrupted: Time function has zero points.\n");
                result.clear();
                return result;
            }

            std::vector <Point> sample_points;
            for (size_t i = 0; i < sample_size; i++) {
                ++check_n_points;

                double x;
                input_tpgr_file >> x;

                double y;
                input_tpgr_file >> y;

                if (x < 0 || x >= period) {
                    KATCH_CONTINUE_STATUS(" ABORT\n");
                    KATCH_ERROR("TPGR file corrupted: x-value not in [0,period).\n");
                    result.clear();
                    return result;
                }

                if (y < 0) {
                    KATCH_CONTINUE_STATUS(" ABORT\n");
                    KATCH_ERROR("TPGR file corrupted: y-value smaller than 0.\n");
                    result.clear();
                    return result;
                }

                if (!sample_points.empty() && x < sample_points.back().x) {
                    KATCH_CONTINUE_STATUS(" ABORT\n");
                    KATCH_ERROR("TPGR file corrupted: x-value smaller than the one before.\n");
                    result.clear();
                    return result;
                }

                sample_points.push_back(Point(x, y));
            }

            assert(sample_points.size() > 0);
//            if(sample_points.size() < 4 && sample_points.size() >1){
//                int point_needed_to_add = 4 - sample_points.size();
//                double x1 = sample_points[0].x;
//                double x2 = sample_points[1].x;
//                double diff = (x2- x1) / (point_needed_to_add + 1);
//                TTF ttf(sample_points.begin(), sample_points.end());
//                std::vector<Point> result_point;
//                result_point.push_back(sample_points[0]);
//                for(int i = 0; i < point_needed_to_add; i ++){
//                    result_point.push_back(Point{x1 + diff*(i+1), ttf.eval(x1 + diff*(i+1))});
//                }
//                for(int i = 1; i < sample_points.size(); i ++){
//                    result_point.push_back(sample_points[i]);
//                }
//                sample_points = result_point;
//            }
            if(sample_points.size() > 1 ){
                TTF ttf(sample_points.begin(), sample_points.end());
                std::vector<Point> result_point;
                for(int hour = 0 ; hour <= 24 ; hour ++){
                    double current_time = hour * 36000.0;
                    double y  = ttf.eval(current_time);
                    result_point.push_back(Point{current_time,y});
                }
                result_point.shrink_to_fit();
                sample_points = result_point;
                if(ttf.get_max()/ttf.get_min() > 10){
                    edge_larger_than_10_times ++;
                }
                total_edges++;
            }

            result.push_back({src,tgt,sample_points });
        }

        if (check_n_edges != n_edges) {

            KATCH_CONTINUE_STATUS(" ABORT\n");
            KATCH_ERROR("TPGR file corrupted: wrong number of edges.\n");
            result.clear();
            return result;
        }

        if (check_n_points != n_points) {
            KATCH_CONTINUE_STATUS(" ABORT\n");
            KATCH_ERROR("TPGR file corrupted: wrong number of bend points\n");
            result.clear();
            return result;
        }

        input_tpgr_file.close();
        KATCH_CONTINUE_STATUS(" OK\n");
        std::cout<<n_points<< std::endl;
        std::cout<<edge_larger_than_10_times<<std::endl;
        std::cout<<total_edges<<std::endl;
        return result;
    }


    std::vector<Edge> add_noise_to_the_graph(const std::vector<Edge>& graph, int mins) {
        double bucket_size = mins*60*10;
        std::vector<Edge> noised_graph;
        for(unsigned i = 0 ; i <graph.size(); i++){
            Edge e ;
            e.src = graph[i].src;
            e.tgt = graph[i].tgt;
            e.sample_points = std::vector<Point>();
            if(graph[i].sample_points.size()==1){
                e.sample_points = graph[i].sample_points;
            }else{
                std::vector<Point> noised_points;
                for(unsigned j = 0 ; j < graph[i].sample_points.size()-1; j++){
                    Point first  = graph[i].sample_points[j];
                    Point second = graph[i].sample_points[j+1];
                    double duration = second.x - first.x;
                    int number_of_points = (duration/bucket_size) -1 ;
                    noised_points.push_back(first);
                    if( second.y > first.y){
                        // ascending order;
                        std::uniform_real_distribution<double> unif(first.y,second.y);
                        std::default_random_engine re;
                        std::vector<double> y;
                        for(unsigned k = 0; k < number_of_points; k++) {
                            double y_v = unif(re);
                            assert(le(y_v,second.y));
                            assert(ge(y_v,first.y));
                            y.push_back(y_v);
                        }
                        sort(y.begin(), y.end());
                        for(auto y_v : y){
                            noised_points.emplace_back(Point {noised_points.back().x+bucket_size, y_v});
                        }
                    }else if (first.y > second.y){
                        std::uniform_real_distribution<double> unif(second.y,first.y);
                        std::default_random_engine re;
                        std::vector<double> y;
                        for(unsigned k = 0; k < number_of_points; k++) {
                            double y_v = unif(re);
                            assert(ge(y_v,second.y));
                            assert(le(y_v,first.y));
                            y.push_back(y_v);
                        }
                        sort(y.begin(), y.end(), std::greater<double>());
                        for(auto y_v : y){
                            noised_points.emplace_back(Point {noised_points.back().x+bucket_size, y_v});
                        }
                    }
                }
                noised_points.push_back(graph[i].sample_points.back());
                assert(noised_points.front() == graph[i].sample_points.front());
                assert(noised_points.back() == graph[i].sample_points.back());
                e.sample_points = noised_points;
            }
            noised_graph.emplace_back(e);
        }
        return noised_graph;
    }


}
#endif //KATCH_TPGR_FORMAT_NOISE_H
