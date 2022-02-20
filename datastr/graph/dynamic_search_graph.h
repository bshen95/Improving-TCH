//
// Created by Bojie Shen on 7/5/21.
//

#ifndef TIME_DEPENDENT_DYNAMIC_SEARCH_GRAPH_H
#define TIME_DEPENDENT_DYNAMIC_SEARCH_GRAPH_H


#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>
#include <stack>

#include "datastr/base/ttf_wrapper.h"
#include "datastr/base/pwl_ttf.h"
#include "datastr/graph/basic.h"
#include <random>
#include "io/vec_io.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
namespace geom = boost::geometry;
typedef geom::model::d2::point_xy<double> point_type;

namespace katch
{

    class DynamicSearchGraph
    {
    private:

        static constexpr double PERIOD = 864000.0;
        static constexpr int BUCKET_SIZE = 8;

        using TTFImpl = PwlTTF<static_cast<int>(PERIOD), BUCKET_SIZE>;

    public:

        using TTFRef = TTFReference<TTFImpl>;
        using TTF = TTFWrapper<TTFImpl>;

    private:

        struct Node
        {
            EdgeIterator _first_edge;
            Node() : _first_edge(INVALID_EDGE_ITERATOR) {}
        };

        static constexpr NodeIterator UNINITIALIZED_NODE_IT = (1 << 28) - 1;

        struct Edge
        {
            NodeIterator _other_node : 28;
            bool _is_constant : 1;
            union
            {
                double _constant_value;
                const TTFImpl* _ttf_impl_ptr;
            };

            Edge(const Edge&) = delete;
            Edge& operator= (const Edge&) = delete;

            double get_free_flow(std::pair<double,double> time_period){
                if(_is_constant)
                    return _constant_value;
                else
                    return _ttf_impl_ptr->get_min(time_period);
            }

            Edge(Edge&& edge)
            {
                _other_node = edge._other_node;
                _is_constant = edge._is_constant;
                if ( _is_constant )
                {
                    _constant_value = edge._constant_value;
                }else{
                    assert( ! _is_constant );
                    _ttf_impl_ptr = edge._ttf_impl_ptr;
                    edge._ttf_impl_ptr = nullptr;
                }

            }

            Edge& operator= (Edge&& edge)
            {
                if ( this != &edge )
                {
                    if ( ! _is_constant )
                        if ( _ttf_impl_ptr )
                            delete _ttf_impl_ptr;

                    _other_node = edge._other_node;
                    _is_constant = edge._is_constant;
                    if ( _is_constant )
                    {
                        _constant_value = edge._constant_value;
                    }else
                    {
                        assert( ! _is_constant );
                        _ttf_impl_ptr = edge._ttf_impl_ptr;
                        edge._ttf_impl_ptr = nullptr;
                    }
                }

                return *this;
            }

            Edge()
                    : _other_node(UNINITIALIZED_NODE_IT),
                      _is_constant(false),
                      _ttf_impl_ptr(nullptr)

            {}

            ~Edge()
            {
                if ( ! _is_constant )
                    if ( _ttf_impl_ptr )
                        delete _ttf_impl_ptr;
                _is_constant = false;
                _ttf_impl_ptr = nullptr;
            }
        };

        std::vector<Node> _nodes;
        std::vector<Edge> _edges;
    public:
        double min_lat, max_lat, min_lon, max_lon;
        std::vector< point_type> _coordinate;
        template <typename EdgeInfo>
        DynamicSearchGraph(std::vector<EdgeInfo>&& edge_list)
                : _nodes(), _edges()
        {

            size_t max_node_it = 0;
            std::for_each( edge_list.begin(), edge_list.end(),
                           [&max_node_it](const EdgeInfo& edge_info) -> void
                           {
                               max_node_it = std::max(max_node_it, edge_info.get_source());
                               max_node_it = std::max(max_node_it, edge_info.get_target());
                           } );

            _nodes.assign(max_node_it + 2, Node());

            std::sort( edge_list.begin(), edge_list.end(),
                       [](const EdgeInfo& lhs, const EdgeInfo& rhs) -> bool
                       {
                           if(lhs.get_source() == rhs.get_source()){
                               return lhs.get_target() < rhs.get_target();
                           }else{
                               return lhs.get_source() < rhs.get_source();
                           }

                       } );

            assert( _nodes[0]._first_edge == INVALID_EDGE_ITERATOR );
            _nodes[0]._first_edge = 0;

            NodeIterator current_source = 0;

            for ( auto it = edge_list.begin() ; it != edge_list.end() ; ++it )
            {
                assert( it->get_source() + 1 < _nodes.size() );

                while ( current_source <= it->get_source() )
                    _nodes[++current_source]._first_edge = _edges.size();

                _edges.push_back(std::move(Edge()));
                _edges.back()._other_node = it->get_target();
                if ( it->get_ttf().is_constant() )
                {
                    _edges.back()._is_constant = true;
                    _edges.back()._constant_value = it->get_ttf().get_constant_value();
                }
                else
                {
                    assert( it->get_ttf().get_ttf_impl_ptr() );

                    _edges.back()._is_constant = false;
                    _edges.back()._ttf_impl_ptr = it->get_ttf().get_ttf_impl_ptr();
                    it->get_ttf().release();
                }

                ++_nodes[current_source]._first_edge;
            }

            while ( current_source <= max_node_it )
                _nodes[++current_source]._first_edge = _edges.size();

            assert( _nodes.back()._first_edge == _edges.size() );

            _edges.emplace_back();

            edge_list.clear();

        }

        size_t get_n_nodes() const noexcept { return _nodes.size() - 1; }
        size_t get_n_edges() const noexcept { return _edges.size() - 1; }

        NodeIterator nodes_begin() const noexcept { return 0; }
        NodeIterator nodes_end() const noexcept { return _nodes.size() - 1; }

        EdgeIterator edges_begin(const NodeIterator u) const noexcept
        {
            assert( u + 1 < _nodes.size() );
            return _nodes[u]._first_edge;
        }

        EdgeIterator edges_end(const NodeIterator u) const noexcept
        {
            assert( u + 1 < _nodes.size() );
            return _nodes[u+1]._first_edge;
        }

        NodeIterator get_other_node(const EdgeIterator& e) const noexcept
        {
            assert( e + 1 < _edges.size() );
            return _edges[e]._other_node;
        }

        double get_free_flow(const EdgeIterator& e, std::pair<double,double> time_period){
            assert( e + 1 < _edges.size() );
            return _edges[e].get_free_flow(time_period);
        }


        bool is_directed_upward(const EdgeIterator e) const noexcept
        {
            assert( e < _edges.size() );
            return false;
        }

        bool is_directed_downward(const EdgeIterator e) const noexcept
        {
            assert( e + 1 < _edges.size() );
            return false;
        }

        double get_upper(const EdgeIterator& e) const noexcept
        {
            assert( e + 1 < _edges.size() );
            if ( _edges[e]._is_constant )
                return _edges[e]._constant_value;

            assert( _edges[e]._ttf_impl_ptr );
            return _edges[e]._ttf_impl_ptr->get_max();
        }

        double get_lower(const EdgeIterator& e) const noexcept
        {
            assert( e + 1 < _edges.size() );
            if ( _edges[e]._is_constant )
                return _edges[e]._constant_value;

            assert( _edges[e]._ttf_impl_ptr );
            return _edges[e]._ttf_impl_ptr->get_min();
        }

        TTFRef get_ttf(const EdgeIterator& e) const noexcept
        {
            assert( e + 1 < _edges.size() );
            if ( _edges[e]._is_constant )
                return TTFRef(_edges[e]._constant_value);
            assert( _edges[e]._ttf_impl_ptr );
            return TTFRef(_edges[e]._ttf_impl_ptr);
        }

        unsigned get_index_size(){
            unsigned size = 0;
            size += _nodes.size()* sizeof(_nodes[0]._first_edge);
            for(auto & _edge : _edges){
                // 4 bytes  for NodeIterator; 1 byte for 4 bool; 4 bytes for uni32_t
                size +=  4;
                if(_edge._is_constant){
                    // each constant value double is 8 bytes
                    size += 8;
                }else{
                    if(_edge._ttf_impl_ptr != nullptr){
                        // each point is (double, double), 8 bytes
                        size += _edge._ttf_impl_ptr->get_n_points()*16;
                        // each bucket is int, 4 bytes
                        size += _edge._ttf_impl_ptr->get_n_bucket()*4;
                        // bucket shift uint32_t
                        size +=  4;
                    }
                };
            }
            return size;
        }
        void print_statistic() {
            unsigned n_constant = 0;
            unsigned n_td = 0;
            unsigned n_points = 0;
            for(const auto& edge: _edges){
                if(edge._is_constant){
                    n_constant++;
                }else{
                    if(edge._ttf_impl_ptr != nullptr){
                        n_points += edge._ttf_impl_ptr->get_n_points();
                        n_td ++;
                    }
                }
            }
            std::cout<<"Total number of nodes: "<<this->get_n_nodes() << std::endl;
            std::cout<<"Total number of edges: "<<this->get_n_edges() << std::endl;
            std::cout<<"Total number of constant edges: "<<n_constant <<" Percentage: "<< (double)n_constant/this->get_n_edges() << std::endl;
            std::cout<<"Total number of td edges: "<<n_td <<" Percentage: "<< (double)n_td/this->get_n_edges()<< std::endl;
            std::cout<<"Total number of turning points : "<<n_points<<std::endl;
            std::cout<<"Average turning points per edges: "<< (double)n_points/this->get_n_edges() << std::endl;
            std::cout<<"Average turning points per td edges: "<< (double)n_points/n_td << std::endl;
        }

        void print_turning_point_info(const std::string& filename) {
            // show turning point every two hours.
            std::vector<unsigned> n_turning_points(24);
            unsigned n_td = 0;
            unsigned n_constant_start = 0;
            for(const auto& edge: _edges){
                if(!edge._is_constant){
                    n_td ++;
                    if(edge._ttf_impl_ptr != nullptr){
                        for(int i = 0; i <24; i++){
                            if( edge._ttf_impl_ptr->get_n_points({i*36000,(i+1)*36000})){
                                n_turning_points[i]++;
                            }
//                            n_turning_points[i] += p;
//                            bool a =0;
                        }
                        if(edge._ttf_impl_ptr->is_constant_start()){
                            n_constant_start++;
                        }
                    }
                }
            }
            std::ofstream myFile(filename);
            myFile<<"time,constant_edges\n";
            for(int i = 0; i < n_turning_points.size(); i++){
                myFile<<std::fixed<<std::setprecision(8)<<i<<","<<(double)n_turning_points[i]/n_td
                      <<"\n";
            }
            myFile.close();
        }

        void print_traffic_index(const std::string& filename) {
            // show traffic index per 20 mins;
            std::vector<double> traffic(24*3);
            unsigned n_td = 0 ;
            int index = 0;
            for(const auto& edge: _edges){
                if(!edge._is_constant){
                    if(edge._ttf_impl_ptr != nullptr){
                        for(int i = 0; i <24*3; i++){
                            traffic[i] += (double)edge._ttf_impl_ptr->eval(i*600*20)/edge._ttf_impl_ptr->get_min();
                        }
                        n_td++;
                    }
                }
                index++;
            }
            std::ofstream myFile(filename);
            myFile<<"time,traffic_index\n";
            for(int i = 0; i < traffic.size(); i++){
                myFile<<std::fixed<<std::setprecision(8)<<i<<","<<(double)traffic[i]/n_td++
                      <<"\n";
            }
            myFile.close();
        }

        void print_edge_traffic_index(const std::string& filename,unsigned id) {
            // show traffic index per 20 mins;
            std::vector<double> traffic(24*3);
            const auto& edge = _edges[(EdgeIterator) id];
            if(!edge._is_constant){
                if(edge._ttf_impl_ptr != nullptr){
                    for(int i = 0; i <24*3; i++){
                        traffic[i] += (double)edge._ttf_impl_ptr->eval(i*600*20);
                    }
                }
            }

            std::ofstream myFile(filename);
            myFile<<"time,traffic_index\n";
            for(int i = 0; i < traffic.size(); i++){
                myFile<<std::fixed<<std::setprecision(8)<<i<<","<<(double)traffic[i]
                      <<"\n";
            }
            myFile.close();
        }


        void print_num_of_directed_edge(){
            unsigned directed_edge = 0;
            unsigned num_edge = 0;
            auto num_nodes = get_n_nodes();

            for(unsigned i  =0 ; i < num_nodes ; i++){
                for(EdgeIterator e = edges_begin(i); e != edges_end(i); e++){
                    num_edge++;
                    NodeIterator other_node = get_other_node(e);
                    bool exist_undirected_edge = false;
                    for(EdgeIterator e1 = edges_begin(other_node); e1 != edges_end(other_node); e1++) {
                        NodeIterator other_node2 = get_other_node(e1);
                        if(other_node2 == i){
                            exist_undirected_edge = true;
                            break;
                        }
                    }
                    if(!exist_undirected_edge){
                        directed_edge ++;
                    }
                }

            }
            std::cout<<"Number of Edges: "<< num_edge<<std::endl;
            std::cout<<"Number of directed edges: "<< directed_edge<<std::endl;





        }


        std::vector<NodeIterator> generate_DFS_ordering2() {

            std::vector<NodeIterator> DFS_ordering(0);
            // Mark all the vertices as not visited

            std::vector<bool> visited(get_n_nodes(), false);
            std::stack<int>stack;
            for(int source_node=0; source_node<get_n_nodes(); ++source_node) {
                if(!visited[source_node]) {
                    stack.push(source_node);
                    while (!stack.empty()) {
                        // Pop a vertex from stack and print it
                        int s = stack.top();
                        stack.pop();

                        // Stack may contain same vertex twice. So
                        // we need to print the popped item only
                        // if it is not visited.
                        if (!visited[s]) {
                            DFS_ordering.push_back(s);
                            visited[s] = true;
                        }

                        // Get all adjacent vertices of the popped vertex s
                        // If a adjacent has not been visited, then push it
                        // to the stack.
                        for(EdgeIterator e = edges_begin(s); e != edges_end(s); e ++ ){
                            if(!visited[get_other_node(e)]){
                                stack.push(get_other_node(e));
                            }
                        }
                    }
                }

            }
            if(is_permutation(DFS_ordering)){
                std::cout<<"ordering correct"<<std::endl;
            }else{
                std::cout<<"ordering incorrect"<<std::endl;
            }
            return DFS_ordering;
        }

        std::vector<NodeIterator> generate_DFS_ordering() {
            std::vector<NodeIterator> DFS_ordering(get_n_nodes(),0);
            // Mark all the vertices as not visited

            std::vector<bool> was_pushed(get_n_nodes(), false);
            std::vector<bool> was_popped(get_n_nodes(), false);
            std::vector<EdgeIterator>next_out = std::vector<EdgeIterator>(get_n_nodes()+1);
            for(int i = 0; i < _nodes.size(); i++){
                next_out[i] = _nodes[i]._first_edge;
            }
            std::stack<NodeIterator>to_visit;
            NodeIterator next_id = 0;
            for(NodeIterator source_node=0; source_node < get_n_nodes(); ++source_node){
                if(!was_pushed[source_node]){

                    to_visit.push(source_node);
                    was_pushed[source_node] = true;

                    while(!to_visit.empty()){
                        NodeIterator x = to_visit.top();
                        to_visit.pop();
                        if(!was_popped[x]){
                            DFS_ordering[x]= next_id++;
                            was_popped[x] = true;
                        }

                        while(next_out[x] != _nodes[x+1]._first_edge){
                            assert(next_out[x]+1 < _edges.size());
                            NodeIterator y = _edges[next_out[x]]._other_node;
                            if(was_pushed[y])
                                ++next_out[x];
                            else{
                                was_pushed[y] = true;
                                to_visit.push(x);
                                to_visit.push(y);
                                ++next_out[x];
                                break;
                            }
                        }
                    }
                }
            }
            DFS_ordering = invert_permutation( DFS_ordering);
            return DFS_ordering;
        }



        void  load_coordinate(const std::string& coordinate_file){
            std::ifstream inNodeFile(coordinate_file);
            if(!inNodeFile)
            {
                std::cout << "Cannot open coordinate file" << coordinate_file<< std::endl;
                return ;
            }
            std::cout << "Reading coordinate file: " << coordinate_file<< std::endl;

            int node_num;
            inNodeFile >> node_num >> min_lat >> max_lat >> min_lon >> max_lon;
            _coordinate.resize(node_num);
            double id;
            double x,y;
            for(int i = 0; i < node_num; i++) {
                inNodeFile >> id >> x >> y;
                _coordinate[i] =  point_type(x,y);
            }
        }

        double get_Euclidean_travel_time(NodeIterator start, NodeIterator target){
            // y is Latitude and x is Longitude;
            //88.5139 is the max speed in the map
            // time unit is 0.1 second, so convert from hours to 0.1 second is 60*60*10;
            double Euclidean_distance_in_km  = util::geo_dist(_coordinate[start].y(), _coordinate[start].x(),
                                                              _coordinate[target].y(), _coordinate[target].x()) /1000;
            return Euclidean_distance_in_km / 88.5139 * 60 * 60 * 10;
//            return 0;
        }


        // may use this to generate query:
        bool is_connected(int s , int t){
            std::vector<bool> visited(this->get_n_nodes(), false);

            // Create a stack for DFS
            std::stack<unsigned> stack;

            // Push the current source node.
            stack.push(s);

            while (!stack.empty())
            {
                // Pop a vertex from stack and print it
                s = stack.top();
                stack.pop();

                // Stack may contain same vertex twice. So
                // we need to print the popped item only
                // if it is not visited.
                if (!visited[s])
                {
                    visited[s] = true;
                    if(s == t){
                        return true;
                    }
                }

                // Get all adjacent vertices of the popped vertex s
                // If a adjacent has not been visited, then push it
                // to the stack.
                for(EdgeIterator e = edges_begin(s);  e != edges_end(s); ++e){
                    if (!visited[get_other_node(e)]){
                        stack.push(get_other_node(e));
                    }

                }
            }
            return false;
        }



        int generate_random_queries( int number_of_query,const std::string& output_file){


            std::cout<<min_lat<<" " << min_lon<<" "<<max_lat<<" " <<max_lon<<std::endl;
            std::cout<<min_lat<<" " << max_lon<<" "<<max_lat<<" " <<min_lon<<std::endl;
            double d1 = util::geo_dist(min_lat,min_lon,max_lat,max_lon);
            double d2 = util::geo_dist(min_lat,max_lon,max_lat,min_lon);



            double max_d = d1>d2 ? d1:d2;
            double bucket_distance = max_d/1024;
            std::vector<std::vector<int>> source= std::vector<std::vector<int>>(10);
            std::vector<std::vector<int>> target= std::vector<std::vector<int>>(10);
            std::cout<<"Max distance: "<< max_d <<std::endl;
            std::default_random_engine e;
            std::uniform_int_distribution<int> u(0, _coordinate.size()-1);
            for(int i = 0 ; i < 10; i ++){
                double current_min = bucket_distance * pow(2,i);
                double current_max = bucket_distance * pow(2,i+1);
                std::cout<<"Distance range from : "<< current_min<<" to " << current_max <<std::endl;
                bool finished = false;
                while(!finished){

                    int  random_source = (random() % static_cast<int>(_coordinate.size()));
                    int  random_target = (random() % static_cast<int>(_coordinate.size()));
//                double distance = geo_dist(lat[random_start],lon[random_start],lat[random_target],lon[random_target]);
                    double distance = util::geo_dist(_coordinate[random_source].x(),_coordinate[random_source].y(),_coordinate[random_target].x(),_coordinate[random_target].y());
                    if(distance <= current_max && distance >= current_min){
                        bool existed = false;
                        for(int j = 0 ; j < source[i].size(); j ++){
                            if((source[i][j] == random_source && target[i][j] == random_target)
                               ||(source[i][j] == random_target && target[i][j] == random_source))
                            {
                                existed  =true;
                                break;
                            }
                        }
                        // make sure it is connected;
                        if(!existed && is_connected(random_source,random_target)){
                            source[i].push_back(random_source);
                            target[i].push_back(random_target);
                        }
                    }

                    if(source[i].size() == number_of_query/10){
                        finished = true;
                    }
                }
                std::cout<<"finished bucket: "<< i <<std::endl;
            }

            std::vector<int>output_source;
            for(const std::vector<int>&s_list : source){
                for(int s : s_list){
                    output_source.push_back(s);
                }
            }

            std::vector<int>output_target;
            for(const std::vector<int>&t_list : target){
                for(int t : t_list){
                    output_target.push_back(t);
                }
            }
            std::vector<double>Euclidean_distance;
            for(int i = 0; i < output_target.size(); i ++){
                unsigned s = output_source[i];
                unsigned t = output_target[i];
                Euclidean_distance.push_back(util::geo_dist(_coordinate[s].x(),_coordinate[s].y(),_coordinate[t].x(),_coordinate[t].y()));
            }
            std::cout<< "finished generated query, size: " << output_source.size()<<std::endl;

            save_vector( output_file + ".source",output_source);
            save_vector( output_file + ".target",output_target);
            save_vector( output_file+ ".Euclidean",Euclidean_distance);
//
//            RT_save_vector( output_file + ".rtsource",output_source);
//            RT_save_vector( output_file + ".rttarget",output_target);
//            RT_save_vector( output_file+ ".rtEuclidean",Euclidean_distance);



        }

        void export_wkt(const std::vector<NodeIterator>& node_id, const std::string& output_file){
            std::ofstream myFile(output_file);
            myFile<<"wkt\n";
            for(auto& id : node_id){
                myFile<< geom::wkt(_coordinate[id]) << "\n";
            }
            myFile.close();
        }

    };



}
#endif //TIME_DEPENDENT_DYNAMIC_SEARCH_GRAPH_H