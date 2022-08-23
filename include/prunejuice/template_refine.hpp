// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string> 
#include <tuple>
#include <unordered_map>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "config.hpp"
#include "constraint_setup.hpp"

namespace prunejuice {

template<typename Vertex, typename Edge, typename VertexData, 
  typename EdgeData = uint64_t, typename Uint = uint64_t>
class pattern_graph_csr {
  public:
    // currently only support one query in a file
    pattern_graph_csr(std::string query_input, bool _directed = true):directed(_directed){
        std::ifstream vertex_input_file(query_input, std::ifstream::in);
        std::string line;
        std::unordered_map<Vertex, int> degree_map;
        bool processed = false;
        while (std::getline(vertex_input_file, line)) {
            Edge edge_id=0;
            if(line[0]=='t'){
                if(line.find("t # -1") != std::string::npos){
                    if (directed == false && !processed){
                        std::sort(edge_list.begin(), edge_list.end());
                        for (auto itr = edge_list.begin(); itr < edge_list.end(); ++itr){
                            auto e = *itr;
                            edges.push_back(std::get<1>(e));
                        }
                        processed = true;
                    }
                    for(auto v : vertices){
                        vertex_degree.push_back(degree_map[v]);
                    }
                    vertex_count = vertices.size();
                    edge_count = edges.size();
                    // calculate the prefix_sum
                    vertices.clear();
                    vertices.push_back(0);
                    for(Vertex i=0;i<vertex_count;i++){
                        vertices.push_back(*(vertices.rbegin())+degree_map[i]);
                    }
#ifdef ENABLE_TDS
                    std::string tds;
                    tds += generate_constraints(query_input);
                    non_local_constraints.push_back(tds);
#endif
                   // calculate the neighbor label count
                    generate_vertex_neighbor_data_count_map();
                }else{
                    continue;
                }
            }else if(line[0]=='v'){
                std::string line_p = line.substr(2, line.size()-2);
                std::istringstream iss(line_p);
                Vertex v_source;
                VertexData v_data;
                iss>>v_source>>v_data;
                degree_map.insert({v_source, 0});
                vertices.push_back(v_source);
                vertex_data.push_back(v_data);
            }else if(line[0]=='e'){
                std::string line_p = line.substr(2, line.size()-2);
                std::istringstream iss(line_p);
                Vertex s(0), t(0);
                EdgeData w(0);
                iss >> s >> t >> w;
                if (directed == true){
                    edges.push_back(t);
                    edge_list.push_back(std::forward_as_tuple(s, t));
                    // to do: add degree
                }else{
                    edge_list.push_back(std::forward_as_tuple(s, t));
                    edge_list.push_back(std::forward_as_tuple(t, s));  
                    ++ degree_map[s];
                    ++ degree_map[t];
                }
                edge_ID.push_back(edge_id);
                edge_data.push_back(w);
                ++ edge_id;
            }else if(line.find("c diameter") != std::string::npos){
                std::string line_p = line.substr(2, line.size()-2);
                std::istringstream iss(line_p);
                std::string e;
                std::vector<std::string> ele;
                while(std::getline(iss, e, ':')){
                    ele.push_back(e);
                }
                diameter = std::stoull(ele[1]);
            }else if(line[0] == 'c'){
                std::string line_p = line.substr(2, line.size()-2);
                non_local_constraints.push_back(line_p);
            }
        }
        vertex_input_file.close();
        // generate adjecency list
        for(auto p : edge_list){
          auto src = std::get<0>(p);
          auto dst = std::get<1>(p);
          auto it = neighbor.find(src);
          if(it == neighbor.end()){
            neighbor.insert({src, {dst}});
          }else{
            it->second.push_back(dst);
          }
        }
    }

    const bool directed; 
    Vertex vertex_count;
    Edge edge_count;
    Edge diameter;  
    std::vector<Edge> vertices;
    std::vector<Edge> vertex_degree;
    std::vector<VertexData> vertex_data;
    std::vector<Vertex> edges;
    std::vector<Edge> edge_ID; 
    std::vector<EdgeData> edge_data; 
    std::vector<std::tuple<Vertex, Vertex>> edge_list;
    std::vector<std::string> non_local_constraints;
    std::vector<std::unordered_map<VertexData, Uint>> vertex_neighbor_data_count_map; 
    std::unordered_map<Vertex, std::vector<Vertex>> neighbor; // only record out-degree

  private:
    Vertex read_vertex_list(std::string vertex_input_filename) { 
      std::ifstream vertex_input_file(vertex_input_filename, std::ifstream::in);
      std::string line;
      while (std::getline(vertex_input_file, line)) {
        std::istringstream iss(line);
        Vertex v_source(0);
        Edge v_degree(0), v_offset(0);
        iss >> v_source >> v_degree >> v_offset;
        vertices.push_back(v_offset);
        vertex_degree.push_back(v_degree);
      }
      vertex_input_file.close();
      vertex_degree.erase(vertex_degree.end() - 1);
      return vertices.size() <= 0 ? 0 : vertices.size() - 1; // vertex count
    }

    void read_vertex_data_list(std::string vertex_data_input_filename) {
      std::ifstream vertex_data_input_file(vertex_data_input_filename, 
        std::ifstream::in);
      std::string line;
      while (std::getline(vertex_data_input_file, line)) {
        std::istringstream iss(line);
        Vertex v_source(0);
        VertexData v_data(0);
        iss >> v_source >> v_data;
        vertex_data.push_back(v_data);
      }
      vertex_data_input_file.close();
    }

    Edge read_edge_list(std::string edge_input_filename) {
      std::ifstream edge_input_file(edge_input_filename, std::ifstream::in);
      std::string line;
      while(std::getline(edge_input_file, line)) {
        std::istringstream iss(line);
        Vertex s(0), t(0);
        iss >> s >> t;
        if (directed == true){
          edges.push_back(t);
          edge_list.push_back(std::forward_as_tuple(s, t));
        }
        else{
          edge_list.push_back(std::forward_as_tuple(s, t));
          edge_list.push_back(std::forward_as_tuple(t, s));
        }
      }
      if (directed == false){
        std::sort(edge_list.begin(), edge_list.end());
        for (auto itr = edge_list.begin(); itr < edge_list.end(); ++itr){
          auto e = *itr;
          edges.push_back(std::get<1>(e));
          }
      }
      
      edge_input_file.close();
      return edges.size(); // edge count
    }

    void read_edge_data_list(std::string edge_data_input_filename) {
      std::ifstream edge_data_input_file(edge_data_input_filename, std::ifstream::in);
      std::string line;
      while(std::getline(edge_data_input_file, line)) {
        std::istringstream iss(line);
        Vertex s(0), t(0);
        Edge e(0);
        EdgeData w(0); 
        iss >> s >> t >> e >> w;
        edge_ID.push_back(e); // TODO: edge IDs should be in a different file
        edge_data.push_back(w); 
      }
      edge_data_input_file.close();   
    }

    Vertex generate_vertex_list() { // assuming graph is undirected
      Vertex vertex_count = 0;
      Vertex max_vertex = std::get<0>(edge_list[edge_list.size() - 1]);
      Vertex l = 0; // edge list index
      Vertex degree = 0;
      Vertex current_vertex = vertex_count;
      Vertex source;
      //Vertex target;       
    
      do {
        auto edge = edge_list[l];
        source = std::get<0>(edge);
        //target = std::get<1>(edge);
        if (source == current_vertex) {
          degree++;
          l++;
        } else {
          //std::cout << current_vertex << std::endl;
          vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
                                                     vertices[vertices.size() - 1]));
          vertex_degree.push_back(degree);

          //VertexData v_data = get_random_uint(rnd_a, rnd_b, rnd_eng);
          //vertex_data.push_back(v_data);

          // write vertex info to file
          //vertex_file << current_vertex << " " << degree << " " << vertices[vertices.size() - 1] << " " << "\n";
   
          //vertex_data_file << current_vertex << " " << v_data << "\n";
      
          // update vertices array
          degree = 0;
          vertex_count++;
          current_vertex = vertex_count;
        }
      } while(current_vertex <= max_vertex); 

      // add the last dummy vertex to the vertex list
      vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
                                                 vertices[vertices.size() - 1]));
      //vertex_file << current_vertex << " " << degree << " " << vertices[vertices.size() - 1] << " "
      //<< "\n";
      //vertex_file.close();   
      //vertex_data_file.close();

      return vertex_count;  
    }

    Vertex generate_vertex_list_2(size_t mpi_rank) { // assuming graph is undirected
      Vertex vertex_count = 0;
      Vertex max_vertex = std::get<0>(edge_list[edge_list.size() - 1]);
      Vertex l = 0; // edge list index
      Vertex degree = 0;
      Vertex current_vertex = vertex_count;
      Vertex source;
      //Vertex target;       
    
      /*do {
        auto edge = edge_list[l];
        source = std::get<0>(edge);
        //target = std::get<1>(edge);
        if (source == current_vertex) {
          degree++;
          l++;
        } else {
          //std::cout << current_vertex << std::endl;
//          vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
//                                                     vertices[vertices.size() - 1]));
//          vertex_degree.push_back(degree);

          //VertexData v_data = get_random_uint(rnd_a, rnd_b, rnd_eng);
          //vertex_data.push_back(v_data);

          // write vertex info to file
          //vertex_file << current_vertex << " " << degree << " " << vertices[vertices.size() - 1] << " " << "\n";
          if (mpi_rank == 0) {
            std::cout << current_vertex << " " << degree << " " << " " << "\n"; 
          }

          //vertex_data_file << current_vertex << " " << v_data << "\n";
      
          // update vertices array
          degree = 0;
          vertex_count++;
          current_vertex = vertex_count;
        }
      } while(current_vertex <= max_vertex); 

      // add the last dummy vertex to the vertex list
//      vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
//                                                 vertices[vertices.size() - 1]));
      //vertex_file << current_vertex << " " << degree << " " << vertices[vertices.size() - 1] << " "
      //<< "\n";
      //vertex_file.close();   
      //vertex_data_file.close();*/
     
      vertices.resize(10);  
      vertices.push_back(0);
      vertices.push_back(2);
      vertices.push_back(8);
      vertices.push_back(10);
      vertices.push_back(12);
      vertices.push_back(13); 		 
      vertices.push_back(14);
      vertices.push_back(15);
      vertices.push_back(16);
 
      vertex_count = vertices.size() - 1; 

      return vertex_count;  
    }


    void read_stat(std::string stat_input_filename) {
      std::ifstream stat_input_file(stat_input_filename, 
        std::ifstream::in);
      std::string line;
      const char delim = ':';
      while(std::getline(stat_input_file, line)) {
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> tokens;
        while(std::getline(iss, token, delim)) {
          tokens.push_back(token);
        }          
        //std::cout << tokens[0] << " " << tokens[1] << std::endl; 
        assert(tokens.size() > 1);
        boost::trim(tokens[0]);
        boost::trim(tokens[1]); 
        if (boost::iequals(tokens[0], "diameter")) {
          diameter = std::stoull(tokens[1]);         
        }
      }
      stat_input_file.close();
    } 

    void generate_vertex_neighbor_data_count_map() {
      vertex_neighbor_data_count_map.resize(vertex_count);
      for (auto v = 0; v < vertex_count; v++) {
        //std::cout << v << " vertex_data " << vertex_data[v]  
        //  << " vertex_degree " << vertex_degree[v] << std::endl;
        //  std::cout << " neighbours : ";
        for (auto e = vertices[v]; e < vertices[v + 1]; e++) {
          auto v_nbr = edges[e];    
          auto v_nbr_vertex_data = vertex_data[v_nbr];
          //std::cout << v_nbr << ", " << v_nbr_vertex_data << " | ";  

          auto find_nbr_vertex_data = vertex_neighbor_data_count_map[v].find(v_nbr_vertex_data);
          if (find_nbr_vertex_data == vertex_neighbor_data_count_map[v].end()) {
            auto insert_status = vertex_neighbor_data_count_map[v].insert({v_nbr_vertex_data, 1});    
            if(!insert_status.second) {
              std::cerr << "Error: failed to add an element to the map." << std::endl; 
              return;
            }
          } else { 
            find_nbr_vertex_data->second++;
          }            
  
        }
        //std::cout << std::endl;
      }  
    }

    void output_graph_stat() {
      std::cout << "Number of vertices: " << vertex_count << std::endl;
      std::cout << "Number of edges: " << edge_count << std::endl;
      std::cout << "Number of vertex data: " << vertex_data.size() << std::endl;
    } 
};

} // end namespace prunejuice 
