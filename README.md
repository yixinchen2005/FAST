# Getting Started

## Required to Build HavoqGT

- GCC 8.1 or later.
- CMake 2.6 or later.
- Boost C++ Libraries 1.64 or later (only headers required).
- Metall (https://github.com/LLNL/metall) 0.5 or later. Metall only serves as a header for inclusion hence it is not necessarily compile and install Metall in advance.

## Build
Several cmake variables should be configured in `CMakeLists.txt`:
- `BOOST_ROOT`
- `METALL_ROOT`
- `METALL_INCLUDE`

An example to build HavoqGT with Spack is:
```bash
cd havoqgt
mkdir build_dir
cd build_dir
cmake ../ \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-std=c++17 -lrt -lstdc++fs -lpthread" \
  -DHAVOQGT_BUILD_TEST=TRUE \
  -DMPIEXEC_NUMPROC_FLAG="-n"
make
make test # option
make install # option

Use `CMAKE_CXX_COMPILER=/path/to/g++` and `MPI_CXX_COMPILER=/path/to/mpic++` CMake options to specify a C++ compiler and a MPI compiler, respectively.
To change the install directory, one can use `CMAKE_INSTALL_PREFIX` CMake option.

## Usage
### Loading the Graph
In directory `havoqgt/data` has several data graphs described by three files:
- `edge_list`, each line `a b` contains an undirected edge $a-b$ of the data graph
- `vertex_label`, each line `a l` describes that vertex $a$ has a label $l$
- `vertex_label_path`, each line `/path/to/vertex_label` is the path of the file `vertex_label`. When the data graph is large, the `vertex_label` can be splited into multiple files for parallel loadding purpose.

Queries in `havoqgt/data/queries` are described by a single file whose formats are as follows:
- `t # 0`, `t # -1` is the start and end of the file.
- `v 0 1` implies vertex $0$ with a label $1$.
- `e 0 1` implies an undirected edge connecting vertex $0$ and $1$.
- `c ...` a constraint for pruning
- Please note that the lines describing vertices and edges are ordered.

Then load the data graph into HavoqGT:
- `-o, -b` is the directory where the graph database resides
- `-u 1` indicates the graph being undirected
- `-np 4` implies using 4 processing and this must be kept consistant with the query execution. This means the number of processes carrying the query execution must be 4 too.

```bash
mpiexec -np 4 ./ingest_edge_list -o ${database}/raw -b ${database}/backup -u 1 havoqgt/data/eu-core/edges.list
```

### Query Execution
The query execution can be carried out as follows:
- `-o` is output directory
- `-p` is the query graph 

```bash
mpiexec -np 4 ./run_pattern_matching_beta_1.1 -i ${database}/raw/ -v havoqgt/data/eu-core/vertex_label_paths.list -p havoqgt/data/queries/path_8_0.graph -o ../../output1/
```

### Execution Configuration
File `havoqgt/include/prunejuice/config.hpp` is the configuration of the project.
- `#define ENABLE_TDS` means using the original executor of PruneJuice. If commented, the executor is switched to the FAST's.
- `#define ENUMERATION_TIME 300` means that the execution time limit is 300s.
- `#define PRINT_RESULT`, `#define PRINT_PRUNED_GRAPH` enables printing the matched results and pruned graph into the `output` directory.

### Query Generation
`gen_query.py` is used to extract some subgraphs from the datagraph randomly and output them in a standard query format.
- `nodes`, the number of nodes in each query
- `num`, number of queries generated
- `type`, types of queries, choices can be [path, tree, cycle, other]

```bash
python gen_query.py --edge_list eu-core/edge_list --label_map eu-core/vertex_label --nodes 8 --num 20 --output_dir queries --type other  
```

After generating the query graphs, one can generate the constraints further.
- `input_label_map`, the vertex label file of the original data graph
- `input`, the directory where the generated queries reside
- `output`, the output directory
- `mode`, policy of the constraint generation. options=[heuristic, random, none]
- `para`, the length of each constraints in ratio, e.g. para=3 means the length of the constraint is 1/3 of the size of the query.
```bash
python generate_c.py --input ./queries --input_label_map ./eu-core/vertex_label --output queries_heuristic --mode heuristic --para 3
```

