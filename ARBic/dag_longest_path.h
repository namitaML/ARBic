#include <algorithm>
#include <iostream>
#include <list>
#include <stack>
#include <vector>

// structure directed acyclic graph(uneighted).
class DagGraph {
 public:
  DagGraph(const short *r1, const short *r2, const std::vector<int> &cols,
           bool reverse = false)
      : row1(r1),
        row2(r2),
        columns(cols),
        reverseMode(reverse),
        V(columns.size()),
        longestPathLen(-1),
        adj(V),
        dist(V, 0),
        prevNodeList(V) {
    make_DAG();
  };
  ~DagGraph()
  {
    for(auto v2v:adj)
  	{
		std::vector<int>(0).swap(v2v);
	}
    //adj.clear();
    std::vector<std::vector<int> >(0).swap(adj);
    std::vector<int>(0).swap(columns);
    std::vector<int>(0).swap(dist);
    std::vector<std::vector<int> >(0).swap(prevNodeList);
  };
  void addEdge(int from, int to) { adj[from].push_back(to); }

  int findLongestPathLength();

  void getAllLongestPaths(std::vector<std::vector<int>> &paths);

 private:
  void make_DAG();
  bool nodeStrictLess(int col1, int col2) const;
  bool nodeEqual(int col1, int col2) const;
  void addLongestPath(std::vector<std::vector<int>> &longestPaths,
                      std::vector<int> &curPath, int curVertex);
  void findlongestpath(std::vector<std::vector<int>> &longestPaths,
                      std::vector<int> &curPath, int curVertex);
  void topoSort(std::stack<int> &resStack);
  void topoSortHelper(int v, std::vector<bool> &visited,
                      std::stack<int> &resStack);
  void permutatePath(std::vector<int> &path, int groupStart, int groupEnd,
                     std::vector<std::vector<int>> &resultPaths);
  void permutateAllPaths(std::vector<std::vector<int>> &paths);

  const short *row1;
  const short *row2;
  std::vector<int> columns;
  bool reverseMode;

  int V;  // total number of vertexes
  int longestPathLen;

  std::vector<std::vector<int> > adj;
  std::vector<int> dist;  // longest path length that ends at this vertex
  std::vector<std::vector<int> > prevNodeList;  // used for backtracking
};
