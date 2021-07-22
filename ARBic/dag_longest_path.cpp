#include "dag_longest_path.h"

using namespace std;

void DagGraph::make_DAG() {
  //cout << "SDFSDF make DAG ...V=" << V << (reverseMode ? " reverse" : "") << endl;
  for (int i = 0; i < V; i++) {
    for (int j = 0; j < V; j++) {
      if (i == j) {
        continue;
      }
//下面的我觉得可以改一下
//      if (nodeStrictLess(i, j)) {
	if (nodeStrictLess(columns[i],columns[j])) {
        //cout << "SDFSDF addEdge " << i << "->" << j << endl;
        addEdge(i, j);
      }
    }
  }
}

bool DagGraph::nodeEqual(int col1, int col2) const {
  return row1[col1] == row1[col2] && row2[col1] == row2[col2];
}
/*
bool DagGraph::nodeStrictLess(int col1, int col2) const {
  int n00 = row1[col1];
  int n01 = row1[col2];
  int n10 = row2[col1];
  int n11 = row2[col2];

  if ((n00 == n01) && (n10 == n11)) {
    // for node with same value, order them by column index
	if(col1<col2)//
		{return true;}//
    //return col1 < col2;
  }

  if (reverseMode) {
    return ((n00 >= n01) && (n10 >= n11));
  } else {
    return ((n00 <= n01) && (n10 <= n11));
  }
}

*/
///修改版
bool DagGraph::nodeStrictLess(int col1, int col2) const {
  int n00 = row1[col1];
  int n01 = row1[col2];
  int n10 = row2[col1];
  int n11 = row2[col2];


  if (reverseMode)
	 {
		if ((n00 == n01) && (n10 == n11))
	 	{
			return col1>col2;
		}
		
		else 
		{	
	   		return false;
	  	}
	 } 
	else 
	{
		if ((n00 == n01) && (n10 == n11))
		 {
			return col1<col2;
		
		}
		else
		{
	 	   return ((n00 >= n01) && (n10 > n11)||(n00 > n01) && (n10 >= n11));
		}
	
	}
	
}







void DagGraph::topoSortHelper(int v, vector<bool> &visited, stack<int> &resStack)
 {
  visited[v] = true;

//  for (auto node : adj[v]) {
    vector<int>::iterator node;
    for(node=adj[v].begin();node!=adj[v].end();++node){
    if (!visited[*node]) {
      topoSortHelper(*node, visited, resStack);
    }
  }

  resStack.push(v);
}

// DFS topo sort algorithm
void DagGraph::topoSort(stack<int> &resStack) {
  vector<bool> visited(V, false);
  for (int i = 0; i < V; i++) {
    if (!visited[i]) {
      topoSortHelper(i, visited, resStack);
    }
  }
}

int DagGraph::findLongestPathLength() {
  // get topo sorted order for the vertexes first
  stack<int> resStack;
  topoSort(resStack);

  while (!resStack.empty()) {
    int u = resStack.top();
    resStack.pop();
    for (auto node : adj[u]) {
      if (dist[node] < dist[u] + 1) {
        dist[node] = dist[u] + 1;
        prevNodeList[node].clear();
        prevNodeList[node].push_back(u);
      } else if (dist[node] == dist[u] + 1) {
        prevNodeList[node].push_back(u);
      }
    }
  }

  longestPathLen = 0;
  for (auto d : dist) {
    if (d > longestPathLen) {
      longestPathLen = d;
    }
  }

  longestPathLen++;

  //cout << "SDFSDF longestPathLen=" << longestPathLen << endl;
  return longestPathLen;
}

void DagGraph::findlongestpath(vector<vector<int> >&longestpaths,vector<int>&longestpath, int curVertex)
{
   longestpath.push_back(curVertex);
   if(prevNodeList[curVertex].empty())
	{
	    longestpaths.push_back(longestpath);
	}
    else
	{
  	    findlongestpath(longestpaths,longestpath,prevNodeList[curVertex][0]);
	}
return;
}




// backtrack to find all the longest paths
void DagGraph::addLongestPath(vector<vector<int>> &longestPaths,
                              vector<int> &curPath, int curVertex) {
  curPath.push_back(curVertex);
  if (prevNodeList[curVertex].empty()) {
    longestPaths.push_back(curPath);
    auto &last = longestPaths.back();
    // need to reverse the newly added path
    reverse(last.begin(), last.end());

    // cout << "SDFSDF PATH=";
    // for (auto i : last) {
    //   cout << "(" << i << "," << row1[i] << "," << row2[i] << ")";
    //}
    // cout << endl;

  } else {
    for (int prev : prevNodeList[curVertex]) {
      addLongestPath(longestPaths, curPath, prev);
    }
  }
  curPath.pop_back();
}
/*
void DagGraph::permutatePath(vector<int> &path, int groupStart, int groupEnd,
                             vector<vector<int>> &resultPaths) {
  if (groupStart == path.size() - 1) {  // last permutation position
    resultPaths.push_back(path);
    return;
  }

  int from = groupStart;
  int end = groupEnd;
  if (groupEnd < 0 || groupStart > groupEnd) {
    // need to find the range with same vertex values
    from = -1;
    for (size_t i = groupStart + 1; i < path.size(); i++) {
      if (nodeEqual(path[i], path[i - 1])) {
        if (from < -1) {
          from = i - 1;
        }
        end = i;
      } else {
        if (from >= 0) {
          break; // found one
        }
      }
    }
  }

  if (from < 0) {
    // no need permutation
    resultPaths.push_back(path);
    return;
  }

  for (int i = from; i <= end; i++) {
    swap(path[i], path[from]);
    permutatePath(path, i + 1, end, resultPaths);
    swap(path[i], path[from]);
  }
}

void DagGraph::permutateAllPaths(vector<vector<int>> &paths) {
  // permutate vetexes with the same value only
  vector<vector<int>> newPaths;
  for (auto& path : paths) {
    permutatePath(path, 0, -1, newPaths);
  }
  paths.swap(newPaths);
}
*/
void DagGraph::getAllLongestPaths(vector<vector<int>> &paths) {
  if (longestPathLen < 0) {
    longestPathLen = findLongestPathLength();
  }

  for (int v = 0; v < dist.size(); v++) {
    if (dist[v] == longestPathLen - 1) {
      vector<int> path;
     // addLongestPath(paths, path, v);
     findlongestpath(paths,path,v);
    }
  }

  //permutateAllPaths(paths);//g
}
