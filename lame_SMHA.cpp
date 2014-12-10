#include <iostream>
#include <array>
#include <cmath>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <cassert>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <algorithm>
//#include <unordered_map>
#include <cstdio>
#include <chrono>
#include <string>

using namespace std;

typedef float Cost;
typedef chrono::high_resolution_clock Clock;
// for a state s, s.first is the board and s.second is its inverse
// permutation: s.first[s.second[i]] = s.second[s.first[i]] = i
typedef pair<vector<int>, vector<int> > State;

// compile-time constants
constexpr int MASK_CLOSED = 1;
constexpr int MASK_CLOSED_ANCHOR = 2;
constexpr int GRID_ROWS = 4;
constexpr int GRID_COLS = 4;
constexpr int NUM_MOVES = 4;
constexpr double TIME_LIMIT = 300;
constexpr Cost INFINITE = 1e30;
constexpr int HEAD_NODE = 0;

#include "lame_SMHA.h"

static int finished = -1;
int SEED = 1;
int NUM_HEURISTICS = 1;
string filename = "test";

// right, up, left, down
constexpr array<int, NUM_MOVES> moves = { 1, -GRID_COLS, -1, GRID_COLS };

Clock::time_point debug_t1, debug_t2;
void debug_time(string str) {
  debug_t1 = debug_t2;
  debug_t2 = Clock::now();
  
  double dt = chrono::duration<double,chrono::seconds::period>(debug_t2-debug_t1).count();
  if (dt > 0.0001)
    cout << str << ' ' << dt << endl;
}

// print a nicely formatted board
void printState(const State& s) {
  for (int i = 0; i < GRID_ROWS; ++i) {
    for (int j = 0; j < GRID_COLS; ++j)
      printf("%3d", s.first[i*GRID_COLS + j]);
    printf("\n");
  }
  printf("\n");
}

Cost manhattanDist(const State& s1, const State& s2) {
  Cost h = 0;
  for (int val = 1; val < GRID_ROWS*GRID_COLS; ++val) {
    // compute distance between val's position in s1 and s2
    int r1 = s1.second[val] / GRID_COLS;
    int c1 = s1.second[val] % GRID_COLS;
    int r2 = s2.second[val] / GRID_COLS;
    int c2 = s2.second[val] % GRID_COLS;
    h += abs(r1 - r2) + abs(c1 - c2);
  }
  return h;
}

Cost misplacedTiles(const State& s1, const State& s2) {
  Cost h = 0;
  for (int val = 1; val < GRID_ROWS*GRID_COLS; ++val) {
    // check if val is in the same position in s1 and s2
    if (s1.second[val] != s2.second[val]) {
      ++h;
    }
  }
  return h;
}

Cost linearConflicts(const State& s1, const State& s2) {
  Cost h = 0;
  // same-row conflicts
  for (int r = 0; r < GRID_ROWS; ++r) {
    vector<int> conflicts(GRID_COLS, 0);
    for (int c1 = 0; c1 < GRID_COLS; ++c1)
      for (int c2 = 0; c2 < c1; ++c2) {
        int val1 = s1.first[r*GRID_COLS + c1];
        int val2 = s2.first[r*GRID_COLS + c2];
        if (s2.second[val1] / GRID_COLS == s2.second[val2] / GRID_COLS &&
          s2.second[val1] < s2.second[val2]) {
          conflicts[c1] |= 1 << c2;
          conflicts[c2] |= 1 << c1;
        }
      }
    while (true) {
      int c1 = 0;
      for (int c = 1; c < GRID_COLS; ++c)
        if (__builtin_popcount(conflicts[c1]) < __builtin_popcount(conflicts[c]))
          c1 = c;
      if (conflicts[c1] <= 0)
        break;
      conflicts[c1] = 0;
      for (int c2 = 0; c2 < GRID_COLS; ++c2)
        conflicts[c2] &= ~(1 << c1);
      h += 2;
    }
  }
  // same-column conflicts
  for (int c = 0; c < GRID_COLS; ++c) {
    vector<int> conflicts(GRID_ROWS, 0);
    for (int r1 = 0; r1 < GRID_ROWS; ++r1)
      for (int r2 = 0; r2 < r1; ++r2) {
        int val1 = s1.first[r1*GRID_COLS + c];
        int val2 = s2.first[r2*GRID_COLS + c];
        if (s2.second[val1] % GRID_COLS == s2.second[val2] % GRID_COLS &&
          s2.second[val1] < s2.second[val2]) {
          conflicts[r1] |= 1 << r2;
          conflicts[r2] |= 1 << r1;
        }
      }
    while (true) {
      int r1 = 0;
      for (int r = 1; r < GRID_ROWS; ++r)
        if (__builtin_popcount(conflicts[r1]) < __builtin_popcount(conflicts[r]))
          r1 = r;
      if (conflicts[r1] <= 0)
        break;
      conflicts[r1] = 0;
      for (int r2 = 0; r2 < GRID_ROWS; ++r2)
        conflicts[r2] &= ~(1 << r1);
      h += 2;
    }
  }
  return h;
}

// get the states directly reachable from s by a single tile move
vector<State> getSuccessors(const State& s)
{
  vector<State> successors;
  // position of the gap or "0-tile"
  int gapPos = s.second[0];
  for (int i = 0; i < NUM_MOVES; ++i) {
    // neighbor of the gap
    int pos = gapPos + moves[i];
    int pr = gapPos / GRID_COLS + moves[i] / GRID_COLS;
    int pc = pos - pr*GRID_COLS;
    if (0 <= pr && pr < GRID_ROWS && 0 <= pc && pc < GRID_COLS)
    {
      int val = s.first[pos];
      // slide the neighbor tile into the gap
      State succ = s;
      succ.first[gapPos] = val;
      succ.first[pos] = 0;
      succ.second[0] = pos;
      succ.second[val] = gapPos;
      successors.push_back(succ);
    }
  }
  return successors;
}

// use s.first to generate its inverse permutation in s.second
void completeHalfState(State& s) {
  s.second.resize(s.first.size());
  for (int i = 0; i < s.first.size(); i++) {
    s.second[s.first[i]] = i;
  }
}

// search parameters: endpoints and weights
State start, goal;
double w1, w2;
// time when the search began
Clock::time_point start_time;
double time_elapsed;
// a dictionary of seen states and their corresponding data
map<State, StateData> data;
// a value that bounds the optimal solution cost
Cost opt_bound;
int total_discovered;
int total_expanded;
vector<Searcher> vec_search;

// determine the heuristics to be used
void computeH(const State& s, StateData& s_data) {
  Cost hMD = manhattanDist(s, goal);
  Cost hLC = linearConflicts(s, goal);
  Cost hMT = misplacedTiles(s, goal);
  for (int i = 0; i < NUM_HEURISTICS; ++i) {
    s_data.h.push_back(vec_search[i].MD*hMD + vec_search[i].LC*hLC + vec_search[i].MT*hMT);
    s_data.iter.push_back(vec_search[i].open.cend());
  }
}

// a process which performs a search, typically assigned to one thread on one machine

Cost Searcher::f(StateData& s_data) { return s_data.g + w1*s_data.h[comm_rank]; }

void Searcher::insert(const State& s, StateData& s_data) {
  // if state is already in open, remove it before re-inserting
  if (s_data.iter[comm_rank] != open.cend())
    open.erase(s_data.iter[comm_rank]);
  // insert state with new key
  s_data.iter[comm_rank] = open.insert(pair<Cost, State>(f(s_data), s));
}

void Searcher::erase(const State& s, StateData& s_data) {
  if (s_data.iter[comm_rank] != open.cend()) {
    open.erase(s_data.iter[comm_rank]);
    s_data.iter[comm_rank] = open.cend();
  }
}

// assumes start, goal, w1 and w2 are already set
void Searcher::init() {
  if (comm_rank == HEAD_NODE) {
    MD = LC = 1;
    MT = 0;
    closing_mask = MASK_CLOSED_ANCHOR;
  }
  else {
    mt19937 gen2(1009 * SEED + comm_rank);
    uniform_real_distribution<> dis(1.0, 5.0);
    MD = dis(gen2);
    LC = dis(gen2);
    MT = dis(gen2);
    closing_mask = MASK_CLOSED;
  }
  num_discovered = 1;
  num_expanded = 0;
}

void Searcher::expand(const State& s, StateData& s_data) {
  assert(!(s_data.mask & closing_mask));
  s_data.mask |= MASK_CLOSED | closing_mask;

  for (int i = 0; i < NUM_HEURISTICS; ++i) {
    vec_search[i].erase(s, s_data);
  }

  for (State& t : getSuccessors(s)) {
    StateData& t_data = data[t];
    if (t_data.g == INFINITE) {
      computeH(t, t_data);
      num_discovered++;
      total_discovered++;
    }

    // critical section for updating g-value and inserting
    if (t_data.g > s_data.g + 1) {
      t_data.g = s_data.g + 1;
      t_data.bp = s.first;
      if (!(t_data.mask & MASK_CLOSED_ANCHOR)) {
        vec_search[0].insert(t, t_data);
        if (!(t_data.mask & MASK_CLOSED)) {
          for (int i = 1; i < NUM_HEURISTICS; ++i) {
            vec_search[i].insert(t, t_data);
          }
        }
      }
    }
  }
}

void Searcher::runSingleIteration() {
  if (comm_rank == HEAD_NODE) {
    if (open.empty()) {
      opt_bound = INFINITE;
    }
    else {
      opt_bound = max(opt_bound, open.cbegin()->first);
    }
  }
  // get a State to expand
  auto it = open.cbegin();

  // search termination condition
  if (data[goal].g <= w2 * opt_bound) {
    finished = comm_rank;
    cout << "Machine " << comm_rank << " found a solution" << endl;
  }
  if (it == open.cend()) {
    return;
  }

  State s = it->second;
  StateData& s_data = data[s];
  num_expanded++;
  total_expanded++;
  expand(s, s_data);
}

// get the path: a sequence of states from start to goal
vector<State> Searcher::getSolution() {
  vector<State> sol;
  State s = goal;
  sol.push_back(s);
  while (s != start) {
    s.first = data[s].bp;
    completeHalfState(s);
    sol.push_back(s);
  }
  reverse(sol.begin(), sol.end());
  return sol;
}


void prepareDistributedSearch() {
  w1 = 1.0;
  w2 = 1.0;
  // make the standard orderly goal state
  for (int i = 1; i < GRID_ROWS*GRID_COLS; i++) {
    goal.first.push_back(i);
  }
  goal.first.push_back(0);
  completeHalfState(goal);
  // make a random start state
  start = goal;
  mt19937 gen(SEED);
  bool parity = false;
  for (int i = 0; i < GRID_ROWS*GRID_COLS - 2; ++i)
  {
    uniform_int_distribution<> dis(i, GRID_ROWS*GRID_COLS - 2);
    int swap_cell = dis(gen);
    if (swap_cell != i) {
      parity = !parity;
    }
    swap(start.first[i], start.first[swap_cell]);
    swap(start.second[start.first[i]],
      start.second[start.first[swap_cell]]);
  }
  // fix the parity to ensure a solution exists
  if (parity) {
    swap(start.first[0], start.first[1]);
    swap(start.second[start.first[0]],
      start.second[start.first[1]]);
  }

  // initialize global data structures
  data.clear();
  opt_bound = 0;
  // the start state is discovered but not yet expanded
  total_discovered = 1;
  total_expanded = 0;

  vec_search.clear();
  vec_search.resize(NUM_HEURISTICS);
  // create data entry for start state
  StateData& start_data = data[start];
  start_data.g = 0;
  computeH(start, start_data);
  start_data.bp = start.first;

  for (int i = 0; i < NUM_HEURISTICS; ++i) {
    vec_search[i].comm_rank = i;
    vec_search[i].open.clear();
    vec_search[i].insert(start, start_data);
    vec_search[i].init();
  }
  start_time = Clock::now();
}

void runDistributedSearch() {
  int benchmark = 0;
  int flag = 0;
  int error = -1;

  cout << "Starting run" << endl;
  int iter = 0;
  Clock::time_point cur_time = Clock::now();
  time_elapsed = 0;
  // repeat until some thread declares the search to be finished
  while (finished == -1 && time_elapsed <= TIME_LIMIT) {
    for (Searcher& searcher : vec_search) {
      searcher.runSingleIteration();
    }
    // timing
    Clock::time_point cur_time = Clock::now();
    time_elapsed = chrono::duration<double,chrono::seconds::period>(cur_time-start_time).count();
    iter++;
    if (total_expanded > benchmark) {
      benchmark += 100000;
      cout << "Total expanded " << total_expanded << endl;
    }
  }
}

int main(int argc, char** argv) {

  for (int i = 1; i < argc; i++) {
    string arg = argv[i];
    if (arg == "-seed" || arg == "-s") {
      SEED = atoi(argv[++i]);
    }
    if (arg == "-filename" || arg == "-f") {
      filename = argv[++i];
    }
    if (arg == "-numh" || arg == "-n") {
      NUM_HEURISTICS = atoi(argv[++i]);
    }
  }

  FILE* fout = fopen((filename + ".csv").c_str(), "a");
  // we can set it up to loop over multiple problem instances
  for (int i = 0; i < 1/*TRIALS*/; i++) {
    // prepare the Searcher processes
    prepareDistributedSearch();

    // run planner and time its execution
    runDistributedSearch();
    
		Cost path_length = data[goal].g;
		// print solution if it was found
		if (path_length <= w2 * opt_bound) {
			/*for (State& s : vec_search[finished].getSolution()) {
				printState(s);
			}*/
      cout << "Found path of length " << path_length << /*" Discovered nodes = " << searcher.num_discovered << ". Expanded nodes: " << searcher.num_expanded << */". Time: " << time_elapsed << endl;
		}
		cout << "Total discovered: " << total_discovered << " total expanded: " << total_expanded << endl;

    // report stats
    //printf("map %d: Path Length=%f Visited Nodes=%d Explored Nodes=%d Planning Time=%f\n", i, path_length, searcher.num_discovered, searcher.num_expanded, time_elapsed);
    fprintf(fout, "%f %f %f %d %f\n", w1, w2, time_elapsed, total_expanded, path_length);
	//cout << "Rank " << comm_rank << " Total discovered: " << total_discovered << " total expanded: " << total_expanded << endl;
  }
  fclose(fout);
}
