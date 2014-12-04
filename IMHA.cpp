#include <array>
#include <cmath>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <cassert>
#include <vector>
#include <map>
#include <algorithm>
//#include <unordered_map>
#include <cstdio>
//#include <chrono>
using namespace std;

typedef double Cost;
typedef chrono::high_resolution_clock Clock;
// for a state s, s.first is the board and s.second is its inverse
// permutation: s.first[s.second[i]] = s.second[s.first[i]] = i
typedef pair<vector<int>,vector<int> > State;

// compile-time constants
constexpr int MASK_CLOSED = 1;
constexpr int GRID_ROWS = 4;
constexpr int GRID_COLS = 4;
constexpr int NUM_THREADS = 4;
constexpr int NUM_MOVES = 4;
constexpr Cost INFINITE = 1e99;
// right, up, left, down
constexpr array<int,NUM_MOVES> moves = {1,-GRID_COLS,-1,GRID_COLS};

// condition variable for ending the search
condition_variable main_cv;

// print a nicely formatted board
void printState(const State& s) {
  for (int i = 0; i < GRID_ROWS; ++i) {
    for (int j = 0; j < GRID_COLS; ++j)
      printf("%3d", s.first[i*GRID_COLS + j]);
    printf("\n");
  }
  printf("\n");
}

// data associated with one state
class StateData {
  public:
    // back-pointer to previous state along the discovered path
    State bp;
    // stores Boolean flags such as whether the state is CLOSED
    int mask;
    // g is the cost of a discovered path, h estimates the remaining cost to goal
    Cost g, h;
    // points to the state's location in the priority queue
    multimap<Cost,State>::const_iterator iter;
    // initialization
    StateData() : mask(0), g(INFINITE), h(-1) {}
};

// get the states directly reachable from s by a single tile move
vector<State> getSuccessors(const State& s)
{
  vector<State> successors;
  // position of the gap or "0-tile"
  int gapPos = s.second[0];
  for (int i = 0; i < NUM_MOVES; ++i)
  {
    // neighbor of the gap
    int pos = gapPos + moves[i];
    int val = s.first[pos];
    if (0 <= pos && pos < GRID_ROWS*GRID_COLS)
    {
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

class Problem {
  public:
    State start;
    State goal;
    double w1, w2;
    vector<multimap<Cost,State> > open;
    vector<map<State,StateData> > data;
    Cost OPEN0_MIN;
    int num_discovered;
    int num_expanded;
    int search_done;
    
    inline Cost manhattan(const State& s1, const State& s2) {
      Cost h = 0;
      for (int val = 1; val < GRID_ROWS*GRID_COLS; ++val) {
        // compute distance between val's position in s1 and s2
        int r1 = s1.second[val] / GRID_COLS;
        int c1 = s1.second[val] % GRID_COLS;
        int r2 = s2.second[val] / GRID_COLS;
        int c2 = s2.second[val] % GRID_COLS;
        h += abs(r1-r2) + abs(c1-c2);
      }
      return h;
    }

    // this function determines the heuristics to be used
    inline Cost pairwiseH(int id, const State& s1, const State& s2) {
      return manhattan(s1, s2);
    }

    inline Cost goalH(int id, const State& s) {
      return pairwiseH(id, s, goal);
    }

    // key for the priority queue
    inline Cost f(int id, StateData& s_data) { return s_data.g + w1*s_data.h; }

    void insert(int id, State& s, StateData& s_data) {
      // if state is already in open[id], remove it before re-inserting
      if (s_data.iter != open[id].cend())
        open[id].erase(s_data.iter);
      // insert state with new key
      s_data.iter = open[id].insert(pair<Cost,State>(f(id, s_data), s));
    }

    // assumes start, goal, w1 and w2 are already set
    void init() {
      open.clear(); open.resize(NUM_THREADS);
      data.clear(); data.resize(NUM_THREADS);
      for (int id = 0; id < NUM_THREADS; ++id) {
        open[id].clear();
        // create entry for start state in data[id]
        StateData& start_data = data[id][start];
        start_data.g = 0;
        start_data.h = goalH(id, start);
        start_data.iter = open[id].cend();
        start_data.bp = start;
        insert(id, start, start_data);
      }
      OPEN0_MIN = f(0, data[0][start]);
      // the start state is discovered but not yet expanded
      num_discovered = 1;
      num_expanded = 0;
      search_done = -1;
    }

    void expand(int id, const State& s, StateData& s_data) {
      assert(!(s_data.mask & MASK_CLOSED));
      s_data.mask |= MASK_CLOSED;

      for (State& t : getSuccessors(s)) {
        StateData& t_data = data[id][t];
        if (t_data.h == -1) {
          t_data.h = goalH(id, t);
          t_data.iter = open[id].cend();
          num_discovered++;
        }

        // critical section for updating g-value and inserting
        if(t_data.g > s_data.g + 1){
          t_data.bp = s;
          t_data.g = s_data.g + 1;
          if (!(t_data.mask & MASK_CLOSED))
            insert(id, t, t_data);
        }
      }
    }

    void astar_thread(int id) {
      // repeat until some thread declares the search to be finished
      while (search_done == -1) {
        // get a State to expand
        auto it = open[id].cbegin();
        if (it == open[id].cend()) {
          // search failed... or did it? maybe OPEN0_MIN will
          // rise enough for existing solution to become valid
          break;
        }
        State s = it->second;
        StateData& s_data = data[id][s];
        open[id].erase(it);
        num_expanded++;
        expand(id, s, s_data);

        // the anchor search sends OPEN0_MIN to all the others
        if (id == 0) {
          if (open[0].empty())
            OPEN0_MIN = INFINITE;
          else
            OPEN0_MIN = min(OPEN0_MIN, open[0].cbegin()->first);
        }
        // search termination condition
        if (data[id][goal].g <= w2 * OPEN0_MIN) {
          if (data[id][goal].g < INFINITE)
            search_done = id;
          main_cv.notify_one();
          break;
        }
      }
    }

    // assumes init() was called
    void solveMaze() {
      mutex cv_mutex;
      unique_lock<mutex> cv_lock(cv_mutex);
      vector<thread> threads(NUM_THREADS);
      // create multiple search threads and run them in parallel
      for(int id = 0; id < NUM_THREADS; ++id) {
        threads[id] = thread([=]{ astar_thread(id); });
      }
      //printf("main thread go to sleep\n");
      main_cv.wait(cv_lock, [=]{ return search_done; });
      //printf("main thread awake\n");
      for(thread& th : threads) {
        th.join();
      }
    
      //if(goal->g >= INFINITE || !(goal->mask & MASK_CLOSED)) {
      //  printf("Queue is empty....failed to find goal!\n");
      //}
    }

    // get the path: a sequence of states from start to goal
    vector<State> getSolution() {
      vector<State> sol;
      State s = goal;
      sol.push_back(s);
      while (s != start) {
        s = data[search_done][s].bp;
        sol.push_back(s);
      }
      reverse(sol.begin(), sol.end());
      return sol;
    }
};

int main(int argc, char** argv) {
  FILE* fout = fopen("stats.csv","w");
  // we can set it up to loop over multiple problem instances
  for(int i=0; i<1/*argc*/; i++) {
    // prepare a problem instance
    Problem problem;
    problem.w1 = 1.4;
    problem.w2 = 1.4;
    problem.goal.first = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0};
    problem.goal.second = {15,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
    // make a random start state
    problem.start = problem.goal;
    random_device rd;
    mt19937 gen(rd());
    bool parity = false;
    for (int i = 0; i < GRID_ROWS*GRID_COLS-2; ++i)
    {
      uniform_int_distribution<> dis(i, GRID_ROWS*GRID_COLS-2);
      int swap_cell = dis(gen);
      if (swap_cell != i) {
        parity = !parity;
      }
      swap(problem.start.first[i], problem.start.first[swap_cell]);
      swap(problem.start.second[problem.start.first[i]],
           problem.start.second[problem.start.first[swap_cell]]);
    }
    // fix the parity to ensure a solution exists
    if (parity) {
      swap(problem.start.first[0], problem.start.first[1]);
      swap(problem.start.second[problem.start.first[0]],
           problem.start.second[problem.start.first[1]]);
    }

    // run planner and time its execution
    Clock::time_point t0 = Clock::now();
    problem.init();
    problem.solveMaze();
    Clock::time_point t1 = Clock::now();

    // print solution if it was found
    Cost path_length = INFINITE;
    if (problem.search_done != -1) {
      path_length = problem.data[problem.search_done][problem.goal].g;
      for (State& s : problem.getSolution()) {
        printState(s);
      }
    }

    // report stats
    double dt = chrono::duration<double, chrono::seconds::period>(t1-t0).count();
    printf("map %d: Path Length=%f Visited Nodes=%d Explored Nodes=%d Planning Time=%f\n",i,path_length,problem.num_discovered,problem.num_expanded,dt);
    fprintf(fout,"%d %f %f %f %d %f\n",NUM_THREADS,problem.w1,problem.w2,dt,problem.num_expanded,path_length);
  }
  fclose(fout);
}
