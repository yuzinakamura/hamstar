#include <array>
#include <cmath>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <cassert>
#include <vector>
#include <queue>
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
constexpr int MASK_CLOSED_ANCHOR = 3;
constexpr int GRID_ROWS = 4;
constexpr int GRID_COLS = 4;
constexpr int NUM_THREADS = 4;
constexpr int NUM_MOVES = 4;
constexpr Cost INFINITE = 1e99;
// right, up, left, down
constexpr array<int,NUM_MOVES> moves = {1,-GRID_COLS,-1,GRID_COLS};

// network communication mutex to prevent simultaneous read and write
mutex netmutex[NUM_THREADS][NUM_THREADS];

// print a nicely formatted board
void printState(const State& s) {
  for (int i = 0; i < GRID_ROWS; ++i) {
    for (int j = 0; j < GRID_COLS; ++j)
      printf("%3d", s.first[i*GRID_COLS + j]);
    printf("\n");
  }
  printf("\n");
}

Cost manhattan(const State& s1, const State& s2) {
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

// a process which performs a search, typically assigned to one thread on one machine
class Searcher {
  public:
    // search parameters: endpoints and weights
    State start, goal;
    double w1, w2;
    // the priority queue or search frontier
    multimap<Cost,State> open;
    // a dictionary of seen states and their corresponding data
    map<State,StateData> data;
    // a value that bounds the optimal solution cost
    Cost OPEN0_MIN;
    // number of states seen and expanded
    int num_discovered;
    int num_expanded;
    // id of the thread which completed the search, or -1 if the search is ongoing
    int search_done;
    // id of the Searcher
    int id;
    // simulate network communication with message queue and corresponding mutexes
    vector<queue<State> > network;

    // this function determines the heuristics to be used
    Cost pairwiseH(const State& s1, const State& s2) {
      return manhattan(s1, s2);
    }

    Cost goalH(const State& s) {
      return pairwiseH(s, goal);
    }

    // key for the priority queue
    Cost f(StateData& s_data) { return s_data.g + w1*s_data.h; }

    void insert(State& s, StateData& s_data) {
      // if state is already in open[id], remove it before re-inserting
      if (s_data.iter != open.cend())
        open.erase(s_data.iter);
      // insert state with new key
      s_data.iter = open.insert(pair<Cost,State>(f(s_data), s));
    }

    // assumes start, goal, w1 and w2 are already set
    void init() {
      open.clear();
      data.clear();
      network.clear(); network.resize(NUM_THREADS);

      // create data entry for start state
      StateData& start_data = data[start];
      start_data.g = 0;
      start_data.h = goalH(start);
      start_data.iter = open.cend();
      start_data.bp = start;
      insert(start, start_data);

      OPEN0_MIN = INFINITE;
      // the start state is discovered but not yet expanded
      num_discovered = 1;
      num_expanded = 0;
      search_done = -1;
    }

    void expand(const State& s, StateData& s_data) {
      assert(!(s_data.mask & MASK_CLOSED));
      s_data.mask |= MASK_CLOSED;

      for (State& t : getSuccessors(s)) {
        StateData& t_data = data[t];
        if (t_data.h == -1) {
          t_data.h = goalH(t);
          t_data.iter = open.cend();
          num_discovered++;
        }

        // critical section for updating g-value and inserting
        if(t_data.g > s_data.g + 1) {
          t_data.bp = s;
          t_data.g = s_data.g + 1;
          if (!(t_data.mask & MASK_CLOSED)) {
            insert(t, t_data);
            unique_lock<mutex> state_lock(netmutex[id][0]);
            network[0].push(t);
          }
        }
      }
    }

    void run() {
      // repeat until some thread declares the search to be finished
      while (search_done == -1) {
        // read messages
        for (int oid = 0; oid < NUM_THREADS; ++oid) {
          unique_lock<mutex> msg_lock(netmutex[id][oid], defer_lock);
          if (msg_lock.try_lock())
            while (!network[oid].empty()) {
              State msg = network[oid].front();
              network[oid].pop();
              // TODO: process the message
            }
        }

        // get a State to expand
        auto it = open.cbegin();
        if (it == open.cend()) {
          // search failed... or did it? maybe OPEN0_MIN will
          // rise enough for existing solution to become valid
          break;
        }
        State s = it->second;
        StateData& s_data = data[s];
        open.erase(it);
        num_expanded++;
        expand(s, s_data);

        // the anchor search sends OPEN0_MIN to all the others
        if (id == 0) {
          if (open.empty())
            OPEN0_MIN = INFINITE;
          else
            OPEN0_MIN = max(OPEN0_MIN, open.cbegin()->first);
        }
        // search termination condition
        if (data[goal].g <= w2 * OPEN0_MIN) {
          if (data[goal].g < INFINITE)
            search_done = id;
          //main_cv.notify_one();
          break;
        }
      }
    }

    // get the path: a sequence of states from start to goal
    vector<State> getSolution() {
      vector<State> sol;
      State s = goal;
      sol.push_back(s);
      while (s != start) {
        s = data[s].bp;
        sol.push_back(s);
      }
      reverse(sol.begin(), sol.end());
      return sol;
    }
};
vector<Searcher> searchers;

void prepareDistributedSearch() {
  // netmutex = vector<vector<mutex> >(NUM_THREADS, vector<mutex>(NUM_THREADS));
  searchers.clear(); searchers.resize(NUM_THREADS);
  // create multiple search threads and run them in parallel
  for(int id = 0; id < NUM_THREADS; ++id) {
  	searchers[id].w1 = 1.4;
  	searchers[id].w2 = 1.4;
  	searchers[id].goal.first = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0};
    searchers[id].goal.second = {15,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
    // make a random start state
    searchers[id].start = searchers[id].goal;
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
      swap(searchers[id].start.first[i], searchers[id].start.first[swap_cell]);
      swap(searchers[id].start.second[searchers[id].start.first[i]],
           searchers[id].start.second[searchers[id].start.first[swap_cell]]);
    }
    // fix the parity to ensure a solution exists
    if (parity) {
      swap(searchers[id].start.first[0], searchers[id].start.first[1]);
      swap(searchers[id].start.second[searchers[id].start.first[0]],
           searchers[id].start.second[searchers[id].start.first[1]]);
    }

  	searchers[id].init();
  }
}

// assumes init() was called
void runDistributedSearch() {
  // condition variable for ending the search
  condition_variable main_cv;
  mutex cv_mutex;
  unique_lock<mutex> cv_lock(cv_mutex);
  vector<thread> threads(NUM_THREADS);
  // create multiple search threads and run them in parallel
  for(int id = 0; id < NUM_THREADS; ++id) {
    threads[id] = thread([=]{ searchers[id].run(); });
  }
  //printf("main thread go to sleep\n");
  main_cv.wait(cv_lock, [=]{ return searchers[0].search_done; });
  //printf("main thread awake\n");
  for(thread& th : threads) {
    th.join();
  }

  //if(goal->g >= INFINITE || !(goal->mask & MASK_CLOSED)) {
  //  printf("Queue is empty....failed to find goal!\n");
  //}
}

int main(int argc, char** argv) {
  FILE* fout = fopen("stats.csv","w");
  // we can set it up to loop over multiple problem instances
  for(int i=0; i<1/*argc*/; i++) {
    // prepare the Searcher processes
    prepareDistributedSearch();

    // run planner and time its execution
    Clock::time_point t0 = Clock::now();
    runDistributedSearch();
    Clock::time_point t1 = Clock::now();

    // print solution if it was found
    Cost path_length = INFINITE;
    if (searchers[0].search_done != -1) {
      path_length = searchers[0].data[searchers[0].goal].g;
      for (State& s : searchers[0].getSolution()) {
        printState(s);
      }
    }

    // report stats
    double dt = chrono::duration<double, chrono::seconds::period>(t1-t0).count();
    printf("map %d: Path Length=%f Visited Nodes=%d Explored Nodes=%d Planning Time=%f\n",i,path_length,searchers[0].num_discovered,searchers[0].num_expanded,dt);
    fprintf(fout,"%d %f %f %f %d %f\n",NUM_THREADS,searchers[0].w1,searchers[0].w2,dt,searchers[0].num_expanded,path_length);
  }
  fclose(fout);
}
