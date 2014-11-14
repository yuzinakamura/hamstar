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
typedef pair<vector<int>,vector<int> > State;
typedef chrono::high_resolution_clock Clock;

constexpr int MASK_CLOSED = 1;
constexpr int GRID_ROWS = 4;
constexpr int GRID_COLS = 4;
constexpr int NUM_THREADS = 4;
constexpr Cost INFINITE = 100000000;
condition_variable main_cv;
array<int,4> moves = {1,-GRID_COLS,-1,GRID_COLS};

void printState(const State& s) {
  for (int i = 0; i < GRID_ROWS; ++i)
  {
    for (int j = 0; j < GRID_COLS; ++j)
      printf("%3d", s.first[i*GRID_COLS + j]);
    printf("\n");
  }
  printf("\n");
}

class StateData {
  public:
    State bp;
    int mask;
    Cost g, h;
    multimap<Cost,State>::const_iterator iter;
    StateData() {
      mask = 0;
      g = INFINITE;
      h = -1;
    }
};

//typedef pair<State,StateData> StateEntry;
vector<State> getSuccessors(const State& s)
{
  vector<State> successors;
  int zeroPos = s.second[0];
  for (int i = 0; i < 4; ++i)
  {
    int pos = zeroPos + moves[i];
    int val = s.first[pos];
    if (0 <= pos && pos < GRID_ROWS*GRID_COLS)
    {
      State succ = s;
      succ.first[zeroPos] = val;
      succ.first[pos] = 0;
      succ.second[0] = pos;
      succ.second[val] = zeroPos;
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
    
    inline Cost pairwiseH(int id, const State& s1, const State& s2) {
      Cost h = 0;
      for (int val = 1; val < GRID_ROWS*GRID_COLS; ++val) {
        int r1 = s1.second[val] / GRID_COLS;
        int c1 = s1.second[val] % GRID_COLS;
        int r2 = s2.second[val] / GRID_COLS;
        int c2 = s2.second[val] % GRID_COLS;
        h += abs(r1-r2) + abs(c1-c2);
      }
      return h;
    }

    inline Cost goalH(int id, const State& s) {
      return pairwiseH(id, s, goal);
    }

    inline Cost f(int id, StateData& s_data) { return s_data.g + w1*s_data.h; }

    void insert(int id, State& s, StateData& s_data) {
      if (s_data.iter != open[id].cend())
        open[id].erase(s_data.iter);
      s_data.iter = open[id].insert(pair<Cost,State>(f(id, s_data), s));
    }

    // assumes start, goal, w1 and w2 are already set
    void init() {
      open.clear(); open.resize(NUM_THREADS);
      data.clear(); data.resize(NUM_THREADS);
      for (int id = 0; id < NUM_THREADS; ++id) {
        open[id].clear();
        StateData& start_data = data[id][start];
        start_data.g = 0;
        start_data.h = goalH(id, start);
        start_data.iter = open[id].cend();
        start_data.bp = start;
        insert(id, start, start_data);
      }
      OPEN0_MIN = f(0, data[0][start]);
      num_discovered = 1;
      num_expanded = 0;
      search_done = -1;
    }

    void expand(int id, const State& s) {
      StateData& s_data = data[id][s];
      assert(!(s_data.mask & MASK_CLOSED));
      s_data.mask |= MASK_CLOSED;

      for (State& t : getSuccessors(s)) {
        StateData& t_data = data[id][t];
        if (t_data.h == -1) {
          t_data.h = goalH(id, t);
          t_data.iter = open[id].cend();
          num_discovered++;
        }

        //critical section for updating g-value and inserting
        if(t_data.g > s_data.g + 1){
          t_data.bp = s;
          t_data.g = s_data.g + 1;
          if (!(t_data.mask & MASK_CLOSED))
            insert(id, t, t_data);
        }
      }
    }

    void astar_thread(int id) {
      while (search_done == -1) {
        //get a State to expand
        auto it = open[id].cbegin();
        if (it == open[id].cend()) {
          // search failed... or did it? maybe OPEN0_MIN will
          // rise enough for existing solution to become valid
          break;
        }
        State s = it->second;
        open[id].erase(it);

        //printState(s);
        num_expanded++;
        expand(id, s);

        if (id == 0) {
          if (open[0].empty())
            OPEN0_MIN = INFINITE;
          else
            OPEN0_MIN = min(OPEN0_MIN, open[0].cbegin()->first);
        }
        if (data[id][goal].g <= w2 * OPEN0_MIN) {
          if (data[id][goal].g < INFINITE)
            search_done = id;
          main_cv.notify_one();
          break;
        }
      }
    }

    void solveMaze() {
      //mutex cv_mutex;
      //unique_lock<mutex> cv_lock(cv_mutex);
      vector<thread> threads(NUM_THREADS);
      for(int id = 0; id < NUM_THREADS; ++id) {
        threads[id] = thread([=]{ astar_thread(id); });
      }
      //printf("main thread go to sleep\n");
      //main_cv.wait(cv_lock, [=]{ return search_done; });
      //printf("main thread awake\n");
      for(thread& th : threads) {
        th.join();
      }
    
      //if(goal->g >= INFINITE || !(goal->mask & MASK_CLOSED)) {
      //  printf("Queue is empty....failed to find goal!\n");
      //}
    }

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
  //loop over maps
  for(int i=0; i<1/*argc*/; i++) {
    for(int j=0; j<1; j++) {
      //load map

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

      //run planner
      Clock::time_point t0 = Clock::now();
      problem.init();
      problem.solveMaze();
      Clock::time_point t1 = Clock::now();

      Cost path_length = INFINITE;
      if (problem.search_done != -1) {
        path_length = problem.data[problem.search_done][problem.goal].g;
        for (State& s : problem.getSolution()) {
          printState(s);
        }
      }

      //report stats
      double dt =  chrono::duration<double, chrono::seconds::period>(t1-t0).count();
      printf("map %d, goal %d: Path Length=%f Visited Nodes=%d Explored Nodes=%d Planning Time=%f\n",i,j,path_length,problem.num_discovered,problem.num_expanded,dt);
      fprintf(fout,"%d %f %f %f %d %f\n",NUM_THREADS,problem.w1,problem.w2,dt,problem.num_expanded,path_length);
    }
  }
  fclose(fout);
}
