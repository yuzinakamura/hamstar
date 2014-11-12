#include <array>
#include <cmath>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <cassert>
#include <vector>
#include <map>
//#include <chrono>
using namespace std;

typedef double Cost;
typedef vector<vector<int> > State;
typedef chrono::high_resolution_clock Clock;

constexpr int MASK_CLOSED = 1;
constexpr int GRID_ROWS = 4;
constexpr int GRID_COLS = 4;
constexpr int NUM_THREADS = 4;
constexpr Cost INFINITE = 1000000000;

class StateData {
  public:
    int mask;
    Cost g, h;
    multimap<Cost,State*>::const_iterator iter;
    StateData() {
      mask = 0;
      g = INFINITE;
      h = -1;
    }
};

//ypedef pair<State,StateData> StateEntry;
getSuccessors(const State& s)
{
  vector<State> successors;
  // insert successors
  return successors;
}

class Problem {
  public:
    State start;
    State goal;
    double w1, w2;
    vector<multimap<Cost,State> > open;
    vector<unordered_map<State,StateData> > data;
    int num_discovered;
    int num_expanded;
    bool search_done;

    
    inline Cost pairwiseH(int id, const State& s1, const State& s2) {
      return 0;
    }

    inline Cost goalH(int id, const State& s) {
      return pairwiseH(id, s, goal);
    }

    inline Cost f(int id, const StateData& s_data) { return s_data.g + w1*s_data.h; }

    void insert(int id, State& s, const StateData& s_data) {
      if (s_data.iter != open[id].end())
        open[id].erase(s_data.iter);
      s_data.iter = open[id].insert(pair<Cost,State>(f(id, s_data), s));
    }

    // assumes start, goal, w1 and w2 are already set
    void init() {
      open.clear(); open.resize(NUM_THREADS);
      data.clear(); data.resize(NUM_THREADS);
      for (int id = 0; id < NUM_THREADS; ++id) {
        StateData start_data& = data[id][start];
        start_data.g = 0;
        start_data.h = goalH(id, start);
        start_data.iter = open[id].cend();
        open[id].clear();
        insert(id, start, start_data);
      }
      num_discovered = 1;
      num_expanded = 0;
      search_done = false;
    }

    void expand(int id, const State* const s) {
      for (State t : getSuccessors(s)) {
        StateData& t_data = data[t];
        if (t_data.h == -1) {
          t_data.h = goalH(t);
          t_data.iter = open.cend();
          num_discovered++;
        }

        //critical section for updating g-value and inserting
        if(t_data.g > s_data.g + 1){
          t_data.g = s_data.g + 1;
          if (!(t->mask & MASK_CLOSED))
            insert(id, t, t_data);
        }
      }
    }

    void astar_thread(int id) {
      while(!search_done) {
        //get a State to expand
        auto it = open[id].cbegin();
        State s = it->second;
        StateData& s_data = data[id][s];
        open[id].erase(it);
        
        num_expands++;
        
        assert(!(s_data.mask & MASK_CLOSED));
        s_data.mask |= MASK_CLOSED;

        if (s == goal) {
          search_done = true;
          //send a message to all processes
          break;
        }

        expand(id, s);
      }
    }

    void solveMaze() {
      Problem problem;
      problem.init();

      mutex cv_mutex;
      unique_lock<mutex> cv_lock(cv_mutex);
      vector<thread> threads(num_threads);
      for(int i=0; i<num_threads; i++) {
        threads[i] = thread(astar_thread);
      }
      //printf("main thread go to sleep\n");
      main_cv.wait(cv_lock, []{return search_done;});
      //printf("main thread awake\n");
      for(thread& th : threads) {
        th.join();
      }
    
      if(goal->g >= INFINITE || !(goal->mask & MASK_CLOSED)) {
        printf("Queue is empty....failed to find goal!\n");
      }
    }
};
    

int main(int argc, char** argv) {
  FILE* fout = fopen("stats.csv","w");
  //loop over maps
  for(int i=1; i<argc; i++) {
    for(int j=0; j<1; j++) {
      //load map

      //run planner
      Clock::time_point t0 = Clock::now();
      solveMaze();
      Clock::time_point t1 = Clock::now();

      //report stats
      double dt =  chrono::duration<double, chrono::seconds::period>(t1-t0).count();
      printf("map %d, goal %d: Path Length=%f Visited Nodes=%d Explored Nodes=%d Planning Time=%f\n",i-1,j,path_length,num_discovered,num_expands,dt);
      fprintf(fout,"%d %f %f %f %d %f\n",NUM_THREADS,w1,w2,dt,num_expands,path_length);
    }
  }
  fclose(fout);
}
