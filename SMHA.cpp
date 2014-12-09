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
#include <mpi.h>
#include <string>
using namespace std;

typedef float Cost;
typedef chrono::high_resolution_clock Clock;
// for a state s, s.first is the board and s.second is its inverse
// permutation: s.first[s.second[i]] = s.second[s.first[i]] = i
typedef pair<vector<int>,vector<int> > State;

// compile-time constants
int SEED = 1;
constexpr int MASK_CLOSED = 1;
constexpr int MASK_CLOSED_ANCHOR = 3;
constexpr int GRID_ROWS = 4;
constexpr int GRID_COLS = 4;
constexpr int NUM_THREADS = 1;
constexpr int NUM_MOVES = 4;
constexpr Cost INFINITE = 1e30;

//communication constants
constexpr int TOHEAD_BUFFER_SIZE = 50; //Number of new g values before a message is sent. The actual buffer is double this, because it needs the node too.
constexpr int FROMHEAD_BUFFER_SIZE = 100; //Number of new g values before a message is sent. The actual buffer is double this, because it needs the node too.
constexpr int DATUM_SIZE = GRID_ROWS * GRID_COLS + 1 + 2; //4x4 state, mask, g, and h. Assumes Cost is same size as int
constexpr int HEAD_NODE = 0;

static int comm_size;
static int comm_rank;

// right, up, left, down
constexpr array<int,NUM_MOVES> moves = {1,-GRID_COLS,-1,GRID_COLS};

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
    h += abs(r1-r2) + abs(c1-c2);
  }
  return h;
}

Cost misplacedTiles(const State& s1, const State& s2) {
  Cost h = 0;
  for (int val = 1; val < GRID_ROWS*GRID_COLS; ++val) {
    // check if val is in the same position in s1 and s2
    if (s1.second[val] != s2.second[val])
      ++h;
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
    // multi-heuristic mixing coefficients
    double MD, LC, MT;
    // the priority queue or search frontier
    multimap<Cost,State> open;
    // a dictionary of seen states and their corresponding data
    map<State,StateData> data;
    // a value that bounds the optimal solution cost
    Cost opt_bound;
    // number of states seen and expanded
    int num_discovered;
    int num_expanded;
    // id of the Searcher
    int id;
    // simulate network communication with message queue and corresponding mutexes
    vector<queue<State> > network;

    //Communication stuff
    int * to_buffer;
    int * from_buffer;
    int buffer_index;
    set<State> updated_states;

    // this function determines the heuristics to be used
    Cost pairwiseH(const State& s1, const State& s2) {
      return MD*manhattanDist(s1,s2) + LC*linearConflicts(s1,s2) + MT*misplacedTiles(s1,s2);
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
      if (comm_rank == HEAD_NODE) {
        MD = LC = 1, MT = 0;
      }
      else {
        mt19937 gen2(comm_rank);
        uniform_real_distribution<> dis(1.0, 5.0);
        MD = dis(gen2);
        LC = dis(gen2);
        MT = dis(gen2);
      }

      open.clear();
      data.clear();
      
      to_buffer = new int[TOHEAD_BUFFER_SIZE * DATUM_SIZE];
      from_buffer = new int[FROMHEAD_BUFFER_SIZE * DATUM_SIZE + 1]; //The +1 is a hack. It's a cost
      buffer_index = 0;

      // create data entry for start state
      StateData& start_data = data[start];
      start_data.g = 0;
      start_data.h = goalH(start);
      start_data.iter = open.cend();
      start_data.bp = start;
      insert(start, start_data);

      opt_bound = 0;
      // the start state is discovered but not yet expanded
      num_discovered = 1;
      num_expanded = 0;
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
            updated_states.insert(s);
          }
        }
      }
    }

    void send_message() {
      ;
    }

    void run() {
      // repeat until some thread declares the search to be finished
      int flag = 0;
      MPI_Status status;
      MPI_Request request;
      int error = -1;

      cout << "Starting run" << endl;
      while (true) {
        if (comm_rank != HEAD_NODE) {
          if (updated_states.size() >= TOHEAD_BUFFER_SIZE) {
            cout << "Non-Anchor sending message" << endl;
            assert(updated_states.size() == TOHEAD_BUFFER_SIZE);
            //Non-anchor sends messages
            int index = 0;
            auto it = updated_states.cbegin();
            for (; index < TOHEAD_BUFFER_SIZE*DATUM_SIZE && it != updated_states.cend(); ++it) {
              const State& s = *it;
              StateData& s_data = data[s];
              for (int pos = 0; pos < GRID_ROWS * GRID_COLS; pos++) {
                to_buffer[index++] = s.first[pos];
              }
              to_buffer[index++] = s_data.mask;
              memcpy(&to_buffer[index++], &s_data.g, sizeof(Cost));
              memcpy(&to_buffer[index++], &s_data.h, sizeof(Cost));
            }
            cout << "Sending message of length " << index << endl;
            for (int x = 0; x < TOHEAD_BUFFER_SIZE*DATUM_SIZE; x += 1) {
              cout << "Sending Message int at part " << x << ": " << to_buffer[x] << endl;
            }
            error = MPI_Isend(to_buffer, TOHEAD_BUFFER_SIZE*DATUM_SIZE, MPI_INT, HEAD_NODE, 0, MPI_COMM_WORLD, &request); //Send data to head node
            if (error != MPI_SUCCESS) {
              cout << "ERROR: " << error << endl;
            }
            updated_states.clear();
          }

          // Non-anchor checks for receiving messages
          error = MPI_Iprobe(HEAD_NODE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
          if (error != MPI_SUCCESS) {
            cout << "ERROR: " << error << endl;
          }
          if (flag && 1 == 0) {
            cout << "Non-Anchor receiving message" << endl;
            error = MPI_Irecv(from_buffer, FROMHEAD_BUFFER_SIZE*DATUM_SIZE + 1, MPI_INT, HEAD_NODE, 0, MPI_COMM_WORLD, &request);
            int index = 0;
            for (int datum = 0; datum < FROMHEAD_BUFFER_SIZE; datum++) {
              State s;
              s.second.resize(GRID_ROWS * GRID_COLS);
              for (int pos = 0; pos < GRID_ROWS * GRID_COLS; pos++) {
                s.first.push_back(from_buffer[index]);
                s.second[from_buffer[index++]] = pos;
              }
              assert(sizeof(Cost) == sizeof(int));
              int maskNew;
              Cost gNew, hNew;
              memcpy(&maskNew, &from_buffer[index++], sizeof(int));
              memcpy(&gNew, &from_buffer[index++], sizeof(Cost));
              memcpy(&hNew, &from_buffer[index++], sizeof(Cost));
              StateData& s_data = data[s];
              s_data.mask |= maskNew;
              if (s_data.g > gNew) {
                s_data.g = gNew;
                s_data.h = hNew;
              }
            }
            memcpy(&opt_bound, &from_buffer[index++], sizeof(Cost)); // open.cbegin()->first
          }
        }
        else {
          // the anchor search sends opt_bound to all the others
          if (open.empty())
            opt_bound = INFINITE;
          else
            opt_bound = max(opt_bound, open.cbegin()->first);

          //Anchor checks for received messages
          for (int i = 0; i < comm_size; i++) {
            if (i == HEAD_NODE) {
              continue;
            }
            error = MPI_Iprobe(i, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
            if (error != MPI_SUCCESS) {
              cout << "ERROR: " << error << endl;
            }
            if (flag) {
              error = MPI_Irecv(to_buffer, TOHEAD_BUFFER_SIZE*DATUM_SIZE, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
              cout << "Anchor receiving message" << endl;
              for (int x = 0; x < TOHEAD_BUFFER_SIZE*DATUM_SIZE; x += 1) {
                cout << "Message int at part " << x << ": " << to_buffer[x] << endl;
              }
              int index = 0;
              for (int datum = 0; datum < TOHEAD_BUFFER_SIZE; datum++) {
                State s;
                s.second.resize(GRID_ROWS * GRID_COLS);
                for (int pos = 0; pos < GRID_ROWS * GRID_COLS; pos++) {
                  s.first.push_back(to_buffer[index]);
                  s.second[to_buffer[index++]] = pos;
                }
                assert(sizeof(Cost) == sizeof(int));
                int maskNew;
                Cost gNew, hNew;
                memcpy(&maskNew, &to_buffer[index++], sizeof(int));
                memcpy(&gNew, &to_buffer[index++], sizeof(Cost));
                memcpy(&hNew, &to_buffer[index++], sizeof(Cost));
                StateData& s_data = data[s];
                s_data.mask |= maskNew;

                cout << "datum " << datum << " Old g/h: " << s_data.g << "," << s_data.h << ". New g/h: " << gNew << "," << hNew << endl;
                if (s_data.g > gNew) {
                  s_data.g = gNew;
                  s_data.h = hNew;
                }
                updated_states.insert(s);
              }
            }
          }

          //Anchor sends messages
          if (updated_states.size() > FROMHEAD_BUFFER_SIZE) {
            cout << "Anchor sending message" << endl;
            //cout << "ABOUT TO COMMUNICATE (HEAD)" << endl;
            int index = 0;
            auto it = updated_states.cbegin();
            for (; index < FROMHEAD_BUFFER_SIZE*DATUM_SIZE && it != updated_states.cend(); ++it) {
              const State& s = *it;
              StateData& s_data = data[s];
              for (int pos = 0; pos < GRID_ROWS * GRID_COLS; pos++) {
                from_buffer[index++] = s.first[pos];
              }
              from_buffer[index++] = s_data.mask;
              memcpy(&from_buffer[index++], &s_data.g, sizeof(Cost));
              memcpy(&from_buffer[index++], &s_data.h, sizeof(Cost));
            }
            updated_states.erase(updated_states.cbegin(), it);
            memcpy(&from_buffer[index++], &opt_bound, sizeof(Cost));
            for (int i = 0; i < comm_size; i++) {
              error = MPI_Isend(from_buffer, FROMHEAD_BUFFER_SIZE*DATUM_SIZE + 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
            }
          }
        }

        // get a State to expand
        auto it = open.cbegin();

        // search termination condition
        if (it == open.cend() || data[goal].g <= w2 * opt_bound) {
          // flush and kill
          // but don't kill the master!!!!!
          break;
        }
        State s = it->second;
        StateData& s_data = data[s];
        open.erase(it);
        num_expanded++;
        expand(s, s_data);
      }


      //Clean-up
      if (comm_rank == HEAD_NODE) {
        while (updated_states.size() > FROMHEAD_BUFFER_SIZE) {
          int index = 0;
          auto it = updated_states.cbegin();
          for (; index < FROMHEAD_BUFFER_SIZE*DATUM_SIZE && it != updated_states.cend(); ++it) {
            const State& s = *it;
            StateData& s_data = data[s];
            for (int pos = 0; pos < GRID_ROWS * GRID_COLS; pos++) {
              from_buffer[index++] = s.first[pos];
            }
            from_buffer[index++] = s_data.mask;
            memcpy(&from_buffer[index++], &s_data.g, sizeof(Cost));
            memcpy(&from_buffer[index++], &s_data.h, sizeof(Cost));
          }
          updated_states.erase(updated_states.cbegin(), it);
          memcpy(&from_buffer[index++], &opt_bound, sizeof(Cost));
          for (int i = 0; i < comm_size; i++) {
            MPI_Isend(from_buffer, FROMHEAD_BUFFER_SIZE*DATUM_SIZE + 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
          }
        }

        if (updated_states.size() > 0) {
          int index = 0;
          int datum_count = 0;
          auto it = updated_states.cbegin();
          for (; index < FROMHEAD_BUFFER_SIZE*DATUM_SIZE && it != updated_states.end(); ++it) {
            const State& s = *it;
            StateData& s_data = data[s];
            for (int pos = 0; pos < GRID_ROWS * GRID_COLS; pos++) {
              from_buffer[index++] = s.first[pos];
            }
            from_buffer[index++] = s_data.mask;
            memcpy(&from_buffer[index++], &s_data.g, sizeof(Cost));
            memcpy(&from_buffer[index++], &s_data.h, sizeof(Cost));
            datum_count++;
          }
          while (datum_count < FROMHEAD_BUFFER_SIZE) {
          memcpy(&from_buffer[index], &from_buffer[index - DATUM_SIZE], DATUM_SIZE * sizeof(int));
          index += DATUM_SIZE;
          datum_count++;
          }
          updated_states.clear();
          memcpy(&from_buffer[index++], &opt_bound, sizeof(Cost));
          for (int i = 0; i < comm_size; i++) {
            MPI_Isend(from_buffer, FROMHEAD_BUFFER_SIZE*DATUM_SIZE + 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
          }
        }
      }
      else {
        while (updated_states.size() >= TOHEAD_BUFFER_SIZE) {
          //Non-anchor sends messages
          int index = 0;
          auto it = updated_states.cbegin();
          for (; index < TOHEAD_BUFFER_SIZE*DATUM_SIZE && it != updated_states.cend(); ++it) {
            const State& s = *it;
            StateData& s_data = data[s];
            for (int pos = 0; pos < GRID_ROWS * GRID_COLS; pos++) {
              to_buffer[index++] = s.first[pos];
            }
            to_buffer[index++] = s_data.mask;
            memcpy(&to_buffer[index++], &s_data.g, sizeof(Cost));
            memcpy(&to_buffer[index++], &s_data.h, sizeof(Cost));
          }
          MPI_Isend(to_buffer, TOHEAD_BUFFER_SIZE*DATUM_SIZE, MPI_INT, HEAD_NODE, 0, MPI_COMM_WORLD, &request); //Send data to head node
          updated_states.erase(updated_states.cbegin(), it);
        }

        if (updated_states.size() > 0) {
          int index = 0;
          int datum_count = 0;
          auto it = updated_states.cbegin();
          for (; index < TOHEAD_BUFFER_SIZE*DATUM_SIZE && it != updated_states.cend(); ++it) {
            const State& s = *it;
            StateData& s_data = data[s];
            for (int pos = 0; pos < GRID_ROWS * GRID_COLS; pos++) {
              to_buffer[index++] = s.first[pos];
            }
            to_buffer[index++] = s_data.mask;
            memcpy(&to_buffer[index++], &s_data.g, sizeof(Cost));
            memcpy(&to_buffer[index++], &s_data.h, sizeof(Cost));
            datum_count++;
          }
          while (datum_count < FROMHEAD_BUFFER_SIZE) {
            memcpy(&to_buffer[index], &to_buffer[index - DATUM_SIZE], DATUM_SIZE * sizeof(int));
            index += DATUM_SIZE;
            datum_count++;
          }
          MPI_Isend(to_buffer, TOHEAD_BUFFER_SIZE*DATUM_SIZE, MPI_INT, HEAD_NODE, 0, MPI_COMM_WORLD, &request); //Send data to head node
          updated_states.clear();
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
Searcher searcher;

void prepareDistributedSearch() {
  searcher.w1 = 2.0;
  searcher.w2 = 2.0;
  // make the standard orderly goal state
  for (int i = 1; i < GRID_ROWS*GRID_COLS; i++) {
    searcher.goal.first.push_back(i);
  }
  searcher.goal.first.push_back(0);
  searcher.goal.second.resize(GRID_ROWS * GRID_COLS);
  for (int i = 1; i < GRID_ROWS*GRID_COLS; i++) {
    searcher.goal.second[searcher.goal.first[i]] = i;
  }
  // make a random start state
  searcher.start = searcher.goal;
  mt19937 gen(SEED);
  bool parity = false;
  for (int i = 0; i < GRID_ROWS*GRID_COLS - 2; ++i)
  {
    uniform_int_distribution<> dis(i, GRID_ROWS*GRID_COLS - 2);
    int swap_cell = dis(gen);
    if (swap_cell != i) {
      parity = !parity;
    }
    swap(searcher.start.first[i], searcher.start.first[swap_cell]);
    swap(searcher.start.second[searcher.start.first[i]],
      searcher.start.second[searcher.start.first[swap_cell]]);
  }
  // fix the parity to ensure a solution exists
  if (parity) {
    swap(searcher.start.first[0], searcher.start.first[1]);
    swap(searcher.start.second[searcher.start.first[0]],
      searcher.start.second[searcher.start.first[1]]);
  }

  searcher.init();
}

int main(int argc, char** argv) {
  MPI_Init(NULL, NULL);

  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  printf("Hello world from processor %s, rank %d"
    " out of %d processors\n",
    processor_name, comm_rank, comm_size);

  if (argc > 1) {
    for (int i = 1; i < argc; i++) {
      std::string arg = argv[i];
      if (arg == "-seed") {
        SEED = atoi(argv[++i]);
      }
    }
  }

  FILE* fout = fopen("stats.csv","w");
  // we can set it up to loop over multiple problem instances
  for(int i=0; i<1/*argc*/; i++) {
    // prepare the Searcher processes
    prepareDistributedSearch();

    // run planner and time its execution
    Clock::time_point t0 = Clock::now();
    searcher.run();
    Clock::time_point t1 = Clock::now();

    // print solution if it was found
    Cost path_length = searcher.data[searcher.goal].g;
    if (path_length < INFINITE) {
      for (State& s : searcher.getSolution()) {
        printState(s);
      }
    }

    // report stats
    double dt = chrono::duration<double, chrono::seconds::period>(t1-t0).count();
    printf("map %d: Path Length=%f Visited Nodes=%d Explored Nodes=%d Planning Time=%f\n",i,path_length,searcher.num_discovered,searcher.num_expanded,dt);
    fprintf(fout,"%d %f %f %f %d %f\n",NUM_THREADS,searcher.w1,searcher.w2,dt,searcher.num_expanded,path_length);
  }
  fclose(fout);

  MPI_Finalize();
}