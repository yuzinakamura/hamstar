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
typedef pair<vector<int>, vector<int> > State;

// compile-time constants
int SEED = 1;
constexpr int MASK_CLOSED = 1;
constexpr int MASK_CLOSED_ANCHOR = 2;
constexpr int GRID_ROWS = 4;
constexpr int GRID_COLS = 4;
constexpr int NUM_THREADS = 1;
constexpr int NUM_MOVES = 4;
constexpr Cost INFINITE = 1e30;

//communication constants
constexpr int COMM_FREQ = 100;
constexpr int BUFFER_SIZE = COMM_FREQ*NUM_MOVES*2; //Number of new g values before a message is sent. The actual buffer is double this, because it needs the node too.
constexpr int DATUM_SIZE = (GRID_ROWS * GRID_COLS)*2 + 1 + 2; //4x4 state, backtrace, mask, g, and h. Assumes Cost is same size as int
constexpr int HEAD_NODE = 0;

static int comm_size;
static int comm_rank;

// right, up, left, down
constexpr array<int, NUM_MOVES> moves = { 1, -GRID_COLS, -1, GRID_COLS };

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
  // g is the cost of a discovered path, h and hAnch estimate the remaining cost to goal
  Cost g, h, hAnch;
  // points to the state's location in the priority queue
  multimap<Cost, State>::const_iterator iter;
  // initialization
  StateData() : mask(0), g(INFINITE) {}
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

// a process which performs a search, typically assigned to one thread on one machine
class Searcher {
public:
  // search parameters: endpoints and weights
  State start, goal;
  double w1, w2;
  // multi-heuristic mixing coefficients
  double MD, LC, MT;
  // the priority queue or search frontier
  multimap<Cost, State> open;
  // a dictionary of seen states and their corresponding data
  map<State, StateData> data;
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
  int * head_buffer;
  int * child_buffer;
  set<State> updated_states;

  // determine the heuristics to be used
  void computeH(const State& s, StateData& s_data) {
    Cost hMD = manhattanDist(s, goal);
    Cost hLC = linearConflicts(s, goal);
    Cost hMT = misplacedTiles(s, goal);
    s_data.hAnch = hMD + hLC;
    s_data.h = MD*hMD + LC*hLC + MT*hMT;
  }

  // key for the priority queue
  Cost f(StateData& s_data) { return s_data.g + w1*s_data.h; }

  void insert(State& s, StateData& s_data) {
    // if state is already in open[id], remove it before re-inserting
    if (s_data.iter != open.cend())
      open.erase(s_data.iter);
    // insert state with new key
    s_data.iter = open.insert(pair<Cost, State>(f(s_data), s));
  }

  // assumes start, goal, w1 and w2 are already set
  void init() {
    if (comm_rank == HEAD_NODE) {
      MD = LC = 1;
      MT = 0;
    }
    else {
      mt19937 gen2(1009 * SEED + comm_rank);
      uniform_real_distribution<> dis(1.0, 5.0);
      MD = dis(gen2);
      LC = dis(gen2);
      MT = dis(gen2);
    }

    open.clear();
    data.clear();

    head_buffer = new int[(BUFFER_SIZE * DATUM_SIZE + 1) * comm_size];
    child_buffer = new int[(BUFFER_SIZE * DATUM_SIZE + 2) * comm_size]; //The +1 is a hack. It's a cost

    // create data entry for start state
    StateData& start_data = data[start];
    start_data.g = 0;
    computeH(start, start_data);
    start_data.iter = open.cend();
    start_data.bp = start;
    insert(start, start_data);

    opt_bound = 0;
    // the start state is discovered but not yet expanded
    num_discovered = 1;
    num_expanded = 0;
  }

  void expand(const State& s, StateData& s_data) {
    // assert(!(s_data.mask & MASK_CLOSED));
    if (comm_rank == HEAD_NODE) {
      if (s_data.mask & MASK_CLOSED_ANCHOR) {
        return;
      }
      s_data.mask |= MASK_CLOSED_ANCHOR;
      s_data.mask |= MASK_CLOSED;
    }
    else {
      if (s_data.mask & MASK_CLOSED) {
        return;
      }
      s_data.mask |= MASK_CLOSED;
    }

    for (State& t : getSuccessors(s)) {
      StateData& t_data = data[t];
      if (t_data.g == INFINITE) {
        computeH(t, t_data);
        t_data.iter = open.cend();
        num_discovered++;
      }

      // critical section for updating g-value and inserting
      if (t_data.g > s_data.g + 1) {
        t_data.g = s_data.g + 1;
        t_data.bp = s;
        if (!(t_data.mask & MASK_CLOSED)) {
          insert(t, t_data);
          updated_states.insert(s);
        }
      }
    }
  }

  State update_state(int * buffer, int& index) { // when passed a buffer, and a starting point in that buffer, interprets a state, updates its stateData, and returns it.
    int position_sum = 0, indexStart = index; // sanity check
    
    State s;
    for (int pos = 0; pos < GRID_ROWS * GRID_COLS; pos++) {
      position_sum += buffer[index];
      assert(0 <= buffer[index] && buffer[index] < GRID_ROWS * GRID_COLS);
      s.first.push_back(buffer[index++]);
    }
    completeHalfState(s);
    if (!(position_sum == (GRID_ROWS*GRID_COLS)*(GRID_ROWS*GRID_COLS - 1) / 2)) {
      cout << "Parse failure from comm. Buffer: ";
      for (int i = index; i < index + DATUM_SIZE; i++) {
        if (i % DATUM_SIZE == 0) {
          cout << endl << "Index " << i << ": ";
        }
        cout << buffer[i] << " ";
      }
    }
    assert(position_sum == (GRID_ROWS*GRID_COLS)*(GRID_ROWS*GRID_COLS-1)/2);

    State bpNew;
    for (int pos = 0; pos < GRID_ROWS * GRID_COLS; pos++) {
      position_sum += buffer[index];
      assert(0 <= buffer[index] && buffer[index] < GRID_ROWS * GRID_COLS);
      bpNew.first.push_back(buffer[index++]);
    }
    completeHalfState(bpNew);

    assert(sizeof(Cost) == sizeof(int));
    int maskNew;
    Cost gNew, hAnchNew;
    memcpy(&maskNew, &buffer[index++], sizeof(int));
    memcpy(&gNew, &buffer[index++], sizeof(Cost));
    memcpy(&hAnchNew, &buffer[index++], sizeof(Cost));

    StateData& s_data = data[s];
    s_data.mask |= maskNew;
    if (s_data.g == INFINITE) {
      if (comm_rank == HEAD_NODE)
        s_data.h = s_data.hAnch = hAnchNew;
      else
        computeH(s, s_data);
    }
    if (s_data.g > gNew) {
      s_data.g = gNew;
      s_data.bp = bpNew;
    }

    assert(index == indexStart + DATUM_SIZE);
    return s;
  }

  void serialize_state(const State * state, StateData * s_data, int * buffer, int& index) { //When passed in a state, it will serialize the state's data and insert it into the buffer starting at index
    int indexStart = index; // sanity check
    
    for (int pos = 0; pos < GRID_ROWS * GRID_COLS; pos++) {
      buffer[index++] = state->first[pos];
    }
    for (int pos = 0; pos < GRID_ROWS * GRID_COLS; pos++) {
      buffer[index++] = s_data->bp.first[pos];
    }
    buffer[index++] = s_data->mask;
    memcpy(buffer+index, &s_data->g, sizeof(Cost));
    index += sizeof(Cost) / sizeof(int);
    memcpy(buffer+index, &s_data->hAnch, sizeof(Cost));
    index += sizeof(Cost) / sizeof(int);

    assert(index == indexStart + DATUM_SIZE);
  }

  void run() {
    // repeat until some thread declares the search to be finished
    int flag = 0;
    MPI_Status status;
    MPI_Request request;
    int error = -1;

    cout << "Starting run" << endl;
    int iter = 0;
    while (true) {
      if (comm_rank == HEAD_NODE) {
        if (open.empty()) {
          opt_bound = INFINITE;
        }
        else {
          opt_bound = max(opt_bound, open.cbegin()->first);
        }
      }
      if (iter % COMM_FREQ == 0) {
        //cout << "Process " << comm_rank << " preparing to communicate " << endl;
        //Handle communication
        //First, Gather to HEAD_NODE
        if (comm_rank == HEAD_NODE) {
          memset(head_buffer, 0, sizeof(int)*(BUFFER_SIZE*DATUM_SIZE + 1)*comm_size);
          //cout << "Process " << comm_rank << " receiving gather" << endl;
          //MPI_Gather(MPI_IN_PLACE, 0, MPI_INT, head_buffer, BUFFER_SIZE*DATUM_SIZE+1, MPI_INT, HEAD_NODE, MPI_COMM_WORLD);
          MPI_Gather(MPI_IN_PLACE, BUFFER_SIZE*DATUM_SIZE + 1, MPI_INT, head_buffer, BUFFER_SIZE*DATUM_SIZE + 1, MPI_INT, HEAD_NODE, MPI_COMM_WORLD);
          //cout << "Process " << comm_rank << " received gather" << endl;
          int index = 0;
          for (int machine = 0; machine < comm_size; machine++) {
            int num_sent = head_buffer[machine*(BUFFER_SIZE*DATUM_SIZE+1)]; //Number of states this child sent in.
            index = machine*(BUFFER_SIZE*DATUM_SIZE + 1) + 1;
            while (num_sent > 0) {
              num_sent--;
              State s = update_state(head_buffer, index);
              updated_states.insert(s);
            }
          }
        }
        else {
          memset(child_buffer, 0, sizeof(int)*(BUFFER_SIZE*DATUM_SIZE + 2)*comm_size);
          child_buffer[0] = updated_states.size();
          //cout << "States updated: " << child_buffer[0] << endl;
          auto it = updated_states.cbegin();
          int index = 1;
          for (; it != updated_states.cend(); it++) {
            const State& s = *it;
            StateData& s_data = data[s];
            serialize_state(&s, &s_data, child_buffer, index);
          }
          updated_states.clear();
          //cout << "Process " << comm_rank << " sending gather" << endl;
          MPI_Gather(child_buffer, BUFFER_SIZE*DATUM_SIZE + 1, MPI_INT, NULL, BUFFER_SIZE*DATUM_SIZE + 1, MPI_INT, HEAD_NODE, MPI_COMM_WORLD);
          //cout << "Process " << comm_rank << " sent gather" << endl;
        }

        //Now, Bcast
        //cout << "Process " << comm_rank << " Beginning bcast" << endl;
        memset(child_buffer, 0, sizeof(int)*(BUFFER_SIZE*DATUM_SIZE + 2)*comm_size);
        if (comm_rank == HEAD_NODE) {
          memcpy(&child_buffer[0], &opt_bound, sizeof(Cost));
          child_buffer[1] = updated_states.size();
          auto it = updated_states.cbegin();
          int index = 2;
          for (; it != updated_states.cend(); it++) {
            const State& s = *it;
            StateData& s_data = data[s];
            serialize_state(&s, &s_data, child_buffer, index);
          }
          updated_states.clear();
        }
        //cout << "Process " << comm_rank << " Sending bcast" << endl;
        MPI_Bcast(child_buffer, (BUFFER_SIZE*DATUM_SIZE + 2)*comm_size, MPI_INT, HEAD_NODE, MPI_COMM_WORLD);
        //cout << "Process " << comm_rank << " sent bcast " << endl;
        if (comm_rank != HEAD_NODE) {
          int index = 0;
          memcpy(&opt_bound, &child_buffer[index++], sizeof(Cost));
          int num_sent = child_buffer[index++];
          assert(num_sent >= 0 && num_sent < 100000000);
          while (num_sent > 0) {
            num_sent--;
            State s = update_state(child_buffer, index);
            //updated_states.insert(s);
          }
        }
        //cout << "Process " << comm_rank << " finished comm" << endl;
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

      iter++;
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
  completeHalfState(searcher.goal);
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

  FILE* fout = fopen("stats.csv", "w");
  // we can set it up to loop over multiple problem instances
  for (int i = 0; i<1/*argc*/; i++) {
    // prepare the Searcher processes
    prepareDistributedSearch();

    // run planner and time its execution
    Clock::time_point t0 = Clock::now();
    searcher.run();
    Clock::time_point t1 = Clock::now();

    double dt = chrono::duration<double, chrono::seconds::period>(t1 - t0).count();
    Cost path_length = searcher.data[searcher.goal].g;
    cout << "Found path of length " << path_length << ". Discovered nodes = " << searcher.num_discovered << ". Expanded nodes: " << searcher.num_expanded << ". Time: " << dt << endl;

    // print solution if it was found
    if (path_length < INFINITE) {
      for (State& s : searcher.getSolution()) {
        printState(s);
      }
    }

    // report stats
    printf("map %d: Path Length=%f Visited Nodes=%d Explored Nodes=%d Planning Time=%f\n", i, path_length, searcher.num_discovered, searcher.num_expanded, dt);
    fprintf(fout, "%d %f %f %f %d %f\n", NUM_THREADS, searcher.w1, searcher.w2, dt, searcher.num_expanded, path_length);
  }
  fclose(fout);

  MPI_Finalize();
}