// TODO: replace final broadcast - just send a normal message back to anchor
#include <array>
#include <cmath>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <cassert>
#include <vector>
#include <queue>
#include <set>
#include <unordered_set>
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
constexpr int MASK_CLOSED = 1;
constexpr int MASK_CLOSED_ANCHOR = 2;
constexpr int GRID_ROWS = 4;
constexpr int GRID_COLS = 4;
constexpr int NUM_MOVES = 4;
constexpr double TIME_LIMIT = 300;
constexpr Cost INFINITE = 1e30;

//communication constants
constexpr bool COMM_ITER = true;
constexpr double COMM_INTERVAL = 0.05;
constexpr int COMM_FREQ = 100000;
constexpr int BUFFER_SIZE = COMM_FREQ*NUM_MOVES*4; //Number of new g values before a message is sent. The actual buffer is double this, because it needs the node too.
constexpr int DATUM_SIZE = (GRID_ROWS * GRID_COLS) + 1 + 2; //4x4 state, backtrace, mask, g, and h. Assumes Cost is same size as int
constexpr int HEAD_NODE = 0;
constexpr int SLOWDOWN = 1;	//Factor to slow down heuristic calculation by. Must be an int and >= 1. Make 1 for no slowdown

static int comm_size;
static int comm_rank;
static int finished = 0;
int SEED = 1;
string filename = "test";

// right, up, left, down
constexpr array<int, NUM_MOVES> moves = { 1, -GRID_COLS, -1, GRID_COLS };

struct state_hash {
	inline std::size_t operator()(const State& v) const {
		std::size_t returnVal = 0;
		for (int i = 1; i < GRID_ROWS*GRID_COLS; i++) {
			returnVal += v.first[i];
			returnVal *= 7;
		}
		return returnVal;
	}
};

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

// data associated with one state
class StateData {
public:
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

// a process which performs a search, typically assigned to one thread on one machine
class Searcher {
public:
  // search parameters: endpoints and weights
  State start, goal;
  double w1, w2;
  // multi-heuristic mixing coefficients
  double MD, LC, MT;
  // Searcher-specific closing mask
  int closing_mask;
  // time when the search began
  Clock::time_point start_time;
  double time_startup; //Time spent on startup
  double time_waiting_gather; //Time idle waiting for comm
  double time_waiting_bcast; //Time idle waiting for comm
  double time_comm; //Time idle waiting for comm
  double time_gather; //Time spent on MPI_Gather-related stuff
  double time_bcast; //Time spend on MPI_Bcast-related stuff
  double time_elapsed; //Time since start of run to finding a solution
  double time_set; //Time spent by HEAD_NODE adding States to updated_states set
  double time_state; //Time spent by HEAD_NODE generating States from incoming data
  // the priority queue or search frontier
  multimap<Cost, State> open;
  // a dictionary of seen states and their corresponding data
  map<State, StateData> data;
  // a value that bounds the optimal solution cost
  Cost opt_bound;
  // number of states seen and expanded
  int num_discovered;
  int num_expanded;
  int total_discovered;
  int total_expanded;

  //Communication stuff
  int * head_buffer;
  int * child_buffer;
  //set<State> updated_states;
  unordered_set< State, state_hash> updated_states;

  // determine the heuristics to be used
  void computeH(const State& s, StateData& s_data) {
	  Cost hMD;
	  Cost hLC;
	  Cost hMT;
	  for (int i = 0; i < SLOWDOWN; i++) {
		  hMD = manhattanDist(s, goal);
		  hLC = linearConflicts(s, goal);
		  hMT = misplacedTiles(s, goal);
	  }
    s_data.hAnch = hMD + hLC;
    s_data.h = MD*hMD + LC*hLC + MT*hMT;
    s_data.iter = open.cend();
  }

  // key for the priority queue
  Cost f(StateData& s_data) { return s_data.g + w1*s_data.h; }

  void insert(const State& s, StateData& s_data) {
    // if state is already in open[id], remove it before re-inserting
    if (s_data.iter != open.cend())
      open.erase(s_data.iter);
    // insert state with new key
    s_data.iter = open.insert(pair<Cost, State>(f(s_data), s));
  }

  void erase(const State& s, StateData& s_data) {
    if (s_data.iter != open.cend()) {
      open.erase(s_data.iter);
      s_data.iter = open.cend();
    }
  }

  // assumes start, goal, w1 and w2 are already set
  void init() {
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

    open.clear();
    data.clear();

    head_buffer = new int[(BUFFER_SIZE * DATUM_SIZE + 4) * comm_size];
    child_buffer = new int[(BUFFER_SIZE * DATUM_SIZE + 5) * comm_size]; //The +1 is a hack. It's a cost

    // create data entry for start state
    StateData& start_data = data[start];
    start_data.g = 0;
    computeH(start, start_data);
    insert(start, start_data);

    opt_bound = 0;
    // the start state is discovered but not yet expanded
    num_discovered = 1;
    num_expanded = 0;
    total_discovered = 1;
    total_expanded = 0;
    start_time = Clock::now();
	time_waiting_gather = 0;
	time_waiting_bcast = 0;
	time_comm = 0;
	time_bcast = 0;
	time_gather = 0;
	time_set = 0;
	time_state = 0;
  }

  void expand(const State& s, StateData& s_data) {
    assert(!(s_data.mask & closing_mask));
    s_data.mask |= MASK_CLOSED | closing_mask;

    erase(s, s_data);

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
        if (!(t_data.mask & closing_mask)) {
          insert(t, t_data);
          updated_states.insert(t);
        }
      }
    }
  }

  // when passed a buffer, and a starting point in that buffer, interprets a state, updates its stateData, and returns it.
  State update_state(int * buffer, int& index) {
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

    assert(sizeof(Cost) == sizeof(int));
    int maskNew;
    Cost gNew, hAnchNew;
    memcpy(&maskNew, &buffer[index++], sizeof(int));
    memcpy(&gNew, &buffer[index++], sizeof(Cost));
    memcpy(&hAnchNew, &buffer[index++], sizeof(Cost));

    StateData& s_data = data[s];
    s_data.mask |= maskNew;
    // if necessary, compute heuristics
    if (s_data.g == INFINITE) {
      if (comm_rank == HEAD_NODE) {
        s_data.h = s_data.hAnch = hAnchNew;
        s_data.iter = open.cend();
      }
      else {
        computeH(s, s_data);
      }
    }
    // if necessary, close this state
    if ((s_data.mask & closing_mask) && s_data.iter != open.cend()) {
      erase(s, s_data);
    }
    // if necessary, update g and bp
    if (s_data.g > gNew) {
      s_data.g = gNew;
      if (!(s_data.mask & closing_mask)) {
        insert(s, s_data);
      }
    }

    assert(index == indexStart + DATUM_SIZE);
    return s;
  }

  //When passed in a state, it will serialize the state's data and insert it into the buffer starting at index
  void serialize_state(const State * state, StateData * s_data, int * buffer, int& index) {
    int indexStart = index; // sanity check
    
    for (int pos = 0; pos < GRID_ROWS * GRID_COLS; pos++) {
      buffer[index++] = state->first[pos];
    }
    buffer[index++] = s_data->mask;
    memcpy(buffer+index, &s_data->g, sizeof(Cost));
    index += sizeof(Cost) / sizeof(int);
    memcpy(buffer+index, &s_data->hAnch, sizeof(Cost));
    index += sizeof(Cost) / sizeof(int);

    assert(index == indexStart + DATUM_SIZE);
  }

  void run() {
	  int benchmark = 0;
    int flag = 0;
    MPI_Status status;
    MPI_Request request;
    int error = -1;

    cout << "Starting run" << endl;
    int iter = 0;
    Clock::time_point last_time = Clock::now();
    // repeat until some thread declares the search to be finished
    while (true) {
		if (num_expanded > benchmark) {
			benchmark += 100000;
			cout << "Machine " << comm_rank << " expanded " << num_expanded << endl;
		}

      if (comm_rank == HEAD_NODE) {
        if (open.empty()) {
          opt_bound = INFINITE;
        }
        else {
          opt_bound = max(opt_bound, open.cbegin()->first);
        }
      }
	  Clock::time_point before_comm = Clock::now();
	  if ((COMM_ITER && iter % COMM_FREQ == 0) || (!COMM_ITER && iter % 10 == 0 && chrono::duration<double, chrono::seconds::period>(Clock::now() - last_time).count() > COMM_INTERVAL)) {
		  last_time = Clock::now();
        //cout << "Process " << comm_rank << " preparing to communicate " << endl;
        //Handle communication
        //First, Gather to HEAD_NODE
		  Clock::time_point before_gather = Clock::now();
        if (comm_rank == HEAD_NODE) {
          memset(head_buffer, 0, sizeof(int)*(BUFFER_SIZE*DATUM_SIZE + 4)*comm_size);
          //cout << "Process " << comm_rank << " receiving gather" << endl;
          //MPI_Gather(MPI_IN_PLACE, 0, MPI_INT, head_buffer, BUFFER_SIZE*DATUM_SIZE+1, MPI_INT, HEAD_NODE, MPI_COMM_WORLD);
		  Clock::time_point before_time = Clock::now();
          MPI_Gather(MPI_IN_PLACE, BUFFER_SIZE*DATUM_SIZE + 4, MPI_INT, head_buffer, BUFFER_SIZE*DATUM_SIZE + 4, MPI_INT, HEAD_NODE, MPI_COMM_WORLD);
		  time_waiting_gather += chrono::duration<double, chrono::seconds::period>(Clock::now() - before_time).count();
          //cout << "Process " << comm_rank << " received gather" << endl;
          int index = 0;
		  total_discovered = num_discovered;
		  total_expanded = num_expanded;
          for (int machine = 0; machine < comm_size; machine++) {
            int num_sent = head_buffer[machine*(BUFFER_SIZE*DATUM_SIZE+4)]; //Number of states this child sent in.
			total_discovered += head_buffer[machine*(BUFFER_SIZE*DATUM_SIZE + 4)+1];
			total_expanded += head_buffer[machine*(BUFFER_SIZE*DATUM_SIZE + 4)+2];
			assert(head_buffer[machine*(BUFFER_SIZE*DATUM_SIZE + 4) + 3] == 1 || head_buffer[machine*(BUFFER_SIZE*DATUM_SIZE + 4) + 3] == 0);
			int in_finished = head_buffer[machine*(BUFFER_SIZE*DATUM_SIZE + 4) + 3];
			finished = max(finished,in_finished);
			
            index = machine*(BUFFER_SIZE*DATUM_SIZE + 4) + 4;
            while (num_sent > 0) {
              num_sent--;
			  Clock::time_point before_state = Clock::now();
              State s = update_state(head_buffer, index);
			  time_state += chrono::duration<double, chrono::seconds::period>(Clock::now() - before_state).count();
			  Clock::time_point before_set = Clock::now();
              updated_states.insert(s);
			  time_set += chrono::duration<double, chrono::seconds::period>(Clock::now() - before_set).count();
            }
          }
        }
        else {
          memset(child_buffer, 0, sizeof(int)*(BUFFER_SIZE*DATUM_SIZE + 4)*comm_size);
          int index = 0;
          child_buffer[index++] = updated_states.size();
		  assert(updated_states.size() <= COMM_FREQ*4);
		  child_buffer[index++] = num_discovered;
		  child_buffer[index++] = num_expanded;
		  child_buffer[index++] = finished;
          //cout << "States updated: " << child_buffer[0] << endl;
          auto it = updated_states.cbegin();
          for (; it != updated_states.cend(); it++) {
            const State& s = *it;
            StateData& s_data = data[s];
            serialize_state(&s, &s_data, child_buffer, index);
          }
          updated_states.clear();
          //cout << "Process " << comm_rank << " sending gather" << endl;
		  Clock::time_point before_time = Clock::now();
          MPI_Gather(child_buffer, BUFFER_SIZE*DATUM_SIZE + 4, MPI_INT, NULL, BUFFER_SIZE*DATUM_SIZE + 4, MPI_INT, HEAD_NODE, MPI_COMM_WORLD);
		  time_waiting_gather += chrono::duration<double, chrono::seconds::period>(Clock::now() - before_time).count();
          //cout << "Process " << comm_rank << " sent gather" << endl;
        }
		time_gather += chrono::duration<double, chrono::seconds::period>(Clock::now() - before_gather).count();

        //Now, Bcast
		Clock::time_point before_bcast = Clock::now();
        //cout << "Process " << comm_rank << " Beginning bcast" << endl;
        memset(child_buffer, 0, sizeof(int)*(BUFFER_SIZE*DATUM_SIZE + 4)*comm_size);
        if (comm_rank == HEAD_NODE) {
          int index = 0;
          memcpy(&child_buffer[index++], &opt_bound, sizeof(Cost));
          child_buffer[index++] = updated_states.size();
		  assert(updated_states.size() <= BUFFER_SIZE*comm_size);
		  child_buffer[index++] = total_discovered;
		  child_buffer[index++] = total_expanded;
		  child_buffer[index++] = finished;
          auto it = updated_states.cbegin();
          for (; it != updated_states.cend(); it++) {
            const State& s = *it;
            StateData& s_data = data[s];
            serialize_state(&s, &s_data, child_buffer, index);
          }
          updated_states.clear();
        }
        //cout << "Process " << comm_rank << " Sending bcast" << endl;
		Clock::time_point before_time = Clock::now();
        MPI_Bcast(child_buffer, (BUFFER_SIZE*DATUM_SIZE + 5)*comm_size, MPI_INT, HEAD_NODE, MPI_COMM_WORLD);
		time_waiting_bcast += chrono::duration<double, chrono::seconds::period>(Clock::now() - before_time).count();
        //cout << "Process " << comm_rank << " sent bcast " << endl;
        if (comm_rank != HEAD_NODE) {
          int index = 0;
          memcpy(&opt_bound, &child_buffer[index++], sizeof(Cost));
          int num_sent = child_buffer[index++];
		  total_discovered = child_buffer[index++];
		  total_expanded = child_buffer[index++];
		  assert(child_buffer[index] == 0 || child_buffer[index] == 1);
		  finished = max(finished, child_buffer[index++]);
          assert(num_sent >= 0 && num_sent < 100000000);
          while (num_sent > 0) {
            num_sent--;
            State s = update_state(child_buffer, index);
            //updated_states.insert(s);
          }
        }
		time_bcast += chrono::duration<double, chrono::seconds::period>(Clock::now() - before_bcast).count();

        //cout << "Process " << comm_rank << " finished comm" << endl;
		if (finished) {
			cout << "Machine " << comm_rank << " received a finish variable " << endl;
			cout << "Time spent wating gather: " << time_waiting_gather << endl;
			cout << "Time spent wating bcast: " << time_waiting_bcast << endl;
			cout << "Time spent comming: " << time_comm << endl;
			cout << "Time spent bcasting: " << time_bcast << endl;
			cout << "Time spent gathering: " << time_gather << endl;
			cout << "Time spent setting: " << time_set << endl;
			cout << "Time spent stating: " << time_state << endl;
			return;
		}
      }
	  time_comm += chrono::duration<double, chrono::seconds::period>(Clock::now() - before_comm).count();

      // get a State to expand
      auto it = open.cbegin();

      // timing
      Clock::time_point cur_time = Clock::now();
      time_elapsed = chrono::duration<double, chrono::seconds::period>(cur_time - start_time).count();

      // search termination condition
      if (data[goal].g <= w2 * opt_bound) {
        // flush and kill, but don't kill the master!!!!!
		  finished = 1;
		  cout << "Machine " << comm_rank << " found a solution" << endl;
		  cout << "Time taken: " << time_elapsed << ". Total expanded: " << total_expanded << " Local expanded: " << num_expanded << endl;
		  cout << "Time waited gathered: " << time_waiting_gather << endl;
		  cout << "Time waited bcasted: " << time_waiting_bcast << endl;
		  cout << "Time commed: " << time_comm << endl;
		  cout << "Time spent bcasted: " << time_bcast << endl;
		  cout << "Time spent gathered: " << time_gather << endl;
		  cout << "Time spent setted: " << time_set << endl;
		  cout << "Time spent stated: " << time_state << endl;

		  //Now do one final comm to inform the others.
		  if (comm_rank == HEAD_NODE) {
			  memset(head_buffer, 0, sizeof(int)*(BUFFER_SIZE*DATUM_SIZE + 4)*comm_size);
			  //cout << "Process " << comm_rank << " receiving gather" << endl;
			  //MPI_Gather(MPI_IN_PLACE, 0, MPI_INT, head_buffer, BUFFER_SIZE*DATUM_SIZE+1, MPI_INT, HEAD_NODE, MPI_COMM_WORLD);
			  MPI_Gather(MPI_IN_PLACE, BUFFER_SIZE*DATUM_SIZE + 4, MPI_INT, head_buffer, BUFFER_SIZE*DATUM_SIZE + 4, MPI_INT, HEAD_NODE, MPI_COMM_WORLD);
			  //cout << "Process " << comm_rank << " received gather" << endl;
			  int index = 0;
			  total_discovered = num_discovered;
			  total_expanded = num_expanded;
			  for (int machine = 0; machine < comm_size; machine++) {
				  int num_sent = head_buffer[machine*(BUFFER_SIZE*DATUM_SIZE + 4)]; //Number of states this child sent in.
				  total_discovered += head_buffer[machine*(BUFFER_SIZE*DATUM_SIZE + 4) + 1];
				  total_expanded += head_buffer[machine*(BUFFER_SIZE*DATUM_SIZE + 4) + 2];
				  assert(head_buffer[machine*(BUFFER_SIZE*DATUM_SIZE + 4) + 3] == 1 || head_buffer[machine*(BUFFER_SIZE*DATUM_SIZE + 4) + 3] == 0);
				  int in_finished = head_buffer[machine*(BUFFER_SIZE*DATUM_SIZE + 4) + 3];
				  finished = max(finished, in_finished);

				  index = machine*(BUFFER_SIZE*DATUM_SIZE + 4) + 4;
				  while (num_sent > 0) {
					  num_sent--;
					  State s = update_state(head_buffer, index);
					  updated_states.insert(s);
				  }
			  }
		  }
		  else {
			  memset(child_buffer, 0, sizeof(int)*(BUFFER_SIZE*DATUM_SIZE + 4)*comm_size);
			  int index = 0;
			  child_buffer[index++] = updated_states.size();
			  child_buffer[index++] = num_discovered;
			  child_buffer[index++] = num_expanded;
			  child_buffer[index++] = finished;
			  //cout << "States updated: " << child_buffer[0] << endl;
			  auto it = updated_states.cbegin();
			  for (; it != updated_states.cend(); it++) {
				  const State& s = *it;
				  StateData& s_data = data[s];
				  serialize_state(&s, &s_data, child_buffer, index);
			  }
			  updated_states.clear();
			  //cout << "Process " << comm_rank << " sending gather" << endl;
			  MPI_Gather(child_buffer, BUFFER_SIZE*DATUM_SIZE + 4, MPI_INT, NULL, BUFFER_SIZE*DATUM_SIZE + 4, MPI_INT, HEAD_NODE, MPI_COMM_WORLD);
			  //cout << "Process " << comm_rank << " sent gather" << endl;
		  }

		  //Now, Bcast
		  //cout << "Process " << comm_rank << " Beginning bcast" << endl;
		  memset(child_buffer, 0, sizeof(int)*(BUFFER_SIZE*DATUM_SIZE + 4)*comm_size);
		  if (comm_rank == HEAD_NODE) {
			  int index = 0;
			  memcpy(&child_buffer[index++], &opt_bound, sizeof(Cost));
			  child_buffer[index++] = updated_states.size();
			  child_buffer[index++] = total_discovered;
			  child_buffer[index++] = total_expanded;
			  child_buffer[index++] = finished;
			  auto it = updated_states.cbegin();
			  for (; it != updated_states.cend(); it++) {
				  const State& s = *it;
				  StateData& s_data = data[s];
				  serialize_state(&s, &s_data, child_buffer, index);
			  }
			  updated_states.clear();
		  }
		  //cout << "Process " << comm_rank << " Sending bcast" << endl;
		  MPI_Bcast(child_buffer, (BUFFER_SIZE*DATUM_SIZE + 5)*comm_size, MPI_INT, HEAD_NODE, MPI_COMM_WORLD);
		  //cout << "Process " << comm_rank << " sent bcast " << endl;
		  if (comm_rank != HEAD_NODE) {
			  int index = 0;
			  memcpy(&opt_bound, &child_buffer[index++], sizeof(Cost));
			  int num_sent = child_buffer[index++];
			  total_discovered = child_buffer[index++];
			  total_expanded = child_buffer[index++];
			  assert(child_buffer[index] == 0 || child_buffer[index] == 1);
			  finished = max(finished, child_buffer[index++]);
			  assert(num_sent >= 0 && num_sent < 100000000);
			  while (num_sent > 0) {
				  num_sent--;
				  State s = update_state(child_buffer, index);
				  //updated_states.insert(s);
			  }
		  }
        break;
      }
      if (it == open.cend()) {
        continue;
      }
      if (time_elapsed > TIME_LIMIT) {
        break;
      }
      State s = it->second;
      StateData& s_data = data[s];
      num_expanded++;
      total_expanded++;
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
      completeHalfState(s);
      sol.push_back(s);
    }
    reverse(sol.begin(), sol.end());
    return sol;
  }
};
Searcher searcher;

void prepareDistributedSearch() {
  searcher.w1 = 1.0;
  searcher.w2 = 1.0;
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

  for (int i = 1; i < argc; i++) {
    string arg = argv[i];
    if (arg == "-seed" || arg == "-s") {
      SEED = atoi(argv[++i]);
    }
    if (arg == "-filename" || arg == "-f") {
      filename = argv[++i];
    }
  }

  FILE* fout = fopen((filename + ".csv").c_str(), "a");
  // we can set it up to loop over multiple problem instances
  for (int i = 0; i < 1/*TRIALS*/; i++) {
    // prepare the Searcher processes
    prepareDistributedSearch();
    searcher.init();

    // run planner and time its execution
    searcher.run();
    
	if (comm_rank == HEAD_NODE) {
		Cost path_length = searcher.data[searcher.goal].g;
		// print solution if it was found
		if (path_length <= searcher.w2 * searcher.opt_bound) {
			/*for (State& s : searcher.getSolution()) {
				printState(s);
			}*/
      cout << "Found path of length " << path_length << " Discovered nodes = " << searcher.num_discovered << ". Expanded nodes: " << searcher.num_expanded << ". Time: " << searcher.time_elapsed << endl;
		}
		cout << "Total discovered: " << searcher.total_discovered << " total expanded: " << searcher.total_expanded << endl;
    fprintf(fout, "%f %f %f %d %f\n", searcher.w1, searcher.w2, searcher.time_elapsed, searcher.total_expanded, path_length);
		//fprintf(fout, "%f %f %f %d %f\n", searcher.w1, searcher.w2, searcher.time_elapsed, searcher.num_expanded, path_length);
	}
  }
  fclose(fout);
  //MPI_Finalize();
}
