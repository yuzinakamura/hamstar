
// data associated with one state
class StateData {
public:
  // back-pointer to previous state along the discovered path
  vector<int> bp;
  // stores Boolean flags such as whether the state is CLOSED
  int mask;
  // g is the cost of a discovered path, h estimates the remaining cost to goal
  Cost g;
  vector<Cost> h;
  // points to the state's location in the priority queue
  vector<multimap<Cost, State>::const_iterator> iter;
  // initialization
  StateData() : mask(0), g(INFINITE) {}
};



// a process which performs a search, typically assigned to one thread on one machine
class Searcher {
public:
  // multi-heuristic mixing coefficients
  double MD, LC, MT;
  // Searcher-specific closing mask
  int closing_mask;
  // number of states seen and expanded
  int num_discovered;
  int num_expanded;
  // id of the Searcher
  int comm_rank;
  // the priority queue or search frontier
  multimap<Cost, State> open;

  // key for the priority queue
  Cost f(StateData& s_data);

  void insert(const State& s, StateData& s_data);

  void erase(const State& s, StateData& s_data);

  // assumes start, goal, w1 and w2 are already set
  void init();

  void expand(const State& s, StateData& s_data);

  void runSingleIteration();

  // get the path: a sequence of states from start to goal
  vector<State> getSolution();
};