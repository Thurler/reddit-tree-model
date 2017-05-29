# reddit-tree-model
Final project for complex networks class at UFRJ. The proposal is to use Kaggle's Reddit dataset to bring forth a model describing Reddit comment trees.

## R(t,p) model
A simple model inspired in thr BA model for the representation of reddit comment threads.

### Parameters:
* t: maximum number of iterations
* p: the probability function that returns the probability that the current comment will receive a reply.

### Initialization:
Starts with a graph containing only one vertex, the tree root (first comment in the thread). A cursor is pointing to this vertex, setting it as the 'current vertex'.

### Step:
For each iteration (from 1 to t):
* Decides with probability p if the current vertex (vertex where the cursor is at) will receive a reply.
  * If it does, then a vertex is added to the graph connecting the current vertex with the new one, directed from the new to the current. Then the cursor goes back to the tree root (first comment in the thread) and the iteration ends.
  * If doesn't receive a reply, then the cursor goes down a level of the tree.
    * If the current vertex has no neighbors, then the cursor goes back to the tree root (first comment in the thread) and the iteration ends.
    * Else, the next vertex is chosen in the next level with probability proportional to it's in-degree, so that the more neighbors the vertex has, the more it is preferred, creating a PA dynamic. The current vertex is set to the selected vertex. Iteration ends.

### Variations
There are a few variations on the basic step that can and have been implemented in this project. They will be documented here (and not only in the code itself) as soon as we have time for it.

### Model simulator

#### Single thread - rp-simulator.py
Creates a single graph file AT --out_dir/allgraphs. Also plot the graph at --out_dir/plot if --draw is set. Other parameters seen using --help.

Example using mostly default values

```python rp-simulator --p 0.01 --out_dir out```

#### Batch process - run-simulation.py
Runs the process (rp-simulator.py) for each configuration --runs times. At the end, gets all the files at --out_dir/allgraphs and, for each detected configuration, concatenates all the graphs in a single file, saving a graph file containing all threads at --out_dir/finalgraphs

Example running 5 times for ten distinct configurations, drawing each thread at --out_dir/plot/:

```python run-simulation.py --out_dir out --runs 5 --draw True --p_min 0.01 --p_max 0.1 --p_step 0.01```
