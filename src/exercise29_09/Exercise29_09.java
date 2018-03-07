/**
 * Class: CSCI2410 Data Structures and Algorithms
 * Instructor: Y. Daniel Liang
 * Description: Uses a text file online to import pieces from a graph. 
 * Once the pieces are imported, information about the graph is printed
 * to the console. 
 * Due: 11/30/2016
 *
 * @author Shaun C. Dobbs
 * @version 1.0
 *
 * I pledge by honor that I have completed the programming assignment
 * independently. I have not copied the code from a student or any source. I
 * have not given my code to any student. * Sign here: Shaun C. Dobbs
 */
//URL https://liveexample.pearsoncmg.com/test/WeightedGraphSample.txt
package exercise29_09;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.*;


public class Exercise29_09 {

    public static void main(String[] args) {
        new Exercise29_09();
    }
    public Exercise29_09() {
        //Create scanner for the URL input
        Scanner s = new Scanner(System.in);

        //Arraylists for various pieces of the graph
        ArrayList<WeightedEdge> edges = new ArrayList<>();
        ArrayList<Integer> vertices = new ArrayList<>();
        ArrayList<String> lines = new ArrayList<>();

        //Get URL and close program if URL is not correct
        System.out.println("Enter a text file URL: ");
        URL url = null;
        try {
            url = new URL(s.nextLine());
        } catch (MalformedURLException e) {
            System.out.println("Invalid URL. Program will close.");
            System.exit(0);
        }

        try {
            URLConnection uc = url.openConnection();
            try (BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(uc.getInputStream()))) {
                String line;
                while ((line = bufferedReader.readLine()) != null) {
                    lines.add(line);
                }
            }
        } catch (IOException e) {
            System.out.println("Error with file I/O. Program will close.");
            System.exit(0);
        }

        //Parsing the integers
        lines.remove(0);
        for (String line : lines) {
            String[] portions = line.split("\\|");

            for (String portion : portions) {
                String[] e = portion.split("\\s*,\\s*");

                int u = Integer.parseInt(e[0].trim());
                int v = Integer.parseInt(e[1].trim());
                double weight = Double.parseDouble(e[2].trim());
                WeightedEdge weightedEdge1 = new WeightedEdge(u, v, weight);
                WeightedEdge weightedEdge2 = new WeightedEdge(v, u, weight);

                if (!vertices.contains(u)) {
                    vertices.add(u);
                }
                if (!vertices.contains(v)) {
                    vertices.add(v);
                }
                edges.add(weightedEdge1);
                edges.add(weightedEdge2);
            }
        }

        Collections.sort(vertices);
        Collections.sort(edges);

        // Create graph and print the appropriate information 
        WeightedGraph<Integer> graph = new WeightedGraph<>(vertices, edges);
        System.out.println("The number of vertices is " + graph.getVertices().size());
        graph.printWeightedEdges();
        WeightedGraph.MST mst = graph.getMinimumSpanningTree();
        System.out.println("Total weight in MST is " + mst.getTotalWeight());
        mst.printTree();
        


    }

    public class WeightedEdge extends Edge
            implements Comparable<WeightedEdge> {

        public double weight; // The weight on edge (u, v)

        /**
         * Create a weighted edge on (u, v)
         */
        public WeightedEdge(int u, int v, double weight) {
            super(u, v);
            this.weight = weight;
        }

        /**
         * Compare two edges on weights
         */
        public int compareTo(WeightedEdge edge) {
            if (weight > edge.weight) {
                return 1;
            } else if (weight == edge.weight) {
                return 0;
            } else {
                return -1;
            }
        }
    }

    public class WeightedGraph<V> extends UnweightedGraph<V> {

        /**
         * Construct an empty
         */
        public WeightedGraph() {
        }

        /**
         * Construct a WeightedGraph from vertices and edged in arrays
         */
        public WeightedGraph(V[] vertices, int[][] edges) {
            createWeightedGraph(java.util.Arrays.asList(vertices), edges);
        }

        /**
         * Construct a WeightedGraph from vertices and edges in list
         */
        public WeightedGraph(int[][] edges, int numberOfVertices) {
            List<V> vertices = new ArrayList<>();
            for (int i = 0; i < numberOfVertices; i++) {
                vertices.add((V) (new Integer(i)));
            }

            createWeightedGraph(vertices, edges);
        }

        /**
         * Construct a WeightedGraph for vertices 0, 1, 2 and edge list
         */
        public WeightedGraph(List<V> vertices, List<WeightedEdge> edges) {
            createWeightedGraph(vertices, edges);
        }

        /**
         * Construct a WeightedGraph from vertices 0, 1, and edge array
         */
        public WeightedGraph(List<WeightedEdge> edges,
                int numberOfVertices) {
            List<V> vertices = new ArrayList<>();
            for (int i = 0; i < numberOfVertices; i++) {
                vertices.add((V) (new Integer(i)));
            }

            createWeightedGraph(vertices, edges);
        }

        /**
         * Create adjacency lists from edge arrays
         */
        private void createWeightedGraph(List<V> vertices, int[][] edges) {
            this.vertices = vertices;

            for (int i = 0; i < vertices.size(); i++) {
                neighbors.add(new ArrayList<Edge>()); // Create a list for vertices
            }

            for (int i = 0; i < edges.length; i++) {
                neighbors.get(edges[i][0]).add(
                        new WeightedEdge(edges[i][0], edges[i][1], edges[i][2]));
            }
        }

        /**
         * Create adjacency lists from edge lists
         */
        private void createWeightedGraph(
                List<V> vertices, List<WeightedEdge> edges) {
            this.vertices = vertices;

            for (int i = 0; i < vertices.size(); i++) {
                neighbors.add(new ArrayList<Edge>()); // Create a list for vertices
            }

            for (WeightedEdge edge : edges) {
                neighbors.get(edge.u).add(edge); // Add an edge into the list
            }
        }

        /**
         * Return the weight on the edge (u, v)
         */
        public double getWeight(int u, int v) throws Exception {
            for (Edge edge : neighbors.get(u)) {
                if (edge.v == v) {
                    return ((WeightedEdge) edge).weight;
                }
            }

            throw new Exception("Edge does not exit");
        }

        /**
         * Display edges with weights
         */
        public void printWeightedEdges() {
            for (int i = 0; i < getSize(); i++) {
                System.out.print(getVertex(i) + " (" + i + "): ");
                for (Edge edge : neighbors.get(i)) {
                    System.out.print("(" + edge.u
                            + ", " + edge.v + ", " + ((WeightedEdge) edge).weight + ") ");
                }
                System.out.println();
            }
        }

        /**
         * Add an edge (u, v, weight) to the graph
         */
        public boolean addEdge(int u, int v, double weight) {
            return addEdge(new WeightedEdge(u, v, weight));
        }

        /**
         * Get a minimum spanning tree rooted at vertex 0
         */
        public MST getMinimumSpanningTree() {
            return getMinimumSpanningTree(0);
        }

        /**
         * Get a minimum spanning tree rooted at a specified vertex
         */
        public MST getMinimumSpanningTree(int startingVertex) {
            // cost[v] stores the cost by adding v to the tree
            double[] cost = new double[getSize()];
            for (int i = 0; i < cost.length; i++) {
                cost[i] = Double.POSITIVE_INFINITY; // Initial cost 
            }
            cost[startingVertex] = 0; // Cost of source is 0

            int[] parent = new int[getSize()]; // Parent of a vertex
            parent[startingVertex] = -1; // startingVertex is the root
            double totalWeight = 0; // Total weight of the tree thus far

            List<Integer> T = new ArrayList<>();

            // Expand T
            while (T.size() < getSize()) {
                // Find smallest cost v in V - T 
                int u = -1; // Vertex to be determined
                double currentMinCost = Double.POSITIVE_INFINITY;
                for (int i = 0; i < getSize(); i++) {
                    if (!T.contains(i) && cost[i] < currentMinCost) {
                        currentMinCost = cost[i];
                        u = i;
                    }
                }

                if (u == -1) {
                    break;
                } else {
                    T.add(u); // Add a new vertex to T
                }
                totalWeight += cost[u]; // Add cost[u] to the tree

                // Adjust cost[v] for v that is adjacent to u and v in V - T
                for (Edge e : neighbors.get(u)) {
                    if (!T.contains(e.v) && cost[e.v] > ((WeightedEdge) e).weight) {
                        cost[e.v] = ((WeightedEdge) e).weight;
                        parent[e.v] = u;
                    }
                }
            } // End of while

            return new MST(startingVertex, parent, T, totalWeight);
        }

        /**
         * MST is an inner class in WeightedGraph
         */
        public class MST extends SearchTree {

            private double totalWeight; // Total weight of all edges in the tree

            public MST(int root, int[] parent, List<Integer> searchOrder,
                    double totalWeight) {
                super(root, parent, searchOrder);
                this.totalWeight = totalWeight;
            }

            public double getTotalWeight() {
                return totalWeight;
            }
        }

        /**
         * Find single source shortest paths
         */
        public ShortestPathTree getShortestPath(int sourceVertex) {
            // cost[v] stores the cost of the path from v to the source
            double[] cost = new double[getSize()];
            for (int i = 0; i < cost.length; i++) {
                cost[i] = Double.POSITIVE_INFINITY; // Initial cost set to infinity
            }
            cost[sourceVertex] = 0; // Cost of source is 0

            // parent[v] stores the previous vertex of v in the path
            int[] parent = new int[getSize()];
            parent[sourceVertex] = -1; // The parent of source is set to -1

            // T stores the vertices whose path found so far
            List<Integer> T = new ArrayList<>();

            // Expand T
            while (T.size() < getSize()) {
                // Find smallest cost v in V - T 
                int u = -1; // Vertex to be determined
                double currentMinCost = Double.POSITIVE_INFINITY;
                for (int i = 0; i < getSize(); i++) {
                    if (!T.contains(i) && cost[i] < currentMinCost) {
                        currentMinCost = cost[i];
                        u = i;
                    }
                }

                if (u == -1) {
                    break;
                } else {
                    T.add(u); // Add a new vertex to T
                }
                // Adjust cost[v] for v that is adjacent to u and v in V - T
                for (Edge e : neighbors.get(u)) {
                    if (!T.contains(e.v)
                            && cost[e.v] > cost[u] + ((WeightedEdge) e).weight) {
                        cost[e.v] = cost[u] + ((WeightedEdge) e).weight;
                        parent[e.v] = u;
                    }
                }
            } // End of while

            // Create a ShortestPathTree
            return new ShortestPathTree(sourceVertex, parent, T, cost);
        }

        /**
         * ShortestPathTree is an inner class in WeightedGraph
         */
        public class ShortestPathTree extends SearchTree {

            private double[] cost; // cost[v] is the cost from v to source

            /**
             * Construct a path
             */
            public ShortestPathTree(int source, int[] parent,
                    List<Integer> searchOrder, double[] cost) {
                super(source, parent, searchOrder);
                this.cost = cost;
            }

            /**
             * Return the cost for a path from the root to vertex v
             */
            public double getCost(int v) {
                return cost[v];
            }

            /**
             * Print paths from all vertices to the source
             */
            public void printAllPaths() {
                System.out.println("All shortest paths from "
                        + vertices.get(getRoot()) + " are:");
                for (int i = 0; i < cost.length; i++) {
                    printPath(i); // Print a path from i to the source
                    System.out.println("(cost: " + cost[i] + ")"); // Path cost
                }
            }
        }
    }

    public interface Graph<V> {

        /**
         * Return the number of vertices in the graph
         */
        public int getSize();

        /**
         * Return the vertices in the graph
         */
        public java.util.List<V> getVertices();

        /**
         * Return the object for the specified vertex index
         */
        public V getVertex(int index);

        /**
         * Return the index for the specified vertex object
         */
        public int getIndex(V v);

        /**
         * Return the neighbors of vertex with the specified index
         */
        public java.util.List<Integer> getNeighbors(int index);

        /**
         * Return the degree for a specified vertex
         */
        public int getDegree(int v);

        /**
         * Print the edges
         */
        public void printEdges();

        /**
         * Clear the graph
         */
        public void clear();

        /**
         * Add a vertex to the graph
         */
        public boolean addVertex(V vertex);

        /**
         * Add an edge (u, v) to the graph
         */
        public boolean addEdge(int u, int v);

        /**
         * Add an edge to the graph
         */
        public boolean addEdge(Edge e);

        /**
         * Remove a vertex v from the graph
         */
        public default boolean remove(V v) {
            return true; // Implementation left as an exercise
        }

        /**
         * Remove an edge (u, v) from the graph
         */
        public default boolean remove(int u, int v) {
            return true; // Implementation left as an exercise
        }

        /**
         * Obtain a depth-first search tree
         */
        public UnweightedGraph<V>.SearchTree dfs(int v);

        /**
         * Obtain a breadth-first search tree
         */
        public UnweightedGraph<V>.SearchTree bfs(int v);
    }

    public class Edge {

        int u;
        int v;

        public Edge(int u, int v) {
            this.u = u;
            this.v = v;
        }

        public boolean equals(Object o) {
            return u == ((Edge) o).u && v == ((Edge) o).v;
        }
    }

    public class UnweightedGraph<V> implements Graph<V> {

        protected List<V> vertices = new ArrayList<>(); // Store vertices
        protected List<List<Edge>> neighbors
                = new ArrayList<>(); // Adjacency lists

        /**
         * Construct an empty graph
         */
        public UnweightedGraph() {
        }

        /**
         * Construct a graph from vertices and edges stored in arrays
         */
        public UnweightedGraph(V[] vertices, int[][] edges) {
            for (int i = 0; i < vertices.length; i++) {
                addVertex(vertices[i]);
            }

            createAdjacencyLists(edges, vertices.length);
        }

        /**
         * Construct a graph from vertices and edges stored in List
         */
        public UnweightedGraph(List<V> vertices, List<Edge> edges) {
            for (int i = 0; i < vertices.size(); i++) {
                addVertex(vertices.get(i));
            }

            createAdjacencyLists(edges, vertices.size());
        }

        /**
         * Construct a graph for integer vertices 0, 1, 2 and edge list
         */
        public UnweightedGraph(List<Edge> edges, int numberOfVertices) {
            for (int i = 0; i < numberOfVertices; i++) {
                addVertex((V) (new Integer(i))); // vertices is {0, 1, ...}
            }
            createAdjacencyLists(edges, numberOfVertices);
        }

        /**
         * Construct a graph from integer vertices 0, 1, and edge array
         */
        public UnweightedGraph(int[][] edges, int numberOfVertices) {
            for (int i = 0; i < numberOfVertices; i++) {
                addVertex((V) (new Integer(i))); // vertices is {0, 1, ...}
            }
            createAdjacencyLists(edges, numberOfVertices);
        }

        /**
         * Create adjacency lists for each vertex
         */
        private void createAdjacencyLists(
                int[][] edges, int numberOfVertices) {
            for (int i = 0; i < edges.length; i++) {
                addEdge(edges[i][0], edges[i][1]);
            }
        }

        /**
         * Create adjacency lists for each vertex
         */
        private void createAdjacencyLists(
                List<Edge> edges, int numberOfVertices) {
            for (Edge edge : edges) {
                addEdge(edge.u, edge.v);
            }
        }

        @Override
        /**
         * Return the number of vertices in the graph
         */
        public int getSize() {
            return vertices.size();
        }

        @Override
        /**
         * Return the vertices in the graph
         */
        public List<V> getVertices() {
            return vertices;
        }

        @Override
        /**
         * Return the object for the specified vertex
         */
        public V getVertex(int index) {
            return vertices.get(index);
        }

        @Override
        /**
         * Return the index for the specified vertex object
         */
        public int getIndex(V v) {
            return vertices.indexOf(v);
        }

        @Override
        /**
         * Return the neighbors of the specified vertex
         */
        public List<Integer> getNeighbors(int index) {
            List<Integer> result = new ArrayList<>();
            for (Edge e : neighbors.get(index)) {
                result.add(e.v);
            }

            return result;
        }

        @Override
        /**
         * Return the degree for a specified vertex
         */
        public int getDegree(int v) {
            return neighbors.get(v).size();
        }

        @Override
        /**
         * Print the edges
         */
        public void printEdges() {
            for (int u = 0; u < neighbors.size(); u++) {
                System.out.print(getVertex(u) + " (" + u + "): ");
                for (Edge e : neighbors.get(u)) {
                    System.out.print("(" + getVertex(e.u) + ", "
                            + getVertex(e.v) + ") ");
                }
                System.out.println();
            }
        }

        @Override
        /**
         * Clear the graph
         */
        public void clear() {
            vertices.clear();
            neighbors.clear();
        }

        @Override
        /**
         * Add a vertex to the graph
         */
        public boolean addVertex(V vertex) {
            if (!vertices.contains(vertex)) {
                vertices.add(vertex);
                neighbors.add(new ArrayList<Edge>());
                return true;
            } else {
                return false;
            }
        }

        @Override
        /**
         * Add an edge to the graph
         */
        public boolean addEdge(Edge e) {
            if (e.u < 0 || e.u > getSize() - 1) {
                throw new IllegalArgumentException("No such index: " + e.u);
            }

            if (e.v < 0 || e.v > getSize() - 1) {
                throw new IllegalArgumentException("No such index: " + e.v);
            }

            if (!neighbors.get(e.u).contains(e)) {
                neighbors.get(e.u).add(e);
                return true;
            } else {
                return false;
            }
        }

        @Override
        /**
         * Add an edge to the graph
         */
        public boolean addEdge(int u, int v) {
            return addEdge(new Edge(u, v));
        }

        @Override
        /**
         * Obtain a DFS tree starting from vertex v
         */
        /**
         * To be discussed in Section 28.7
         */
        public SearchTree dfs(int v) {
            List<Integer> searchOrder = new ArrayList<>();
            int[] parent = new int[vertices.size()];
            for (int i = 0; i < parent.length; i++) {
                parent[i] = -1; // Initialize parent[i] to -1
            }
            // Mark visited vertices
            boolean[] isVisited = new boolean[vertices.size()];

            // Recursively search
            dfs(v, parent, searchOrder, isVisited);

            // Return a search tree
            return new SearchTree(v, parent, searchOrder);
        }

        /**
         * Recursive method for DFS search
         */
        private void dfs(int u, int[] parent, List<Integer> searchOrder,
                boolean[] isVisited) {
            // Store the visited vertex
            searchOrder.add(u);
            isVisited[u] = true; // Vertex v visited

            for (Edge e : neighbors.get(u)) {
                if (!isVisited[e.v]) {
                    parent[e.v] = u; // The parent of vertex e.v is u
                    dfs(e.v, parent, searchOrder, isVisited); // Recursive search
                }
            }
        }

        @Override
        /**
         * Starting bfs search from vertex v
         */
        /**
         * To be discussed in Section 28.9
         */
        public SearchTree bfs(int v) {
            List<Integer> searchOrder = new ArrayList<>();
            int[] parent = new int[vertices.size()];
            for (int i = 0; i < parent.length; i++) {
                parent[i] = -1; // Initialize parent[i] to -1
            }
            java.util.LinkedList<Integer> queue
                    = new java.util.LinkedList<>(); // list used as a queue
            boolean[] isVisited = new boolean[vertices.size()];
            queue.offer(v); // Enqueue v
            isVisited[v] = true; // Mark it visited

            while (!queue.isEmpty()) {
                int u = queue.poll(); // Dequeue to u
                searchOrder.add(u); // u searched
                for (Edge e : neighbors.get(u)) {
                    if (!isVisited[e.v]) {
                        queue.offer(e.v); // Enqueue w
                        parent[e.v] = u; // The parent of w is u
                        isVisited[e.v] = true; // Mark it visited
                    }
                }
            }

            return new SearchTree(v, parent, searchOrder);
        }

        /**
         * Tree inner class inside the AbstractGraph class
         */
        /**
         * To be discussed in Section 28.6
         */
        public class SearchTree {

            private int root; // The root of the tree
            private int[] parent; // Store the parent of each vertex
            private List<Integer> searchOrder; // Store the search order

            /**
             * Construct a tree with root, parent, and searchOrder
             */
            public SearchTree(int root, int[] parent,
                    List<Integer> searchOrder) {
                this.root = root;
                this.parent = parent;
                this.searchOrder = searchOrder;
            }

            /**
             * Return the root of the tree
             */
            public int getRoot() {
                return root;
            }

            /**
             * Return the parent of vertex v
             */
            public int getParent(int v) {
                return parent[v];
            }

            /**
             * Return an array representing search order
             */
            public List<Integer> getSearchOrder() {
                return searchOrder;
            }

            /**
             * Return number of vertices found
             */
            public int getNumberOfVerticesFound() {
                return searchOrder.size();
            }

            /**
             * Return the path of vertices from a vertex to the root
             */
            public List<V> getPath(int index) {
                ArrayList<V> path = new ArrayList<>();

                do {
                    path.add(vertices.get(index));
                    index = parent[index];
                } while (index != -1);

                return path;
            }

            /**
             * Print a path from the root to vertex v
             */
            public void printPath(int index) {
                List<V> path = getPath(index);
                System.out.print("A path from " + vertices.get(root) + " to "
                        + vertices.get(index) + ": ");
                for (int i = path.size() - 1; i >= 0; i--) {
                    System.out.print(path.get(i) + " ");
                }
            }

            /**
             * Print the whole tree
             */
            public void printTree() {
                System.out.println("Root is: " + vertices.get(root));
                System.out.print("Edges: ");
                for (int i = 0; i < parent.length; i++) {
                    if (parent[i] != -1) {
                        // Display an edge
                        System.out.print("(" + vertices.get(parent[i]) + ", "
                                + vertices.get(i) + ") ");
                    }
                }
                System.out.println();
            }
        }
    }

}
