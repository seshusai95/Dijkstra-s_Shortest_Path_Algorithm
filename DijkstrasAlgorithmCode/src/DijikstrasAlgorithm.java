import java.io.*;
import java.util.*;
import java.util.Map.Entry;

public class DijikstrasAlgorithm {
 //Variable to store start vertex
 private static String START = "";
 //Variable to store end vertex
 private static String END = "";
 //Variable to store vertex 1 from the given input
 public static String node1 = "";
 //Variable to store vertex 2 from the given input
 public static String node2 = "";
 //Variable to store edge type from the given input
 public static int type = 0;
 //Variable to store edge's alpha from the given input
 public static double alpha = 0;
 //Variable to store edge's beta from the given input
 public static double beta = 0;
 //Variable to store edge's mean calculated from the given input
 public static double mean = 0;
 //Variable to store edge's variance calculated  from the given input
 public static double variance = 0;
 //Variable to store edge's cSquare calculated from the given input
 public static double cSquare = 0;
 public static double dist = 0;
 public static int inputStartEndCount = 0;
 //Variable to store number of hops for each criteria
 public static int hopCount = 0;
 public static int edgesCount = 0;
 public static int graphEdgesCount = 0;
 //Variable used to loop in the program for each of the 6 criteria's namely mean value, optimist, pessimist, double pessi, stable and own
 public static int criteriaCount = 0;
 public static int valuesCount = 0;
 //Map to store given input edges and their details
 static Map < String, ArrayList < String >> inputMap = new HashMap < String, ArrayList < String >> ();
 //Map to store given input edges and assign them to Graph for each iteration of criteria's
 static Map < String, ArrayList < String >> edgesMap = new HashMap < String, ArrayList < String >> ();
 //Map to store mean, optimis, pessimist, stable and own path values
 static Map < String, ArrayList < String >> valuesMap = new HashMap < String, ArrayList < String >> ();
 //Variable to store shortest path for Mean values criteria
 public static String shortestPathVar1="";
 //Variable to store shortest path for optimist criteria
 public static String shortestPathVar2="";
 //Variable to store shortest path for pessimist criteria		 
 public static String shortestPathVar3="";
 //Variable to store shortest path for double pessimist criteria
 public static String shortestPathVar4="";
 //Variable to store shortest path for stable criteria
 public static String shortestPathVar5="";
 //Variable to store shortest path for own criteria
 public static String shortestPathVar6="";
 

 public static void main(String[] args) throws NumberFormatException, IOException {
  acceptInputs(); //Stores input in a map
  for (int i = 0; i < 6; i++) {
   criteriaCount++;
   hopCount = 0;
   calculatevalues();  //Calculates mean, variance and cSquare for edges
   Graph.Edge[] GRAPH = new Graph.Edge[edgesMap.size()];
   int graphCount = 0;
   //Assign values to graph
   for (String key: edgesMap.keySet()) {
    GRAPH[graphCount] = new Graph.Edge(edgesMap.get(key).get(0), edgesMap.get(key).get(1), Double.parseDouble(edgesMap.get(key).get(2)));
    graphCount++;
   }  
   edgesMap.clear();
   Graph g = new Graph(GRAPH);
   //To display output
   switch(i){
   case 0:
	   System.out.println("\nMean Value Path:\nDijikstra's shortest path for the given input with least total 'expected value' is ");
	   break;
   case 1:
	   System.out.println("\nOptimist Path:\nDijikstra's shortest path for the given input with least 'expected value - standard deviation' is ");
	   break;
   case 2:
	   System.out.println("\nPessimist Path:\nDijikstra's shortest path for the given input with least 'expected value + standard deviation' is ");
	   break;
   case 3:
	   System.out.println("\nDouble Pessimist Path:\nDijikstra's shortest path for the given input with least 'expected value + 2 * standard deviation' is ");
	   break;
   case 4:
	   System.out.println("\nStable Path:\nDijikstra's shortest path for the given input with least total 'squared coefficient of variation' is ");
	   break;
   case 5:
	   System.out.println("\nMean+CSquare Path:\nDijikstra's shortest path for the given input with least 'expected value + squared coefficient' is ");
	   break;
   default:
	   break;
   }
   g.dijkstra(START.trim()); //Passing start vertex
   g.printPath(END.trim()); //Passing end vertex

   System.out.println("Number of hops for the path traversed above is " + DijikstrasAlgorithm.hopCount);
  }
  linksUsedByPaths(); //Displays edges and criterias using those edges
  System.out.println("\n \t\t \u03BC \t\t\t \u03BC-\u03C3 \t\t\t \u03BC+\u03C3 \t\t\t \u03BC+2*\u03C3 \t\t\t CSq \t\t\t Mean+CSq");
  comparePathValues();
 }
 private static void comparePathValues() {
	 String tempCompare="";
	 String shortestPathVar="";
	 float criteria1=0;
	 float criteria2=0;
	 float criteria3=0;
	 float criteria4=0;
	 float criteria5=0;
	 float criteria6=0;
 for(int i=1;i<7;i++){
 
	 switch (i) {
	case 1:
		shortestPathVar="";
		shortestPathVar=shortestPathVar1;
		System.out.println("\nMean\t");
		break;
	case 2:
		shortestPathVar="";
		shortestPathVar=shortestPathVar2;
		System.out.println("Opt\t");
		break;
	case 3:
		shortestPathVar="";
		shortestPathVar=shortestPathVar3;
		System.out.println("Pmst\t");
	break;
	case 4:
		shortestPathVar="";
		shortestPathVar=shortestPathVar4;
		System.out.println("DPess\t");
	break;
	case 5:
		shortestPathVar="";
		shortestPathVar=shortestPathVar5;
		System.out.println("Stbl\t");

	break;
	case 6:
		shortestPathVar="";
		shortestPathVar=shortestPathVar6;
		System.out.println("M+SQq\t");

	break;
	default:
		break;
	}
	 for (String key: valuesMap.keySet()) {
			tempCompare=valuesMap.get(key).get(0) + "," + valuesMap.get(key).get(1);
			if(shortestPathVar.contains(tempCompare)){
				criteria1= (float) (criteria1+ Double.parseDouble(valuesMap.get(key).get(2))) ;
				criteria2= (float) (criteria2+ Double.parseDouble(valuesMap.get(key).get(3))) ;
				criteria3= (float) (criteria3+ Double.parseDouble(valuesMap.get(key).get(4))) ;
				criteria4= (float) (criteria4+ Double.parseDouble(valuesMap.get(key).get(5))) ;
				criteria5= (float) (criteria5+ Double.parseDouble(valuesMap.get(key).get(6))) ;
				criteria6= (float) (criteria6+ Double.parseDouble(valuesMap.get(key).get(7))) ;
			}
			
		 }System.out.println("\t\t" + criteria1 + "\t\t" + criteria2+ "\t\t" + criteria3+ "\t\t" + criteria4+ "\t\t" + criteria5
					+ "\t\t" + criteria6 );
		  tempCompare="";
		  criteria1=0;
		  criteria2=0;
		  criteria3=0;
		  criteria4=0;
		  criteria5=0;
		  criteria6=0;	 
 }	
}
private static void acceptInputs() throws IOException {
  int inputCount = 0;
  System.out.println("Enter graph inputs");
  BufferedReader stdin = new BufferedReader(new InputStreamReader(System.in));
  String line;
  while ((line = stdin.readLine()) != null && line.length() != 0) {
   String[] input = line.split(",");
   // Input's first line where start and end nodes are given
   if (input.length == 3) {
    inputStartEndCount++;
    if (inputStartEndCount == 1) {
     START = input[1];
     END = input[2];					    	
    } else {
     System.err.printf("\n Start and end nodes are given more than once. Please verify your input");
     System.exit(1);
    }
    // Input's other lines where edges with type, alpha, beta are provided
   } else if (input.length == 6) {
    char inputFirstElement = input[0].charAt(0);
    if (inputFirstElement == 'E') {
     node1 = input[1];
     node2 = input[2];
     type = Integer.parseInt(input[3]);
     alpha = Double.parseDouble(input[4]);
     beta = Double.parseDouble(input[5]);
     inputCount++;
     ArrayList < String > inputEdges = new ArrayList < String > ();
     inputEdges.add(node1);
     inputEdges.add(node2);
     inputEdges.add(Integer.toString(type));
     inputEdges.add(Double.toString(alpha));
     inputEdges.add(Double.toString(beta));
     inputMap.put(Integer.toString(inputCount), inputEdges);
    } else {
     System.err.printf("\n Invalid input graph edge. First letter of input has to be 'E' indicating Edge");
     System.exit(1);
    }
   } else {
    System.err.printf("\n Improper graph input. Please verify and try again.");
    System.exit(1);
   }
  }

 }

 private static void linksUsedByPaths() {
	 String inputEdges="";
	 String criteria1="";
	 String criteria2="";
	 String criteria3="";
	 String criteria4="";
	 String criteria5="";
	 String criteria6="";
	 System.out.println("\nLinks used by paths:");
	 System.out.println("\nEdges \t\t MV \t Op \t Ps \t DP \t St \t Own");
	 for (String key: inputMap.keySet()) {
		 inputEdges="";
		 inputEdges=inputMap.get(key).get(0) + "," + inputMap.get(key).get(1);
		 criteria1="";
		 criteria2="";
		 criteria3="";
		 criteria4="";
		 criteria5="";
		 criteria6="";
		 if(shortestPathVar1.contains(inputEdges)){
			 criteria1="*"; 
		 }
		 if(shortestPathVar2.contains(inputEdges)){
			 criteria2="*"; 
		 }
		 if(shortestPathVar3.contains(inputEdges)){
			 criteria3="*"; 
		 }
		 if(shortestPathVar4.contains(inputEdges)){
			 criteria4="*"; 
		 }
		 if(shortestPathVar5.contains(inputEdges)){
			 criteria5="*"; 
		 }
		 if(shortestPathVar6.contains(inputEdges)){
			 criteria6="*"; 
		 }
		 System.out.println("Edge(" + inputMap.get(key).get(0) + "," + inputMap.get(key).get(1) + ")"
				 + "\t" + criteria1
				 + "\t" + criteria2
				 + "\t" + criteria3
				 + "\t" + criteria4
				 + "\t" + criteria5
				 + "\t" + criteria6);
	}
 }

 private static void calculatevalues() throws NumberFormatException, IOException {
  if (edgesMap.isEmpty()) {
   graphEdgesCount = 0;
   for (Entry < String, ArrayList < String >> entry: inputMap.entrySet()) {
    ArrayList < String > value = entry.getValue();
    double edgeWeight = 0;
    double standDeviation=0;
    int inputType = Integer.parseInt(value.get(2));
    switch (inputType) {
     case 1:
      node1 = value.get(0);
      node2 = value.get(1);
      alpha = Double.parseDouble(value.get(3));
      beta = Double.parseDouble(value.get(4));
      mean = alpha;
      variance = 0;
      cSquare = 0;
      standDeviation= (double) Math.sqrt(variance);
      edgeWeight=criteriaBasedEdgeWeightsCal(node1, node2, criteriaCount,mean,cSquare,standDeviation);
      
      ArrayList < String > alEdges1 = new ArrayList < String > ();
      alEdges1.add(node1);
      alEdges1.add(node2);
      alEdges1.add(Double.toString(edgeWeight));
      graphEdgesCount++;
      edgesMap.put(Integer.toString(graphEdgesCount), alEdges1);
      break;
     case 2:
      node1 = value.get(0);
      node2 = value.get(1);
      alpha = Double.parseDouble(value.get(3));
      beta = Double.parseDouble(value.get(4));
      mean = (alpha + beta) / 2;
      variance = (double)(Math.pow((beta - alpha), 2) / 12);
      cSquare = (double)(0.333333 * Math.pow((beta - alpha) / (alpha + beta), 2));
      standDeviation= (double) Math.sqrt(variance);
      edgeWeight=criteriaBasedEdgeWeightsCal(node1, node2, criteriaCount,mean,cSquare,standDeviation);
     
      ArrayList < String > alEdges2 = new ArrayList < String > ();
      alEdges2.add(node1);
      alEdges2.add(node2);
      alEdges2.add(Double.toString(edgeWeight));
      graphEdgesCount++;
      edgesMap.put(Integer.toString(graphEdgesCount), alEdges2);
      break;
     case 3:
      node1 = value.get(0);
      node2 = value.get(1);
      alpha = Double.parseDouble(value.get(3));
      beta = Double.parseDouble(value.get(4));
      mean = 1 / alpha;
      variance = (double)(1 / Math.pow(alpha, 2));
      cSquare = 1;
      standDeviation= (double) Math.sqrt(variance);
      edgeWeight=criteriaBasedEdgeWeightsCal(node1, node2, criteriaCount,mean,cSquare,standDeviation);
      
      ArrayList < String > alEdges3 = new ArrayList < String > ();
      alEdges3.add(node1);
      alEdges3.add(node2);
      alEdges3.add(Double.toString(edgeWeight));
      graphEdgesCount++;
      edgesMap.put(Integer.toString(graphEdgesCount), alEdges3);
      break;
     case 4:
      node1 = value.get(0);
      node2 = value.get(1);
      alpha = Double.parseDouble(value.get(3));
      beta = Double.parseDouble(value.get(4));
      mean = beta + (1 / alpha);
      variance = (double)(1 / Math.pow(alpha, 2));
      cSquare = (double) Math.pow(1 / (1 + (beta * alpha)), 2);
      standDeviation= (double) Math.sqrt(variance);
      edgeWeight=criteriaBasedEdgeWeightsCal(node1, node2, criteriaCount,mean,cSquare,standDeviation);
      
      ArrayList < String > alEdges4 = new ArrayList < String > ();
      alEdges4.add(node1);
      alEdges4.add(node2);
      alEdges4.add(Double.toString(edgeWeight));
      graphEdgesCount++;
      edgesMap.put(Integer.toString(graphEdgesCount), alEdges4);
      break;
     case 5:
      node1 = value.get(0);
      node2 = value.get(1);
      alpha = Double.parseDouble(value.get(3));
      beta = Double.parseDouble(value.get(4));
      mean = alpha;
      variance = beta;
      cSquare = (double)(beta / Math.pow(alpha, 2));
      standDeviation= (double) Math.sqrt(variance);
      edgeWeight=criteriaBasedEdgeWeightsCal(node1, node2, criteriaCount,mean,cSquare,standDeviation);
      
      ArrayList < String > alEdges5 = new ArrayList < String > ();
      alEdges5.add(node1);
      alEdges5.add(node2);
      alEdges5.add(Double.toString(edgeWeight));
      graphEdgesCount++;
      edgesMap.put(Integer.toString(graphEdgesCount), alEdges5);
      break;
     case 6:
      node1 = value.get(0);
      node2 = value.get(1);
      alpha = Double.parseDouble(value.get(3));
      beta = Double.parseDouble(value.get(4));
      mean = alpha;
      variance = (double)((Math.pow(alpha, 2)) * beta);
      cSquare = beta;
      standDeviation= (double) Math.sqrt(variance);
      edgeWeight=criteriaBasedEdgeWeightsCal(node1, node2, criteriaCount,mean,cSquare,standDeviation);
      
      ArrayList < String > alEdges6 = new ArrayList < String > ();
      
      alEdges6.add(node1);
      alEdges6.add(node2);
      alEdges6.add(Double.toString(edgeWeight));
      graphEdgesCount++;
      edgesMap.put(Integer.toString(graphEdgesCount), alEdges6);
      break;
     default:
      break;
    }
   }
  }
 }
private static double criteriaBasedEdgeWeightsCal(String node1, String node2, int criteriaCount, double mean, double cSquare2, double standDeviation) {
	double tempEdgeWeight=0;
 switch (criteriaCount) {
 case 1:
	  tempEdgeWeight = mean;
    break;
   case 2:
	   tempEdgeWeight = mean-standDeviation;
    break;
   case 3:
	   tempEdgeWeight = mean + standDeviation;
    break;
   case 4:
	   tempEdgeWeight = mean + 2 * standDeviation;
    break;
   case 5:
	   tempEdgeWeight = cSquare;
    break;
   case 6:
	   tempEdgeWeight = mean + cSquare;
    break;
   default:
    break;
  }
 double tempValueStore=0;
 if(DijikstrasAlgorithm.criteriaCount==6){
	 ArrayList < String > alValues = new ArrayList < String > ();
     valuesCount++; 
	 alValues.add(node1);
	 alValues.add(node2);
	 tempValueStore=mean;
	 alValues.add(Double.toString(tempValueStore));
	 tempValueStore=mean-standDeviation;
	 alValues.add(Double.toString(tempValueStore));
	 tempValueStore=mean+standDeviation;
	 alValues.add(Double.toString(tempValueStore));
	 tempValueStore=mean+ 2* standDeviation;
	 alValues.add(Double.toString(tempValueStore));
	 tempValueStore= cSquare ;
	 alValues.add(Double.toString(tempValueStore));
	 tempValueStore=mean + cSquare ;
	 alValues.add(Double.toString(tempValueStore));
     valuesMap.put(Integer.toString(valuesCount), alValues); 
 } 
return tempEdgeWeight;
}
}
class Graph {
 private static Map < String, Vertex > graph; // mapping of vertex names to Vertex objects, built from a set of Edges

 /** One edge of the graph (only used by Graph constructor) */
 public static class Edge {
  public final String v1, v2;
  public final double dist;
  public Edge(String v1, String v2, double d) {
   this.v1 = v1;
   this.v2 = v2;
   this.dist = d;
  }
 }

 /** One vertex of the graph, complete with mappings to neighbouring vertices */
 public static class Vertex implements Comparable < Vertex > {
  public final String name;
  public double dist = Double.MAX_VALUE; // MAX_VALUE assumed to be infinity
  public Vertex previous = null;
  public final Map < Vertex,
  Double > neighbours = new HashMap < > ();

  public Vertex(String name) {
   this.name = name;
  }

  private void printPath() {
   if (this == this.previous) {
	   switch(DijikstrasAlgorithm.criteriaCount){
	   case 1:
		   DijikstrasAlgorithm.shortestPathVar1 = DijikstrasAlgorithm.shortestPathVar1.concat(this.name).concat(",");
		   break;
	   case 2:
		   DijikstrasAlgorithm.shortestPathVar2 = DijikstrasAlgorithm.shortestPathVar2.concat(this.name).concat(",");
		   break;
	   case 3:
		   DijikstrasAlgorithm.shortestPathVar3 = DijikstrasAlgorithm.shortestPathVar3.concat(this.name).concat(",");
		   break;
	   case 4:
		   DijikstrasAlgorithm.shortestPathVar4 = DijikstrasAlgorithm.shortestPathVar4.concat(this.name).concat(",");
		   break;
	   case 5:
		   DijikstrasAlgorithm.shortestPathVar5 = DijikstrasAlgorithm.shortestPathVar5.concat(this.name).concat(",");
		   break;
	   case 6:
		   DijikstrasAlgorithm.shortestPathVar6 = DijikstrasAlgorithm.shortestPathVar6.concat(this.name).concat(",");
		   break;
	   default:
		   break;		   
	   }
    System.out.printf("%s", this.name);
   } else if (this.previous == null) {
    System.out.printf("%s(unreached)", this.name);
   } else {
    this.previous.printPath();
    DijikstrasAlgorithm.hopCount++;
    switch(DijikstrasAlgorithm.criteriaCount){
	   case 1:
		   DijikstrasAlgorithm.shortestPathVar1 = DijikstrasAlgorithm.shortestPathVar1.concat(this.name).concat(",");
		   break;
	   case 2:
		   DijikstrasAlgorithm.shortestPathVar2 = DijikstrasAlgorithm.shortestPathVar2.concat(this.name).concat(",");
		   break;
	   case 3:
		   DijikstrasAlgorithm.shortestPathVar3 = DijikstrasAlgorithm.shortestPathVar3.concat(this.name).concat(",");
		   break;
	   case 4:
		   DijikstrasAlgorithm.shortestPathVar4 = DijikstrasAlgorithm.shortestPathVar4.concat(this.name).concat(",");
		   break;
	   case 5:
		   DijikstrasAlgorithm.shortestPathVar5 = DijikstrasAlgorithm.shortestPathVar5.concat(this.name).concat(",");
		   break;
	   case 6:
		   DijikstrasAlgorithm.shortestPathVar6 = DijikstrasAlgorithm.shortestPathVar6.concat(this.name).concat(",");
		   break;
	   default:
		   break;		   
	   }
    System.out.printf(" -> %s(%f)", this.name, this.dist);
   }

  }

  public int compareTo(Vertex other) {
   if (dist == other.dist)
    return name.compareTo(other.name);

   return Double.compare(dist, other.dist);
  }

  @Override public String toString() {
   return "(" + name + ", " + dist + ")";
  }
 }

 /** Builds a graph from a set of edges */
 public Graph(Edge[] edges) {
  graph = new HashMap < > (edges.length);

  //one pass to find all vertices
  for (Edge e: edges) {
   if (!graph.containsKey(e.v1)) graph.put(e.v1, new Vertex(e.v1));
   if (!graph.containsKey(e.v2)) graph.put(e.v2, new Vertex(e.v2));
  }

  //another pass to set neighbouring vertices
  for (Edge e: edges) {
	  //graph.get(e.v1).neighbours.put(graph.get(e.v2), e.dist);
   graph.get(e.v1).neighbours.put(graph.get(e.v2), e.dist);
   graph.get(e.v2).neighbours.put(graph.get(e.v1), e.dist); // for undirected graph
  }
 }

 /** Runs dijkstra using a specified source vertex */
 public void dijkstra(String startName) {
  if (!graph.containsKey(startName)) {
   System.err.printf("\n Graph doesn't contain start vertex " + startName);
   System.exit(1);
  }
  final Vertex source = graph.get(startName);
  NavigableSet < Vertex > q = new TreeSet < > ();

  // set-up vertices
  for (Vertex v: graph.values()) {
   v.previous = v == source ? source : null;
   v.dist = v == source ? 0 : Double.MAX_VALUE;
   q.add(v);
  }

  dijkstra(q);
 }

 /** Implementation of dijkstra's algorithm using a binary heap. */
 private void dijkstra(final NavigableSet < Vertex > q) {
  Vertex u, v;
  while (!q.isEmpty()) {

   u = q.pollFirst(); // vertex with shortest distance (first iteration will return source)
   if (u.dist == Double.MAX_VALUE) break; // we can ignore u (and any other remaining vertices) since they are unreachable

   //look at distances to each neighbour
   for (Map.Entry < Vertex, Double > a: u.neighbours.entrySet()) {
    v = a.getKey(); //the neighbour in this iteration

    final double alternateDist = u.dist + a.getValue();
    if (alternateDist < v.dist) { // shorter path to neighbour found
     q.remove(v);
     v.dist = alternateDist;
     v.previous = u;
     q.add(v);
    }
   }
  }
 }

 /** Prints a path from the source to the specified vertex */
 public void printPath(String endName) {
   if (!graph.containsKey(endName)) {
    System.err.printf("\n Graph doesn't contain end vertex " + endName);
    System.exit(1);
   }

   graph.get(endName).printPath();
   System.out.println();
  }
  /** Prints the path from the source to every vertex (output order is not guaranteed) */
 public static void printAllPaths() {
  for (Vertex v: graph.values()) {
   v.printPath();
   System.out.println();
  }
 }


}