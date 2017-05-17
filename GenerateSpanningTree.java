import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;
import javax.imageio.ImageIO;
import javax.swing.JApplet;
import javax.swing.JFrame;

/*
 * This class implements Comparator class to manually override compare function.
 * Compare method takes two points and uses orientation method to compute their orientation. 
 */
class PrivateComparator implements Comparator<Point> {

	public int compare(Point p1, Point p2) {

		int o = GenerateSpanningTree.orientation(GenerateSpanningTree.p0, p1, p2);
		if (o == 0)
			return (GenerateSpanningTree.distSq(GenerateSpanningTree.p0, p2) >= GenerateSpanningTree
					.distSq(GenerateSpanningTree.p0, p1)) ? -1 : 1;

		return (o == 2) ? -1 : 1;
	}

}

/*
 * This is the main class in which the whole process takes place. It takes
 * inputs as a form of vertices and their points and compute all crossing free
 * spanning trees. After calculating total number of crossing free spanning
 * trees, it generates a crossing free spanning tree from a novel algorithm and
 * inputs this tree to a novel markov chain algorithm which randomly chooses an
 * edge and generates new crossing free spanning tree.
 * 
 * Markov Chain algorithm is used as a procedure to find the number of
 * iterations after which we are able to reach all the crossing free spanning
 * trees.
 * 
 * In this, I will be writting Crossing free Spanning Trees as CFST.
 */
public class GenerateSpanningTree extends JApplet {

	// Applet colors
	final static Color bg = Color.white;
	final static Color fg = Color.black;
	final static Color red = Color.red;
	final static Color white = Color.white;
	// Basic strokes which changes the width of line based on the value used.
	final static BasicStroke stroke = new BasicStroke(2.0f);
	final static BasicStroke wideStroke = new BasicStroke(80.0f);

	final static float dash1[] = { 10.0f };
	final static BasicStroke dashed = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, dash1,
			0.0f);
	// This hashmap is used to store vertices and their points.
	static Map<String, Point> hmapForPoints = new HashMap<String, Point>();
	// This arraylist is used to store the lines for the CFST generated from
	// novel algorithm.
	static List<Line> listOfSpanningTreeLines = new ArrayList<Line>();
	// This structure(ArrayList of Lists) is used to store all CFSTs while
	// enumerating through them.
	List<List<String>> resultantCrossingFreeSpanningTrees = new ArrayList<List<String>>();
	// This HashMap keeps count of all CFSTs.
	static Map<List<String>, Integer> hmapForMarkovChainCount = new HashMap<List<String>, Integer>();
	static Point p0;
	static int couponCollector;
	static double MConstant = 0.5772156649;
	static int count = 0;
	Graphics g;

	/*
	 * This method calculates the distance between two points.
	 */
	static int distSq(Point p1, Point p2) {
		return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
	}

	/*
	 * This method calculates the orientation of one line combining p and q with
	 * respect to r.
	 */
	static int orientation(Point p, Point q, Point r) {
		int val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);

		if (val == 0)
			return 0; // colinear
		return (val > 0) ? 1 : 2; // clock or counterclock wise
	}

	/*
	 * This method calculates if all points are collinear or not. If yes, it
	 * will return TrueA if they are parallel to y-axis. TrueB if they are
	 * parallel to x-axis and TrueC "HighPoint" "LowPoint" if they are at an
	 * angle and HighPoint is the topmost point and LowPoint is the lowest
	 * point.
	 */
	private static String pointsCollinear() {

		Set<String> hset = hmapForPoints.keySet();
		Iterator it = hset.iterator();
		String tempPoint = (String) it.next();

		// All x-coordinates are same. (Parallel to y-axis)
		boolean isCollinear = true;
		for (Map.Entry<String, Point> h : hmapForPoints.entrySet()) {
			if (h.getValue().x != hmapForPoints.get(tempPoint).x) {
				isCollinear = false;
				break;
			}
		}
		if (isCollinear)
			return "TrueA";

		// All y-coordinates are same. (Parallel to x-axis)
		isCollinear = true;
		for (Map.Entry<String, Point> h : hmapForPoints.entrySet()) {
			if (h.getValue().y != hmapForPoints.get(tempPoint).y) {
				isCollinear = false;
				break;
			}
		}
		if (isCollinear)
			return "TrueB";

		// All points are at an angle.
		int low = Integer.MAX_VALUE, high = Integer.MIN_VALUE;
		String lowPoint = "", highPoint = "";
		for (Map.Entry<String, Point> h : hmapForPoints.entrySet()) {
			Point p = h.getValue();
			if (p.x > high) {
				high = p.x;
				highPoint = h.getKey();
			}
			if (p.x < low) {
				low = p.x;
				lowPoint = h.getKey();
			}
		}

		// Generate a slope from low and high points.
		// see if all the points have same slope
		// if yes return TrueC
		float slope;
		Point h = hmapForPoints.get(highPoint);
		Point l = hmapForPoints.get(lowPoint);

		slope = (h.y - l.y) / (h.x - l.x);

		for (Map.Entry<String, Point> hm : hmapForPoints.entrySet()) {
			if (hm.getKey().equals(lowPoint) || hm.getKey().equals(highPoint))
				continue;
			if (slope != (hm.getValue().y - l.y) / (hm.getValue().x - l.x))
				return "False";
		}
		if (!isCollinear)
			return "TrueC " + highPoint + " " + lowPoint;
		return "False";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.Container#paint(java.awt.Graphics) This function prints
	 * Lines and points for each JFrame.
	 */
	public void paint(Graphics g) {

		Graphics2D g2 = (Graphics2D) g;
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g2.setStroke(new BasicStroke(1.4f));
		// public Line2D.Double(double x1,double y1,double x2,double y2)
		// x1 - the X coordinate of the start point
		// y1 - the Y coordinate of the start point
		// x2 - the X coordinate of the end point
		// y2 - the Y coordinate of the end point

		// Point x,y
		// x=x-r/2, y=y-r/2
		for (Map.Entry<String, Point> h : hmapForPoints.entrySet()) {
			g2.drawLine(h.getValue().x, h.getValue().y, h.getValue().x, h.getValue().y);
			g2.drawOval(h.getValue().x - 5, h.getValue().y - 5, 10, 10);
		}

		List<String> list = resultantCrossingFreeSpanningTrees.get(count);
		for (int i = 0; i < list.size(); i++) {
			g2.drawLine(hmapForPoints.get(list.get(i).charAt(0) + "").x,
					hmapForPoints.get(list.get(i).charAt(0) + "").y, hmapForPoints.get(list.get(i).charAt(1) + "").x,
					hmapForPoints.get(list.get(i).charAt(1) + "").y);
		}

	}

	/*
	 * This method enumerates all CFSTs from a graph. The algorithm is written
	 * inside the function.
	 */
	private void generateAllSpanningTrees() {

		// This method generates all the spanning trees using divide and conquer
		// algorithm.
		// Algorithm:
		// 1: Choose a random vertex and call it vref. Remove the pendant edges
		// in the graph starting with vref.
		// Since we are dealing with complete graph, choose any vertex as every
		// one will be same.
		// Also, there will be no pendant edge in complete graph.
		// 2: Generate all the combinations of edges from vref. Also, for each
		// combination,
		// generate secondary components.
		// 3: Calculate the secondary connectors and main connectors for each
		// combination.
		// 4: Create trees based on these connectors and partitions.

		// Step 0: Construct a complete graph containing the vertices all
		// connected with each other.
		Map<String, List<String>> hmapOfCompleteGraph = new HashMap<String, List<String>>();
		hmapOfCompleteGraph = generateCompleteGraph();

		// Step 1:
		Iterator it = hmapForPoints.keySet().iterator();
		String vref = (String) (it.next());

		// Step 2:
		Set<String> hset = new HashSet<String>(hmapForPoints.keySet());
		hset.remove(vref);
		Set<String> numberOfEdges = generateEdgesForVref(vref, hset);
		// remove edges containing vref from completegraph
		hmapOfCompleteGraph = removeVrefEdges(hmapOfCompleteGraph, vref);
		Map<Set<String>, List<SecondarypartitionClass>> hmapOfSecondaryPartition = new HashMap<Set<String>, List<SecondarypartitionClass>>();
		Set<Set<String>> hsetForCombinations = generateCombinations(numberOfEdges);
		Iterator iterate = hsetForCombinations.iterator();
		while (iterate.hasNext()) {
			Set<String> tempSet = (Set<String>) iterate.next();
			hmapOfSecondaryPartition.put(tempSet, generateSecondaryPartition(tempSet));
		}

		// Step 3:
		// Generate secondary connectors.
		hmapOfSecondaryPartition = generateSecondaryConnectors(hmapOfSecondaryPartition);

		// Generate Main Connectors
		Map<Set<String>, List<String>> hmapForMainConnectors = new HashMap<Set<String>, List<String>>();
		hmapForMainConnectors = generateMainConnectors(vref, hmapOfSecondaryPartition);

		// Step 4: Generate trees
		// a: Connect Primary partition + Secondary Partition + secondary
		// connector + main connector
		List<List<String>> resultantSpanningTrees = new ArrayList<List<String>>();
		resultantSpanningTrees = generateTrees1(hmapOfSecondaryPartition, hmapForMainConnectors);

		// b:
		// For every primary partition, we know the main connectors. We have to
		// make sets of 2,3.... n from main connectors.
		// Add primary partition in each set and check if no of edges in the set
		// is less than or equal to n-1.
		// Compute unvisited vertices for each and generate secondary
		// connectors, main connectors for every set and combine trees as part
		// a.
		generateTreesB(hmapForMainConnectors, resultantSpanningTrees, vref);

		resultantCrossingFreeSpanningTrees = checkCrossingFreeForAllSpanningTrees(resultantSpanningTrees);
		rearrangeSpanningTrees();
		resultantCrossingFreeSpanningTrees = new ArrayList<List<String>>(hmapForMarkovChainCount.keySet());
	}

	/*
	 * This method checks the crossing nature of each spanning tree. It removes
	 * all the spanning trees which are crossing.
	 */
	private static List<List<String>> checkCrossingFreeForAllSpanningTrees(List<List<String>> resultantSpanningTrees) {

		List<List<String>> resultList = new ArrayList<List<String>>();
		int count = 0;
		for (int i = 0; i < resultantSpanningTrees.size(); i++) {
			List<String> tempList = resultantSpanningTrees.get(i);
			if (checkSpanningTree(tempList)) {
				boolean intersect = false;
				// Check if crossing free or not if yes then add it to
				// resultList
				for (int j = 0; j < tempList.size() - 1; j++) {
					String str1 = tempList.get(j);
					String u = str1.charAt(0) + "";
					String v = str1.charAt(1) + "";
					Point p1 = hmapForPoints.get(u);
					Point q1 = hmapForPoints.get(v);
					for (int k = j + 1; k < tempList.size(); k++) {
						String str2 = tempList.get(k);
						String w = str2.charAt(0) + "";
						String x = str2.charAt(1) + "";
						Point p2 = hmapForPoints.get(w);
						Point q2 = hmapForPoints.get(x);
						if (u.equals(w) || u.equals(x) || v.equals(w) || v.equals(x)) {

						} else {
							if (checkIfIntersect(p1, q1, p2, q2))// return true;
							{
								intersect = true;
							}
						}
						if (intersect)
							break;

					}
					if (intersect)
						break;
				}
				if (!intersect)
					resultList.add(tempList);
			}

		}

		return resultList;
	}

	/*
	 * This method will check if the tree present in the ArrayList is Spanning
	 * tree or not. It checks two properties of Spanning trees: 1: The number of
	 * edges should be n-1 2: Every vertex must be connected with each other.
	 */
	private static boolean checkSpanningTree(List<String> tempList) {

		// If the no of edges are not n-1.
		if (tempList.size() != hmapForPoints.size() - 1)
			return false;

		// Generate undirected hmap to counter the situation where one edge
		// couldn't be reached
		// because it was stored in the data structure in the opposite order
		// such as for edge AB,
		// if the data structure contains A -> B but we are only trying to
		// access through B.
		Map<String, List<String>> hmap = new HashMap<String, List<String>>();
		for (int i = 0; i < tempList.size(); i++) {
			String str = tempList.get(i);
			String firstV = str.charAt(0) + "";
			String secondV = str.charAt(1) + "";
			if (hmap.containsKey(firstV)) {
				List<String> temp = hmap.get(firstV);
				temp.add(secondV);
				hmap.put(firstV, temp);
			} else {
				List<String> temp = new ArrayList<String>();
				temp.add(secondV);
				hmap.put(firstV, temp);
			}
			if (hmap.containsKey(secondV)) {
				List<String> temp = hmap.get(secondV);
				temp.add(firstV);
				hmap.put(secondV, temp);
			} else {
				List<String> temp = new ArrayList<String>();
				temp.add(firstV);
				hmap.put(secondV, temp);
			}
		}

		// Breadth First Implementation
		Set<String> hset = new HashSet<String>();
		Iterator it = hmap.keySet().iterator();
		String root = (String) it.next();
		Queue<String> queue = new LinkedList<String>();
		queue.add(root);
		hset.add(root);
		while (!queue.isEmpty()) {
			String tempString = queue.poll();
			List<String> t = hmap.get(tempString);
			for (int i = 0; i < t.size(); i++) {
				if (!hset.contains(t.get(i))) {
					queue.add(t.get(i));
					hset.add(t.get(i));
				}
			}
		}

		if (hset.size() != hmapForPoints.size())
			return false;
		else
			return true;

	}

	/*
	 * This method is the phase B of Step 4 enumeration algorithm.
	 */
	private void generateTreesB(Map<Set<String>, List<String>> hmapForMainConnectors,
			List<List<String>> resultantSpanningTrees, String vref) {

		int count = 0;
		Set<Set<String>> tempHset = new HashSet<Set<String>>();
		Map<Set<String>, Set<String>> hmapContainingNewPrimaryAndPrimary = new HashMap<Set<String>, Set<String>>();
		for (Map.Entry<Set<String>, List<String>> h : hmapForMainConnectors.entrySet()) {

			Set<String> tempSet = h.getKey();
			List<String> tempValue = h.getValue();
			if (tempValue.size() != 0 && hmapForPoints.size() - tempSet.size() > 2) {

				Set<String> numberOfEdges = new HashSet<String>(tempValue);
				Set<Set<String>> hsetForCombinations = generateCombinations(numberOfEdges);
				hsetForCombinations = removeSingleCombinations(hsetForCombinations);
				hsetForCombinations = addPrimaryPartition(hsetForCombinations, tempSet);
				count += hsetForCombinations.size();

				Iterator it = hsetForCombinations.iterator();
				while (it.hasNext()) {

					Set<String> tSet = (Set<String>) it.next();
					if (tSet.size() == hmapForPoints.size() - 1) {
						List<String> newList = new ArrayList<String>(tSet);
						resultantSpanningTrees.add(newList);
					} else if (tSet.size() < hmapForPoints.size() - 1) {
						tempHset.add(tSet);
						hmapContainingNewPrimaryAndPrimary.put(tSet, tempSet);
					}
				}
			}
		}

		// Now tempHset contains primary partitions. All the partitions which
		// are of size 5 are added to the resultant.
		Map<Set<String>, List<SecondarypartitionClass>> hmapOfSecondaryPartition = new HashMap<Set<String>, List<SecondarypartitionClass>>();
		Iterator iterate = tempHset.iterator();
		while (iterate.hasNext()) {
			Set<String> tempSet = (Set<String>) iterate.next();
			hmapOfSecondaryPartition.put(tempSet, generateSecondaryPartitionB(tempSet, vref));
		}
		hmapOfSecondaryPartition = generateSecondaryConnectors(hmapOfSecondaryPartition);

		Map<Set<String>, List<String>> hmapForMainConnectorsPartB = new HashMap<Set<String>, List<String>>();
		hmapForMainConnectorsPartB = generateMainConnectors(vref, hmapOfSecondaryPartition);

		// Remove the main connectors which are already there.
		hmapForMainConnectorsPartB = filterPreviousMainConnectors(hmapContainingNewPrimaryAndPrimary,
				hmapForMainConnectorsPartB, hmapForMainConnectors);
		List<List<String>> resultantSpanningTreesPartB = new ArrayList<List<String>>();
		resultantSpanningTreesPartB = generateTrees1(hmapOfSecondaryPartition, hmapForMainConnectorsPartB);

		for (int i = 0; i < resultantSpanningTreesPartB.size(); i++) {
			List<String> l = resultantSpanningTreesPartB.get(i);
			if (l.size() == hmapForPoints.size() - 1)
				resultantSpanningTrees.add(l);
		}

	}

	/*
	 * This is a subfunction which removes the previous main connectors used.
	 */
	private Map<Set<String>, List<String>> filterPreviousMainConnectors(
			Map<Set<String>, Set<String>> hmapContainingNewPrimaryAndPrimary,
			Map<Set<String>, List<String>> hmapForMainConnectorsPartB,
			Map<Set<String>, List<String>> hmapForMainConnectors) {

		Map<Set<String>, List<String>> tempCopy = new HashMap<Set<String>, List<String>>(hmapForMainConnectorsPartB);
		for (Map.Entry<Set<String>, List<String>> h : tempCopy.entrySet()) {

			Set<String> key = h.getKey();
			List<String> value = new ArrayList<String>(h.getValue());

			Set<String> parentPrimaryPartition = hmapContainingNewPrimaryAndPrimary.get(key);
			List<String> primaryMainConnectorList = hmapForMainConnectors.get(parentPrimaryPartition);

			for (int i = 0; i < primaryMainConnectorList.size(); i++) {
				String temp = primaryMainConnectorList.get(i);
				value.remove(temp);
			}
			hmapForMainConnectorsPartB.put(key, value);
		}
		return hmapForMainConnectorsPartB;
	}

	/*
	 * This function generates secondary partitions for the phase B of Step 4
	 * enumeration algorithm.
	 */
	private List<SecondarypartitionClass> generateSecondaryPartitionB(Set<String> tempSet, String vref) {

		Iterator it = tempSet.iterator();
		List<String> unVisitedVertexlist = new ArrayList<String>(hmapForPoints.keySet());
		List<SecondarypartitionClass> l = new ArrayList<SecondarypartitionClass>();
		while (it.hasNext()) {
			String temp = (String) it.next();
			for (int i = 0; i < temp.length(); i++) {
				String t = temp.charAt(i) + "";
				if (unVisitedVertexlist.contains(t))
					unVisitedVertexlist.remove(t);
			}
		}

		if (unVisitedVertexlist.size() == 1) {
			l.add(new SecondarypartitionClass(null, unVisitedVertexlist.get(0)));
		} else if (unVisitedVertexlist.size() == 0) {
			l.add(new SecondarypartitionClass(null, ""));
		} else if (unVisitedVertexlist.size() % 2 == 0) {
			// unvisited vertex will not be there.
			// Wrote a helper function which will generate an edge of two.
			List<List<String>> unvisitedEdges = generateSecondaryEdges(unVisitedVertexlist);
			for (int i = 0; i < unvisitedEdges.size(); i++) {
				l.add(new SecondarypartitionClass(unvisitedEdges.get(i), ""));
			}
		} else {
			// unvisited vertex will be there.
			for (int i = 0; i < unVisitedVertexlist.size(); i++) {
				String unvisitedVert = unVisitedVertexlist.get(i);
				List<String> newList = new ArrayList<String>(unVisitedVertexlist);
				newList.remove(unvisitedVert);
				List<List<String>> unvisitedEdges = generateSecondaryEdges(newList);
				for (int j = 0; j < unvisitedEdges.size(); j++) {
					l.add(new SecondarypartitionClass(unvisitedEdges.get(j), unvisitedVert));
				}
			}
		}
		return l;

	}

	/*
	 * This is a helper function for phase B while generating primary
	 * partitions. It adds primary partition to the combinations of main
	 * connectors.
	 */
	private Set<Set<String>> addPrimaryPartition(Set<Set<String>> hsetForCombinations, Set<String> tempSet) {

		Set<Set<String>> temp = new HashSet<Set<String>>();
		Iterator it = hsetForCombinations.iterator();
		int targetSize = hmapForPoints.size() - tempSet.size();
		while (it.hasNext()) {
			Set<String> t = (Set<String>) it.next();
			if (t.size() < targetSize) {
				Set<String> t1 = new HashSet<String>(t);
				t1.addAll(tempSet);
				temp.add(t1);
			}

		}

		return temp;
	}

	/*
	 * This is a helper function for phase B while generating primary
	 * partitions. It removes all the combinations which are of size 1.
	 */
	private Set<Set<String>> removeSingleCombinations(Set<Set<String>> hsetForCombinations) {

		Set<Set<String>> temp = new HashSet<Set<String>>();
		Iterator it = hsetForCombinations.iterator();
		while (it.hasNext()) {
			Set<String> t = (Set<String>) it.next();
			if (t.size() >= 2)
				temp.add(t);
		}

		return temp;
	}

	/*
	 * This is a helper function for phase B while generating primary
	 * partitions. It joins components together to form trees.
	 */
	private List<List<String>> generateTrees1(Map<Set<String>, List<SecondarypartitionClass>> hmapOfSecondaryPartition,
			Map<Set<String>, List<String>> hmapForMainConnectors) {

		List<List<String>> resultantSpanningTrees = new ArrayList<List<String>>();
		for (Map.Entry<Set<String>, List<SecondarypartitionClass>> h : hmapOfSecondaryPartition.entrySet()) {
			Set<String> tempSet = h.getKey();
			List<String> setToList = new ArrayList<String>(tempSet);
			List<SecondarypartitionClass> tempValue = h.getValue();
			for (int k = 0; k < tempValue.size(); k++) {
				SecondarypartitionClass temp = tempValue.get(k);
				List<List<String>> tempList = addSeondaryConnectorWithSecondaryPartition(temp.unvisitedEdges,
						temp.listOfSecondaryConnectors);

				tempList = addSeondaryConnectorWithSecondaryPartition(setToList, tempList);

				tempList = addAllWithMainConnectors(tempList, hmapForMainConnectors.get(tempSet));
				resultantSpanningTrees.addAll(tempList);
			}
		}

		return resultantSpanningTrees;
	}

	/*
	 * This is also a helper function which adds all the combinations with the
	 * primary partitions.
	 */
	private List<List<String>> addAllWithMainConnectors(List<List<String>> tempList, List<String> list) {

		List<List<String>> resultList = new ArrayList<List<String>>();
		if (list.size() == 0)
			return tempList;
		for (int i = 0; i < tempList.size(); i++) {
			for (int j = 0; j < list.size(); j++) {
				List<String> list1 = new ArrayList<String>(tempList.get(i));
				String str = list.get(j);
				list1.add(str);
				resultList.add(list1);
			}
		}

		return resultList;
	}

	/*
	 * This is also a helper function which adds all the Secondary partitions
	 * with the secondary connectors.
	 */
	private List<List<String>> addSeondaryConnectorWithSecondaryPartition(List<String> unvisitedEdges,
			List<List<String>> listOfSecondaryConnectors) {

		List<List<String>> tempList = new ArrayList<List<String>>();
		if (unvisitedEdges == null) {
			return listOfSecondaryConnectors;
		}
		if (listOfSecondaryConnectors == null) {
			tempList.add(unvisitedEdges);
			return tempList;
		}
		for (int i = 0; i < listOfSecondaryConnectors.size(); i++) {
			List<String> l = new ArrayList<String>(listOfSecondaryConnectors.get(i));
			for (int j = 0; j < unvisitedEdges.size(); j++) {
				l.add(unvisitedEdges.get(j));
			}
			tempList.add(l);
		}

		return tempList;
	}

	/*
	 * This method is used to calculate the main connectors by computing the
	 * unvisited vertices and visited vertices by primary partition and
	 * calculating the possible edges barring all edges including vref.
	 */
	private Map<Set<String>, List<String>> generateMainConnectors(String vref,
			Map<Set<String>, List<SecondarypartitionClass>> hmapOfSecondaryPartition) {

		Map<Set<String>, List<String>> hmapForMainConnectors = new HashMap<Set<String>, List<String>>();
		Set<Set<String>> hset = hmapOfSecondaryPartition.keySet();
		Iterator it = hset.iterator();
		while (it.hasNext()) {
			Set<String> tempSet = (Set<String>) it.next();
			List<String> listOfMainConnectors = new ArrayList<String>();
			Set<String> differentPrimaryVertices = findDifferentVertices(tempSet);
			Set<String> differentSecondaryVertices = findSecondaryVertices(hmapForPoints, differentPrimaryVertices);
			differentPrimaryVertices.remove(vref);
			// Find all pairs for primary and secondary and store it into list
			// and put it into map
			Iterator it1 = differentPrimaryVertices.iterator();
			while (it1.hasNext()) {
				String key1 = (String) it1.next();
				Iterator it2 = differentSecondaryVertices.iterator();
				while (it2.hasNext()) {
					String key2 = (String) it2.next();
					listOfMainConnectors.add(key1 + key2);
				}
			}

			hmapForMainConnectors.put(tempSet, listOfMainConnectors);
		}

		return hmapForMainConnectors;
	}

	/*
	 * This function calculates the vertices which are not present in primary
	 * set.
	 */
	private Set<String> findSecondaryVertices(Map<String, Point> hmapForPoints, Set<String> differentPrimaryVertices) {

		Set<String> hset = new HashSet<String>(hmapForPoints.keySet());
		Iterator it = differentPrimaryVertices.iterator();
		while (it.hasNext()) {
			String str = (String) it.next();
			hset.remove(str);
		}

		return hset;
	}

	/*
	 * This helper method will return set of different vertices present in the
	 * tempSet.
	 */
	private Set<String> findDifferentVertices(Set<String> tempSet) {

		Set<String> temp = new HashSet<String>();
		Iterator it = tempSet.iterator();
		while (it.hasNext()) {
			String str = (String) it.next();

			for (int i = 0; i < str.length(); i++)
				temp.add(str.charAt(i) + "");
		}
		return temp;
	}

	/*
	 * This method generates secondary connectors for each primary partition. It
	 * takes each primary partition, get all the unvisited edges and vertices,
	 * and find out the combinations that can connect firstly all the unvisited
	 * edges together and then add set of all edges which will connect unvisited
	 * vertex together(if present).
	 */
	private Map<Set<String>, List<SecondarypartitionClass>> generateSecondaryConnectors(
			Map<Set<String>, List<SecondarypartitionClass>> hmapOfSecondaryPartition) {

		Map<Set<String>, List<SecondarypartitionClass>> temp = new HashMap<Set<String>, List<SecondarypartitionClass>>(
				hmapOfSecondaryPartition);
		for (Map.Entry<Set<String>, List<SecondarypartitionClass>> h : temp.entrySet()) {
			Set<String> set = h.getKey();
			List<SecondarypartitionClass> resultList = new ArrayList<SecondarypartitionClass>();
			List<SecondarypartitionClass> tempList = h.getValue();
			for (int i = 0; i < tempList.size(); i++) {
				SecondarypartitionClass spc = tempList.get(i);
				if (spc.unvisitedEdges == null) {
					resultList.add(new SecondarypartitionClass(spc.unvisitedEdges, spc.unvisitedVertex, null));
				} else if (spc.unvisitedVertex.equals("")) {
					if (spc.unvisitedEdges.size() == 1) {
						resultList.add(new SecondarypartitionClass(spc.unvisitedEdges, spc.unvisitedVertex, null));
					} else {
						List<List<String>> list = secondaryConnectors(spc.unvisitedEdges);
						resultList.add(new SecondarypartitionClass(spc.unvisitedEdges, spc.unvisitedVertex, list));
					}

				} else {
					List<List<String>> list = secondaryConnectors(spc.unvisitedEdges);
					list = updateWithUnvisitedVertex(list, spc.unvisitedVertex);
					resultList.add(new SecondarypartitionClass(spc.unvisitedEdges, spc.unvisitedVertex, list));
				}
			}
			hmapOfSecondaryPartition.put(set, resultList);
		}

		return hmapOfSecondaryPartition;
	}

	/*
	 * This function adds every possible edge formed by the unvisited vertex
	 * with the vertices connected by the edges present in the list.
	 */
	private List<List<String>> updateWithUnvisitedVertex(List<List<String>> list, String unvisitedVertex) {

		List<List<String>> list1 = new ArrayList<List<String>>();
		for (int i = 0; i < list.size(); i++) {
			List<String> lStr = list.get(i);
			Set<String> hset = findVertices(lStr);
			Iterator it = hset.iterator();
			while (it.hasNext()) {
				String str = (String) it.next();
				List<String> tempList = new ArrayList<String>();
				tempList.add(str + unvisitedVertex);
				list1.add(tempList);
			}
		}

		return list1;
	}

	/*
	 * This function will return all the different vertices present in the
	 * ArrayList.
	 */
	private Set<String> findVertices(List<String> lStr) {

		Set<String> hset = new HashSet<String>();
		for (int i = 0; i < lStr.size(); i++) {
			String temp = lStr.get(i);
			for (int j = 0; j < temp.length(); j++) {
				hset.add(temp.charAt(j) + "");
			}
		}

		return hset;
	}

	/*
	 * This method finds secondary connectors for list of unvisited edges.
	 */
	private static List<List<String>> secondaryConnectors(List<String> arr) {

		List<List<String>> temp1 = new ArrayList<List<String>>();
		if (arr.size() == 1) {
			temp1.add(arr);
			return temp1;
		}
		Queue<String> queue = new LinkedList<String>();
		queue.add(arr.get(0));

		for (int i = 1; i < arr.size(); i++) {
			List<String> temp = new ArrayList<String>();
			List<String> l1 = new ArrayList<String>();
			while (!queue.isEmpty()) {
				String t = queue.poll();
				l1 = generateDifferentVertices(t, l1);
			}
			List<String> l2 = new ArrayList<String>();
			l2 = generateDifferentVertices(arr.get(i), l2);
			for (int j = 0; j < l1.size(); j++) {
				String first = l1.get(j);
				for (int k = 0; k < l2.size(); k++) {
					String second = l2.get(k);
					queue.add(first + second);
					temp.add(first + second);
				}
			}
			temp1.add(temp);

		}
		List<List<String>> t = findCombinations(temp1);

		return t;
	}

	private static List<List<String>> findCombinations(List<List<String>> temp1) {

		List<List<String>> temp = new ArrayList<List<String>>();
		Queue<List<String>> queue = new LinkedList<List<String>>();
		List<String> l = temp1.get(0);
		for (int i = 0; i < l.size(); i++) {
			String str = l.get(i);
			List<String> listString = new ArrayList<String>();
			listString.add(str);
			queue.add(listString);
		}

		for (int i = 1; i < temp1.size(); i++) {
			Queue<List<String>> queue2 = new LinkedList<List<String>>();
			List<String> list2 = temp1.get(i);
			while (!queue.isEmpty()) {
				List<String> list = queue.poll();
				for (int k = 0; k < list2.size(); k++) {
					List<String> tempList = new ArrayList<String>(list);
					tempList.add(list2.get(k));
					queue2.add(tempList);
				}
			}
			queue = queue2;
		}

		while (!queue.isEmpty()) {
			List<String> lS = queue.poll();
			temp.add(lS);
		}

		return temp;
	}

	/*
	 * This method returns different vertices present in the arraylist.
	 */
	private static List<String> generateDifferentVertices(String t, List<String> l1) {

		for (int i = 0; i < t.length(); i++) {
			String temp = t.charAt(i) + "";
			if (!l1.contains(temp))
				l1.add(temp);
		}
		return l1;
	}

	/*
	 * This method generates secondary partitions for every primary partition.
	 * Here, we will calculate unvisited vertices for each primary vertex. After
	 * calculating unvisited vertices, we will make all different possible
	 * combinations of two for the in a set such that for 5 vertices A, B, C, D,
	 * E the combinations will be {AB,CD,E}, {AC, BD, E},{AD, BC, E} {AB, CE, D}
	 * , {AC, BE, D}, {AE, BC, D} and so on. Here single vertex shows that it
	 * can't be joined together with any other vertex, so it will be counted as
	 * unvisited vertex only.
	 */
	private List<SecondarypartitionClass> generateSecondaryPartition(Set<String> tempSet) {

		Iterator it = tempSet.iterator();
		List<String> unVisitedVertexlist = new ArrayList<String>(hmapForPoints.keySet());
		List<SecondarypartitionClass> l = new ArrayList<SecondarypartitionClass>();
		while (it.hasNext()) {
			String temp = (String) it.next();
			for (int i = 0; i < temp.length(); i++) {
				String t = temp.charAt(i) + "";
				if (unVisitedVertexlist.contains(t))
					unVisitedVertexlist.remove(t);
			}
		}

		if (unVisitedVertexlist.size() == 1) {
			l.add(new SecondarypartitionClass(null, unVisitedVertexlist.get(0)));
		} else if (unVisitedVertexlist.size() == 0) {
			l.add(new SecondarypartitionClass(null, ""));
		} else if (unVisitedVertexlist.size() % 2 == 0) {
			// unvisited vertex will not be there.
			// Write a helper function which will generate an edge of two.
			List<List<String>> unvisitedEdges = generateSecondaryEdges(unVisitedVertexlist);
			for (int i = 0; i < unvisitedEdges.size(); i++) {
				l.add(new SecondarypartitionClass(unvisitedEdges.get(i), ""));
			}
		} else {
			// unvisited vertex will be there.
			for (int i = 0; i < unVisitedVertexlist.size(); i++) {
				String unvisitedVert = unVisitedVertexlist.get(i);
				List<String> newList = new ArrayList<String>(unVisitedVertexlist);
				newList.remove(unvisitedVert);
				List<List<String>> unvisitedEdges = generateSecondaryEdges(newList);
				for (int j = 0; j < unvisitedEdges.size(); j++) {
					l.add(new SecondarypartitionClass(unvisitedEdges.get(j), unvisitedVert));
				}
			}
		}

		return l;
	}

	/*
	 * This method will generate edges that combines the set of unvisited
	 * vertices of two.
	 */
	private List<List<String>> generateSecondaryEdges(List<String> newList) {

		List<String> str = new ArrayList<String>();
		combinations2(newList, str, "");
		List<List<String>> unvisitedEdges = new ArrayList<List<String>>();
		for (int i = 0; i < str.size(); i++) {
			unvisitedEdges.add(convertStringToList(str.get(i)));
		}
		return unvisitedEdges;
	}

	/*
	 * Helper function for upper function
	 */
	private static void combinations2(List<String> arr, List<String> finalList, String s) {

		if (arr.size() == 2) {
			s = s + arr.get(0) + arr.get(1);
			finalList.add(s);
		} else {
			String temp = arr.get(0);

			for (int i = 1; i < arr.size(); i++) {
				List<String> tempList = new ArrayList<String>(arr);
				tempList.remove(temp);
				tempList.remove(arr.get(i));
				combinations2(tempList, finalList, s + temp + arr.get(i) + " ");
			}
		}
	}

	/*
	 * This function takes in String containing edges followed by spaces and
	 * returns list of strings removing spaces.
	 */
	private List<String> convertStringToList(String str) {

		String[] arr = str.split(" ");
		return Arrays.asList(arr);
	}

	/*
	 * This method removes all the edges in the complete graph which contains
	 * vref.
	 */
	private Map<String, List<String>> removeVrefEdges(Map<String, List<String>> hmapOfCompleteGraph, String vref) {

		hmapOfCompleteGraph.remove(vref);
		Iterator it = hmapOfCompleteGraph.keySet().iterator();
		while (it.hasNext()) {
			String tempKey = (String) it.next();
			List<String> tempList = hmapOfCompleteGraph.get(tempKey);
			tempList.remove(vref);
			hmapOfCompleteGraph.put(tempKey, tempList);
		}

		return hmapOfCompleteGraph;
	}

	/*
	 * This method basically takes in a set of strings and generates all the
	 * possible combinations that can be formed from them. Let us say we have 3
	 * numbers in a set. Then we can generate 7 combinations from it.
	 */
	private Set<Set<String>> generateCombinations(Set<String> numberOfEdges) {

		Set<Set<String>> hsetForCombinations = new HashSet<Set<String>>();
		int lengthOfSet = numberOfEdges.size();
		int numberOfSets = (int) Math.pow(2, lengthOfSet);
		List<String> al = new ArrayList<String>(numberOfEdges);
		for (int i = 1; i < numberOfSets; i++) {
			String bin = convertInBinary(i, lengthOfSet);
			Set<String> s = getSubset(bin, al);
			hsetForCombinations.add(s);
		}

		return hsetForCombinations;
	}

	/*
	 * This method converts the binary form of string into values representing
	 * them.
	 */
	private Set<String> getSubset(String bin, List<String> al) {

		Set<String> s = new HashSet<String>();
		for (int i = bin.length() - 1; i >= 0; i--) {
			if (bin.charAt(i) == '1') {
				String val = al.get(i);
				s.add(val);
			}
		}

		return s;
	}

	/*
	 * This method converts a number into its binary form.
	 */
	private String convertInBinary(int i, int lengthOfSet) {

		String bin = Integer.toBinaryString(i);
		bin = String.format("%0" + lengthOfSet + "d", Integer.parseInt(bin));

		return bin;
	}

	/*
	 * This function returns all the edges present in the graph containing vref.
	 * Since this is a complete graph, then the number of edges will be n-1.
	 */
	private Set<String> generateEdgesForVref(String vref, Set<String> hset) {

		Set<String> numberOfEdges = new HashSet<String>();

		Iterator it = hset.iterator();
		while (it.hasNext()) {
			String tempKey = (String) it.next();
			numberOfEdges.add(vref + tempKey);
		}

		return numberOfEdges;
	}

	/*
	 * This function will return a hashmap containing all the vertices connected
	 * to each other.
	 */
	private Map<String, List<String>> generateCompleteGraph() {

		Map<String, List<String>> hmapOfCompleteGraph = new HashMap<String, List<String>>();
		Iterator it1 = hmapForPoints.keySet().iterator();
		while (it1.hasNext()) {
			String temp1 = (String) it1.next();
			Iterator it2 = hmapForPoints.keySet().iterator();
			List<String> arrList = new ArrayList<String>();
			while (it2.hasNext()) {
				String temp2 = (String) it2.next();
				if (!temp2.equals(temp1))
					arrList.add(temp2);
			}
			hmapOfCompleteGraph.put(temp1, arrList);
		}
		return hmapOfCompleteGraph;
	}

	/*
	 * This is the implementation of randomized Markov Chain principle.
	 */
	public static void markovChainRule() {

		// Choose a random vertex and call it root.
		// Create a directed graph from the spanning tree directing towards the
		// root.
		// Implement the Markov chain algorithm.
		Iterator it = hmapForPoints.keySet().iterator();
		String root = (String) (it.next());

		// Generate a directed graph directing towards root.
		Map<String, List<String>> hmapOfDirectedGraph = generateDirectedGraph(root);
		for (int z = 1; z <= couponCollector; z++) {
			boolean isContain = false;
			while (!isContain) {
				int c = 0;
				while (c < 9000) {
					// Implement markov chain algorithm on hmapOfDiretedGraphs.
					// Choose a random edge and a random direction (u,v)
					// if edge is not part of previous spanning tree and u is
					// not the root and it is crossing free
					// T= TUe-parent edge of u.
					// else stay at T.
					String randomEdge = generateRandomEdge();
					String[] arr = randomEdge.split(" ");
					String u = arr[0];
					String v = arr[1];
					// Direction is from u to v by default
					boolean isSpanningTree = notPreviousSpanningTree(arr, hmapOfDirectedGraph);
					if (!isSpanningTree) {
						continue;
					}

					if (u.equals(root)) {
						continue;
					}

					boolean isCrossingFree = checkCrossingFree(arr, hmapOfDirectedGraph);

					if (isCrossingFree) {
						continue;
					}

					if (isSpanningTree && !u.equals(root) && !isCrossingFree) {
						c++;
						String parentOfU = "";
						boolean isFalse = false;
						for (Map.Entry<String, List<String>> h : hmapOfDirectedGraph.entrySet()) {
							List<String> tempList = h.getValue();
							for (int i = 0; i < tempList.size(); i++) {
								if (tempList.get(i).equals(u)) {
									parentOfU = h.getKey();
									isFalse = true;
									break;
								}
							}
							if (isFalse)
								break;
						}

						List<String> temp = hmapOfDirectedGraph.get(parentOfU);
						if (temp.size() == 1)
							hmapOfDirectedGraph.remove(parentOfU);
						else {
							temp.remove(u);
							hmapOfDirectedGraph.put(parentOfU, temp);
						}

						// Add e to the previous spanning tree
						if (hmapOfDirectedGraph.containsKey(v)) {
							List<String> tempList = hmapOfDirectedGraph.get(v);
							tempList.add(u);
							hmapOfDirectedGraph.put(v, tempList);
						} else {
							List<String> tempList = new ArrayList<String>();
							tempList.add(u);
							hmapOfDirectedGraph.put(v, tempList);
						}

					}
				}
				Set<String> treeSet = new TreeSet<String>();
				for (Map.Entry<String, List<String>> h : hmapOfDirectedGraph.entrySet()) {
					String key = h.getKey();
					List<String> value = h.getValue();
					for (int x = 0; x < value.size(); x++) {
						treeSet.add(key + value.get(x));
					}
				}
				List<String> tempList = new ArrayList<String>(treeSet);
				if (hmapForMarkovChainCount.containsKey(tempList)) {
					isContain = true;
					hmapForMarkovChainCount.put(tempList, hmapForMarkovChainCount.get(tempList) + 1);
				}
			}
		}
	}

	/*
	 * Given three colinear points p, q, r, the function checks if point q lies
	 * on line segment 'pr'
	 */
	public static boolean onSegment(Point p, Point q, Point r) {

		if (q.x <= Math.max(p.x, r.x) && q.x >= Math.min(p.x, r.x) && q.y <= Math.max(p.y, r.y)
				&& q.y >= Math.min(p.y, r.y))
			return true;

		return false;
	}

	/*
	 * The function that returns true if line segment 'p1q1' and 'p2q2'
	 * intersect.
	 */
	public static boolean checkIfIntersect(Point p1, Point q1, Point p2, Point q2) {

		// Find the four orientations needed for general and
		// special cases
		int o1 = orientation(p1, q1, p2);
		int o2 = orientation(p1, q1, q2);
		int o3 = orientation(p2, q2, p1);
		int o4 = orientation(p2, q2, q1);

		// General case
		if (o1 != o2 && o3 != o4)
			return true;

		// Special Cases
		// p1, q1 and p2 are colinear and p2 lies on segment p1q1
		if (o1 == 0 && onSegment(p1, p2, q1))
			return true;

		// p1, q1 and p2 are colinear and q2 lies on segment p1q1
		if (o2 == 0 && onSegment(p1, q2, q1))
			return true;

		// p2, q2 and p1 are colinear and p1 lies on segment p2q2
		if (o3 == 0 && onSegment(p2, p1, q2))
			return true;

		// p2, q2 and q1 are colinear and q1 lies on segment p2q2
		if (o4 == 0 && onSegment(p2, q1, q2))
			return true;

		return false; // Doesn't fall in any of the above cases
	}

	/*
	 * This method returns if the edge represented by an array is present in the
	 * hashmap, representing spanning tree or not.
	 */
	private static boolean notPreviousSpanningTree(String[] arr, Map<String, List<String>> hmapOfDirctedGraph) {

		String u = arr[0];
		String v = arr[1];
		if (hmapOfDirctedGraph.containsKey(u)) {
			List<String> listOfU = hmapOfDirctedGraph.get(u);
			for (int i = 0; i < listOfU.size(); i++) {
				if (listOfU.get(i).equals(v))
					return false;
			}
		}
		if (hmapOfDirctedGraph.containsKey(v)) {
			List<String> listOfV = hmapOfDirctedGraph.get(v);
			for (int i = 0; i < listOfV.size(); i++) {
				if (listOfV.get(i).equals(u))
					return false;
			}
		}

		return true;
	}

	/*
	 * This method checks if the edge represented by an array is crossing any
	 * edge present in the hashmap. If yes, then it will return true.
	 */
	private static boolean checkCrossingFree(String[] arr, Map<String, List<String>> hmapOfDirctedGraph) {

		String u = arr[0];
		String v = arr[1];
		Point p1 = hmapForPoints.get(u);
		Point q1 = hmapForPoints.get(v);
		Iterator it = hmapOfDirctedGraph.keySet().iterator();
		while (it.hasNext()) {
			String key = (String) it.next();
			Point p2 = hmapForPoints.get(key);
			List<String> tempList = new ArrayList<String>();
			for (int i = 0; i < tempList.size(); i++) {
				Point q2 = hmapForPoints.get(tempList.get(i));
				if (u.equals(key) || u.equals(tempList.get(i)) || v.equals(key) || v.equals(tempList.get(i))) {

				} else {
					// This function will return true if they intersect.
					if (checkIfIntersect(p1, q1, p2, q2))
						return true;
				}
			}
		}

		return false;
	}

	/*
	 * This method generates random edge by generating random numbers and
	 * finding the vertex representing that number.
	 */
	private static String generateRandomEdge() {

		Random rand = new Random();
		int randomUValue = rand.nextInt(hmapForPoints.size());
		Iterator it = hmapForPoints.keySet().iterator();
		int k = -1;
		String u = "";
		while (k != randomUValue) {
			u = (String) it.next();
			k++;
		}
		int randomVValue = rand.nextInt(hmapForPoints.size());
		while (randomUValue == randomVValue) {
			randomVValue = rand.nextInt(hmapForPoints.size());
		}
		Iterator it1 = hmapForPoints.keySet().iterator();
		k = -1;
		String v = "";
		while (k != randomVValue) {

			v = (String) it1.next();
			k++;
		}

		return u + " " + v;
	}

	/*
	 * This method will take CFST generated from the novel CFST algorithm and
	 * will convert it into directed graph directing towards random root.
	 */
	private static Map<String, List<String>> generateDirectedGraph(String root) {

		Set<Integer> hset = new HashSet<Integer>();
		Map<String, List<String>> hmapForDirectedGraph = new HashMap<String, List<String>>();
		Queue<String> queue = new LinkedList<String>();
		queue.add(root);

		while (!queue.isEmpty()) {
			String tempString = queue.poll();
			for (int i = 0; i < listOfSpanningTreeLines.size(); i++) {
				Line l = listOfSpanningTreeLines.get(i);
				if (!hset.contains(i) && l.firstPoint.equals(tempString)) {
					hset.add(i);
					queue.add(l.secondPoint);
					if (hmapForDirectedGraph.containsKey(tempString)) {
						List<String> tempList = hmapForDirectedGraph.get(tempString);
						tempList.add(l.secondPoint);
						hmapForDirectedGraph.put(tempString, tempList);
					} else {
						List<String> tempList = new ArrayList<String>();
						tempList.add(l.secondPoint);
						hmapForDirectedGraph.put(tempString, tempList);
					}
				} else if (!hset.contains(i) && l.secondPoint.equals(tempString)) {
					hset.add(i);
					queue.add(l.firstPoint);
					if (hmapForDirectedGraph.containsKey(tempString)) {
						List<String> tempList = hmapForDirectedGraph.get(tempString);
						tempList.add(l.firstPoint);
						hmapForDirectedGraph.put(tempString, tempList);
					} else {
						List<String> tempList = new ArrayList<String>();
						tempList.add(l.firstPoint);
						hmapForDirectedGraph.put(tempString, tempList);
					}
				}
			}
		}

		return hmapForDirectedGraph;
	}

	/*
	 * This is the novel CFST algorithm which generates CFST by forming a convex
	 * hull around it and joining the convex hull with the rest of the inner
	 * vertices.
	 */
	public static void generateSpanningTree() {

		String returnString;
		if (hmapForPoints.size() <= 1)
			System.out.println("Not Possible");
		else if (hmapForPoints.size() == 2) {
			// Brute force algorithm
			String firstPoint, secondPoint;
			Set<String> hset = hmapForPoints.keySet();
			Iterator it = hset.iterator();
			firstPoint = (String) it.next();
			secondPoint = (String) it.next();
			// Draw line between them.
		} else if ((returnString = pointsCollinear()).contains("True")) {
			// Now we know that the lines are collinear.
			// If it is A then the line will be parallel to y-axis.
			// If it is B then it is parallel to x-axis.
			// If it is C then it is at an angle.
			String lowPoint = "", highPoint = "";
			if (returnString.charAt(4) == 'A') {
				int low = Integer.MAX_VALUE, high = Integer.MIN_VALUE;
				for (Map.Entry<String, Point> h : hmapForPoints.entrySet()) {
					Point p = h.getValue();
					if (p.y > high) {
						high = p.y;
						highPoint = h.getKey();
					}
					if (p.y < low) {
						low = p.y;
						lowPoint = h.getKey();
					}
				}
				// High point and low point are there.
				// just draw the line.
			} else if (returnString.charAt(4) == 'B') {
				int low = Integer.MAX_VALUE, high = Integer.MIN_VALUE;
				for (Map.Entry<String, Point> h : hmapForPoints.entrySet()) {
					Point p = h.getValue();
					if (p.x > high) {
						high = p.x;
						highPoint = h.getKey();
					}
					if (p.x < low) {
						low = p.x;
						lowPoint = h.getKey();
					}
				}
				// High point and low point are there.
				// just draw the line.
			} else {
				String tempSt = returnString.substring(6);
				String[] arr = tempSt.split(" ");
				highPoint = arr[0];
				lowPoint = arr[1];
				// Draw the line.
			}
		} else if (hmapForPoints.size() == 3) {
			String firstPoint, secondPoint, thirdPoint;
			Set<String> hset = hmapForPoints.keySet();
			Iterator it = hset.iterator();
			firstPoint = (String) it.next();
			secondPoint = (String) it.next();
			thirdPoint = (String) it.next();
			// Draw two lines
		} else {
			// First generate convex hull
			// get the inner points and calculate minimum distance with the
			// outer vertices.
			// get the min distance and draw the line.
			Collection<Point> c = hmapForPoints.values();
			List<Point> list = new ArrayList<Point>(c);
			listOfSpanningTreeLines = convexHull(list, list.size());
			Set<String> outerPoints = new HashSet<String>();
			for (int i = 0; i < listOfSpanningTreeLines.size(); i++) {
				outerPoints.add(listOfSpanningTreeLines.get(i).firstPoint);
				outerPoints.add(listOfSpanningTreeLines.get(i).secondPoint);
			}
			Set<String> innerPoints = new HashSet<String>();
			Iterator it = hmapForPoints.keySet().iterator();
			while (it.hasNext()) {
				String str = (String) it.next();
				if (!outerPoints.contains(str)) {
					innerPoints.add(str);
				}
			}
			// Lines Inside circle
			generateInnerLines(innerPoints, outerPoints);
		}
	}

	/*
	 * This method will generate inner lines for the convex hull to complete
	 * CFST. It will compute the closest outer vertex for each inner vertex and
	 * join them.
	 */
	private static void generateInnerLines(Set<String> innerPoints, Set<String> outerPoints) {

		Iterator it = innerPoints.iterator();
		while (it.hasNext()) {
			String innerPoint = (String) it.next();
			String closestOuterPoint = getClosestPoint(innerPoint, outerPoints);
			listOfSpanningTreeLines.add(new Line(innerPoint, closestOuterPoint));
		}

	}

	/*
	 * This function will find the closest outer point for inner vertex.
	 */
	private static String getClosestPoint(String innerPoint, Set<String> outerPoints) {

		int distance = Integer.MAX_VALUE;
		String closestPoint = "";
		Iterator it = outerPoints.iterator();
		while (it.hasNext()) {
			String outerPoint = (String) it.next();
			int temp = distSq(hmapForPoints.get(outerPoint), hmapForPoints.get(innerPoint));
			if (distance > temp) {
				distance = temp;
				closestPoint = outerPoint;
			}
		}

		return closestPoint;
	}

	/*
	 * This function will generate a convex hull for the set of points.
	 */
	public static List<Line> convexHull(List<Point> points, int n) {

		// Find the bottommost point
		int ymin = points.get(0).y, min = 0;
		for (int i = 1; i < n; i++) {
			int y = points.get(i).y;

			// Pick the bottom-most or chose the left
			// most point in case of tie
			if ((y < ymin) || (ymin == y && points.get(i).x < points.get(min).x)) {
				ymin = points.get(i).y;
				min = i;
			}
		}

		// Place the bottom-most point at first position
		Collections.swap(points, 0, min);

		// Sort n-1 points with respect to the first point.
		// A point p1 comes before p2 in sorted ouput if p2
		// has larger polar angle (in counterclockwise
		// direction) than p1
		p0 = points.get(0);
		Collections.sort(points, new PrivateComparator());

		// If two or more points make same angle with p0,
		// Remove all but the one that is farthest from p0
		// Remember that, in above sorting, our criteria was
		// to keep the farthest point at the end when more than
		// one points have same angle.
		int m = 1;
		for (int i = 1; i < n; i++) {
			// Keep removing i while angle of i and i+1 is same
			// with respect to p0
			while (i < n - 1 && orientation(p0, points.get(i), points.get(i + 1)) == 0)
				i++;

			points.set(m, points.get(i));
			m++;
		}

		// If modified array of points has less than 3 points,
		// convex hull is not possible
		if (m < 3)
			return new ArrayList<Line>();

		Stack<Point> S = new Stack<Point>();
		S.push(points.get(0));
		S.push(points.get(1));
		S.push(points.get(2));

		// Process remaining n-3 points
		for (int i = 3; i < m; i++) {
			// Keep removing top while the angle formed by
			// points next-to-top, top, and points[i] makes
			// a non-left turn
			while (orientation(secondElementFromTop(S), S.peek(), points.get(i)) != 2)
				S.pop();
			S.push(points.get(i));
		}
		Stack<String> tempStack = new Stack<String>();
		while (!S.isEmpty()) {
			Point p = S.pop();
			for (Map.Entry<String, Point> th : hmapForPoints.entrySet()) {
				if (th.getValue().x == p.x && th.getValue().y == p.y) {
					tempStack.push(th.getKey());
					break;
				}
			}
		}
		List<Line> listOfPoints = new ArrayList<Line>();
		String firstPoint = "", secondPoint = "";
		while (!tempStack.empty()) {
			String str = tempStack.pop();
			if (firstPoint.equals("")) {
				firstPoint = str;
				continue;
			} else if (!firstPoint.equals("") && secondPoint.equals("")) {
				secondPoint = str;
			}
			if (!firstPoint.equals("") && !secondPoint.equals("")) {
				listOfPoints.add(new Line(firstPoint, secondPoint));
				firstPoint = secondPoint;
				secondPoint = "";
			}

		}
		return listOfPoints;
	}

	/*
	 * This function returns second topmost point in the stack.
	 */
	private static Point secondElementFromTop(Stack<Point> s) {

		Point top = s.pop();
		Point secondFromTop = s.peek();
		s.push(top);

		return secondFromTop;
	}

	/*
	 * This is the main function from where the execution starts. It takes input
	 * from a file containing vertices followed by x-axis and y-axis(seperated
	 * by spaces). It will generate all CFSTs based on points and their
	 * locations. After calculating all CFSTs, we will generate novel CFST and
	 * use Markov Chain to randomly generate new CFST and will run till we don't
	 * get all CFSTs. Following this, we will print all CFSTs.
	 */
	public static void main(String[] args) throws NumberFormatException, IOException {

		GenerateSpanningTree fp = new GenerateSpanningTree();
		FileReader fr = new FileReader(new File("src/TestFile"));
		BufferedReader br = new BufferedReader(fr);
		String str = "";
		while ((str = br.readLine()) != null && !str.equals("TestCases")) {
			String[] arr = str.split(" ");
			hmapForPoints.put(arr[0].toUpperCase(), new Point(Integer.parseInt(arr[1]), Integer.parseInt(arr[2])));
		}

		fp.generateAllSpanningTrees();
		fp.generateSpanningTree();
		int nValue = hmapForMarkovChainCount.size();
		couponCollector = (int) ((nValue * Math.log(nValue)) + MConstant * nValue);
		fp.markovChainRule();

		int size = fp.resultantCrossingFreeSpanningTrees.size();
		do {
			JFrame f = new JFrame("Final_Project");
			f.addWindowListener(new WindowAdapter() {
				public void windowClosing(WindowEvent e) {
					System.exit(0);
				}
			});
			JApplet applet = fp;
			f.getContentPane().add("Center", applet);
			applet.init();
			f.pack();
			f.setSize(new Dimension(1000, 800));
			f.setVisible(true);
			fp.saveImage(applet, count);

			++count;
		} while (count < size);
		System.out.println("Program Finished");
	}

	/*
	 * This function rearranges the spanning in the sorting order. Firstly we
	 * will sort the characters in the string and then the strings. Ex: for AB,
	 * XF, DA, CA the result will be AB, AC, AD, FX
	 */
	private void rearrangeSpanningTrees() {

		for (int i = 0; i < resultantCrossingFreeSpanningTrees.size(); i++) {
			List<String> t = resultantCrossingFreeSpanningTrees.get(i);
			Set<String> treeSet = new TreeSet<String>();
			for (int j = 0; j < t.size(); j++) {
				String str = t.get(j);
				char[] chars = str.toCharArray();
				Arrays.sort(chars);
				String sorted = new String(chars);
				treeSet.add(sorted);
			}
			List<String> list = new ArrayList<String>(treeSet);
			hmapForMarkovChainCount.put(list, 0);

		}

	}

	/*
	 * Below is the function to save the JFrame generated to save it at a
	 * specified location.
	 */
	private void saveImage(JApplet panel, int i) {

		BufferedImage img = new BufferedImage(panel.getWidth(), panel.getHeight(), BufferedImage.TYPE_INT_RGB);
		panel.paint(img.getGraphics());
		try {
			ImageIO.write(img, "png", new File("C:/Users/sawan/Desktop/Results/CFST" + i + ".png"));

		} catch (Exception e) {
			System.out.println("panel not saved" + e.getMessage());
		}
	}
}