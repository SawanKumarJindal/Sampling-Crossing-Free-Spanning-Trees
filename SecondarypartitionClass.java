import java.util.List;
import java.util.Map;
/*
 * This class is used as a basic entity whose object will be stored in the main class.
 * It will store unvisited edges, unvisited vertices, secondary connectors for the primary partitions
 * 
 */
public class SecondarypartitionClass {
	public List<String> unvisitedEdges;
	public String unvisitedVertex;
	public List<String> secondaryConnectors;
	public List<List<String>> listOfSecondaryConnectors;
	public SecondarypartitionClass(List<String> unvisitedEdges, String unvisitedVertex)
	{
		this.unvisitedEdges = unvisitedEdges;
		this.unvisitedVertex = unvisitedVertex;
	}
	public SecondarypartitionClass(List<String> unvisitedEdges, String unvisitedVertex, List<List<String>> listOfSecondaryConnectors)
	{
		this.unvisitedEdges = unvisitedEdges;
		this.unvisitedVertex = unvisitedVertex;
		this.listOfSecondaryConnectors=listOfSecondaryConnectors;
	}
}
