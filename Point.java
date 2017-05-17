/*
 * This class is used as a basic entity whose object will be stored in the main class.
 * This class basically contains x and y values of the vertices of the graph. 
 */
public class Point {
	
	int x,y;
	public Point(int x,int y)
	{
		this.x=x;
		this.y=y;
	}
	public void display()
	{
		System.out.println("X:"+x+" Y:"+y);
	}

}
