/*
 * This class is used as a basic entity whose object will be stored in the main class.
 * It contains first and second points which joins together to form a line.
 */
public class Line {
	String firstPoint,secondPoint;
	public Line(String firstPoint, String secondPoint)
	{
		this.firstPoint = firstPoint;
		this.secondPoint = secondPoint;
	}
	public void display()
	{
		System.out.println("FirstPoint:"+firstPoint+" SecondPoint:"+secondPoint);
	}
}
