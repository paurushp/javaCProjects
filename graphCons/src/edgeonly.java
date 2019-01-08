import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
public class Graph {
	 public static void main (String args[]) throws FileNotFoundException{
	 
int row;//row and columns read from console
System.out.println("Please enter the number of S genes and hit return key");
Scanner in = new Scanner(System.in);
//Input from console
row=in.nextInt();
System.out.println("Please enter the number of time steps and hit return key");
Scanner step = new Scanner(System.in);
int time=step.nextInt();
//Matrix initialization
int col=row;
int array[][] = new int[row][col];
//read matrix data from file
Scanner myScanner = new Scanner (new File("/home/bit/praveen/Desktop/adj_mat.txt")); 
	//Outer loop for each row                                                                                                                                                                                                                                                                                                                                           
for (int i = 0; i < row; i++) {
	
	//Inner loop for columns
	for (int j = 0; j < col; j++) {
		//Read from file
		 array[i][j] = myScanner.nextInt(); 
		 //decision braces for edges and time delays based on matrix entries
		 if (array[i][j]!=0)
			{
			 System.out.println(time+"\t"+"Gene"+i+"\t"+"Gene"+j+"\t"+"Edge=True"+"\t"+"t="+array[i][j]);
		
			}
		else 
			{
			System.out.println(time+"\t"+"Gene"+i+"\t"+"Gene"+j+"\t"+"Edge=False"+"\t"+"t="+array[i][j]);
			}
			 System.out.println();
	 }
           /* Error reported if number of genes is more than that in the data matrix 
              at java.util.Scanner.throwFor(Scanner.java:838)
	          at java.util.Scanner.next(Scanner.java:1461)
	          at java.util.Scanner.nextInt(Scanner.java:2091)
	          at java.util.Scanner.nextInt(Scanner.java:2050)
	          at Graph.main(Graph.java:21)*/
                                                                                                                                                                                                                                            
    }
  }}


