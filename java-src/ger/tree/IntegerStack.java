package ger.tree;

public class IntegerStack
{

    private int[] array;
    private int head;

    public IntegerStack(int capacity)
    {
	array = new int[capacity];
	head = 0;
    }

    public void clean()
    {
	head = 0;
    }

    public void push(int x)
    {
	array[head] = x;
	head++;
    }

    public int pop()
    {
	head--;
	int x = array[head];
	return x;
    }

    public boolean isEmpty()
    {
	return head == 0;
    }

    @Override
    public String toString()
    {
	StringBuilder sb = new StringBuilder();
	for (int i = head - 1; i >= 0; i--)
	{
	    sb.append("->").append(array[i]);
	}
	return sb.toString();
    }

    /** ha vannak gyerekei az x-nek akkor beteszi a verembe */
    public void pushChildren(DoubleArrayTree tree, int x)
    {

	int c1 = tree.getChildA(x);
	if (c1 != -1)
	{
	    this.push(c1);
	}
	int c2 = tree.getChildB(x);
	if (c2 != -1)
	{
	    this.push(c2);
	}
    }

}
