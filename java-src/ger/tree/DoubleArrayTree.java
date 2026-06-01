package ger.tree;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ger.maddison.Node;
import ger.maddison.NodeUtils;
import ger.maddison.NwkNode;

/**
 * Double-array tree representation for phylogenetic trees.
 *
 * <p>
 * This representation can be initialized from a tree represented by
 * {@link NwkNode} / {@link Node}. Each tree node is assigned a compact integer
 * index in the range {@code 0..N-1}. This makes it possible to store per-node
 * values in primitive arrays, such as {@code boolean[]}, {@code int[]}, or
 * {@code String[]}, instead of using {@link Set} or {@link Map}-based
 * structures.
 * </p>
 *
 * <p>
 * The representation is useful in the Maddison simulation pipeline because many
 * operations repeatedly attach binary states, event labels, or masks to all
 * nodes of the phylogenetic tree. Array-based access is more compact and faster
 * than map-based node lookup.
 * </p>
 */
public class DoubleArrayTree
{
    static final int[] emptyChildArray = new int[0]; // Reused empty child array for leaf nodes.

    int root;
    int[] child;
    int[] parent;
    String[] names;

    private int firstFreeArrayPosition;

    private int alloc(int size)
    {
	if (size == 0)
	{
	    return -1;
	}
	// Double the backing-array sizes if more capacity is needed.
	if (child.length < firstFreeArrayPosition + size)
	{
	    int N = child.length * 2;
	    int[] tmp = new int[N];
	    System.arraycopy(child, 0, tmp, 0, child.length);
	    child = tmp;

	    tmp = new int[N];
	    System.arraycopy(parent, 0, tmp, 0, parent.length);
	    parent = tmp;

	    String[] tmpS = new String[N];
	    System.arraycopy(names, 0, tmpS, 0, names.length);
	    names = tmpS;
	}

	int tmp = firstFreeArrayPosition;
	firstFreeArrayPosition += size;
	return tmp;

    }

    private void trimArraySizes()
    {
   	// Trim the backing arrays to the number of used positions.
	int N = firstFreeArrayPosition;
	int[] tmp = new int[N];
	System.arraycopy(child, 0, tmp, 0, N);
	child = tmp;

	tmp = new int[N];
	System.arraycopy(parent, 0, tmp, 0, N);
	parent = tmp;

	String[] tmpS = new String[N];
	System.arraycopy(names, 0, tmpS, 0, N);
	names = tmpS;

    }

    public int size()
    {
	return firstFreeArrayPosition;
    }

    public int getRoot()
    {
	return root;
    }

    public String getName(int x)
    {
	return names[x];
    }

    public int getParent(int x)
    {
	return parent[x];
    }

    public boolean isLeaf(int x)
    {
	return child[x] == -1; // TODO
    }

    public int[] getChildren(int x)
    {

	if (isLeaf(x))
	{
	    return emptyChildArray;
	}

	int c0 = child[x];
	if (parent[c0] == x && parent[c0 + 1] == x)
	{
	    return new int[] { c0, c0 + 1 };
	} else if (parent[c0] == x)
	{
	    return new int[] { c0 };
	} else
	{
		// This branch should never be reached for a valid double-array tree.
	    return emptyChildArray;
	}
    }

    public int getChildA(int x)
    {
	int c = child[x];
	if (c == -1)
	{
	    return -1;
	}
	int res = (parent[c] == x) ? c : -1;
	return res;
    }

    public int getChildB(int x)
    {
	int c = child[x] + 1;
	if (c == 0 || c >= parent.length)  // This tests whether child[x] was -1 before adding 1.
	{
	    return -1;
	}
	int res = (parent[c] == x) ? c : -1;
	return res;
    }

    /**
     * Return the internal node-name array.
     *
     * <p>
     * The returned array is the internal backing array and <b>must not be modified</b> by
     * callers.
     * </p>
     */
    public String[] getNames()
    {
	return names;
    }

    // --------------------------------------------
    // --------------------------------------------
    // CONSTRUCTOR

    public DoubleArrayTree(File file) throws IOException
    {
	this(NwkNode.parseNWK(file));

    }

    public DoubleArrayTree(File file, String subroot) throws IOException
    {
	this(NodeUtils.findFirst(NwkNode.parseNWK(file), subroot));

    }
    
    public DoubleArrayTree(Node treeRoot)
    {
	int N = 10; // initial array size
	this.child = new int[N];
	this.parent = new int[N];
	this.names = new String[N];

	this.root = 0;
	this.firstFreeArrayPosition = 1;

	build(treeRoot, -1, 0);

	trimArraySizes();
    }

    /**
     * Recursively build the double-array representation.
     *
     * @param n
     *            node to be represented at the current position
     * @param parent
     *            array index of the parent node, or {@code -1} for the root
     * @param position
     *            array index assigned to the current node
     */
    private void build(Node n, int parent, int position)
    {

	this.parent[position] = parent;
	this.names[position] = n.getName();

	List<Node> children = n.getChildren();
	int chilPosition = alloc(children.size());
	this.child[position] = chilPosition;

	int offset = 0;
	for (Node c : children)
	{
	    build(c, position, chilPosition + offset);
	    offset++;
	}

    }

    // --------------------------------------------------------------

    /***
     * Wrapper that simulates the Node interface.
     * 
     * It contains only an integer as position and a reference to the
     * {@link DoubleArrayTree}
     * 
     */
    class DANode implements Node
    {

	final int pos;

	public DANode(int pos)
	{
	    this.pos = pos;
	}

	@Override
	public boolean isTerminal()
	{
	    return DoubleArrayTree.this.isLeaf(pos);
	}

	@Override
	public List<Node> getChildren()
	{
	    int[] children = DoubleArrayTree.this.getChildren(pos);
	    ArrayList<Node> list = new ArrayList<Node>(children.length);

	    for (int c : children)
	    {
		list.add(new DANode(c));
	    }
	    return list;
	}

	@Override
	public String getName()
	{
	    return DoubleArrayTree.this.getName(pos);

	}

	@Override
	public String toString()
	{
	    StringBuilder sb = new StringBuilder();
	    List<Node> children = getChildren();
	    if (!children.isEmpty())
	    {
		sb.append('(');
		for (Iterator<Node> iterator = children.iterator(); iterator.hasNext();)
		{
		    Node node = (Node) iterator.next();
		    sb.append(node.toString());
		    sb.append(",");
		}
		sb.setLength(sb.length() - 1); // Remove the trailing comma.
		sb.append(')');
	    }

	    sb.append(getName());
	    // sb.append(':');
	    // sb.append(length);
	    return sb.toString();
	}
    }

    // ------------------------

   
}
