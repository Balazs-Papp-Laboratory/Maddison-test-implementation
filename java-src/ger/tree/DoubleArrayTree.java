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

/***
 * Dupla tombos fa reprezentacio fiolgenetikai fáknak.
 * 
 * Csak {@link NwkNode}-ban reprezentalt fából lehet inicializálni.
 * 
 * Azért jó, mart minden node-nak van egy poziciós indexe, ami azt a node-ot
 * reprezentálja. ezek szépen 0-tol N-ig mennek. Lehet tombokkel true/false
 * értékeket vagy Stringeket rendelni minden nodejahoz a fának, ami gyorsabb és
 * tömörebb mit {@link Set}-et vagy {@link Map}-et használni. Mivel gyorsítani
 * szeretnem a szamitast erre szukseg lesz.
 */
public class DoubleArrayTree
{
    static final int[] emptyChildArray = new int[0]; // sokszor kell ures tombot visszaadni. Olyankor ezt a konstanst
						     // adom

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
	// ha kell megnovelem a tomb mereteket duplajara
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
	// A tomb mereteket lecsokkentem, hogy minimalisra
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
	    // elvileg ide soha nem jön a vezerles
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
	if (c == 0 | c >= parent.length) // igazablol -1 re tesztelem, csak már hozzaadtam 1-et
	{
	    return -1;
	}
	int res = (parent[c] == x) ? c : -1;
	return res;
    }

    /**
     * fontos, h nem szabad modositani
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
	this(NodeUtils.findFirst(NwkNode.parseNWK(file),subroot));

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

    /***
     * @param n
     *            ez a node kepzodik le
     * @param parent
     *            ez a node apjanak a pozicioja a tombos reprezentacioba
     * @param position
     *            az az aktualis n node pozicioja, ahova teszi. rekurziv fa epito
     *            fuggveny, amit a konstruktor hiv
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
    // --------------------------------------------

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
		sb.setLength(sb.length() - 1); // remove the last column
		sb.append(')');
	    }

	    sb.append(getName());
	    // sb.append(':');
	    // sb.append(length);
	    return sb.toString();
	}
    }

    // ------------------------

    public static void main(String[] args) throws IOException
    {
	
	DoubleArrayTree dat = new DoubleArrayTree(new File("data/tree/ST_tree_renamedtips_mrca.nwk"),"Node16060");

//	// Node tree = NwkNode.parseNWK(new File("data/tree/mini.nwk"));
//	Node tree = NwkNode.parseNWK(new File("data/tree/ST_tree_renamedtips_mrca.nwk"));
//	// Node tree = NwkNode.parseNWK(new
//	// File("data/tree/rapidnj_tree_nonegative_outgroupMRCARooted.nwk"));
//
//	DoubleArrayTree dat = new DoubleArrayTree(tree);

	DANode root = dat.new DANode(0);

//	System.out.println(tree);
	System.out.println(root);

	System.out.println(dat.size());
//	System.out.println(NodeUtils.collectAllDescendants(tree).size());
	System.out.println(NodeUtils.collectAllDescendants(root).size());

	System.out.println("Done");

    }
}
