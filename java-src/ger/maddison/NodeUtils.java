package ger.maddison;

import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;

public class NodeUtils
{

    /**
     * rekurzívan begyűjti a gyerekeit. A gyökeret és az összes alatta lévő node-ot
     * bepakolja egy Set-be.
     * 
     * @param r
     * @return
     */
    public static Set<Node> collectAllDescendants(Node r)
    {
	LinkedList<Node> gray = new LinkedList<Node>();
	gray.push(r);
	LinkedHashSet<Node> set1 = new LinkedHashSet<Node>();

	while (!gray.isEmpty())
	{
	    Node x = gray.removeLast();
	    boolean notAlreadyContaind = set1.add(x);
	    if (notAlreadyContaind)
	    {
		gray.addAll(x.getChildren());
	    }

	}
	return set1;
    }

    public static Node findFirst(Node r, String name)
    {

	LinkedList<Node> gray = new LinkedList<Node>();
	gray.push(r);

	while (!gray.isEmpty())
	{
	    Node x = gray.removeLast();
	    if (name.equals(x.getName()) )
	    {
		return x;
	    }
	    gray.addAll(x.getChildren());
	}

	return null;
    }

    public static String[] collectNamesOfAllDescendants(Node r)
    {
	LinkedHashSet<String> nameSet1 = new LinkedHashSet<String>();

	Set<Node> set1 = collectAllDescendants(r);

	for (Node n : set1)
	{
	    nameSet1.add(n.getName());
	}
	return nameSet1.toArray(new String[nameSet1.size()]);
    }

    public static Set<Node> filterNodeSetByNames(Set<Node> nodeSet, Set<String> nameSet)
    {
	LinkedHashSet<Node> set1 = new LinkedHashSet<Node>();

	for (Node n : nodeSet)
	{
	    if (nameSet.contains(n.getName()))
	    {
		set1.add(n);
	    }
	}
	return set1;
    }

    ////////////////////////////////

    /**
     * Itt nem is tudom mit akarok. Valamiathivas R-bol
     * 
     * @param keys
     * @param values
     * @return
     */
    public static Map<String, Integer> createMap(String[] keys, int[] values)
    {
	if (keys.length != values.length)
	{
	    throw new IllegalArgumentException("Different length of keys and values parameter");
	}

	LinkedHashMap<String, Integer> map = new LinkedHashMap<String, Integer>(keys.length);
	for (int i = 0; i < keys.length; i++)
	{
	    Integer previousValue = map.put(keys[i], values[i]);
	    if (previousValue != null)
	    {
		throw new IllegalArgumentException(
			"keys parameter contains replicates. keys[" + i + "]=\"" + keys[i] + "\"");
	    }
	}

	return map;
    }

    /**
     * Itt nem is tudom mit akarok. Valamiathivas R-bol
     * 
     * @param elements
     * @param values
     * @return
     */
    public static Set<String> createMap(String[] elements)
    {

	LinkedHashSet<String> map = new LinkedHashSet<String>(elements.length);
	for (int i = 0; i < elements.length; i++)
	{

	    boolean isNewInTheSet = map.add(elements[i]);
	    if (!isNewInTheSet)
	    {
		throw new IllegalArgumentException(
			"keys parameter contains replicates. elemets[" + i + "]=\"" + elements[i] + "\"");
	    }
	}

	return map;
    }

    public static int[][] ize(Map<String, Integer> map1, Map<String, Integer> map2, Node r)
    {

	System.out.println("IZE invoked");
	int[][] result = new int[2][2];

	LinkedList<Node> gray = new LinkedList<Node>();
	gray.push(r);
	LinkedHashSet<Node> set1 = new LinkedHashSet<Node>();

	while (!gray.isEmpty())
	{
	    Node x = gray.removeLast();
	    boolean notAlreadyContaind = set1.add(x);
	    if (notAlreadyContaind)
	    {
		gray.addAll(x.getChildren());

		// Integer x_v1_obj = map1.get(x.getName());
		// Integer x_v2_obj = map2.get(x.getName());
		//
		// if (x_v1_obj == null || x_v2_obj == null) {
		// System.out.println("ERROR: '" + x.getName() + "' nincs a map-ben.");
		// } else {

		int x_v1 = map1.get(x.getName());
		int x_v2 = map2.get(x.getName());

		int present1 = (x_v1 != 0 ? 1 : 0);
		int present2 = (x_v2 != 0 ? 1 : 0);
		result[present1][present2]++;
		// for(Node y :x.getChildren())
		// {
		// int y_v1 = map1.get(y.getName());
		// int y_v2 = map2.get(y.getName());
		//
		//
		// boolean gain1= x_v1==0 & y_v1!=0;
		// boolean gain2= x_v2==0 & y_v2!=0;
		// boolean loss1= x_v1!=0 & y_v1==0;
		// boolean loss2= x_v2!=0 & y_v2==0;
		//
		// }
	    }

	}

	return result;
    }

    public static int[] count_gene_presence_absence(Map<String, Integer> map1, Node r)
    {

	Set<Node> allNodes = collectAllDescendants(r);

	int[] result = new int[2];

	for (Node n : allNodes)
	{

	    int gene1 = map1.get(n.getName()); // itt tobbnyire 0/1 érték van, de lehet 2,3... is
	    int present = (gene1 != 0 ? 1 : 0);
	    result[present]++;
	}

	return result;
    }

    /// ---------------
    public static int[] countLossAndGainEvents(Map<String, Integer> map1, Node r)
    {

	Set<Node> allNodes = collectAllDescendants(r);

	int[] result = new int[4];

	for (Node parent : allNodes)
	{
	    int parent_gene1 = map1.get(parent.getName()); // itt tobbnyire 0/1 érték van, de lehet 2,3... is
	    int present_parent_gene1 = (parent_gene1 != 0 ? 1 : 0);

	    for (Node child : parent.getChildren())
	    {

		int child_gene1 = map1.get(child.getName());
		int present_child_gene1 = (child_gene1 != 0 ? 1 : 0);

		int idx_gene1 = 2 * present_parent_gene1 + present_child_gene1;
		result[idx_gene1]++;
	    }

	}

	return result;
    }

    // ----

    public static int[][] count_loss_and_gain_eventses(Map<String, Integer> map1, Map<String, Integer> map2, Node r)
    {

	Set<Node> allNodes = collectAllDescendants(r);

	int[][] result = new int[4][4];

	for (Node parent : allNodes)
	{

	    for (Node child : parent.getChildren())
	    {

		int parent_gene1 = map1.get(parent.getName()); // itt tobbnyire 0/1 érték van, de lehet 2,3... is
		int parent_gene2 = map2.get(parent.getName());
		int child_gene1 = map1.get(child.getName());
		int child_gene2 = map2.get(child.getName());

		int present_parent_gene1 = (parent_gene1 != 0 ? 1 : 0);
		int present_parent_gene2 = (parent_gene2 != 0 ? 1 : 0);
		int present_child_gene1 = (child_gene1 != 0 ? 1 : 0);
		int present_child_gene2 = (child_gene2 != 0 ? 1 : 0);

		int idx_gene1 = 2 * present_parent_gene1 + present_child_gene1;
		int idx_gene2 = 2 * present_parent_gene2 + present_child_gene2;
		result[idx_gene1][idx_gene2]++;
	    }

	}

	return result;
    }

    public static int[][] count_loss_and_gain_eventses2(Map<Node, Integer> map1, Map<Node, Integer> map2, Node r)
    {

	Set<Node> allNodes = collectAllDescendants(r);

	int[][] result = new int[4][4];

	for (Node parent : allNodes)
	{

	    for (Node child : parent.getChildren())
	    {

		int parent_gene1 = map1.get(parent); // itt tobbnyire 0/1 érték van, de lehet 2,3... is
		int parent_gene2 = map2.get(parent);
		int child_gene1 = map1.get(child);
		int child_gene2 = map2.get(child);

		int present_parent_gene1 = (parent_gene1 != 0 ? 1 : 0);
		int present_parent_gene2 = (parent_gene2 != 0 ? 1 : 0);
		int present_child_gene1 = (child_gene1 != 0 ? 1 : 0);
		int present_child_gene2 = (child_gene2 != 0 ? 1 : 0);

		int idx_gene1 = 2 * present_parent_gene1 + present_child_gene1;
		int idx_gene2 = 2 * present_parent_gene2 + present_child_gene2;
		result[idx_gene1][idx_gene2]++;
	    }

	}

	return result;
    }

    /***
     * Kivalogatja azokat a node-okat ahol presence/absence/gain vagy loss esemeny
     * van.
     * 
     * (Az esemenyek valojaban két node kozotti elen vannak, de a gyerek nodehoz
     * rendeljuk. A fuggveny sem a beadott node-okbol valogat, hanem a gyerekeikbol,
     * de ugyis mindig a fa osszes node-jat adjuk oda, vagy legalabb egy lefele zart
     * reszfat)
     * 
     * @param allNodes
     * @param map1
     * @param presentParent
     * @param presentChild
     * @return
     */
    public static Set<Node> seletNodes(Set<Node> allNodes, Map<String, Integer> map1, boolean presentParent,
	    boolean presentChild)
    {

	Set<Node> resultSet = new LinkedHashSet<Node>();
	for (Node parent : allNodes)
	{

	    for (Node child : parent.getChildren())
	    {

		int parent_gene1 = map1.get(parent.getName()); // itt tobbnyire 0/1 érték van, de lehet 2,3... is
		int child_gene1 = map1.get(child.getName());

		boolean present_parent_gene1 = (parent_gene1 != 0);
		boolean present_child_gene1 = (child_gene1 != 0);

		if (present_child_gene1 == presentChild && present_parent_gene1 == presentParent)
		{
		    resultSet.add(child);
		}
	    }

	}

	return resultSet;
    }

    public static Set<String> getNameSetOfNodeSet(Set<Node> nodeSet)
    {

	Set<String> resultSet = new LinkedHashSet<String>(nodeSet.size());
	for (Node node : nodeSet)
	{
	    resultSet.add(node.getName());
	}

	return resultSet;
    }

    public static String[] seletNamesOfNodes(Set<Node> allNodes, Map<String, Integer> map1, boolean presentParent,
	    boolean presentChild)
    {
	Set<Node> s = seletNodes(allNodes, map1, presentParent, presentChild);
	Set<String> s2 = getNameSetOfNodeSet(s);
	return s2.toArray(new String[s2.size()]);
    }

    // csak kiprobalni, hogy az R-nek visszaadott tomboknel mi a sorrend.
    public static int[][] test1()
    {
	int[][] result = new int[4][4];
	result[0][0] = 33;
	result[0][3] = 2;
	result[2][0] = 1;
	return result;

    }

}
