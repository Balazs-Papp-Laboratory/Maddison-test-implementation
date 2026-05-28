package ger.tree;

/***
 * Utilities for {@link DoubleArrayTree}
 */

public class DATUtils {

	// /**
	// * rekurzívan begyűjti a gyerekeit. A gyökeret és az összes alatta lévő
	// node-ot
	// * bepakolja egy Set-be.
	// *
	// * @param r
	// * @return
	// */
	// public static Set<Node> collectAllDescendants(Node r)
	// {
	// LinkedList<Node> gray = new LinkedList<Node>();
	// gray.push(r);
	// LinkedHashSet<Node> set1 = new LinkedHashSet<Node>();
	//
	// while (!gray.isEmpty())
	// {
	// Node x = gray.removeLast();
	// boolean notAlreadyContaind = set1.add(x);
	// if (notAlreadyContaind)
	// {
	// gray.addAll(x.getChildren());
	// }
	//
	// }
	// return set1;
	// }
	//
	// public static String[] collectNamesOfAllDescendants(Node r)
	// {
	// LinkedHashSet<String> nameSet1 = new LinkedHashSet<String>();
	//
	// Set<Node> set1 = collectAllDescendants(r);
	//
	// for (Node n : set1)
	// {
	// nameSet1.add(n.getName());
	// }
	// return nameSet1.toArray(new String[nameSet1.size()]);
	// }
	//
	// public static Set<Node> filterNodeSetByNames(Set<Node> nodeSet, Set<String>
	// nameSet)
	// {
	// LinkedHashSet<Node> set1 = new LinkedHashSet<Node>();
	//
	// for (Node n : nodeSet)
	// {
	// if (nameSet.contains(n.getName()))
	// {
	// set1.add(n);
	// }
	// }
	// return set1;
	// }

	////////////////////////////////
	//
	// public static int[][] ize(Map<String, Integer> map1, Map<String, Integer>
	// map2, Node r)
	// {
	//
	// System.out.println("IZE invoked");
	// int[][] result = new int[2][2];
	//
	// LinkedList<Node> gray = new LinkedList<Node>();
	// gray.push(r);
	// LinkedHashSet<Node> set1 = new LinkedHashSet<Node>();
	//
	// while (!gray.isEmpty())
	// {
	// Node x = gray.removeLast();
	// boolean notAlreadyContaind = set1.add(x);
	// if (notAlreadyContaind)
	// {
	// gray.addAll(x.getChildren());
	//
	// // Integer x_v1_obj = map1.get(x.getName());
	// // Integer x_v2_obj = map2.get(x.getName());
	// //
	// // if (x_v1_obj == null || x_v2_obj == null) {
	// // System.out.println("ERROR: '" + x.getName() + "' nincs a map-ben.");
	// // } else {
	//
	// int x_v1 = map1.get(x.getName());
	// int x_v2 = map2.get(x.getName());
	//
	// int present1 = (x_v1 != 0 ? 1 : 0);
	// int present2 = (x_v2 != 0 ? 1 : 0);
	// result[present1][present2]++;
	// // for(Node y :x.getChildren())
	// // {
	// // int y_v1 = map1.get(y.getName());
	// // int y_v2 = map2.get(y.getName());
	// //
	// //
	// // boolean gain1= x_v1==0 & y_v1!=0;
	// // boolean gain2= x_v2==0 & y_v2!=0;
	// // boolean loss1= x_v1!=0 & y_v1==0;
	// // boolean loss2= x_v2!=0 & y_v2==0;
	// //
	// // }
	// }
	//
	// }
	//
	// return result;
	// }

	public static int[] count_gene_presence_absence(int[] geneValues, DoubleArrayTree tree, int root) {

		int[] result = new int[2];

		IntegerStack gray = new IntegerStack(tree.size());
		gray.push(root);

		while (!gray.isEmpty()) {

			int x = gray.pop();

			int present = (geneValues[x] != 0 ? 1 : 0);
			result[present]++;

			gray.pushChildren(tree, x);

		}

		return result;
	}

	/// ---------------

	/**
	 * @return int[4] the order of counts: absent, gain,loss, present
	 */
	public static int[] count_loss_and_gain_eventses(int[] map1, DoubleArrayTree tree, int root) {

		int[] result = new int[4];

		IntegerStack gray = new IntegerStack(tree.size());
		gray.pushChildren(tree, root);

		while (!gray.isEmpty()) {

			int x = gray.pop();
			int p = tree.getParent(x);

			int parent_gene1 = map1[p]; // itt tobbnyire 0/1 érték van, de lehet 2,3... is
			int child_gene1 = map1[x];

			int present_parent_gene1 = (parent_gene1 != 0 ? 1 : 0);
			int present_child_gene1 = (child_gene1 != 0 ? 1 : 0);

			// order: absent, gain,loss, present
			int idx_gene1 = 2 * present_parent_gene1 + present_child_gene1;
			result[idx_gene1]++;

			gray.pushChildren(tree, x);
		}

		return result;
	}

	/**
	 * @return int[4] the order of counts: absent, gain,loss, present
	 */
	public static int[] count_loss_and_gain_eventses_in_black_area(int[] map1, DoubleArrayTree tree,
			boolean[] blackNodeIdx) {

		int[] result = new int[4];

		assert (tree.getRoot() == 0);
		for (int x = 1; x < map1.length; x++) {
			if (blackNodeIdx[x]) {
				int p = tree.getParent(x);

				// itt tobbnyire 0/1 érték van, de lehet 2,3... is
				int present_parent_gene1 = (map1[p] != 0 ? 1 : 0);
				int present_child_gene1 = (map1[x] != 0 ? 1 : 0);

				// order: absent, gain,loss, present
				int idx_gene1 = 2 * present_parent_gene1 + present_child_gene1;
				result[idx_gene1]++;
			}
		}

		return result;
	}

	// ----

	public static int[] count_loss_and_gain_eventses(int[] map1, int[] map2, DoubleArrayTree tree, int root) {

		int[] result = new int[4 * 4];

		IntegerStack gray = new IntegerStack(tree.size());
		gray.pushChildren(tree, root);

		while (!gray.isEmpty()) {

			int x = gray.pop();
			int p = tree.getParent(x);

			int idx_gene1 = eventIndex(map1, p, x);
			int idx_gene2 = eventIndex(map2, p, x);
			int idxCombined = idx_gene1 + 4 * idx_gene2;
			result[idxCombined]++;

			gray.pushChildren(tree, x);
		}

		return result;
	}

	static final int ABSENT = 0;
	static final int GAIN = 1;
	static final int LOSS = 2;
	static final int PRESENT = 3;

	/** az idx lehetseges ertekei ABSENT ,GAIN, LOSS, PRESENT [0,1,2,3] */
	private static final int eventIndex(int[] map, int parent, int child) {
		int present_parent_gene1 = (map[parent] != 0 ? 1 : 0);
		int present_child_gene1 = (map[child] != 0 ? 1 : 0);

		int idx = 2 * present_parent_gene1 + present_child_gene1;
		return idx;

	}

	/***
	 * Kivalogatja azokat a node-okat ahol presence/absence/gain vagy loss esemeny
	 * van.
	 * 
	 * Visszaad egy boolean[]-ot ahol mindennodehoz hozzárendeli, hogy a kért
	 * eseméynt tortent-e rajta.
	 * 
	 * Az esemenyek valojaban két node kozotti elen vannak, de a gyerek nodehoz
	 * rendeljuk. A root node index poziciojan mindig false lesz.
	 * 
	 * @return boolean[]
	 */
	public static boolean[] selectNodes(int[] geneValue, DoubleArrayTree tree, int parentValues, int childValue) {
		assert geneValue.length == tree.size();

		boolean[] result = new boolean[geneValue.length];

		IntegerStack gray = new IntegerStack(tree.size());
		gray.pushChildren(tree, tree.getRoot());

		while (!gray.isEmpty()) {

			int x = gray.pop();
			int p = tree.getParent(x);

			result[x] = geneValue[x] == childValue && geneValue[p] == parentValues;

			gray.pushChildren(tree, x);
		}

		return result;
	}

	public static boolean[] nagateBooleanArray(boolean[] arr) {

		boolean[] result = new boolean[arr.length];

		for (int i = 0; i < arr.length; i++) {
			result[i] = !arr[i];
		}
		return result;
	}

	public static boolean[] orOperatorOfBooleanArrays(boolean[] arrA, boolean[] arrB) {
		if (arrA.length != arrB.length) {
			throw new RuntimeException();
		}

		boolean[] result = new boolean[arrA.length];

		for (int i = 0; i < arrA.length; i++) {
			result[i] = arrA[i] || arrB[i];
		}
		return result;
	}

	public static int countTrueValues(boolean[] arr) {
		int cnt = 0;
		for (int i = 0; i < arr.length; i++) {
			if (arr[i]) {
				cnt++;
			}
		}
		return cnt;
	}

	// csak kiprobalni, hogy az R-nek visszaadott tomboknel mi a sorrend.
	public static int[][] test1() {
		int[][] result = new int[4][4];
		result[0][0] = 33;
		result[0][3] = 2;
		result[2][0] = 1;
		return result;

	}

}
