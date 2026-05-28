package ger.maddison;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class NwkNode implements Node {

	private ArrayList<Node> children = new ArrayList<Node>(2);
	public final String name;
	public final double length;

	public NwkNode(List<Node> children, String name, double length) {
		this.children.addAll(children);
		this.name = name;
		this.length = length;
	}

	@Override
	public boolean isTerminal() {
		return children.isEmpty();
	}

	@Override
	public List<Node> getChildren() {

		return children;
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		if (!children.isEmpty()) {
			sb.append('(');
			for (Iterator iterator = children.iterator(); iterator.hasNext();) {
				Node node = (Node) iterator.next();
				sb.append(node.toString());
				sb.append(",");
			}
			sb.setLength(sb.length() - 1); // remove the last column
			sb.append(')');
		}

		sb.append(name);
//		sb.append(':');
//		sb.append(length);
		return sb.toString();
	}
	//-------------------------------------------------------------------------------------

	private static class Node_and_idx {
		NwkNode n;
		int idx;

		public Node_and_idx(NwkNode n, int idx) {
			super();
			this.n = n;
			this.idx = idx;
		}

	}

	
	public static NwkNode parseNWK(File f) throws IOException {

		BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(f)));
		try {
			String nwkString = in.readLine();
			NwkNode node1 = parseNWK(nwkString);

			// I test if the rest of the file is empty
			String s = in.readLine();
			while ("".equals(s)) {
				s = in.readLine();
			}
			if (s != null) {
				throw new IllegalArgumentException(
						"NWK file shall contains only 1 row. It contines more. " + f.getPath());
			}

			return node1;
		} finally {

			in.close();
		}

	}

	public static NwkNode parseNWK(String pattern) {
		pattern = pattern.subSequence(0, pattern.length()-1) + ")"; //TODO ez kis trukkozes, hogy a sor vegen levő pontosvesszöt lecserelem csuko zarojelre
		Node_and_idx parse1 = parseNWK(pattern, 0);
		if (parse1.idx != pattern.length() - 1) {
			throw new IllegalArgumentException();
		}
		return parse1.n;
	}

	public static Node_and_idx parseNWK(String pattern, int idx) {

		ArrayList<Node> childList = new ArrayList<Node>();

		char c1 = pattern.charAt(idx);

		if (c1 == '(') {
			while (true) {
				idx++;
				Node_and_idx sub1 = parseNWK(pattern, idx);
				childList.add(sub1.n);
				idx = sub1.idx;

				// idx++;
				c1 = pattern.charAt(idx);

				if (c1 == ')') {
					idx++;
					break;
				} else if (c1 == ',') {
					// do nothong, the while loop reads the next brother node
				} else {
					throw new IllegalStateException();
				}

			}

		}
		// elolvastuk a gyerkeket, most kell jöjjön a neve

		// if(c1!='(' & c1!=')' & c1!=',')
		StringBuilder nameAndLength = new StringBuilder();
		while (true) {
			c1 = pattern.charAt(idx);

			if (c1 == '(') {
				throw new IllegalStateException();
			} else if (c1 == ')' | c1 == ',') {
				break;
			}
			nameAndLength.append(c1);
			idx++;

		}

		String name;
		double length;
		int t1 = nameAndLength.toString().indexOf(':');
		if (t1 == -1) {
			name = nameAndLength.toString();
			length = Double.NaN;
		} else {
			name = nameAndLength.substring(0, t1);
			String lengthString = nameAndLength.substring(t1 + 1);
			length = Double.parseDouble(lengthString);
		}

		return new Node_and_idx(new NwkNode(childList, name, length), idx);

	}
	//-------------------------------------------------------------------------------------
	public static void main(String[] args) throws IOException {
	
//		 BufferedReader in = new BufferedReader(new InputStreamReader(new	 FileInputStream("data/mini.nwk")));
		BufferedReader in = new BufferedReader(
				new InputStreamReader(new FileInputStream("data/rapidnj_tree_nonegative_outgroupMRCARooted.nwk")));

		String nwkString = in.readLine();
		in.close();

		System.out.println(parseNWK(nwkString));
	}

}
